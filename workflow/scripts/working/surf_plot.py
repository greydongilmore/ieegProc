#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 21:18:35 2022

@author: greydon
"""

from nilearn import plotting, surface
import nibabel as nb
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from matplotlib.collections import PolyCollection
from pathlib import Path
from nipype.interfaces.base import File, traits
from mpl_toolkits import mplot3d
from collections import namedtuple
import matplotlib.colors as mcolors

def normal_vectors(vertices, faces):
	tris = vertices[faces]
	n = np.cross(tris[::, 1] - tris[::, 0], tris[::, 2] - tris[::, 0])
	n = normalize_v3(n)
	return n

def normalize_v3(arr):
	"""Normalize a numpy array of 3 component vectors shape=(n,3)"""
	lens = np.sqrt(arr[:, 0] ** 2 + arr[:, 1] ** 2 + arr[:, 2] ** 2)
	arr[:, 0] /= lens +  1e-12 
	arr[:, 1] /= lens +  1e-12 
	arr[:, 2] /= lens +  1e-12 
	return arr

def frustum(left, right, bottom, top, znear, zfar):
	M = np.zeros((4, 4), dtype=np.float32)
	M[0, 0] = +2.0 * znear / (right - left)
	M[1, 1] = +2.0 * znear / (top - bottom)
	M[2, 2] = -(zfar + znear) / (zfar - znear)
	M[0, 2] = (right + left) / (right - left)
	M[2, 1] = (top + bottom) / (top - bottom)
	M[2, 3] = -2.0 * znear * zfar / (zfar - znear)
	M[3, 2] = -1.0
	return M


def perspective(fovy, aspect, znear, zfar):
	h = np.tan(0.5 * np.radians(fovy)) * znear
	w = h * aspect
	return frustum(-w, w, -h, h, znear, zfar)

def translate(x, y, z):
	return np.array(
		[[1, 0, 0, x], [0, 1, 0, y], [0, 0, 1, z], [0, 0, 0, 1]], dtype=float
	)

def xrotate(theta):
	t = np.pi * theta / 180
	c, s = np.cos(t), np.sin(t)
	return np.array(
		[[1, 0, 0, 0], [0, c, -s, 0], [0, s, c, 0], [0, 0, 0, 1]], dtype=float
	)

def yrotate(theta):
	t = np.pi * theta / 180
	c, s = np.cos(t), np.sin(t)
	return np.array(
		[[c, 0, s, 0], [0, 1, 0, 0], [-s, 0, c, 0], [0, 0, 0, 1]], dtype=float
	)


def surf_from_cifti(ts_cifti,struct_names=['CIFTI_STRUCTURE_CORTEX_LEFT','CIFTI_STRUCTURE_CORTEX_RIGHT']):
	"""
	Gets the time series of cortical surface vertices (Left and Right)
	Args:
		ts_cifti (cifti obj) - cifti object of time series
	Returns:
		cii (cifti object) - contains the time series for the cortex
	"""
	# get brain axis models
	bmf = ts_cifti.header.get_axis(1)
	# print(dir(bmf))
	# get the data array with all the time points, all the structures
	ts_array = ts_cifti.get_fdata()
	ts_list = []
	for idx, (nam,slc,bm) in enumerate(bmf.iter_structures()):
		# just get the cortical surfaces
		if nam in struct_names: 
			values = np.full((ts_array.shape[0],bmf.nvertices[nam]),np.nan)
			# get the values corresponding to the brain model
			values[:,bm.vertex] = ts_array[:, slc]
			ts_list.append(values)
		else:
			break
	return ts_list


views = traits.List(
	[{
		"view": "lateral",
		"hemi": "left"
	}, {
		"view": "medial",
		"hemi": "left"
	}, {
		"view": "lateral",
		"hemi": "right"
	}, {
		"view": "medial",
		"hemi": "right"
	}],
	usedefault=True,
	desc="List of dictionaries describing views "
	" to display per map",
	inner_traits=traits.Dict(
		key_trait=traits.Enum(values=["view", "hemi"]),
		value_trait=traits.Enum(
			values=["lateral", "medial", "dorsal", "ventral"])))

Hemispheres = namedtuple("Hemispheres", ["left", "right"])

darkness = traits.Float(
	0.3,
	usedefault=True,
	desc="Multiplicative factor of bg_img onto foreground map",
	mandatory=False)

visualize_all_maps = traits.Bool(False,
				  usedefault=True,
				  desc="Visualize all mappings in "
				  "mapping file, if false will visualize "
				  "only the first mapping")

zero_nan = traits.Bool(False,
		      usedefault=True,
		      desc="Display NaNs as zeros")

colormap = traits.String("magma",
			 usedefault=True,
			 desc="Colormap to use to plot mapping",
			 mandatory=False)

#%%


debug = True
if debug:
    class dotdict(dict):
        """dot.notation access to dictionary attributes"""
        __getattr__ = dict.get
        __setattr__ = dict.__setitem__
        __delattr__ = dict.__delitem__

    class Namespace:
        def __init__(self, **kwargs):
            self.__dict__.update(kwargs)

    isub = 'sub-P093'
    data_dir = r'/home/greydon/Documents/data/SEEG/derivatives'

    input = dotdict({'t1': f'{data_dir}/atlasreg/{isub}/{isub}_desc-n4_T1w.nii.gz',
                    'pet': f'{data_dir}/atlasreg/{isub}/{isub}_space-T1w_desc-rigid_pet.nii.gz',
                     'gm': f'{data_dir}/atlasreg/{isub}/{isub}_label-GM_desc-atropos3seg_probseg.nii.gz',
                     'wm': f'{data_dir}/atlasreg/{isub}/{isub}_label-GM_desc-atropos3seg_probseg.nii.gz',
                     'mask': f'{data_dir}/atlasreg/{isub}/{isub}_desc-brain_from-MNI152NLin2009cSym_reg-affine_mask.nii.gz',
                     'seg': f'{data_dir}/atlasreg/{isub}/{isub}_desc-atroposKseg_dseg.nii.gz',
                     })

    output = dotdict({
        't1_grad': f'{data_dir}/atlasreg/{isub}/{isub}_desc-magnitude_T1w.nii.gz',
        't1_intensity': f'{data_dir}/atlasreg/{isub}/{isub}_desc-intensity_T1w.nii.gz',
        'png_mask': f'{data_dir}/atlasreg/{isub}/qc/{isub}_desc-brain_maskqc.png',
        'png_seg': f'{data_dir}/atlasreg/{isub}/qc/{isub}_desc-segmentation_segqc.png'
    })

    snakemake = Namespace(output=output, input=input)


lh_pial_file = f'{data_dir}/fastsurfer/{isub}/fastsurfer/surf/lh.pial'
rh_pial_file = f'{data_dir}/fastsurfer/{isub}/fastsurfer/surf/rh.pial'
lh_sulc_file = f'{data_dir}/fastsurfer/{isub}/fastsurfer/surf/lh.sulc'
rh_sulc_file = f'{data_dir}/fastsurfer/{isub}/fastsurfer/surf/rh.sulc'
lh_infl_file = f'{data_dir}/fastsurfer/{isub}/fastsurfer/surf/lh.inflated'
rh_infl_file = f'{data_dir}/fastsurfer/{isub}/fastsurfer/surf/rh.inflated'


lh_data = nb.freesurfer.read_geometry(lh_pial_file)
rh_data = nb.freesurfer.read_geometry(rh_pial_file)

lh_sulc_data = nb.freesurfer.read_morph_data(lh_sulc_file)
rh_sulc_data = nb.freesurfer.read_morph_data(rh_sulc_file)

lh_infl_data = nb.freesurfer.read_geometry(lh_infl_file)
rh_infl_data = nb.freesurfer.read_geometry(rh_infl_file)



gii_meshes=[]
for isurf in [lh_pial_file,rh_pial_file]:
	surf = Path(isurf)
	out = surf.with_suffix('.surf.gii')
	
	vertices, faces = nb.freesurfer.read_geometry(isurf)
	vertices = nb.gifti.GiftiDataArray(vertices, 'NIFTI_INTENT_POINTSET',
										datatype='NIFTI_TYPE_FLOAT32')
	faces = nb.gifti.GiftiDataArray(faces, 'NIFTI_INTENT_TRIANGLE',
									 datatype='NIFTI_TYPE_INT32')
	gii = nb.gifti.GiftiImage(darrays=[vertices, faces])
	
	nb.save(gii, out)
	
	coords, faces = nb.load(out).get_arrays_from_intent(nb.nifti1.intent_codes['NIFTI_INTENT_POINTSET'])[0].data, \
		 nb.load(out).get_arrays_from_intent(nb.nifti1.intent_codes['NIFTI_INTENT_TRIANGLE'])[0].data
	
	gii_meshes.append([coords, faces, None])

map_hemi = Hemispheres(left=(gii_meshes[0][0],gii_meshes[0][1],None), right=(gii_meshes[1][0],gii_meshes[1][1],None))
bg_hemi = Hemispheres(left=None, right=None)


num_views =[("lateral","left"),("medial","right"),("lateral","left"),("medial","right")]
num_maps = 1
vmin, vmax = None, None


surf_mesh=[np.array(coords), np.array(faces)]


w, h = plt.figaspect(num_maps / (len(num_views)))
fig, axs = plt.subplots(num_maps,
						len(num_views),
						subplot_kw={'projection': '3d'},
						figsize=(w, h))
fig.set_facecolor("black")
fig.tight_layout()


for i, a in enumerate(axs.flat):
	a.set_facecolor("black")

	view_ind = i % len(num_views)
	map_ind = i // len(num_views)

	view =num_views[view_ind][0]
	hemi =num_views[view_ind][1]

	display_map = getattr(map_hemi, hemi)
	display_bg = getattr(bg_hemi, hemi)

	v, t, m = display_map
	m = v[map_ind]
	if zero_nan:
		m[np.isnan(m)] = 0

	# Plot
	plotting.plot_surf([v, t],
					surf_map=None,
					#bg_map=None,
					cmap=colormap,
					axes=a,
					hemi=hemi,
					bg_on_data=True,
					vmin=vmin,
					vmax=vmax)

plt.draw()
plt.savefig(self._out_report)




plotting.plot_surf(lh_infl_file,surf_mesh,
							hemi='left',view='medial', colorbar=True, cmap='hot', title='gifti left')


all_ver = np.concatenate([lh_data[0], rh_data[0]], axis=0)
all_face = np.concatenate([lh_data[1], rh_data[1]+lh_data[0].shape[0]], axis=0)
surf_mesh = [all_ver, all_face]
bg_map = np.concatenate((lh_sulc_data, rh_sulc_data))


fig = plotting.view_surf(surf_mesh, bg_map, symmetric_cmap=False)
fig.open_in_browser()

#plot PET surface
surf_over = surface.vol_to_surf(snakemake.input.pet, surf_mesh)
surf_over = np.nan_to_num(surf_over)
fig = plotting.view_surf(surf_mesh, surf_over, bg_map=bg_map,symmetric_cmap=False)
fig.open_in_browser()


fig = plotting.plot_surf(surf_mesh, surf_over, bg_map=bg_map,
                         view='anterior', engine='plotly',
                         ouput_file="/home/greydon/Downloads/text.html", alpha=1,
                         bg_on_data=True, darkness=0.8, title=None, symmetric_cmap=True)

fig.show()
fig.figure.write_html("/home/greydon/Downloads/sample_surface.html")

if len(rh_sulc_data) < len(lh_sulc_data):
    lh_sulc_data = lh_sulc_data[:len(rh_sulc_data)]


surf_mesh_dict = {}
surf_mesh_dict['pial_left'] = lh_data
surf_mesh_dict['pial_right'] = rh_data
surf_mesh_dict['sulc_left'] = lh_sulc_data
surf_mesh_dict['sulc_right'] = rh_sulc_data
surf_mesh_dict['infl_left'] = lh_infl_data
surf_mesh_dict['infl_right'] = rh_infl_data

from brain4views import plot_surf4

lh_over = surface.vol_to_surf(snakemake.input.pet, lh_pial_file, radius=1, kind="ball")
rh_over = surface.vol_to_surf(snakemake.input.pet, rh_pial_file, radius=1, kind="ball")

plot_surf4(
    [lh_pial_file, rh_pial_file],
    #overlays=[lh_over, rh_over],
    avg_method="mean",
    colorbar=True,
    output_file="/home/greydon/Downloads/test.png",
	dpi=350
)

plot_surf4(
	[lh_pial_file, rh_pial_file],
	#overlays=[lh_over, rh_over],
	vmin=-1.2,
	vmax=1.2,
	threshold=0.3,
	cmap="RdBu_r",
	avg_method="mean",
	title="Correlation (Z)",
	colorbar=True,
	output_file="/home/greydon/Downloads.test.png",
)



# initiate figure
fig = plt.figure(figsize=(8, 6))


overlays=[lh_over,rh_over]
sulc_maps=[lh_sulc_file,rh_sulc_file]

rotations_both = [[90, 270], [270, 90]]
bg_map = np.concatenate((overlays))

overlay_faces_final=[]
label_colors=None

for m, mesh in enumerate([lh_pial_file,rh_pial_file]):
	vertices, faces = surface.load_surf_mesh(mesh)
	vertices = vertices.astype(np.float64)
	faces = faces.astype(int)
	
	# Set up lighting, intensity and shading
	rotations = rotations_both[m]
	vert_range = max(vertices.max(0) - vertices.min(0))
	vertices = (vertices - (vertices.max(0) + vertices.min(0)) / 2) / vert_range
	face_normals = normal_vectors(vertices, faces)
	light = np.array([0, 0, 1])
	intensity = np.dot(face_normals, light)
	shading = 0.7  # shading 0-1. 0=none. 1=full
	# top 20% all become fully colored
	denom = np.percentile(intensity, 80) - np.min(intensity)
	intensity = (1 - shading) + shading * (intensity - np.min(intensity)) / denom
	intensity[intensity > 1] = 1
	
	# initiate array for face colors
	face_colors = np.ones((faces.shape[0], 4))
	
	mask = np.zeros(vertices.shape[0]).astype(bool)
	
	sulc = surface.load_surf_data(sulc_maps[m])
	assert sulc.shape[0] == vertices.shape[0]
	
	sulc_faces = np.mean(sulc[faces], axis=1)
	
	if sulc_faces.min() != sulc_faces.max():
		neg_sulc = np.where(sulc_faces <= 0)
		pos_sulc = np.where(sulc_faces > 0)
		sulc_faces[neg_sulc] = 0
		sulc_faces[pos_sulc] = 1
	
	label_masks = []
	label_mask_faces = [np.median(L[faces], axis=1) for L in label_masks]
	
	
	overlay = surface.load_surf_data(overlays[m])
	
	cmap = plt.cm.get_cmap(plt.rcParamsDefault["image.cmap"])
	greys = plt.get_cmap("Greys", 512)
	greys_narrow = ListedColormap(greys(np.linspace(0.42, 0.58, 256)))
	face_colors = greys_narrow(sulc_faces)
	
	
	kept_indices = np.arange(sulc_faces.shape[0])

	
	overlay_faces = np.mean(lh_over[faces.astype(int)], axis=1)
	overlay_faces_final.append(overlay_faces)
	
	vmin = np.nanmin(overlay_faces)
	vmax = np.nanmax(overlay_faces)
	overlay_faces = overlay_faces - vmin
	overlay_faces = overlay_faces / (vmax - vmin)
	face_colors[kept_indices] = cmap(overlay_faces[kept_indices])
	
	# assign label faces to appropriate color
	for i, L in enumerate(label_mask_faces):
		L_idx = np.where(L >= 0.5)
		# blend (multiply) label color with underlying color
		face_colors[L_idx] = face_colors[L_idx] * label_colors[i]

	face_colors[:, 0] *= intensity
	face_colors[:, 1] *= intensity
	face_colors[:, 2] *= intensity

	for i, view in enumerate(rotations):
		MVP = (
			perspective(25, 1, 1, 100)
			@ translate(0, 0, -3)
			@ yrotate(view)
			@ xrotate(270)
		)
		# translate coordinates based on viewing position
		V = np.c_[vertices, np.ones(len(vertices))] @ MVP.T
		V /= V[:, 3].reshape(-1, 1)
		V = V[faces]
		# triangle coordinates
		T = V[:, :, :2]
		# get Z values for ordering triangle plotting
		Z = -V[:, :, 2].mean(axis=1)
		# sort the triangles based on their z coordinate
		Zorder = np.argsort(Z)
		T, C = T[Zorder, :], face_colors[Zorder, :]
		# add subplot and plot PolyCollection
		ax = fig.add_subplot(
			2,
			2,
			m + 2 * i + 1,
			xlim=[-1, +1],
			ylim=[-0.6, +0.6],
			frameon=False,
			aspect=1,
			xticks=[],
			yticks=[],
		)
		collection = PolyCollection(
			T, closed=True, antialiased=False, facecolor=C, edgecolor=C, linewidth=0
		)
		collection.set_alpha(1)
		ax.add_collection(collection)

from matplotlib.colors import Normalize, to_rgba_array
from matplotlib.cm import ScalarMappable
our_cmap = plt.get_cmap(cmap)
norm = Normalize(vmin=vmin, vmax=vmax)
bounds = np.linspace(vmin, vmax, our_cmap.N)

ticks = [vmin, vmax]


# we need to create a proxy mappable
proxy_mappable = ScalarMappable(cmap=our_cmap, norm=norm)
proxy_mappable.set_array(overlay_faces)
cax = plt.axes([0.38, 0.466, 0.24, 0.024])
plt.colorbar(
	proxy_mappable,
	cax=cax,
	boundaries=bounds,
	ticks=ticks,
	orientation="horizontal",
)

fig.text(0.25, 0.975, "Left", ha="center", va="top")
fig.text(0.75, 0.975, "Right", ha="center", va="top")
fig.text(0.025, 0.75, "Lateral", ha="left", va="center", rotation=90)
fig.text(0.025, 0.25, "Medial", ha="left", va="center", rotation=90)
fig.subplots_adjust(left=0, right=1, top=1, bottom=0, wspace=0, hspace=0)


pet_surf = surface.vol_to_surf(snakemake.input.pet, lh_pial_file)


fig = plotting.view_img_on_surf(snakemake.input.pet,surf_mesh_dict)
fig.open_in_browser()




plotting.view_surf(snakemake.input.pet, surf_mesh_dict, bg_map=pet_surf)

fig = plotting.plot_surf(surf_mesh, surf_map=bg_map, bg_map=test,
                         view='anterior', cmap=plt.cm.Greys_r, engine='plotly',
                         ouput_file="/home/greydon/Downloads/text.html", alpha=1,
                         bg_on_data=True, darkness=0.8, title=None, symmetric_cmap=True)

fig.show()
