#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 21:57:30 2022

@author: greydon
"""

import os
import nibabel as nb
import numpy as np
import pandas as pd
import regex as re
import matplotlib.pyplot as plt
from nilearn import plotting,surface
from nilearn.plotting.displays import PlotlySurfaceFigure
import plotly.graph_objs as go

AXIS_CONFIG = {
    "showgrid": False,
    "showline": False,
    "ticks": "",
    "title": "",
    "showticklabels": False,
    "zeroline": False,
    "showspikes": False,
    "spikesides": False,
    "showbackground": False,
}

LAYOUT = {
	"scene": {f"{dim}axis": AXIS_CONFIG for dim in ("x", "y", "z")},
	"paper_bgcolor": "#fff",
	"hovermode": False,
	"showlegend":True,
	"legend":{
		"itemsizing": "constant",
		"groupclick":"togglegroup",
		"yanchor":"top",
		"y":0.8,
		"xanchor":"left",
		"x":0.05,
		"title_font_family":"Times New Roman",
		"font":{
			"size":20
		},
		"bordercolor":"Black",
		"borderwidth":1
	},
	"margin": {"l": 0, "r": 0, "b": 0, "t": 0, "pad": 0},
}

CAMERAS = {
    "left": {
        "eye": {"x": -1.5, "y": 0, "z": 0},
        "up": {"x": 0, "y": 0, "z": 1},
        "center": {"x": 0, "y": 0, "z": 0},
    },
    "right": {
        "eye": {"x": 1.5, "y": 0, "z": 0},
        "up": {"x": 0, "y": 0, "z": 1},
        "center": {"x": 0, "y": 0, "z": 0},
    },
    "dorsal": {
        "eye": {"x": 0, "y": 0, "z": 1.5},
        "up": {"x": 0, "y": 1, "z": 0},
        "center": {"x": 0, "y": 0, "z": 0},
    },
    "ventral": {
        "eye": {"x": 0, "y": 0, "z": -1.5},
        "up": {"x": 0, "y": 1, "z": 0},
        "center": {"x": 0, "y": 0, "z": 0},
    },
    "anterior": {
        "eye": {"x": 0, "y": 1.5, "z": 0},
        "up": {"x": 0, "y": 0, "z": 1},
        "center": {"x": 0, "y": 0, "z": 0},
    },
    "posterior": {
        "eye": {"x": 0, "y": -1.5, "z": 0},
        "up": {"x": 0, "y": 0, "z": 1},
        "center": {"x": 0, "y": 0, "z": 0},
    },
}

lighting_effects = dict(ambient=0.4, diffuse=0.5, roughness = 0.9, specular=0.6, fresnel=0.2)

def determine_groups(iterable, numbered_labels=False):
	values = []
	for item in iterable:
		temp=None
		if re.findall(r"([a-zA-Z]+)([0-9]+)([a-zA-Z]+)", item):
			temp = "".join(list(re.findall(r"([a-zA-Z]+)([0-9]+)([a-zA-Z]+)", item)[0]))
		elif '-' in item:
			temp=item.split('-')[0]
		else:
			if numbered_labels:
				temp=''.join([x for x in item if not x.isdigit()])
				for sub in ("T1","T2"):
					if sub in item:
						temp=item.split(sub)[0] + sub
			else:
				temp=item
		if temp is None:
			temp=item
		
		values.append(temp)
	
	vals,indexes,count = np.unique(values, return_index=True, return_counts=True)
	values_unique = [values[index] for index in sorted(indexes)]
	
	return values_unique,count



def remove_nan(nifti_in):
	"""remove_nan is a method needed after a registration performed by spmregister : insted of filling space with 0, nan
	are used to extend the PET space. We propose to replace them with 0s.
	Args:
		(string) volname : path to the Nifti volume where NaNs need to be replaced by 0s
	Returns:
		(string) Path to the volume in Nifti that does not contain any NaNs
	"""
	
	# Load the volume and get the data
	data = np.nan_to_num(nifti_in.get_fdata(dtype="float32"))
	
	# Now create final image (using header of original image), and save it in current directory
	nifti_out = nb.Nifti1Image(data, nifti_in.affine, header=nifti_in.header)
	
	return nifti_out

hemi = ["lh", "rh"]
surf_suffix = ["pial", "white", "inflated"]

def readRegMatrix(trsfPath):
	with open(trsfPath) as (f):
		return np.loadtxt(f.readlines())

def resample_img(img,voxel_dims=[1,1,1],target_shape=None):
	
	target_affine = img.affine.copy()
	
	# Decompose the image affine to allow scaling
	u,s,v = np.linalg.svd(target_affine[:3,:3],full_matrices=False)
	
	# Calculate the translation part of the affine
	spatial_dimensions = (img.header['dim'] * img.header['pixdim'])[1:4]
	resize_factor = spatial_dimensions / voxel_dims
	real_resize_factor = np.round(img.shape * resize_factor) / img.shape
	new_spacing = spatial_dimensions / real_resize_factor
	
	# Reconstruct the affine
	target_affine[:3,:3] = u @ np.diag(voxel_dims) @ v
	resampled_img = image.resample_img(img, target_affine=target_affine, target_shape=target_shape,interpolation='nearest')
	return resampled_img

#%%

from nilearn import image
from mne.transforms import apply_trans

isub = 'sub-P022'
data_dir = r'/home/greydon/Documents/data/clinical'
out_dir = f'{data_dir}/derivatives/petsurfer/{isub}'
gifti_out = os.path.join(data_dir,'derivatives','gifti',isub,"surf")
html_out = f"{data_dir}/derivatives/atlasreg/{isub}/qc/{isub}_space-native_desc-affine_electrodes.html"

if not os.path.exists(out_dir):
	os.makedirs(out_dir)


coords_fcsv=f'{data_dir}/derivatives/seega_coordinates/{isub}/{isub}_space-native_SEEGA.tsv'
xfm_noncontrast=f'{data_dir}/derivatives/atlasreg/{isub}/{isub}_acq-noncontrast_desc-rigid_from-noncontrast_to-contrast_type-ras_xfm.txt'


t1_fname = f'{data_dir}/derivatives/fastsurfer/{isub}/mri/orig.mgz'
lh_pial_file = f'{data_dir}/derivatives/fastsurfer/{isub}/surf/lh.pial'
rh_pial_file = f'{data_dir}/derivatives/fastsurfer/{isub}/surf/rh.pial'
lh_sulc_file = f'{data_dir}/derivatives/fastsurfer/{isub}/surf/lh.sulc'
rh_sulc_file = f'{data_dir}/derivatives/fastsurfer/{isub}/surf/rh.sulc'

t1_obj = nb.load(t1_fname)
Torig = t1_obj.header.get_vox2ras_tkr()
fs_transform=(t1_obj.affine-Torig)+np.eye(4)
t1_transform=readRegMatrix(xfm_noncontrast)

verl,facel=nb.freesurfer.read_geometry(lh_pial_file)
verr,facer=nb.freesurfer.read_geometry(rh_pial_file)

all_ver = np.concatenate([verl, verr], axis=0)
all_face = np.concatenate([facel, facer+verl.shape[0]], axis=0)
surf_mesh = [all_ver, all_face]

all_ver_shift=(apply_trans(fs_transform, all_ver))
all_ver_shift=(apply_trans(np.linalg.inv(t1_transform), all_ver_shift))


lh_sulc_data = nb.freesurfer.read_morph_data(lh_sulc_file)
rh_sulc_data = nb.freesurfer.read_morph_data(rh_sulc_file)
bg_map = np.concatenate((lh_sulc_data, rh_sulc_data))


mesh_3d = go.Mesh3d(x=all_ver_shift[:,0], y=all_ver_shift[:,1], z=all_ver_shift[:,2], i=all_face[:,0], j=all_face[:,1], k=all_face[:,2],opacity=.1,color='grey')
					lighting=lighting_effects)


value=np.arange(.1,.6,.05)

if os.path.exists(coords_fcsv):
	df = pd.read_table(coords_fcsv,sep='\t',header=0)
	if coords_fcsv.endswith("SEEGA.fcsv"):
		groups,n_members = determine_groups(df['label'].tolist(), True)
	else:
		groups,n_members = determine_groups(df['label'].tolist())
	
	df['group']=np.repeat(groups,n_members)
	
	cmap = plt.get_cmap('rainbow')
	colors=np.repeat(cmap(np.linspace(0, 1, len(groups))), n_members, axis=0)
	
	data=[mesh_3d]
	for igroup in groups:
		idx = [i for i,x in enumerate(df['label'].tolist()) if igroup in x]
		data.append(go.Scatter3d(
			x = df['x'][idx].values,
			y = df['y'][idx].values,
			z = df['z'][idx].values,
			name=igroup,
			mode = "markers+text",
			text=df['label'][idx].tolist(),
			textfont=dict(
				family="sans serif",
				size=16,
				color="black"
			),
			textposition = "top center",
			marker=dict(
				size=5,
				line=dict(
					width=1,
				),
				color=['rgb({},{},{})'.format(int(r*256),int(g*256),int(b*256)) for r,g,b,h in colors[idx]],
				opacity=1
				)))
	
	fig = go.Figure(data=data)
	fig.update_layout(scene_camera=CAMERAS['left'],
					  legend_title_text="Electrodes",
					  **LAYOUT)
	
	steps = []
	for i in range(len(value)):
		step = dict(
			label = str(f"{value[i]:.2f}"),
			method="restyle",
			args=['opacity', [value[i]]+(len(data)-1)*[1]],
		)
		steps.append(step)
	
	sliders = [dict(
		currentvalue={"visible": True,"prefix": "Opacity: ","font":{"size":16}},
		active=0,
		steps=steps,
		x=.35,y=.1,len=.3,
		pad={"t": 8},
	)]
	
	fig.update_layout(sliders=sliders)
	fig.show('browser')
	
	fig.write_html(os.path.splitext(html_out)[0]+".png")
	


#%%

pet_path = f'{data_dir}/derivatives/atlasreg/{isub}/{isub}_pet.nii.gz'
pet_reg = f'{data_dir}/derivatives/atlasreg/{isub}/{isub}_space-T1w_desc-rigid_pet.nii.gz'
pet_xfm = f'{data_dir}/derivatives/atlasreg/{isub}/{isub}_desc-rigid_from-pet_to-T1w_type-ras_xfm.txt'
suvr_filename = os.path.join(out_dir,f"{isub}_desc-suvr_pet.nii.gz")
mask = f'{data_dir}/derivatives/atlasreg/{isub}/{isub}_desc-brain_from-MNI152NLin2009cSym_reg-affine_mask.nii.gz'

if not os.path.exists(suvr_filename):
	# Load PET data (they must be in gtmsegspace, or same space as label file)
	pet = nb.load(pet_path)
	pet = remove_nan(pet)
	
	pet_reg_nifti = nb.load(pet_reg)
	pet_reg_nifti = remove_nan(pet_reg_nifti)
	pet_reg_data=pet_reg_nifti.get_fdata(dtype="float32")
	
	# Load mask
	eroded_mask_nifti = nb.load(mask)
	eroded_mask = eroded_mask_nifti.get_fdata(dtype="float32").copy()
	eroded_mask = eroded_mask > 0
	
	# check that eroded mask is not null
	mask_size = sum(sum(sum(eroded_mask)))
	
	# Mask unwanted values to determine mean uptake value
	pet_activity = eroded_mask * pet_reg_data
	mean_pet_activity = sum(sum(sum(pet_activity))) / mask_size
	
	# Then normalize PET data by this mean activity
	suvr_pet_data = pet_reg_data / mean_pet_activity
	suvr = nb.Nifti1Image(suvr_pet_data, eroded_mask_nifti.affine, header=pet.header)
	
	voxsize = (pet.header["pixdim"])[1:4]
	resampled_image=resample_img(suvr, voxsize, target_shape=pet.shape[:3])
	nb.save(resampled_image, suvr_filename)


surf_over = surface.vol_to_surf(suvr_filename, surf_mesh)
