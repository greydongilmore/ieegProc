import os
import nibabel as nb
import numpy as np
import pandas as pd
import regex as re
import matplotlib.pyplot as plt
from nilearn.plotting.displays import PlotlySurfaceFigure
import plotly.graph_objs as go
from mne.transforms import apply_trans


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

def determine_groups(iterable):
	values = []
	for item in iterable:
		if re.findall(r"([a-zA-Z]+)([0-9]+)([a-zA-Z]+)", item):
			temp = "".join(list(re.findall(r"([a-zA-Z]+)([0-9]+)([a-zA-Z]+)", item)[0]))
		elif '-' in item:
			temp=item.split('-')[0]
		else:
			temp="".join(x for x in item if not x.isdigit())
		
		values.append(temp)
	
	vals,indexes,count = np.unique(values, return_index=True, return_counts=True)
	values = [values[index] for index in sorted(indexes)]
	
	return values,count

hemi = ["lh", "rh"]
surf_suffix = ["pial", "white", "inflated"]

def readRegMatrix(trsfPath):
	with open(trsfPath) as (f):
		return np.loadtxt(f.readlines())

#%%


debug = False
if debug:
	class dotdict(dict):
		"""dot.notation access to dictionary attributes"""
		__getattr__ = dict.get
		__setattr__ = dict.__setitem__
		__delattr__ = dict.__delitem__
	
	class Namespace:
		def __init__(self, **kwargs):
			self.__dict__.update(kwargs)
	
	isub="006"
	datap=r'/home/greydon/Documents/data/SEEG_peds/derivatives'
	
	input=dotdict({
		't1_fname':datap+f'/fastsurfer/sub-P{isub}/mri/orig.mgz',
		'fcsv':datap+ f'/seega_coordinates/sub-P{isub}/sub-P{isub}_space-native_SEEGA.tsv',
		'xfm_noncontrast':datap+f'/atlasreg/sub-P{isub}/sub-P{isub}_acq-noncontrast_desc-rigid_from-noncontrast_to-contrast_type-ras_xfm.txt',
	})
	
	output=dotdict({
		'html':datap+f'/atlasreg/sub-P{isub}/sub-P{isub}_space-native_electrodes.html',
	})
	
	params=dotdict({
		'lh_pial':datap+f'/fastsurfer/sub-P{isub}/surf/lh_pial',
		'rh_pial':datap+f'/fastsurfer/sub-P{isub}/surf/rh_pial',
		'lh_sulc':datap+f'/fastsurfer/sub-P{isub}/surf/lh_sulc',
		'rh_sulc':datap+f'/fastsurfer/sub-P{isub}/surf/rh_sulc',
	})

	snakemake = Namespace(output=output, input=input,params=params)

t1_obj = nb.load(snakemake.input.t1_fname)
Torig = t1_obj.header.get_vox2ras_tkr()
fs_transform=(t1_obj.affine-Torig)+np.eye(4)

verl,facel=nb.freesurfer.read_geometry(snakemake.params.lh_pial)
verr,facer=nb.freesurfer.read_geometry(snakemake.params.rh_pial)

all_ver = np.concatenate([verl, verr], axis=0)
all_face = np.concatenate([facel, facer+verl.shape[0]], axis=0)
surf_mesh = [all_ver, all_face]

all_ver_shift=(apply_trans(fs_transform, all_ver))

if snakemake.input.xfm_noncontrast:
	t1_transform=readRegMatrix(snakemake.input.xfm_noncontrast)
	all_ver_shift=(apply_trans(np.linalg.inv(t1_transform), all_ver_shift))


lh_sulc_data = nb.freesurfer.read_morph_data(snakemake.params.lh_sulc)
rh_sulc_data = nb.freesurfer.read_morph_data(snakemake.params.rh_sulc)
bg_map = np.concatenate((lh_sulc_data, rh_sulc_data))


mesh_3d = go.Mesh3d(x=all_ver_shift[:,0], y=all_ver_shift[:,1], z=all_ver_shift[:,2], i=all_face[:,0], j=all_face[:,1], k=all_face[:,2],opacity=.1,color='grey',alphahull=.1,
					lighting=lighting_effects)

value=np.arange(.1,.6,.05)

df = pd.read_table(os.path.splitext(snakemake.input.fcsv)[0]+".tsv",sep='\t',header=0)
groups,n_members=determine_groups(df['label'].tolist())
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
			size=14,
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
fig.write_html(snakemake.output.html)



