#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 24 16:02:30 2021

@author: greydon
"""

import nrrd
import re
import pandas as pd
import nibabel as nb
import numpy as np
import matplotlib.pyplot as plt


def bbox2(img):
	rows = np.any(img, axis=(1, 2))
	cols = np.any(img, axis=(0, 2))
	z = np.any(img, axis=(0, 1))
	
	ymin, ymax = np.where(rows)[0][[0, -1]]
	xmin, xmax = np.where(cols)[0][[0, -1]]
	zmin, zmax = np.where(z)[0][[0, -1]]
	return img[ymin:ymax+1, xmin:xmax+1, zmin:zmax+1]

def sparsify(a):
	ncols = int(a.max()) + 1
	if ncols < 3:
		ncols = 3
	out = np.zeros( (a.size,ncols), dtype=np.uint8)
	out[np.arange(a.size),a.ravel()] = 1
	out.shape = a.shape + (ncols,)
	out = np.transpose(out, axes=(3,0,1,2))
	return out

def bounding_box(seg):
	x = np.any(np.any(seg, axis=0), axis=1)
	y = np.any(np.any(seg, axis=1), axis=1)
	z = np.any(np.any(seg, axis=1), axis=0)
	ymin, ymax = np.where(y)[0][[0, -1]]
	xmin, xmax = np.where(x)[0][[0, -1]]
	zmin, zmax = np.where(z)[0][[0, -1]]
	bbox = np.array([ymin,ymax,xmin,xmax,zmin,zmax])
	return bbox

def get_shape_origin(img_data):
	bbox = bounding_box(img_data)
	ymin, ymax, xmin, xmax, zmin, zmax = bbox
	shape = list(np.array([ymax-ymin, xmax-xmin, zmax-zmin]) + 1)
	origin = [ymin, xmin, zmin]
	return shape, origin

def rgbToHex(color):
	""" Converts RGB colour to HEX.
	"""
	r = int(color[0]*255)
	g = int(color[1]*255)
	b = int(color[2]*255)
	rgb2hex = "#{:02x}{:02x}{:02x}".format(r,g,b)
	return rgb2hex

def write_nrrd(data_obj, out_file,atlas_labels):
	
	data=data_obj.get_fdata()
	
	keyvaluepairs = {}
	keyvaluepairs['dimension'] = 3
	keyvaluepairs['encoding'] = 'gzip'
	keyvaluepairs['kinds'] = ['domain', 'domain', 'domain']
	keyvaluepairs['space'] = 'right-anterior-superior'
	keyvaluepairs['space directions'] = data_obj.affine[:3,:3].T
	keyvaluepairs['type'] = 'double'
	
	box = bounding_box(data)
	seg_cut = data[box[0]:box[1]+1,box[2]:box[3]+1,box[4]:box[5]+1]
	shape, origin = get_shape_origin(data)
	origin = nb.affines.apply_affine(data_obj.affine, np.array([origin]))

	keyvaluepairs['sizes'] = np.array([*shape])
	keyvaluepairs['space origin'] = origin[0]
	
	
	cmap = plt.get_cmap('hsv')
	color_maps=cmap(np.linspace(0, 1, int(np.max(data)))).tolist()
	
	for i in range(int(np.max(data))):
		col_lut=np.array(atlas_labels[atlas_labels['label']==i+1]['lut'].values[0]+[255])/255
		name = 'Segment{}'.format(i)
		keyvaluepairs[name + '_Color'] = ' '.join([f"{a:10.3f}" for a in col_lut])
		keyvaluepairs[name + '_ColorAutoGenerated'] = '1'
		keyvaluepairs[name + '_Extent'] = f'0 {shape[0]-1} 0 {shape[1]-1} 0 {shape[2]-1}'
		keyvaluepairs[name + '_ID'] = 'Segment_{}'.format(i+1)
		keyvaluepairs[name + '_LabelValue'] = '{}'.format(i+1)
		keyvaluepairs[name + '_Layer'] = '0'
		if 'hemi' in list(atlas_labels):
			keyvaluepairs[name + '_Name'] = '_'.join([atlas_labels[atlas_labels['label']==i+1]['hemi'].values[0], atlas_labels[atlas_labels['label']==i+1]['name'].values[0]])
		else:
			keyvaluepairs[name + '_Name'] = atlas_labels[atlas_labels['label']==i+1]['name'].values[0]
		keyvaluepairs[name + '_NameAutoGenerated'] = 1
		keyvaluepairs[name + '_Tags'] = 'TerminologyEntry:Segmentation category' +\
			' and type - 3D Slicer General Anatomy list~SRT^T-D0050^Tissue~SRT^' +\
			'T-D0050^Tissue~^^~Anatomic codes - DICOM master list~^^~^^|'

	keyvaluepairs['Segmentation_ContainedRepresentationNames'] = 'Binary labelmap|'
	keyvaluepairs['Segmentation_ConversionParameters'] = 'placeholder'
	keyvaluepairs['Segmentation_MasterRepresentation'] = 'Binary labelmap'
	
	nrrd.write(out_file, seg_cut, keyvaluepairs)


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
	
	isub="P031"
	data_dir=r'/home/greydon/Documents/data/fNIRS-fmri/derivatives/'
	repo_path = r'/home/greydon/Downloads/tpl-MNI152NLin2009aAsym'
	
	input=dotdict({
				'segs':'/home/greydon/Downloads/tpl-MNI152NLin2009aAsym/tpl-MNI152NLin2009aAsym_atlas-Glasser_dseg.nii.gz',
				})
	
	params=dotdict({
				'atlas_labels':repo_path + r'/tpl-MNI152NLin2009aAsym_atlas-Glasser_dseg.tsv',
				'atlas_colors': '/home/greydon/Documents/GitHub/seeg2bids-pipeline/resources/generic_colors.txt',
				})
	
	output=dotdict({
				'seg_nrrd':"/home/greydon/Downloads/tpl-MNI152NLin2009aAsym/tpl-MNI152NLin2009aAsym_atlas-Glasser_dseg.seg.nrrd"
				})
	
	snakemake = Namespace(output=output, input=input,params=params)


atlas_labels = pd.read_table(snakemake.params.atlas_labels)

if not all(x in list(atlas_labels) for x in {'r','g','b'}):
	atlas_colors = pd.read_csv(snakemake.params.atlas_colors, sep="\t",header=0)

	col_lut_out=[]
	for ilabel in list(atlas_labels['label']):
		col_lut=sum([re.findall(r'\d+', str(x)) for x in atlas_colors[atlas_colors['index'].values==int(ilabel)][['r','g','b']].to_numpy()[0]],[])
		col_lut_out.append([int(x) for x in col_lut])

	atlas_labels['lut']=col_lut_out
else:
	atlas_labels['lut']=atlas_labels[['r','g','b']].values.tolist()

data_obj=nb.load(snakemake.input.segs)

write_nrrd(data_obj, snakemake.output.seg_nrrd, atlas_labels)


#%%


# atlas_labels = pd.read_table(r'/home/greydon/Documents/GitHub/trajectoryGuideModules/resources/ext_libs/space/tpl-MNI152NLin2009cAsym/atlases/Tian_Subcortex_S4_3T_label.txt',header=None)
# atlas_labels['label']=atlas_labels.index+1
# atlas_colors = pd.read_csv(snakemake.params.atlas_colors, sep="\t",header=0)
# atlas_labels.rename(columns={0:'name'}, inplace=True)

# col_lut_out=[]
# for ilabel in list(atlas_labels['label']):
# 	col_lut=sum([re.findall(r'\d+', str(x)) for x in atlas_colors[atlas_colors['index'].values==int(ilabel)][['r','g','b']].to_numpy()[0]],[])
# 	col_lut_out.append([int(x) for x in col_lut])

# atlas_labels['lut']=col_lut_out
# atlas_labels['hemi']=[x.split('-')[-1] for x in atlas_labels['name']]
# atlas_labels['name']=['-'.join(x.split('-')[:-1]).replace('-','') for x in atlas_labels['name']]

# data_obj=nb.load('/home/greydon/Downloads/Tian_Subcortex_S4_3T_2009cAsym.nii.gz')
# out_file='/home/greydon/Documents/GitHub/trajectoryGuideModules/resources/ext_libs/space/tpl-MNI152NLin2009cAsym/atlases/Tian_Subcortex_S4_3T_2009cAsym.seg.nrrd'

# write_nrrd(data_obj, out_file, atlas_labels)

# with open(r'/home/greydon/Documents/GitHub/trajectoryGuideModules/resources/ext_libs/space/template_model_dictionary.json') as (file):
# 	template_model_dictionary = json.load(file)


# melbourne_dir=r'/home/greydon/Documents/GitHub/trajectoryGuideModules/resources/ext_libs/space/tpl-MNI152NLin2009cAsym/atlases/melbourne'

# for imodel in [x for x in os.listdir(melbourne_dir) if x.endswith('.vtk')]:
# 	modelName=imodel.split('desc-')[-1].split('.vtk')[0]
# 	if modelName not in list(template_model_dictionary['melbourne']):
# 		template_model_dictionary['melbourne'][modelName]={}
# 		template_model_dictionary['melbourne'][modelName]['main']=""
# 		template_model_dictionary['melbourne'][modelName]['sub']=""
# 		template_model_dictionary['melbourne'][modelName]['label']=atlas_labels[atlas_labels['name']==modelName]['name'].values[0]
# 		template_model_dictionary['melbourne'][modelName]['color']=rgbToHex(np.array(atlas_labels[atlas_labels['name']==modelName]['lut'].values[0])/255)
# 		template_model_dictionary['melbourne'][modelName]['visible']=True


# json_output = json.dumps(template_model_dictionary, indent=4)
# with open(r'/home/greydon/Documents/GitHub/trajectoryGuideModules/resources/ext_libs/space/template_model_dictionary.json', 'w') as (fid):
# 	fid.write(json_output)
# 	fid.write('\n')
	





