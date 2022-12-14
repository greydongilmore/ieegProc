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
import itertools
import matplotlib.pyplot as plt


def bbox2(img):
	rows = np.any(img, axis=(1, 2))
	cols = np.any(img, axis=(0, 2))
	z = np.any(img, axis=(0, 1))
	
	ymin, ymax = np.where(rows)[0][[0, -1]]
	xmin, xmax = np.where(cols)[0][[0, -1]]
	zmin, zmax = np.where(z)[0][[0, -1]]
	return img[ymin:ymax+1, xmin:xmax+1, zmin:zmax+1]

def sorted_nicely(data, reverse = False):
	convert = lambda text: int(text) if text.isdigit() else text
	alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
	
	return sorted(data, key = alphanum_key, reverse=reverse)

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
		name = 'Segment{}'.format(i)
		keyvaluepairs[name + '_Color'] = ' '.join([f"{a:10.3f}" for a in color_maps[i]])
		keyvaluepairs[name + '_ColorAutoGenerated'] = '1'
		keyvaluepairs[name + '_Extent'] = f'0 {shape[0]-1} 0 {shape[1]-1} 0 {shape[2]-1}'
		keyvaluepairs[name + '_ID'] = 'Segment_{}'.format(i+1)
		keyvaluepairs[name + '_LabelValue'] = '{}'.format(i+1)
		keyvaluepairs[name + '_Layer'] = '0'
		keyvaluepairs[name + '_Name'] = '_'.join([atlas_labels[atlas_labels['label']==i+1]['hemi'].values[0], atlas_labels[atlas_labels['label']==i+1]['name'].values[0]])
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
	
	isub="P101"
	data_dir=r'/home/greydon/Documents/data/SEEG/derivatives/'
	repo_path = r'/home/greydon/Documents/GitHub'
	
	input=dotdict({
				'segs':data_dir +  f'atlasreg/sub-{isub}/sub-{isub}_desc-dilated_atlas-CerebrA_from-MNI152NLin2009cSym_reg-SyN_dseg.nii.gz',
				})
	
	params=dotdict({
				'atlas_labels':repo_path + r'/seeg2bids-pipeline/resources/tpl-MNI152NLin2009cSym/tpl-MNI152NLin2009cSym_atlas-CerebrA_dseg.tsv',
				})
	
	output=dotdict({
				'seg_nrrd':data_dir + f"seega_scenes/sub-{isub}/sub-{isub}_desc-segmentations.seg.nrrd"
				})
	
	snakemake = Namespace(output=output, input=input,params=params)


atlas_labels = pd.read_table(snakemake.params.atlas_labels)

data_obj=nb.load(snakemake.input.segs)

write_nrrd(data_obj, snakemake.output.seg_nrrd, atlas_labels)


#%%


