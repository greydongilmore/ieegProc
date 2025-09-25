#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 23 02:40:37 2025

@author: greydon
"""

import ants
import numpy as np

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
	
	isub="D151"
	data_dir=r'/home/greydon/Documents/data/lhsc_seeg/derivatives/atlasreg'
	
	input=dotdict({
				't1_path':f'{data_dir}/sub-{isub}/sub-{isub}_desc-n4_T1w.nii.gz',
				'mask':f'{data_dir}/sub-{isub}/sub-{isub}_label-brain_desc-affine_from-MNI152NLin2009cSym_mask.nii.gz',
				})
	params=dotdict({
				'k':3,
				'm':'[0.2,1x1x1]',
				'c':'[3,0]',
				})
	output=dotdict({
				'tissue_priors':[f'{data_dir}/sub-{isub}/sub-{isub}_label-CSF.nii.gz',
					 f'{data_dir}/sub-{isub}/sub-{isub}_label-GM.nii.gz',
					 f'{data_dir}/sub-{isub}/sub-{isub}_label-WM.nii.gz'],
				'seg':f'{data_dir}/sub-{isub}/sub-{isub}_desc-atroposKseg_dseg.nii.gz',
				})
	snakemake = Namespace(output=output, input=input,params=params)


t1 = ants.image_read(snakemake.input.t1_path)
mask = ants.image_read(snakemake.input.mask)
seg = ants.atropos(
    a=t1,
    x=mask,
    i=f'kmeans[{snakemake.params.k}]',
    m=f'{snakemake.params.m}',
    c=f'{snakemake.params.c}'
)
seg_img = seg['segmentation']
prob_imgs = seg['probabilityimages']
means = []
for k in range(1, 4):
    cls_mask = seg_img == k
    vals = t1[cls_mask > 0]
    means.append(np.mean(vals) if vals.size > 0 else -np.inf)
order = np.argsort(means)  # ascending
class_by_intensity = {int(order[0])+1: 1,  # CSF -> 1
                  int(order[1])+1: 2,  # GM  -> 2
                  int(order[2])+1: 3}  # WM  -> 3
relabeled = seg_img.clone()
relabeled.set_direction(seg_img.direction)
relabeled.set_spacing(seg_img.spacing)
relabeled.set_origin(seg_img.origin)
relabeled_np = relabeled.numpy().copy()
for src, dst in class_by_intensity.items():
    relabeled_np[seg_img.numpy() == src] = dst
relabeled = ants.from_numpy(relabeled_np, origin=seg_img.origin, spacing=seg_img.spacing, direction=seg_img.direction)
prob_by_label = [None, None, None]
for src_label, dst_label in class_by_intensity.items():
    prob_by_label[dst_label - 1] = prob_imgs[src_label - 1]

ants.image_write(relabeled, snakemake.output.seg)
ants.image_write(prob_by_label[0], snakemake.output.tissue_priors[0])
ants.image_write(prob_by_label[1], snakemake.output.tissue_priors[1])
ants.image_write(prob_by_label[2], snakemake.output.tissue_priors[2])