#!/usr/bin/env python3

import numpy as np
import nibabel as nib
import SimpleITK as sitk


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
	
	isub='sub-P109'
	data_dir=r'/home/greydon/Documents/data/SEEG/derivatives/atlasreg'
	
	input=dotdict({
				'xfm': f'{data_dir}/{isub}/{isub}_desc-rigid_from-ct_to-T1w_type-ras_xfm.txt',
				'flo': f'{data_dir}/{isub}/{isub}_ct.nii.gz',
				})
	
	output=dotdict({
		'warped_subj':f'{data_dir}/{isub}/{isub}_space-T1w_desc-rigid_ct.nii.gz',
	})
	
	snakemake = Namespace(output=output, input=input)

#load transform matrix
transformMatrix = np.loadtxt(snakemake.input.xfm)

#load source nifti and apply transform matrix
flo_img=nib.load(snakemake.input.flo)
transform = np.dot(transformMatrix, flo_img.affine)

flo_img_trans = nib.Nifti1Image(flo_img.get_fdata(), header=flo_img.header, affine=transform)
flo_img_trans.set_qform(flo_img_trans.affine, 1)
nib.save(flo_img_trans,snakemake.output.warped_subj)


