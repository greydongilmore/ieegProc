#!/usr/bin/env python3

import numpy as np
import nibabel as nib


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
	
	isub='sub-P108'
	data_dir=r'/home/greydon/Documents/data/SEEG/derivatives/atlasreg'
	
	input=dotdict({
				'xfm': f'{data_dir}/{isub}/{isub}_desc-rigid_from-noncontrast_to-contrast_type-ras_xfm.txt',
				'flo': f'{data_dir}/{isub}/{isub}_acq-noncontrast_T1w.nii.gz',
				})
	
	output=dotdict({
		'warped_subj':f'{data_dir}/{isub}/{isub}_acq-noncontrast_space-T1w_desc-rigid_T1w.nii.gz',
	})
	
	snakemake = Namespace(output=output, input=input)


transformMatrix = np.loadtxt(snakemake.input.xfm)

flo_img=nib.load(snakemake.input.flo)
transform = np.dot( np.linalg.inv(transformMatrix), flo_img.affine)
flo_img_trans = nib.Nifti1Image(flo_img.get_fdata(), header=flo_img.header, affine=transform)
flo_img_trans.header.set_data_dtype('float32')
nib.save(flo_img_trans,snakemake.output.warped_subj)

