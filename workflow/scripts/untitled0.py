
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 19:22:44 2022

@author: greydon
"""
import ants
import numpy as np
from collections import Counter
from nilearn import plotting
from svgutils.transform import SVGFigure, GroupElement,fromstring
from svgutils.compose import Unit
from tempfile import TemporaryDirectory
from pathlib import Path
from uuid import uuid4
import re
import matplotlib.pyplot as plt
import matplotlib as mpl
from io import BytesIO
import base64
import os
import nibabel as nb
from skfuzzy import cmeans
from scipy.ndimage import gaussian_gradient_magnitude as gradient
import scipy
import SimpleITK as sitk
import glob
import h5py
import vtk

from file_loader.xltek_loader import XltekLoader

loader = XltekLoader(r'/home/greydon/Documents/data/ieeg/Pellinen~ Tyso_0a8e74db-39e7-42d9-95bf-211c618510f4')
ret = loader.load()


def getVTKMatrixFromNumpyMatrix(numpyMatrix):
	dimensions = len(numpyMatrix) - 1
	if dimensions == 2:
		vtkMatrix = vtk.vtkMatrix3x3()
	elif dimensions == 3:
		vtkMatrix = vtk.vtkMatrix4x4()
	else:
		raise ValueError('Unknown matrix dimensions.')
	for row in range(dimensions + 1):
		for col in range(dimensions + 1):
			vtkMatrix.SetElement(row, col, numpyMatrix[row, col])
	return vtkMatrix

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
	
	isub='sub-P232'
	data_dir=r'/home/greydon/Documents/data/DBS/derivatives/trajGuide'
	
	input=dotdict({
				't1': f'{data_dir}/derivatives/{isub}/{isub}*_T1w.nii.gz',
				't2':f'{data_dir}/derivatives/{isub}/{isub}_space-T1w*_T2w.nii.gz',
				'mni_warp':f'{data_dir}/derivatives/{isub}/space/{isub}_from-subject_to-*_type-composite_xfm.h5',
				
				})
	
	output=dotdict({
		'png_out':f'{data_dir}/derivatives/{isub}/summaries/{isub}_desc-contacts_qc.png',
	})
	
	snakemake = Namespace(output=output, input=input)

endPoint=np.array([-34.63,66.05,58.3])
startPoint=np.array([-14.13,44.46,-18.7])

orientation = nb.aff2axcodes(t1_img.affine)

if glob.glob(snakemake.input.t1):
	#t1_img = nb.load(glob.glob(snakemake.input.t1)[0])
	t1_img = sitk.ReadImage(glob.glob(snakemake.input.t1)[0])
	t1_img_data = t1_img.get_fdata()

if glob.glob(snakemake.input.t2):
	t2_img = nb.load(glob.glob(snakemake.input.t2)[0])
	t2_img_data = t2_img.get_fdata()

if glob.glob(snakemake.input.mni_warp):
	warp_fname=glob.glob(snakemake.input.mni_warp)[0]
	mni_warp = h5py.File(warp_fname)
	
	transParam= mni_warp["TransformGroup"]['0']["TransformParameters"][()]
	transParamFixed= mni_warp["TransformGroup"]['0']["TransformFixedParameters"][()].astype(np.int16)
	datasetr = np.reshape(transParam,(transParamFixed[0],transParamFixed[1],transParamFixed[2],3))
	
	affine = np.eye(4)
	affine[0,3]=-1*transParamFixed[3]
	affine[1,3]=-1*transParamFixed[4]
	affine[2,3]=transParamFixed[5]
	
	lps2ras=np.diag([-1, -1, 1, 1])
	transform_lps=np.dot(lps2ras,np.dot(affine,lps2ras))
	affine_vtk = getVTKMatrixFromNumpyMatrix(transform_lps)
	
	img = nb.Nifti1Image(datasetr,affine)
	img.header.set_xyzt_units("mm","sec")
	NEW_NAME=os.path.join(os.path.dirname(warp_fname), os.path.basename(warp_fname).split('.h5')[0] + '_out.nii.gz')
	nb.save(img,NEW_NAME)


	t2_img_data = t2_img.get_fdata()


