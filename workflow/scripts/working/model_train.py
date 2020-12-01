#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 19:50:47 2020

@author: greydon
"""
import numpy as np
import pandas as pd
import pydicom
import os
import scipy.ndimage
import nibabel as nib
from bids.layout import BIDSLayout
import matplotlib.pyplot as plt
import glob
import shutil
import base64
from io import BytesIO
from mne.transforms import apply_trans

from skimage import measure, morphology #scikit-image
from mpl_toolkits.mplot3d.art3d import Poly3DCollection #3d-plotting

def sorted_nicely(lst):
	convert = lambda text: int(text) if text.isdigit() else text
	alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
	sorted_lst = sorted(lst, key = alphanum_key)
	
	return sorted_lst

def show_slice(img_data, inds):

	unique_z=np.unique([x[2] for x in inds])
	
	cols=7
	rows=len(unique_z)//cols
	if len(unique_z)%cols:
		rows=rows+1
		
	fig,ax = plt.subplots(rows,cols,figsize=[20,20])
	row_cnt=0
	col_cnt=0
	for islice in unique_z:
		ax[row_cnt,col_cnt].set_title('slice %d' % islice)
		ax[row_cnt,col_cnt].imshow(img_data[:,:,islice].T, cmap="gray")
		ax[row_cnt,col_cnt].axis('off')
		plot_annot=[x for x in inds if x[2] == islice]
		for iannot in plot_annot:
			ax[row_cnt,col_cnt].plot(iannot[0],iannot[1], marker='.', markersize=2)
		if col_cnt==(cols-1):
			row_cnt+=1
			col_cnt=0
		else:
			col_cnt+=1
	
	
			
def get_vox2ras_tkr(img):
	'''Get the vox2ras-tkr transform. Inspired
	by get_vox2ras_tkr in
	'''
	ds = img.header.get_zooms()[:3]
	ns = np.array(img.shape[:3]) * ds / 2.0
	v2rtkr = np.array([[-ds[0], 0, 0, ns[0]],
					   [0, 0, ds[2], -ns[2]],
					   [0, -ds[1], 0, ns[1]],
					   [0, 0, 0, 1]], dtype=np.float32)
	return v2rtkr

def resample_img(img):
	
	voxel_dims=[1,1,1]
	target_affine = img.affine.copy()
	
	# Decompose the image affine to allow scaling
	u,s,v = np.linalg.svd(target_affine[:3,:3],full_matrices=False)
	
	# Calculate the translation part of the affine
	spatial_dimensions = (img.header['dim'] * img.header['pixdim'])[1:4]
	resize_factor = spatial_dimensions / voxel_dims
	real_resize_factor = np.round(img.shape * resize_factor) / img.shape
	new_spacing = spatial_dimensions / real_resize_factor
	
	# Reconstruct the affine
	target_affine[:3,:3] = u @ np.diag(new_spacing) @ v
	
	

# 	# Calculate the translation affine as a proportion of the real world
# 	# spatial dimensions
# 	image_center_as_prop = img.affine[0:3,3] / spatial_dimensions
# 	
# 	# Calculate the equivalent center coordinates in the target image
# 	dimensions_of_target_image = (np.array(voxel_dims) * np.array(target_shape))
# 	target_center_coords =  dimensions_of_target_image * image_center_as_prop
# 	
# 	target_affine[:3,3] = target_center_coords
	
	resampled_img = image.resample_img(img, target_affine=target_affine, interpolation='nearest')
	resampled_img.header.set_zooms((np.absolute(new_spacing)))


bids_dir=r'/media/veracrypt6/projects/iEEG/imaging/clinical/bids'
seega_dir=r'/media/veracrypt6/projects/iEEG/imaging/clinical/deriv/seega_scenes'
seega_coords_dir=r'/media/veracrypt6/projects/iEEG/imaging/clinical/deriv/seega_coordinates_old'
out_dir=r'/media/veracrypt6/projects/iEEG/imaging/clinical/deriv/model_pred'

layout = BIDSLayout(bids_dir)

if not os.path.exists(os.path.join(out_dir, 'ground_truth')):
	os.makedirs(os.path.join(out_dir, 'ground_truth'))
	
for isub in layout.get_subjects():
	
	infile=glob.glob(seega_dir+f'/sub-{isub}/ct_to_mri.nii.gz')
	in_tsv=glob.glob(seega_coords_dir+f'/sub-{isub}/sub-{isub}_space-native_SEEGA.tsv')
	if infile and in_tsv:
		outfile = f"sub-{isub}_"+ os.path.splitext(os.path.splitext(os.path.basename(infile[0]))[0])[0]+'_resample.nii'
		if not os.path.exists(os.path.join(out_dir, 'preprocessed', outfile)):
		
			# Console output
			print("\nPreprocessing data from patient " + str(i) + "\n");
		
			img = nib.load(os.path.join())
		
			# Resample data to be consistent / fixed at 1 mm x 1 mm x 1 mm
			img_resampled = resample(img)
			
			print("Shape before resampling\t", img.shape)
			print("Shape after resampling\t", img_resampled.shape)
		
			nib.save(os.path.join(out_dir, 'preprocessed', outfile), patient_pixels)
		
		if not os.path.exists(os.path.join(out_dir, 'ground_truth', os.path.basename(in_tsv[0]))):
			shutil.copyfile(in_tsv[0], os.path.join(out_dir, 'ground_truth', os.path.basename(in_tsv[0])))

subjects = sorted_nicely([x.split('sub-')[-1].split('_')[0] for x in os.listdir(os.path.join(out_dir,'ground_truth'))])
for isub in subjects:
	coord_tsv = pd.read_table(glob.glob(out_dir+f'/ground_truth/sub-{isub}_space-native_SEEGA.tsv')[0])
	img = nib.load(glob.glob(out_dir+f'/preprocessed/sub-{isub}_ct_to_mri_resample.nii')[0])
	img_data=img.get_fdata()
	
	vox2ras_tkr = get_vox2ras_tkr(img)
	ras2vox_tkr = np.linalg.inv(vox2ras_tkr)
	
	coords = coord_tsv[['x','y','z']].to_numpy()

	inds = []
	for i in range(len(coords)):
		vec = np.hstack([coords[i,:],1])
	
		#dseg_affine is used to xfm indices to RAS coords, 
		# so we use the inverse to go the other way
		tvec = np.linalg.inv(img.affine) @ vec.T
		inds.append(np.round(tvec[:3]).astype('int'))
	
	show_slice(img_data, inds)
	