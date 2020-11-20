#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 24 01:19:14 2020

@author: greydon
"""
import pandas as pd
import numpy as np
import nibabel as nib
import os

#read fcsv electrodes file
df_elec = pd.read_table(snakemake.input.fcsv,sep=',',header=2)
df_elec
df_atlas = pd.read_table(snakemake.input.dseg_tsv)
df_atlas

#load up tissue probability, warped from template
tissue_prob_vol = dict()
tissue_prob_elec = dict()

for label,nii in zip(snakemake.config['tissue_labels'], snakemake.input.tissue_seg):
	print(label)
	print(nii)
	tissue_prob_vol[label] = nib.load(nii).get_fdata()
	tissue_prob_elec[label] = list()

#load dseg nii (as integer)
dseg_nii = nib.load(snakemake.input.dseg_nii)
dseg_vol = dseg_nii.get_fdata().astype('int')

#get affine from image, so we can go from RAS coords to array indices
dseg_affine = dseg_nii.affine
dseg_affine

#get coords from fcsv
coords = df_elec[['x','y','z']].to_numpy()


labelnames = []

for i in range(len(coords)):

	vec = np.hstack([coords[i,:],1])

	#dseg_affine is used to xfm indices to RAS coords, 
	# so we use the inverse to go the other way
	tvec = np.linalg.inv(dseg_affine) @ vec.T   
	inds = np.round(tvec[:3]).astype('int')
	labelnum = dseg_vol[inds[0],inds[1],inds[2]]
	
	
	if labelnum >0:
		labelnames.append(df_atlas.loc[df_atlas['label']==labelnum,'name'].to_list()[0])
	else:
		labelnames.append('None')
		
	for label in snakemake.config['tissue_labels']:
		tissue_prob_elec[label].append(tissue_prob_vol[label][inds[0],inds[1],inds[2]])
	
#add new columns to existing dataframe
df_elec['atlas_label'] = labelnames
for label in snakemake.config['tissue_labels']:
	df_elec[label] = tissue_prob_elec[label]

#create new dataframe with selected variables and save it
out_df = df_elec[['label','atlas_label'] + snakemake.config['tissue_labels'] + ['x','y','z']]
out_df.to_csv(snakemake.output.tsv,sep='\t',float_format='%.3f',index=False)
out_df.to_excel(os.path.splitext(snakemake.output.tsv)[0]+'.xlsx',float_format='%.3f',index=False)

out_df