#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 24 01:19:14 2020

@author: greydon
"""
import pandas as pd
import numpy as np
import nibabel as nib
import os,glob
import regex as re

def lookup_atlas_label(df_template, coords_columns, dseg_nii, df_atlas, fuzzy_dist=2):
	dseg_vol = dseg_nii.get_fdata().astype('int')
	dseg_affine = dseg_nii.affine
	coords = df_template[coords_columns].to_numpy()
	labelnames = []
	fuzzy_list=[]
	voxel_list=[]
	for i in range(len(coords)):
		coords_vx = np.round(nib.affines.apply_affine(np.linalg.inv(dseg_affine),coords[i,:])).astype(int)
		used_fuzzy = 0
		try:
			voxel_value = dseg_vol[coords_vx[0], coords_vx[1],coords_vx[2]]
		except:
			voxel_value = 0
		
		if (fuzzy_dist is not None) and (voxel_value == 0):
			fuzzy_diameter = fuzzy_dist * 2 + 1
			distances_mat = np.zeros((fuzzy_diameter, fuzzy_diameter, fuzzy_diameter))
			for x in range(fuzzy_diameter):
				for y in range(fuzzy_diameter):
					for z in range(fuzzy_diameter):
						distances_mat[x, y, z] = np.linalg.norm(np.array([fuzzy_dist, fuzzy_dist, fuzzy_dist])-np.array([x, y, z]))
			
			# check if the distances box will exceed image boundaries (super edgy case)
			trim = np.zeros((3, 2), dtype=int)
			for i in range(3):
				if coords_vx[i] - fuzzy_dist < 0:
					trim[i, 0] = fuzzy_dist - coords_vx[i]
				if coords_vx[i] + fuzzy_dist + 1 > dseg_vol.shape[i]:
					trim[i, 1] = coords_vx[i] + fuzzy_dist + 1 - dseg_vol.shape[i]
			
			assert np.all(trim >= 0), 'Trim (' + str(trim) + ') should be non-negative'
			
			# get nearest voxel that is not zero, but less than specified voxels away
			selected_atlasdata = dseg_vol[
									(coords_vx[0] - fuzzy_dist + trim[0, 0]):(coords_vx[0] + fuzzy_dist + 1 - trim[0, 1]),
									(coords_vx[1] - fuzzy_dist + trim[1, 0]):(coords_vx[1] + fuzzy_dist + 1 - trim[1, 1]),
									(coords_vx[2] - fuzzy_dist + trim[2, 0]):(coords_vx[2] + fuzzy_dist + 1 - trim[2, 1])
								]
			trimmed_distances = distances_mat[
				trim[0, 0]:(-1 * trim[0, 1]) if trim[0, 1] != 0 else None,
				trim[1, 0]:(-1 * trim[1, 1]) if trim[1, 1] != 0 else None,
				trim[2, 0]:(-1 * trim[2, 1]) if trim[2, 1] != 0 else None]
			
			distances = np.ma.masked_where((selected_atlasdata == 0) | (trimmed_distances > fuzzy_dist),trimmed_distances)
			nearest_voxel = np.unravel_index(np.argmin(distances),distances.shape)
			voxel_value = selected_atlasdata[nearest_voxel]
			used_fuzzy = 1
		
		if voxel_value > 0:
			region_name = df_atlas.loc[df_atlas['index']==voxel_value,'name'].to_list()[0]
		else:
			region_name = np.nan
			voxel_value = np.nan
		
		labelnames.append(region_name)
		fuzzy_list.append(used_fuzzy)
		voxel_list.append(voxel_value)
	
	out_df=df_template.copy()
	out_df.insert(1,'atlas_label',labelnames)
	out_df.insert(2,'fuzzy',fuzzy_list)
	out_df.insert(3,'vox_val',voxel_list)
	return out_df

def sorted_nicely(lst):
	convert = lambda text: int(text) if text.isdigit() else text
	alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
	sorted_lst = sorted(lst, key = alphanum_key)
	
	return sorted_lst


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
	
	isub="P177"
	data_dir=r'/home/greydon/Documents/GitHub/seeg_accuracy_2025/data/fcsv_new'
	repo_path = r'/home/greydon/Documents/GitHub/seeg_accuracy_2025/resources'

	input=dotdict({
				'fcsv_template':f'{data_dir}/sub-{isub}/sub-{isub}_space-MNI152NLin2009cSym_SEEGA.fcsv',
				'dseg_tsv':f'{repo_path}/tpl-MNI152NLin2009cSym_res-1_atlas-USCL_dseg.tsv',
				'dseg_nii':f'{repo_path}/tpl-MNI152NLin2009cSym_res-1_atlas-USCL_dseg_dilated3.nii.gz',
				})
	
	output=dotdict({
				'tsv':f'{data_dir}/sub-{isub}/sub-{isub}_space-MNI152NLin2009cSym_electrodes.tsv',
				'xlsx':f'{data_dir}/sub-{isub}/sub-{isub}_space-MNI152NLin2009cSym_electrodes.xlsx',
				})
	config=dotdict({'fuzzy_dist':2,
				})
	snakemake = Namespace(output=output, input=input,config=config)

fuzzy_dist=snakemake.config.fuzzy_dist

#read fcsv electrodes file
df_template = pd.read_table(snakemake.input.fcsv_template,sep=',',header=2)

df_atlas = pd.read_table(snakemake.input.dseg_tsv)

#load dseg nii (as integer)
dseg_nii = nib.load(snakemake.input.dseg_nii)

out_df = lookup_atlas_label(df_template, dseg_nii, df_atlas)

out_df.to_csv(snakemake.output.tsv,sep='\t',float_format='%.3f',index=False)
out_df.to_excel(snakemake.output.xlsx,index=False)

#%%

isub='sub-F011'
fuzzy_dist=3


data_dir=r'/home/greydon/Documents/GitHub/seeg_accuracy_2025/data/fcsv_new'
repo_path = r'/home/greydon/Documents/GitHub/seeg_accuracy_2025/resources'

df_atlas = pd.read_table(f'{repo_path}/tpl-MNI152NLin2009cSym_res-1_atlas-USCL_dseg.tsv')
dseg_nii = nib.load(f'{repo_path}/tpl-MNI152NLin2009cSym_res-1_atlas-USCL_dseg_dilated3.nii.gz')

atlas_df=[]
for isub in sorted_nicely([x for x in os.listdir( data_dir)]):
	df_template = pd.read_table(f'{data_dir}/{isub}/{isub}_space-MNI152NLin2009cSym_actual.fcsv',sep=',',header=2)
	out_df = lookup_atlas_label(df_template, dseg_nii, df_atlas)
	atlas_df.append(out_df)


atlas_df_out=pd.concat(atlas_df)

atlas_df_out.to_csv(f'{os.path.dirname(os.path.dirname(data_dir))}/output/space-MNI152NLin2009cSym_atlas-USC_electrodes.tsv',sep='\t',float_format='%.3f',index=False)
atlas_df_out.to_excel(f'{os.path.dirname(os.path.dirname(data_dir))}/output/space-MNI152NLin2009cSym_atlas-USC_electrodes.xlsx', na_rep='NA',index=False)

