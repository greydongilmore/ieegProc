import os
import shutil
import pandas as pd
import numpy as np
import re
import csv
import glob
import math


#Euclidean distance calculation
def euclidianDistanceCalc(xyz_planned, xyz_actual):
	if xyz_planned.shape[0]>1:
		euc_dist=[]
		for ipoint in range(xyz_planned.shape[0]):
			plan_act_diff = xyz_planned[ipoint] - xyz_actual[ipoint]
			euc_dist.append(math.sqrt(sum(plan_act_diff**2)))
	else:
		plan_act_diff = xyz_planned - xyz_actual
		euc_dist = math.sqrt(sum(plan_act_diff**2))
	return euc_dist

#Radial distance calculation
def radialDistanceCalc(pt, xyz_entry, xyz_target):
	if xyz_entry.shape[0]>1:
		dist3d=[]
		for ipoint in range(xyz_entry.shape[0]):
			x1_minus_pt = pt[ipoint] - xyz_entry[ipoint]
			x2_minus_x1 = xyz_target[ipoint] - xyz_entry[ipoint]
		
			sumsq_x1_minus_pt = sum(x1_minus_pt * x1_minus_pt)
			sumsq_x2_minus_x1 = sum(x2_minus_x1 * x2_minus_x1)
		
			mydotprod = np.dot(x1_minus_pt, x2_minus_x1)
		
			dist3d.append(np.sqrt((sumsq_x1_minus_pt * sumsq_x2_minus_x1 - (mydotprod * mydotprod))/sumsq_x2_minus_x1))
	else:
		x1_minus_pt = pt - xyz_entry
		x2_minus_x1 = xyz_target - xyz_entry
	
		sumsq_x1_minus_pt = sum(x1_minus_pt * x1_minus_pt)
		sumsq_x2_minus_x1 = sum(x2_minus_x1 * x2_minus_x1)
	
		mydotprod = np.dot(x1_minus_pt, x2_minus_x1)
	
		dist3d = np.sqrt((sumsq_x1_minus_pt * sumsq_x2_minus_x1 - (mydotprod * mydotprod))/sumsq_x2_minus_x1)
	return dist3d

#Radial angle calculation
def ptLineAngleCalc(pt, x_entry, x_target):
	if x_entry.shape[0]>1:
		deg_angle=[]
		for ipoint in range(x_entry.shape[0]):
			try:
				x1_minus_pt = pt[ipoint] - x_entry[ipoint]
				x2_minus_x1 = x_target[ipoint] - x_entry[ipoint]
			
				sumsq_x1_minus_pt = sum(x1_minus_pt**2)
				sumsq_x2_minus_x1 = sum(x2_minus_x1**2)
			
				mydotprod = np.dot(x1_minus_pt, x2_minus_x1) # sum of products of elements
			
				rad_angle = math.acos(mydotprod/(np.sqrt(sumsq_x1_minus_pt)*np.sqrt(sumsq_x2_minus_x1)))
				deg_angle.append(math.degrees(rad_angle))
			except:
				deg_angle.append(np.nan)
				print(f"Check point {ipoint}")
	else:
		x1_minus_pt = pt - x_entry
		x2_minus_x1 = x_target - x_entry
	
		sumsq_x1_minus_pt = sum(x1_minus_pt**2)
		sumsq_x2_minus_x1 = sum(x2_minus_x1**2)
	
		mydotprod = np.dot(x1_minus_pt, x2_minus_x1) # sum of products of elements
	
		rad_angle = math.acos(mydotprod/(np.sqrt(sumsq_x1_minus_pt)*np.sqrt(sumsq_x2_minus_x1)))
		deg_angle = math.degrees(rad_angle)
	return deg_angle

#Line angle calculation
def lineLineAngleCalc(a_entry, a_target, b_entry, b_target):
	if a_entry.shape[0]>1:
		deg_angle=[]
		for ipoint in range(a_entry.shape[0]):
			try:
				vectorA = a_target[ipoint] - a_entry[ipoint]
				vectorB = b_target[ipoint] - b_entry[ipoint]
			
				sumsq_vectorA = sum(vectorA**2)
				sumsq_vectorB = sum(vectorB**2)
			
				mydotprod = sum(vectorA*vectorB)
			
				rad_angle = math.acos(mydotprod/(np.sqrt(sumsq_vectorA)*np.sqrt(sumsq_vectorB)))
				deg_angle.append(math.degrees(rad_angle))
			except:
				deg_angle.append(np.nan)
				print(f"Check point {ipoint}")
	else:
		vectorA = a_target - a_entry
		vectorB = b_target - b_entry
	
		sumsq_vectorA = sum(vectorA**2)
		sumsq_vectorB = sum(vectorB**2)
	
		mydotprod = sum(vectorA*vectorB)
	
		rad_angle = math.acos(mydotprod/(np.sqrt(sumsq_vectorA)*np.sqrt(sumsq_vectorB)))
		deg_angle = math.degrees(rad_angle)
	return deg_angle

def sorted_nicely(lst):
	convert = lambda text: int(text) if text.isdigit() else text
	alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
	sorted_lst = sorted(lst, key = alphanum_key)
	
	return sorted_lst

def determine_groups(iterable):
	values = []
	for item in iterable:
		if re.findall(r"([a-zA-Z]+)([0-9]+)([a-zA-Z]+)", item):
			temp = "".join(list(re.findall(r"([a-zA-Z]+)([0-9]+)([a-zA-Z]+)", item)[0]))
		else:
			temp="".join(x for x in item if not x.isdigit())
		
		values.append(temp)
	
	indexes = np.unique(values, return_index=True)[1]
	values_unique = [values[index] for index in sorted(indexes)]
	
	return values_unique,values

chan_label_dic ={
	3:"RD10R-SP03X",
	4:"RD10R-SP04X",
	5:"RD10R-SP05X",
	6:"RD10R-SP06X",
	7:"RD10R-SP07X"
}

subject_id=f"sub-{snakemake.config['subject_prefix']}" + snakemake.wildcards.subject
output_dir=os.path.join(snakemake.config['out_dir'], 'derivatives', 'slicer_scene', subject_id)

planned=False
actual=False

if len(snakemake.input.fcsv_files)>0:
	for ifcsv in snakemake.input.fcsv_files:
		output_filen=os.path.join(output_dir, '_'.join([subject_id, os.path.basename(str(ifcsv)).split('_')[-1]]))
		shutil.copy(str(ifcsv), output_filen)
		if 'planned' in os.path.basename(output_filen):
			planned=True
		elif 'actual' in os.path.basename(output_filen):
			actual=True

if len(snakemake.input.fcsv_acpc)>0:
	if os.path.exists(str(snakemake.input.fcsv_acpc)):
		output_filen=os.path.join(output_dir, '_'.join([subject_id, 'acpc.fcsv']))
		shutil.copy(str(snakemake.input.fcsv_acpc), output_filen)

if len(snakemake.input.noncontrast_t1w)>0:
	if os.path.exists(str(snakemake.input.noncontrast_t1w)):
		output_filen=os.path.join(output_dir, '_'.join([subject_id, 'acq-noncontrast_T1w.nii.gz']))
		shutil.copy(str(snakemake.input.noncontrast_t1w), output_filen)

if len(snakemake.input.contrast_t1w)>0:
	if os.path.exists(str(snakemake.input.contrast_t1w)):
		output_filen=os.path.join(output_dir, '_'.join([subject_id, 'acq-contrast_T1w.nii.gz']))
		shutil.copy(str(snakemake.input.contrast_t1w), output_filen)

if len(snakemake.input.segs)>0:
	for iseg in snakemake.input.segs:
		output_filen=os.path.join(output_dir, '_'.join([subject_id, iseg.split('_')[1] +'.nii.gz']))
		shutil.copy(str(iseg), output_filen)

if len(snakemake.input.atlas_segs)>0:
	output_filen=os.path.join(output_dir, '_'.join([subject_id, 'desc-segmentations.nii.gz']))
	shutil.copy(str(snakemake.input.atlas_segs), output_filen)

if len(snakemake.input.ct_file)>0:
	output_filen=os.path.join(output_dir, '_'.join([subject_id, '_ct.nii.gz']))
	shutil.copy(str(snakemake.input.ct_file), output_filen)

if len(snakemake.input.pet_file)>0:
	output_filen=os.path.join(output_dir, '_'.join([subject_id, '_pet.nii.gz']))
	shutil.copy(str(snakemake.input.pet_file), output_filen)