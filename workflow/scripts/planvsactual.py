#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 19 20:15:44 2022

@author: greydon
"""

import os
import pandas as pd
import numpy as np
import re
import csv
import shutil
import glob
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D
import math

##################################
## Metrics Function Definitions ##
##################################

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

#%%

data_dir=[r'/media/veracrypt6/projects/SEEG/old_database/derivatives/seega_coordinates',r'/media/veracrypt6/projects/SEEG/derivatives/seega_coordinates']

output_filename='/media/veracrypt6/projects/SEEG/derivatives/seega_coordinates'

elec_data=[]
error_files=[]
for idir in data_dir:
	patients = sorted_nicely([x for x in os.listdir(idir) if os.path.isdir(os.path.join(idir, x))])
	for isub in patients:
		patient_files=glob.glob(f"{os.path.join(idir,isub)}/*.fcsv")
		spaces = np.unique([os.path.basename(x).split('_')[1] for x in patient_files if 'space' in x])
		
		for ispace in spaces:
			file_data={}
			for ifile in [x for x in patient_files if ispace in os.path.basename(x)]:
				fcsv_data = pd.read_csv(ifile, skiprows=3, header=None)
				fcsv_data.rename(columns={0:'node_id', 1:'x', 2:'y', 3:'z', 4:'ow', 5:'ox',
						6:'oy', 7:'oz', 8:'vis', 9:'sel', 10:'lock',
						11:'label', 12:'description', 13:'associatedNodeID'}, inplace=True)
				if ifile.lower().endswith('actual.fcsv'):
					file_data['actual']=fcsv_data
				elif ifile.lower().endswith('planned.fcsv'):
					file_data['planned']=fcsv_data
				elif ifile.lower().endswith('seega.fcsv'):
					file_data['seega']=fcsv_data
			
			groupsActual,actual_all = determine_groups(np.array(file_data['actual']['label'].values))
			groupsPlanned, planned_all = determine_groups(np.array(file_data['planned']['label'].values))
			if 'seega' in list(file_data):
				groupsSeega, seega_all = determine_groups(np.array(file_data['seega']['label'].values))
			
			for igroup in set(groupsActual).intersection(groupsActual):
				try:
					elec_temp={}
					elec_temp['subject']=isub
					elec_temp['space']=ispace.split('-')[-1]
					elec_temp['electrode']=igroup
					
					if 'seega' in list(file_data):
						elec_data_temp = file_data['seega'][pd.DataFrame(seega_all)[0]==igroup][['x','y','z']].to_numpy()
						dist = np.mean(np.linalg.norm(elec_data_temp[:-1,:] - elec_data_temp[1:,:], axis=1))
						idx_ielec,val_ielec = min(enumerate(list(chan_label_dic)), key=lambda x: abs(x[1]-dist))
						elec_temp['electrodeType']=chan_label_dic[val_ielec]
						elec_temp['numContacts']=file_data['seega'][pd.DataFrame(seega_all)[0]==igroup].shape[0]
					
					elec_temp['actualTipX']=file_data['actual'][pd.DataFrame(actual_all)[0]==igroup]['x'].values[0]
					elec_temp['actualTipY']=file_data['actual'][pd.DataFrame(actual_all)[0]==igroup]['y'].values[0]
					elec_temp['actualTipZ']=file_data['actual'][pd.DataFrame(actual_all)[0]==igroup]['z'].values[0]
					elec_temp['actualEntryX']=file_data['actual'][pd.DataFrame(actual_all)[0]==igroup]['x'].values[1]
					elec_temp['actualEntryY']=file_data['actual'][pd.DataFrame(actual_all)[0]==igroup]['y'].values[1]
					elec_temp['actualEntryZ']=file_data['actual'][pd.DataFrame(actual_all)[0]==igroup]['z'].values[1]
					elec_temp['plannedTipX']=file_data['planned'][pd.DataFrame(actual_all)[0]==igroup]['x'].values[0]
					elec_temp['plannedTipY']=file_data['planned'][pd.DataFrame(actual_all)[0]==igroup]['y'].values[0]
					elec_temp['plannedTipZ']=file_data['planned'][pd.DataFrame(actual_all)[0]==igroup]['z'].values[0]
					elec_temp['plannedEntryX']=file_data['planned'][pd.DataFrame(actual_all)[0]==igroup]['x'].values[1]
					elec_temp['plannedEntryY']=file_data['planned'][pd.DataFrame(actual_all)[0]==igroup]['y'].values[1]
					elec_temp['plannedEntryZ']=file_data['planned'][pd.DataFrame(actual_all)[0]==igroup]['z'].values[1]
					elec_data.append(elec_temp)
				except:
					error_files.append([isub,igroup])


elec_data_raw=pd.DataFrame(elec_data)

elec_data_final=[]
for ispace in ('native','acpc'):
	
	elec_data_space=elec_data_raw[elec_data_raw['space']==ispace].copy()
	
	xyz_planned_entry = np.array([elec_data_space['plannedEntryX'], elec_data_space['plannedEntryY'], elec_data_space['plannedEntryZ']]).T
	xyz_actual_entry = np.array([elec_data_space['actualEntryX'], elec_data_space['actualEntryY'], elec_data_space['actualEntryZ']]).T
	xyz_planned_target = np.array([elec_data_space['plannedTipX'], elec_data_space['plannedTipY'], elec_data_space['plannedTipZ']]).T
	xyz_actual_target = np.array([elec_data_space['actualTipX'], elec_data_space['actualTipY'], elec_data_space['actualTipZ']]).T
	
	# calculate euclidean distance
	elec_data_space['euclid_dist'] = euclidianDistanceCalc(xyz_planned_target, xyz_actual_target)
	elec_data_space['radial_dist_target'] = radialDistanceCalc(xyz_actual_target, xyz_planned_entry, xyz_planned_target)
	elec_data_space['radial_dist_entry'] = radialDistanceCalc(xyz_actual_entry, xyz_planned_entry, xyz_planned_target)
	elec_data_space['radial_angle'] = ptLineAngleCalc(xyz_actual_target, xyz_planned_entry, xyz_planned_target)
	elec_data_space['line_angle'] = lineLineAngleCalc(xyz_actual_entry, xyz_actual_target, xyz_planned_entry, xyz_planned_target)
	
	elec_data_out=elec_data_space[['subject','plannedTipX','plannedTipY','plannedTipZ',
								 'actualTipX','actualTipY','actualTipZ','electrode','electrodeType',
								 'euclid_dist','radial_dist_target','radial_dist_entry','radial_angle','line_angle']]
	
	# remove outliers
	large_depth_error = elec_data_out[np.abs(elec_data_out['euclid_dist']) > 10]
	elec_data_out=elec_data_out.drop(large_depth_error.index)
	
	elec_data_out = elec_data_out.rename(columns={'subject': 'subjid', 
							'plannedTipX': 'X_planned',
							'plannedTipY': 'Y_planned',
							'plannedTipZ': 'Z_planned',
							'actualTipX': 'X_actual',
							'actualTipY': 'Y_actual',
							'actualTipZ': 'Z_actual',
							'electrode': 'label',
							'electrodeType': 'description'
	})
	
	elec_data_out['type']=np.repeat(ispace,elec_data_out.shape[0])
	elec_data_out.to_csv(os.path.join(output_filename,f'output_space-{ispace}.csv'), index=False , float_format='%.3f')
	
	elec_data_final.append(elec_data_out)

elec_data_final=pd.concat(elec_data_final)

elec_data_native=elec_data_final[elec_data_final['type']=='native']
elec_data_acpc=elec_data_final[elec_data_final['type']=='acpc']

elec_data_native['euclid_dist'].describe().round(2)
elec_data_acpc['euclid_dist'].describe().round(2)


vals,cnts=np.unique(elec_data_native['label'],return_counts=True)
label_counts=pd.DataFrame({'label':vals,'count':cnts})

label_counts_mask = label_counts[np.abs(label_counts['count']) < 10]
label_counts=label_counts.drop(label_counts_mask.index).reset_index(drop=True)

errors=[]
for ilabel in label_counts['label']:
	errors.append(elec_data_native[elec_data_native['label']==ilabel]['euclid_dist'].mean())

label_counts['error']=errors


fig = plt.figure(figsize=(16,14))
ax = fig.add_subplot(111, projection='3d')
scatter = ax.scatter(elec_data_acpc['X_actual'].values, elec_data_acpc['Y_actual'].values, elec_data_acpc['Z_actual'].values,s=2)

fig = plt.figure(figsize=(16,14))
ax = fig.add_subplot(111, projection='3d')
scatter = ax.scatter(elec_data_acpc['X_actual'].values, elec_data_acpc['Y_actual'].values, elec_data_acpc['Z_actual'].values,s=2)

ax.tick_params(axis='x', labelrotation=45)
ax.tick_params(axis='y', labelrotation=45)
ax.set_xlabel('X axis: M-L (mm)', fontsize=14, fontweight='bold', labelpad=14)
ax.set_ylabel('Y axis: A-P (mm)', fontsize=14, fontweight='bold', labelpad=14)
ax.set_zlabel('Z axis: I-S (mm)', fontsize=14, fontweight='bold', labelpad=14)
ax.xaxis._axinfo['label']['space_factor'] = 3.8
ax.view_init(elev=20, azim=50)
ax.figure.canvas.draw()
fig.tight_layout()
plt.legend(handles=scatter.legend_elements()[0], labels=[str(x) for x in np.unique(voxel_info[:,3])])




















