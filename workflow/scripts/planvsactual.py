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


data_dir=[r'/media/veracrypt6/projects/SEEG/old_database/derivatives/seega_coordinates',r'/media/veracrypt6/projects/SEEG/derivatives/seega_coordinates']


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

# calculate TPLE and remove entry data from dataframe
elec_data_raw['TPLE'] = np.sqrt(np.square(elec_data_raw['actualTipX'] - elec_data_raw['plannedTipX']) + 
							 np.square(elec_data_raw['actualTipY'] - elec_data_raw['plannedTipY']) + 
							 np.square(elec_data_raw['actualTipZ'] - elec_data_raw['plannedTipZ'])).round(2)

# check for outliers in Z axis
large_depth_error = elec_data_raw[np.abs(elec_data_raw['TPLE']) > 10]
print('Outliers in Z axis (> 10mm):\n\n', large_depth_error['TPLE'])

# remove outliers
elec_data=elec_data_raw.drop(large_depth_error.index)
print('\nNew dataframe shape:', elec_data.shape)

elec_data_native=elec_data[elec_data['space']=='native']
elec_data_acpc=elec_data[elec_data['space']=='acpc']

elec_data_native['TPLE'].describe().round(2)
elec_data_acpc['TPLE'].describe().round(2)


vals,cnts=np.unique(elec_data_native['electrode'],return_counts=True)
label_counts=pd.DataFrame({'label':vals,'count':cnts})

label_counts_mask = label_counts[np.abs(label_counts['count']) < 10]
label_counts=label_counts.drop(label_counts_mask.index).reset_index(drop=True)

errors=[]
for ilabel in label_counts['label']:
	errors.append(elec_data_native[elec_data_native['electrode']==ilabel]['TPLE'].mean())

label_counts['error']=errors

fig = plt.figure(figsize=(16,14))
ax = fig.add_subplot(111, projection='3d')
scatter = ax.scatter(elec_data_acpc['actualTipX'].values, elec_data_acpc['actualTipY'].values, elec_data_acpc['actualTipZ'].values,s=2)
scatter = ax.scatter(elec_data_acpc['actualEntryX'].values, elec_data_acpc['actualEntryY'].values, elec_data_acpc['actualEntryZ'].values,s=2)

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




















