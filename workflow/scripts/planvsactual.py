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
			temp=item
		
		values.append(temp)
	
	indexes = np.unique(values, return_index=True)[1]
	values = [values[index] for index in sorted(indexes)]
	
	return values

data_dir=r'/media/veracrypt6/projects/SEEG/derivatives/seega_coordinates'

patients = sorted_nicely([x for x in os.listdir(data_dir) if os.path.isdir(os.path.join(data_dir, x))])


elec_data=[]
error_files=[]
for isub in patients:
	patient_files=glob.glob(f"{os.path.join(data_dir,isub)}/*.fcsv")
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
			
		groupsActual = determine_groups(np.array(file_data['actual']['label'].values))
		groupsPlanned = determine_groups(np.array(file_data['planned']['label'].values))
		
		for igroup in set(groupsActual).intersection(groupsActual):
			try:
				elec_temp={}
				elec_temp['subject']=isub
				elec_temp['space']=ispace.split('-')[-1]
				elec_temp['electrode']=igroup
				elec_temp['electrodeType']=file_data['seega'][file_data['seega']['label'].str.contains(igroup)]['description'].values[0]
				elec_temp['numContacts']=file_data['seega'][file_data['seega']['label'].str.contains(igroup)].shape[0]
				elec_temp['actualTipX']=file_data['actual'].iloc[::2, :][file_data['actual'].iloc[::2, :]['label']==igroup]['x'].values[0]
				elec_temp['actualTipY']=file_data['actual'].iloc[::2, :][file_data['actual'].iloc[::2, :]['label']==igroup]['y'].values[0]
				elec_temp['actualTipZ']=file_data['actual'].iloc[::2, :][file_data['actual'].iloc[::2, :]['label']==igroup]['z'].values[0]
				elec_temp['actualEntryX']=file_data['actual'].iloc[1::2, :][file_data['actual'].iloc[1::2, :]['label']==igroup]['x'].values[0]
				elec_temp['actualEntryY']=file_data['actual'].iloc[1::2, :][file_data['actual'].iloc[1::2, :]['label']==igroup]['y'].values[0]
				elec_temp['actualEntryZ']=file_data['actual'].iloc[1::2, :][file_data['actual'].iloc[1::2, :]['label']==igroup]['z'].values[0]
				elec_temp['plannedTipX']=file_data['planned'].iloc[::2, :][file_data['actual'].iloc[::2, :]['label']==igroup]['x'].values[0]
				elec_temp['plannedTipY']=file_data['planned'].iloc[::2, :][file_data['actual'].iloc[::2, :]['label']==igroup]['y'].values[0]
				elec_temp['plannedTipZ']=file_data['planned'].iloc[::2, :][file_data['actual'].iloc[::2, :]['label']==igroup]['z'].values[0]
				elec_temp['plannedEntryX']=file_data['planned'].iloc[1::2, :][file_data['actual'].iloc[1::2, :]['label']==igroup]['x'].values[0]
				elec_temp['plannedEntryY']=file_data['planned'].iloc[1::2, :][file_data['actual'].iloc[1::2, :]['label']==igroup]['y'].values[0]
				elec_temp['plannedEntryZ']=file_data['planned'].iloc[1::2, :][file_data['actual'].iloc[1::2, :]['label']==igroup]['z'].values[0]
				elec_data.append(elec_temp)
			except:
				error_files.append([isub,igroup])


print(np.unique(error_files))

elec_data_raw=pd.DataFrame(elec_data)

# calculate TPLE and remove entry data from dataframe
elec_data_raw['TPLE'] = np.sqrt(np.square(elec_data_raw['actualTipX'] - elec_data_raw['plannedTipX']) + 
							 np.square(elec_data_raw['actualTipY'] - elec_data_raw['plannedTipY']) + 
							 np.square(elec_data_raw['actualTipZ'] - elec_data_raw['plannedTipZ'])).round(1)

# check for outliers in Z axis
large_depth_error = elec_data_raw[np.abs(elec_data_raw['TPLE']) > 10]
print('Outliers in Z axis (> 10mm):\n\n', large_depth_error['TPLE'])

# remove outliers
elec_data=elec_data_raw.drop(large_depth_error.index)
print('\nNew dataframe shape:', elec_data.shape)

tple = elec_data['TPLE']
tple.describe().round(1)
