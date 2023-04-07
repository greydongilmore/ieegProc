#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import csv
from collections import ChainMap
import pandas as pd
import numpy as np
import os


def determineFCSVCoordSystem(input_fcsv):
	# need to determine if file is in RAS or LPS
	# loop through header to find coordinate system
	coordFlag = re.compile('# CoordinateSystem')
	verFlag = re.compile('# Markups fiducial file version')
	headFlag = re.compile('# columns')
	coord_sys=None
	headFin=None
	ver_fin=None
	with open(input_fcsv, 'r+') as fid:
		rdr = csv.DictReader(filter(lambda row: row[0]=='#', fid))
		row_cnt=0
		for row in rdr:
			cleaned_dict={k:v for k,v in row.items()}
			if None in list(cleaned_dict):
				cleaned_dict['# columns'] = cleaned_dict.pop(None)
			
			if any(coordFlag.match(x) for x in list(cleaned_dict.values()) if not isinstance(x,list)):
				coordString = list(filter(coordFlag.match,  list(cleaned_dict.values())))
				assert len(coordString)==1
				coord_sys = coordString[0].split('=')[-1].strip()
			if any(verFlag.match(x) for x in list(cleaned_dict)):
				verString = list(filter(verFlag.match,  list(cleaned_dict)))
				assert len(verString)==1
				ver_fin = verString[0].split('=')[-1].strip()
			if any(headFlag.match(x) for x in list(cleaned_dict)):
				headString = list(filter(headFlag.match, list(cleaned_dict)))
				headFin=['id']+cleaned_dict[headString[0]]
			row_cnt +=1
	
	if headFin is not None:
		headFin=dict(ChainMap(*[{i:x} for i,x in enumerate(headFin)]))
	
	if any(x in coord_sys for x in {'LPS','1'}):
		df = pd.read_csv(input_fcsv, skiprows=3, header=None)
		df[1] = -1 * df[1] # flip orientation in x
		df[2] = -1 * df[2] # flip orientation in y
		
		with open(input_fcsv, 'w') as fid:
			fid.write("# Markups fiducial file version = 4.11\n")
			fid.write("# CoordinateSystem = 0\n")
			fid.write("# columns = id,x,y,z,ow,ox,oy,oz,vis,sel,lock,label,desc,associatedNodeID\n")
		
		df.rename(columns={0:'node_id', 1:'x', 2:'y', 3:'z', 4:'ow', 5:'ox',
							6:'oy', 7:'oz', 8:'vis', 9:'sel', 10:'lock',
							11:'label', 12:'description', 13:'associatedNodeID'}, inplace=True)
		
		df['associatedNodeID']= pd.Series(np.repeat('',df.shape[0]))
		df.round(3).to_csv(input_fcsv, sep=',', index=False, line_terminator="", mode='a', header=False)
		
		print(f"Converted LPS to RAS: {os.path.dirname(input_fcsv)}/{os.path.basename(input_fcsv)}")
	return ver_fin,headFin

def determine_groups(iterable, numbered_labels=False):
	values = []
	for item in iterable:
		temp=None
		if re.findall(r"([a-zA-Z]+)([0-9]+)([a-zA-Z]+)", item):
			temp = "".join(list(re.findall(r"([a-zA-Z]+)([0-9]+)([a-zA-Z]+)", item)[0]))
		elif '-' in item:
			temp=item.split('-')[0]
		else:
			if numbered_labels:
				temp=''.join([x for x in item if not x.isdigit()])
				for sub in ("T1","T2"):
					if sub in item:
						temp=item.split(sub)[0] + sub
			else:
				temp=item
		if temp is None:
			temp=item
		
		values.append(temp)
	
	vals,indexes,count = np.unique(values, return_index=True, return_counts=True)
	values_unique = [values[index] for index in sorted(indexes)]
	
	return values_unique,vals

def sorted_nicely(data, reverse = False):
	convert = lambda text: int(text) if text.isdigit() else text
	alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
	
	return sorted(data, key = alphanum_key, reverse=reverse)

def normDir(direction):
	return np.array(direction) / np.linalg.norm(np.array(direction))

def mag_vec(P1, P2):
	if isinstance(P1, list):
		P1 = np.array(P1)
	if isinstance(P1, list):
		P2 = np.array(P2)
	DirVec = P2-P1
	MagVec = np.sqrt([np.square(DirVec[0]) + np.square(DirVec[1]) + np.square(DirVec[2])])
	return MagVec

def norm_vec(P1, P2):
	if isinstance(P1, list):
		P1 = np.array(P1)
	if isinstance(P2, list):
		P2 = np.array(P2)
	DirVec = P2-P1
	MagVec = np.sqrt([np.square(DirVec[0]) + np.square(DirVec[1]) + np.square(DirVec[2])])
	NormVec = np.array([float(DirVec[0] / MagVec), float(DirVec[1] / MagVec), float(DirVec[2] / MagVec)])
	return NormVec

#%%

electrodeModels = {}

RD10RSP03 = {
		'num_groups': 10,
		'num_contacts': 10,
		'encapsultation': 1.0,
		'lead_shift':1.5,
		'lead_tail': 15,
		'contact_size': 2.29,
		'contact_spacing': [0.71,0.71,0.71,0.71,0.71,0.71,0.71,0.71,0.71],
		'diameter': 0.86,
		'electrode_1': [1,2,3,4,5,6,7,8,9,10],
		'electrode_2': [1,2,3,4,5,6,7,8,9,10],
		'contact_label':['','','','','','','', ''],
		'lead_type': 'linear',
		'filename': 'rd10rsp03'
	 }
	 
electrodeModels['3 mm'] = RD10RSP03

RD10RSP04 = {
		'num_groups': 10,
		'num_contacts': 10,
		'encapsultation': 1.0,
		'lead_shift':1.5,
		'lead_tail': 15,
		'contact_size': 2.29,
		'contact_spacing': [1.71,1.71,1.71,1.71,1.71,1.71,1.71,1.71,1.71],
		'diameter': 0.86,
		'electrode_1': [1,2,3,4,5,6,7,8,9,10],
		'electrode_2': [1,2,3,4,5,6,7,8,9,10],
		'contact_label':['','','','','','','', ''],
		'lead_type': 'linear',
		'filename': 'rd10rsp04'
	 }
	 
electrodeModels['4 mm'] = RD10RSP04

RD10RSP05 = {
		'num_groups': 10,
		'num_contacts': 10,
		'encapsultation': 1.0,
		'lead_shift':1.5,
		'lead_tail': 15,
		'contact_size': 2.29,
		'contact_spacing': [2.71,2.71,2.71,2.71,2.71,2.71,2.71,2.71,2.71],
		'diameter': 0.86,
		'electrode_1': [1,2,3,4,5,6,7,8,9,10],
		'electrode_2': [1,2,3,4,5,6,7,8,9,10],
		'contact_label':['','','','','','','', ''],
		'lead_type': 'linear',
		'filename': 'rd10rsp05'
	 }
	 
electrodeModels['5 mm'] = RD10RSP05

RD10RSP06 = {
		'num_groups': 10,
		'num_contacts': 10,
		'encapsultation': 1.0,
		'lead_shift':1.5,
		'lead_tail': 15,
		'contact_size': 2.29,
		'contact_spacing':[3.71,3.71,3.71,3.71,3.71,3.71,3.71,3.71,3.71],
		'diameter': 0.86,
		'electrode_1': [1,2,3,4,5,6,7,8,9,10],
		'electrode_2': [1,2,3,4,5,6,7,8,9,10],
		'contact_label':['','','','','','','', ''],
		'lead_type': 'linear',
		'filename': 'rd10rsp06'
	 }
	 
electrodeModels['6 mm'] = RD10RSP06

RD10RSP07 = {
		'num_groups': 10,
		'num_contacts': 10,
		'encapsultation': 1.0,
		'lead_shift':1.5,
		'lead_tail': 15,
		'contact_size': 2.29,
		'contact_spacing': [4.71,4.71,4.71,4.71,4.71,4.71,4.71,4.71,4.71],
		'diameter': 0.86,
		'electrode_1': [1,2,3,4,5,6,7,8,9,10],
		'electrode_2': [1,2,3,4,5,6,7,8,9,10],
		'contact_label':['','','','','','','', ''],
		'lead_type': 'linear',
		'filename': 'rd10rsp07'
	 }
	 
electrodeModels['7 mm'] = RD10RSP07

MM16DSP05 = {
		'num_groups': 8,
		'num_contacts': 8,
		'encapsultation': 1.5,
		'lead_shift':1.5,
		'lead_tail': 15,
		'contact_size': 2.29,
		'contact_spacing': [2.71,2.71,2.71,2.71,2.71,2.71,2.71,2.71,2.71],
		'diameter': 1.3,
		'electrode_1': [1,2,3,4,5,6,7,8],
		'electrode_2': [1,2,3,4,5,6,7,8],
		'contact_label':['','','','','','','', ''],
		'lead_type': 'linear',
		'filename': 'mm16dsp05'
	 }
	 
electrodeModels['MM16D-SP05X'] = MM16DSP05

BF09RSP51X = {
		'num_groups': 9,
		'num_contacts': 9,
		'encapsultation': 1.5,
		'lead_shift':1.5,
		'lead_tail': 15,
		'contact_size': 1.6,
		'contact_spacing': [1.4, 3.4, 3.4, 3.4, 3.4, 3.4, 3.4, 3.4],
		'diameter': 1.3,
		'electrode_1': [1,2,3,4,5,6,7,8],
		'electrode_2': [1,2,3,4,5,6,7,8],
		'contact_label':['','','','','','','', ''],
		'lead_type': 'linear',
		'filename': 'bf09rsp51x'
	 }
	 
electrodeModels['BF09R-SP51X'] = BF09RSP51X

D0805AM = {
		'num_groups': 5,
		'num_contacts': 5,
		'encapsultation': 0,
		'lead_shift':0,
		'lead_tail': 15,
		'contact_size': 2.0,
		'contact_spacing': [1.5, 1.5, 1.5, 1.5],
		'diameter': 0.8,
		'electrode_1': [1,2,3,4,5],
		'electrode_2': [1,2,3,4,5],
		'contact_label':['','','','',''],
		'lead_type': 'linear',
		'filename': 'd0805am'
	 }

electrodeModels['5 contact'] = D0805AM

D0808AM = {
		'num_groups': 8,
		'num_contacts': 8,
		'encapsultation': 0,
		'lead_shift':0,
		'lead_tail': 15,
		'contact_size': 2.0,
		'contact_spacing': [1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5],
		'diameter': 0.8,
		'electrode_1': [1,2,3,4,5,6,7,8],
		'electrode_2': [1,2,3,4,5,6,7,8],
		'contact_label':['','','','','','','',''],
		'lead_type': 'linear',
		'filename': 'd0808am'
	 }

electrodeModels['8 contact'] = D0808AM

D0810AM = {
		'num_groups': 10,
		'num_contacts': 10,
		'encapsultation': 0,
		'lead_shift':0,
		'lead_tail': 15,
		'contact_size': 2.0,
		'contact_spacing': [1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5],
		'diameter': 0.8,
		'electrode_1': [1,2,3,4,5,6,7,8,9,10],
		'electrode_2': [1,2,3,4,5,6,7,8,9,10],
		'contact_label':['','','','','','','','','',''],
		'lead_type': 'linear',
		'filename': 'd0810am'
	 }

electrodeModels['10 contact'] = D0810AM

D0812AM = {
		'num_groups': 12,
		'num_contacts': 12,
		'encapsultation': 0,
		'lead_shift':0,
		'lead_tail': 15,
		'contact_size': 2.0,
		'contact_spacing': [1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5],
		'diameter': 0.8,
		'electrode_1': [1,2,3,4,5,6,7,8,9,10,11,12],
		'electrode_2': [1,2,3,4,5,6,7,8,9,10,11,12],
		'contact_label':['','','','','','','','','','','',''],
		'lead_type': 'linear',
		'filename': 'd0812am'
	 }

electrodeModels['12 contact'] = D0812AM

D0815AM = {
		'num_groups': 15,
		'num_contacts': 15,
		'encapsultation': 0,
		'lead_shift':0,
		'lead_tail': 15,
		'contact_size': 2.0,
		'contact_spacing': [1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5],
		'diameter': 0.8,
		'electrode_1': [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15],
		'electrode_2': [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15],
		'contact_label':['','','','','','','','','','','','','','',''],
		'lead_type': 'linear',
		'filename': 'd0815am'
	 }

electrodeModels['15 contact'] = D0815AM

D0818AM = {
		'num_groups': 18,
		'num_contacts': 18,
		'encapsultation': 0,
		'lead_shift':0,
		'lead_tail': 18,
		'contact_size': 2.0,
		'contact_spacing': [1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5],
		'diameter': 0.8,
		'electrode_1': [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18],
		'electrode_2': [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18],
		'contact_label':['','','','','','','','','','','','','','','','','',''],
		'lead_type': 'linear',
		'filename': 'd0818am'
	 }

electrodeModels['18 contact'] = D0818AM
