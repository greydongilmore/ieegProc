#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
import pandas as pd
import numpy as np
import glob
import math
from statistics import NormalDist
import json
import dataframe_image as dfi
os.chdir('/home/greydon/Documents/GitHub/seeg2bids-pipeline/workflow/scripts/working')
from helpers import determineFCSVCoordSystem,determine_groups,norm_vec,mag_vec


##################################
## Metrics Function Definitions ##
##################################


#Euclidean distance calculation
def euclidianDistanceCalc(xyz_planned, xyz_actual):
	if xyz_planned.ndim>1:
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
	if xyz_entry.ndim>1:
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
	if x_entry.ndim>1:
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
	if a_entry.ndim>1:
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

def mean_confidence_interval(data, confidence=0.95):
	dist = NormalDist.from_samples(data[~np.isnan(data)])
	z = NormalDist().inv_cdf((1 + confidence) / 2.)
	h = dist.stdev * z / ((len(data) - 1) ** .5)
	return dist.mean - h, dist.mean + h


chan_label_dic ={
	3:"RD10R-SP03X",
	4:"RD10R-SP04X",
	5:"RD10R-SP05X",
	6:"RD10R-SP06X",
	7:"RD10R-SP07X"
}


controlpoints_dict={
	"id": "", 	#"vtkMRMLMarkupsFiducialNode_1",
	"label": "",
	"description": "",
	"associatedNodeID": "", #"vtkMRMLScalarVolumeNode1",
	"position": [],
	"orientation": [-1.0, -0.0, -0.0, -0.0, -1.0, -0.0, 0.0, 0.0, 1.0],
	"selected": True,
	"locked": True,
	"visibility": True,
	"positionStatus": "defined"
}





#%%

debug = True
write_lines = False


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
	#data_dir=r'/media/data/data/SEEG/derivatives/seega_scenes'
	data_dir=r'/home/greydon/Documents/data/SEEG/derivatives/seega_scenes'
	
	input=dotdict({
				'isub': isub,
				'data_dir':data_dir,
				})
	params=dotdict({
				'sample_line': r'/home/greydon/Documents/GitHub/seeg2bids-pipeline/resources/sample_line.mrk.json',
				})
	output=dotdict({
		'out_svg':f'{data_dir}/{isub}/{isub}_errors.svg',
		'out_excel':f'{data_dir}/{isub}/{isub}_error_metrics.xlsx',
	})
	
	snakemake = Namespace(output=output, params=params,input=input)

if write_lines:
	with open(snakemake.params.sample_line) as (file):
		sample_line = json.load(file)



isub = snakemake.input.isub
data_dir = snakemake.input.data_dir

patient_files = glob.glob(f"{os.path.join(data_dir,isub)}/*csv")


file_data={}
for ifile in [x for x in patient_files if not x.endswith('empty.csv')]:
	determineFCSVCoordSystem(ifile)
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
	elif ifile.lower().endswith('acpc.fcsv'):
		file_data['acpc']=fcsv_data

groupsPlanned, planned_all = determine_groups(np.array(file_data['planned']['label'].values))
label_set=set(groupsPlanned)

if 'actual' in list(file_data):
	groupsActual, actual_all = determine_groups(np.array(file_data['actual']['label'].values))
	label_set=set(groupsActual).intersection(groupsPlanned)

if 'seega' in list(file_data):
	groupsSeega, seega_all = determine_groups(np.array(file_data['seega']['label'].values), True)

mcp_point=None

if 'acpc' in list(file_data):
	ac_point = file_data['acpc'][file_data['acpc']['label'].str.lower() == 'ac'][['x','y','z']].values
	pc_point = file_data['acpc'][file_data['acpc']['label'].str.lower() == 'pc'][['x','y','z']].values
	mcp_point = ((ac_point+pc_point)/2)[0]


elec_data=[]
vtk_cnt=1
for igroup in label_set:
	elec_temp={}
	elec_temp['subject']=isub
	elec_temp['electrode']=igroup
	elec_temp['group']=igroup[1:]
	elec_temp['side']='L' if igroup.startswith('L') else 'R'
	
	if 'seega' in list(file_data):
		seeg_idx=[i for i,x in enumerate(file_data['seega']['label'].values) if x.startswith(igroup)]
		elec_data_temp = file_data['seega'].loc[seeg_idx,['x','y','z']].to_numpy()
		dist = np.mean(np.linalg.norm(elec_data_temp[:-1,:] - elec_data_temp[1:,:], axis=1))
		idx_ielec,val_ielec = min(enumerate(list(chan_label_dic)), key=lambda x: abs(x[1]-dist))
		elec_temp['electrodeType']=chan_label_dic[val_ielec]
		elec_temp['numContacts']=file_data['seega'].loc[seeg_idx].shape[0]
	
	if 'planned' in list(file_data):
		planned_idx=[i for i,x in enumerate(file_data['planned']['label'].values) if x.startswith(igroup)]
		elec_temp['plannedTipX']=file_data['planned'].loc[planned_idx,'x'].values[0]
		elec_temp['plannedTipY']=file_data['planned'].loc[planned_idx,'y'].values[0]
		elec_temp['plannedTipZ']=file_data['planned'].loc[planned_idx,'z'].values[0]
		elec_temp['plannedEntryX']=file_data['planned'].loc[planned_idx,'x'].values[1]
		elec_temp['plannedEntryY']=file_data['planned'].loc[planned_idx,'y'].values[1]
		elec_temp['plannedEntryZ']=file_data['planned'].loc[planned_idx,'z'].values[1]
	
	if 'actual' in list(file_data):
		actual_idx=[i for i,x in enumerate(file_data['actual']['label'].values) if x.startswith(igroup)]
		elec_temp['actualTipX']=file_data['actual'].loc[actual_idx,'x'].values[0]
		elec_temp['actualTipY']=file_data['actual'].loc[actual_idx,'y'].values[0]
		elec_temp['actualTipZ']=file_data['actual'].loc[actual_idx,'z'].values[0]
		elec_temp['actualEntryX']=file_data['actual'].loc[actual_idx,'x'].values[1]
		elec_temp['actualEntryY']=file_data['actual'].loc[actual_idx,'y'].values[1]
		elec_temp['actualEntryZ']=file_data['actual'].loc[actual_idx,'z'].values[1]
		
		#need to account for Ad-Tech electrodes encapsulation at the tip. Planned target is electrode tip but
		#actual tip is the edge of the first contact
		mag = mag_vec(file_data['planned'].loc[planned_idx,['x','y','z']].values[0],
		  file_data['planned'].loc[planned_idx,['x','y','z']].values[1])
		norm = norm_vec(file_data['planned'].loc[planned_idx,['x','y','z']].values[0],
		  file_data['planned'].loc[planned_idx,['x','y','z']].values[1])
		plannedTipOffset=file_data['planned'].loc[planned_idx,['x','y','z']].values[1]-(norm*(mag-1))
		
		elec_temp['plannedOffsetX']=elec_temp['plannedTipX']
		elec_temp['plannedOffsetY']=elec_temp['plannedTipY']
		elec_temp['plannedOffsetZ']=elec_temp['plannedTipZ']
		
		xyz_planned_entry = np.array([elec_temp['plannedEntryX'], elec_temp['plannedEntryY'], elec_temp['plannedEntryZ']])
		xyz_actual_entry = np.array([elec_temp['actualEntryX'], elec_temp['actualEntryY'], elec_temp['actualEntryZ']]).T
		xyz_planned_target = np.array([elec_temp['plannedOffsetX'], elec_temp['plannedOffsetY'], elec_temp['plannedOffsetZ']]).T
		xyz_actual_target = np.array([elec_temp['actualTipX'], elec_temp['actualTipY'], elec_temp['actualTipZ']]).T
		
		elec_temp['euclid_dist_target'] = euclidianDistanceCalc(xyz_planned_target, xyz_actual_target)
		elec_temp['euclid_dist_entry'] = euclidianDistanceCalc(xyz_planned_entry, xyz_actual_entry)
		elec_temp['radial_dist_target'] = radialDistanceCalc(xyz_actual_target, xyz_planned_entry, xyz_planned_target)
		elec_temp['radial_dist_entry'] = radialDistanceCalc(xyz_actual_entry, xyz_planned_entry, xyz_planned_target)
		
		if not np.array_equal(np.round(xyz_actual_target,2), np.round(xyz_planned_target,2)):
			try:
				elec_temp['radial_angle'] = ptLineAngleCalc(xyz_actual_target, xyz_planned_entry, xyz_planned_target)
				elec_temp['line_angle'] = lineLineAngleCalc(xyz_actual_entry, xyz_actual_target, xyz_planned_entry, xyz_planned_target)
			except:
				pass
	
	if mcp_point is not None:
		if 'actual' in list(file_data):
			elec_temp['actualTipX_mcp'] = elec_temp['actualTipX']-mcp_point[0]
			elec_temp['actualTipY_mcp']= elec_temp['actualTipY']-mcp_point[1]
			elec_temp['actualTipZ_mcp']= elec_temp['actualTipZ']-mcp_point[2]
			elec_temp['actualEntryX_mcp']= elec_temp['actualEntryX']-mcp_point[0]
			elec_temp['actualEntryY_mcp']= elec_temp['actualEntryY']-mcp_point[1]
			elec_temp['actualEntryZ_mcp']= elec_temp['actualEntryZ']-mcp_point[2]
		if 'planned' in list(file_data):
			elec_temp['plannedTipX_mcp']= elec_temp['plannedOffsetX']-mcp_point[0]
			elec_temp['plannedTipY_mcp']= elec_temp['plannedOffsetY']-mcp_point[1]
			elec_temp['plannedTipZ_mcp']= elec_temp['plannedOffsetZ']-mcp_point[2]
			elec_temp['plannedEntryX_mcp']= elec_temp['plannedEntryX']-mcp_point[0]
			elec_temp['plannedEntryY_mcp']= elec_temp['plannedEntryY']-mcp_point[1]
			elec_temp['plannedEntryZ_mcp']= elec_temp['plannedEntryZ']-mcp_point[2]
	
	if write_lines:
		for itype in ['actual','planned']:
			if itype in list(file_data):
				line_out=sample_line.copy()
				entryp=controlpoints_dict.copy()
				targetp=controlpoints_dict.copy()
				
				entryp["id"]= f"vtkMRMLMarkupsFiducialNode{vtk_cnt}"
				entryp["associatedNodeID"]= f"vtkMRMLScalarVolumeNode{vtk_cnt}"
				entryp["label"]=igroup
				entryp["position"]=[elec_temp[f'{itype}EntryX'],elec_temp[f'{itype}EntryY'],elec_temp[f'{itype}EntryZ']]
				
				targetp["id"]= f"vtkMRMLMarkupsFiducialNode{vtk_cnt+1}"
				targetp["associatedNodeID"]= f"vtkMRMLScalarVolumeNode{vtk_cnt+1}"
				targetp["label"]=igroup
				targetp["position"]=[elec_temp[f'{itype}TipX'],elec_temp[f'{itype}TipY'],elec_temp[f'{itype}TipZ']]
				
				line_out['markups'][0]["controlPoints"]=[targetp,entryp]
				
				json_output = json.dumps(line_out, indent=4)
				with open(os.path.join(data_dir, isub,f'{igroup}_{itype}.mrk.json'), 'w') as (fid):
					fid.write(json_output)
					fid.write('\n')
				
				vtk_cnt+=2
	
	elec_data.append(elec_temp)

elec_data_raw=pd.DataFrame(elec_data)



elec_table=elec_data_raw[['electrode','euclid_dist_target', 'radial_dist_target', 'euclid_dist_entry','radial_dist_entry','radial_angle','line_angle']].round(2)
for item in list(elec_table)[1:]:
	elec_table[item]=elec_table[item].astype(float)

elec_table = elec_table.sort_values('electrode').set_index('electrode').round(2)
elec_table_styled=elec_table.style.applymap(lambda x: "background-color:#ccffcc;" if x<2 else 'background-color:#ffff00;' if x>=2 and x<3 else "background-color:#ffcccc;")

writer = pd.ExcelWriter(snakemake.output.out_excel, engine='openpyxl')
elec_table_styled.to_excel(writer,sheet_name='Sheet1', float_format='%.2f')
book = writer.book
book._named_styles['Normal'].number_format = '#,##0.00'
writer.save()

dfi.export(elec_table_styled, snakemake.output.out_svg)


