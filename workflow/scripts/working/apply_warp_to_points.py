#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess
import re
import pandas as pd
import csv


def run_command(cmdLineArguments):
	subprocess.run(cmdLineArguments, stdout=subprocess.PIPE,stderr=subprocess.STDOUT, shell=True)

def determineFCSVCoordSystem(input_fcsv):
	# need to determine if file is in RAS or LPS
	# loop through header to find coordinate system
	coordFlag = re.compile('# CoordinateSystem')
	verFlag = re.compile('# Markups fiducial file version')
	headFlag = re.compile('# columns')
	coord_sys=None
	headFin=None
	ver_fin=None
	
	with open(input_fcsv, 'r') as myfile:
		firstNlines=myfile.readlines()[0:3]
	
	for row in firstNlines:
		row=re.sub("[\s\,]+[\,]","",row).replace("\n","")
		cleaned_dict={row.split('=')[0].strip():row.split('=')[1].strip()}
		if None in list(cleaned_dict):
			cleaned_dict['# columns'] = cleaned_dict.pop(None)
		if any(coordFlag.match(x) for x in list(cleaned_dict)):
			coord_sys = list(cleaned_dict.values())[0]
		if any(verFlag.match(x) for x in list(cleaned_dict)):
			verString = list(filter(verFlag.match,  list(cleaned_dict)))
			assert len(verString)==1
			ver_fin = verString[0].split('=')[-1].strip()
		if any(headFlag.match(x) for x in list(cleaned_dict)):
			headFin=list(cleaned_dict.values())[0].split(',')
	
	return coord_sys,headFin

def convertSlicerRASFCSVtoAntsLPSCSV( input_fcsv, output_csv, coord_system):
	df = pd.read_csv(input_fcsv, skiprows=2, usecols=['x','y','z']) # first 2 rows of fcsv not necessary for header
	if any(x in coord_system for x in {'RAS','0'}):
		df['x'] = -1 * df['x'] # flip orientation in x
		df['y'] = -1 * df['y'] # flip orientation in y
	
	# need to add extra column 't' for ANTS
	df['t'] = [0] * df['x'].shape[0]
	df.to_csv( output_csv, index=False )

def convertAntsLPSCSVtoSlicerRASFCSV( input_csv, output_fcsv, ref_fcsv, coord_system):
	f = open( ref_fcsv, 'r' )
	lines = f.readlines()
	f.close()

	# orienting the image image back to RAS from LPS
	input_df = pd.read_csv( input_csv, usecols=['x','y','z'] ) # use reference fcsv as template
	df = pd.read_csv( ref_fcsv, skiprows=2 ) # use reference fcsv as template
	
	if any(x in coord_system for x in {'RAS','0'}):
		df['x'] = -1 * input_df['x'] # flip orientation in x
		df['y'] = -1 * input_df['y'] # flip orientation in y
	else:
		df['x'] = input_df['x']
		df['y'] = input_df['y']
		
	df['z'] = input_df['z'] # normal orientation in z
	df.to_csv( output_fcsv, index=False )

	# add in extracted Slicer header
	with open( output_fcsv, 'r+' ) as f:
		old = f.read() # read all the old csv file info
		f.seek(0) # rewind, start at the top
		f.write( lines[0] + lines[1] + old ) # add expected Slicer header

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
	
	input=dotdict({
		'fcsv': '/home/greydon/Documents/data/SEEG_peds/derivatives/seeg_coordinates/sub-P016/sub-P016_space-native_SEEGA.fcsv',
		'xfm_composite':'/home/greydon/Documents/data/SEEG_peds/derivatives/atlasreg/sub-P016/sub-P016_from-subject_to-MNIPediatricAsymCohort6_InverseComposite.h5',
	})
	
	output=dotdict({
		'fcsv_fname_warped':'/home/greydon/Documents/data/SEEG_peds/derivatives/seeg_coordinates/sub-P016/sub-P016_space-MNIPediatricAsymCohort6_SEEGA.fcsv',
	})
	
	snakemake = Namespace(output=output, input=input)

fcsv_type=os.path.splitext(os.path.basename(snakemake.input.fcsv))[0].split("_")[-1]
tmp_slicer_to_LPS_csv = os.path.join(os.path.dirname(snakemake.output.fcsv_fname_warped), f"tmp_slicer_to_LPS_{fcsv_type}.csv")
tmp_slicer_to_LPS_transformed_csv = os.path.join(os.path.dirname(snakemake.output.fcsv_fname_warped),  f"tmp_slicer_to_LPS_transformed-warp_{fcsv_type}.csv")

coordSys,headFin=determineFCSVCoordSystem(snakemake.input.fcsv)
convertSlicerRASFCSVtoAntsLPSCSV(snakemake.input.fcsv, tmp_slicer_to_LPS_csv,coordSys)

cmd = ' '.join([f'{snakemake.params.ants+snakemake.params.ants_ext}',
	  '-d', str(3),
	  '-i', '"'+tmp_slicer_to_LPS_csv+'"',
	  '-o', '"'+tmp_slicer_to_LPS_transformed_csv+'"',
	  '-t','"'+snakemake.input.xfm_composite+'"'])

run_command(cmd)

convertAntsLPSCSVtoSlicerRASFCSV(tmp_slicer_to_LPS_transformed_csv, snakemake.output.fcsv_fname_warped, snakemake.input.fcsv, coordSys)

os.remove(tmp_slicer_to_LPS_csv)
os.remove(tmp_slicer_to_LPS_transformed_csv)


