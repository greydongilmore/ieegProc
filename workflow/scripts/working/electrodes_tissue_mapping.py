#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 05:22:25 2025

@author: greydon
"""

import pandas as pd
import numpy as np
import re,os
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils.dataframe import dataframe_to_rows
from openpyxl.utils import get_column_letter
pd.set_option("future.no_silent_downcasting", True)

debug = False

# Define styles
blue_font = Font(color="0000FF")
red_font = Font(color="FF0000")
black_font = Font(color="000000")
grey_fill = PatternFill(start_color="A9A9A9", end_color="A9A9A9", fill_type="solid")
alignment = Alignment(horizontal="center", vertical="center")
border = Border(
	left=Side(style="thin"), right=Side(style="thin"),
	top=Side(style="thin"), bottom=Side(style="thin")
)

# Helper function to create ranges
def create_ranges(labels):
	leads = [re.match(r'^(.*?)-', label).group(1) for label in labels]
	contacts = [int(re.match(r'.*-(\d+)$', label).group(1)) for label in labels]
	unique_leads = set(leads)
	ranges = []
	for lead in unique_leads:
		lead_contacts = [contacts[i] for i in range(len(contacts)) if leads[i] == lead]
		sorted_contacts = sorted(lead_contacts)
		start_idx = 0
		for i in range(1, len(sorted_contacts)):
			if sorted_contacts[i] != sorted_contacts[i-1] + 1:
				ranges.append(f"{lead}-{sorted_contacts[start_idx]:02d}" if start_idx == i-1 else f"{lead}-{sorted_contacts[start_idx]:02d}-{sorted_contacts[i-1]:02d}")
				start_idx = i
		# Handle the last segment
		ranges.append(f"{lead}-{sorted_contacts[start_idx]:02d}-{sorted_contacts[-1]:02d}" if start_idx != len(sorted_contacts)-1 else f"{lead}-{sorted_contacts[start_idx]:02d}")
	return ', '.join(ranges)


if debug:
	class dotdict(dict):
		"""dot.notation access to dictionary attributes"""
		__getattr__ = dict.get
		__setattr__ = dict.__setitem__
		__delattr__ = dict.__delitem__
	
	class Namespace:
		def __init__(self, **kwargs):
			self.__dict__.update(kwargs)
	
	input_file=r"/home/greydon/Documents/data/emory_seeg/derivatives/atlasreg/sub-EMOP0249/sub-EMOP0249_desc-nonlin_atlas-CerebrA_from-MNI152NLin2009cSym_electrodes.xlsx"
	
	input=dotdict({
				'input_file': input_file,
				})
	output=dotdict({
		'out_excel': os.path.join(os.path.dirname(input_file), os.path.basename(input_file).replace("_electrodes","_tissue_map")),
	})
	
	snakemake = Namespace(output=output, input=input)


#%%


# path to the *_electrodes.xlsx file
#input_file=r"/home/greydon/Documents/data/emory_seeg/derivatives/atlasreg/sub-EMOP0249/sub-EMOP0249_desc-nonlin_atlas-CerebrA_from-MNI152NLin2009cSym_electrodes.xlsx"

input_file=snakemake.input.input_file

#out_excel=os.path.join(os.path.dirname(input_file), os.path.basename(input_file).replace("_electrodes","_tissue_map"))

out_excel=snakemake.output.out_excel

if not os.path.exists(out_excel):
	# Read the first sheet of the uploaded Excel file
	data = pd.read_excel(input_file, sheet_name=0)
	
	# Ensure required columns exist
	required_columns = ["label", "atlas_label", "GM", "WM", "CSF"]
	if not all(col in data.columns for col in required_columns):
		raise ValueError("The input file is missing required columns: label, atlas_label, GM, WM, CSF.")
	
	# Update atlas_label for CSF and WM
	data['atlas_label'] = np.where(data['CSF'] > 0.7, "CSF", np.where(data['WM'] > 0.7, "WM", data['atlas_label']))
	
	# Extract leads and contacts
	data['contact_number'] = data['label'].str.extract(r'-(\d+)$').astype(float)
	max_contacts = int(data['contact_number'].max())
	leads = data['label'].str.extract(r'^(.*?)-')[0].unique()
	
	reformatted_table_tmp={}
	# Populate reformatted table
	for lead_idx,lead in enumerate(leads):
		lead_data = data[data['label'].str.startswith(lead + '-')]
		reformatted_table_tmp[lead]=[]
		for _, row in lead_data.iterrows():
			reformatted_table_tmp[lead].append(str(row['atlas_label']))
		
		while len(reformatted_table_tmp[lead]) < max_contacts:
			reformatted_table_tmp[lead].append(np.nan)
	
	reformatted_table = pd.DataFrame.from_dict(reformatted_table_tmp,orient='index')
	reformatted_table = reformatted_table.rename(columns=dict(zip(reformatted_table.columns, [str(i) for i in range(1, max_contacts + 1)])))
	reformatted_table=reformatted_table.reset_index(names='lead')
	reformatted_table = reformatted_table.fillna('NaN')
	reformatted_table = reformatted_table.replace('nan', np.nan, regex=True)
	reformatted_table = reformatted_table.replace('NaN', np.nan, regex=True)
	reformatted_table = reformatted_table.fillna('nan')
	
	# Extract CSF and WM electrodes
	csf_data = data[data['atlas_label'] == "CSF"]
	wm_data = data[data['atlas_label'] == "WM"]
	csf_list = create_ranges(csf_data['label'])
	wm_list = create_ranges(wm_data['label'])
	
	# Generate summary statistics
	total_electrodes = len(reformatted_table)
	total_contacts = reformatted_table.notna().sum().sum() - total_electrodes
	
	# Right side is marked a few different mays over time (R*, *Rd, *')
	# need to search for all markers
	right_electrodes = sum(reformatted_table['lead'].str.startswith('R'))
	if right_electrodes <=1:
		right_electrodes = sum(reformatted_table['lead'].str.endswith("'"))
		if right_electrodes <=1:
			right_electrodes = sum(reformatted_table['lead'].str.endswith("Rd"))
			r_idx=[i for i,x in enumerate(reformatted_table['lead'].str.endswith("Rd")) if x ==True]
			right_contacts = reformatted_table[reformatted_table['lead'].str.endswith("Rd")].notna().sum().sum() - right_electrodes
		else:
			r_idx=[i for i,x in enumerate(reformatted_table['lead'].str.endswith("'")) if x ==True]
			right_contacts = reformatted_table[reformatted_table['lead'].str.endswith("'")].notna().sum().sum() - right_electrodes
	else:
		r_idx=[i for i,x in enumerate(reformatted_table['lead'].str.startswith('R')) if x ==True]
		right_contacts = reformatted_table[reformatted_table['lead'].str.startswith('R')].notna().sum().sum() - right_electrodes
	
	left_electrodes = total_electrodes-right_electrodes
	if right_electrodes==0:
		left_contacts = reformatted_table.notna().sum().sum() - left_electrodes
		l_idx=list(np.arange(0,left_electrodes,1))
	else:
		l_idx=list(np.arange(0,total_electrodes,1))
		l_idx=list(set(l_idx).difference(set(r_idx)))
	
	csf_contacts = len(csf_data)
	wm_contacts = len(wm_data)
	
	
	wb = Workbook()
	ws = wb.active
	
	for row in dataframe_to_rows(reformatted_table, index=False, header=True):
		ws.append(row)
	
	# Merge contiguous identical cells
	for row_idx in range(2, len(reformatted_table) + 1):
		start_col = None
		for col_idx in range(2, max_contacts + 2):
			current_value = ws.cell(row=row_idx+1, column=col_idx).value
			prev_value = ws.cell(row=row_idx+1, column=col_idx - 1).value if col_idx > 2 else None
			
			if prev_value is None or current_value != prev_value:
				if start_col is not None and col_idx > start_col + 1:
					ws.merge_cells(start_row=row_idx+1, start_column=start_col, end_row=row_idx+1, end_column=col_idx - 1)
				start_col = col_idx
		if start_col is not None and start_col <= max_contacts + 1:
			ws.merge_cells(start_row=row_idx+1, start_column=start_col, end_row=row_idx+1, end_column=max_contacts + 1)
	
	for col_idx in range(1, reformatted_table.shape[1] + 1):
		cell = ws.cell(row=1, column=col_idx)
		cell.font = black_font
		cell.alignment = alignment
		cell.border = border
	
	for row_idx in range(1,len(reformatted_table)+2):
		cell = ws.cell(row=row_idx, column=1)
		cell.font = black_font
		cell.alignment = alignment
		cell.border = border
		if row_idx-1 in l_idx:
				for col_idx in range(2, reformatted_table.shape[1] + 1):
						cell = ws.cell(row=row_idx+1, column=col_idx)
						cell.font = blue_font
						cell.alignment = alignment
						cell.border = border
		elif row_idx-1 in r_idx:
				for col_idx in range(2, reformatted_table.shape[1] + 1):
						cell = ws.cell(row=row_idx+1, column=col_idx)
						cell.font = red_font
						cell.alignment = alignment
						cell.border = border
		# Grey fill for NA cells
		for col_idx in range(2, reformatted_table.shape[1] + 1):
			if ws.cell(row=row_idx + 1, column=col_idx).value=='nan':
					cell = ws.cell(row=row_idx+1, column=col_idx)
					cell.fill = grey_fill
					cell.alignment = alignment
					cell.border = border
					cell.font = black_font
	
	column_widths = []
	for row in ws.iter_rows():
		for i, cell in enumerate(row):
			try:
				column_widths[i] = max(column_widths[i], len(str(cell.value)))
			except IndexError:
				column_widths.append(len(str(cell.value)))
	
	# set min column width to 6
	column_widths=[x if x > 6 else 6 for x in column_widths]
	
	for col, value in enumerate(column_widths):
		ws.column_dimensions[get_column_letter(col+1)].width = value
	
	# Save the workbook to a temporary file
	wb.save(out_excel)
