#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 22:50:23 2020

@author: greydon
"""

import struct
import numpy as np
import pandas as pd
import re
import glob
import os
import shutil



def determine_groups(iterable, numbered_labels=False):
	values = []
	pattern = re.compile(r'^(\d+[A-Za-z]+)')
	
	for item in iterable:
		item=str(item).strip()
		temp=None
		if pattern.match(item):
			temp = "".join(pattern.match(item)[0])
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
	vals=vals[indexes.argsort()]
	count=count[indexes.argsort()]
	
	return vals,count

def get_montage(ifile):
	
	ignore_keys={'IChannelId','IInputId','ISiteId','ITypeId','OChannelId','OTypeId','GroupId'}
	
	mtg_file = np.fromfile(ifile, dtype='uint8')
	mtg_file_tmp = "".join([struct.unpack('s', x)[0].decode('ISO-8859-1') for x in mtg_file])
	mtg_file_tmp = re.findall(r'\(.\(..*?\)\)', mtg_file_tmp)
	
	#chan info
	chans_info = [x.replace('(.','(') for x in mtg_file_tmp if x.startswith('(.(."ChanIndex"')]
	
	chan_info=[]
	for ichan in range(len(chans_info)):
		chan_info_tmp=[eval(re.findall(r'\(.*?\)',chans_info[ichan])[0].replace('((','('))]+[eval(x) for x in re.findall(r'\(.*?\)',chans_info[ichan])[2:] if not any( y in x for y in ignore_keys)]
		chan_info.append({key: value for (key, value) in chan_info_tmp})
		
	chan_info_df=pd.DataFrame(chan_info)
	chan_info_df=chan_info_df.loc[chan_info_df["From_Name"].isin([0])].reset_index(drop=True)
	chan_info_df["From_Name"]=[str(x) for x in chan_info_df["From_Name"].values]
	chan_info_df["To_Name"]=[str(x) for x in chan_info_df["To_Name"].values]
	
	groups, n_members = determine_groups(np.array(chan_info_df['To_Name'].values),numbered_labels=True)
	
	group_lbl=[]
	zero_idx=None
	for igroup,imember in zip(groups,n_members):
		if igroup == '0' or igroup == '':
			zero_idx=[i for i,x in enumerate(chan_info_df['To_Name'].values) if x=='0']
		else:
			group_lbl.append(np.repeat(igroup,imember))
	
	group_lbl=[str(x) for x in  np.concatenate(group_lbl).tolist()]
	
	if zero_idx is not None:
		[group_lbl.insert(x, '0') for x in zero_idx]
	
	chan_info_df['Group']=group_lbl
	chan_info_df=chan_info_df[~(chan_info_df["Group"].isin(['0']))].reset_index(drop=True)
	chan_info_df['ChanIndex']=list(chan_info_df.index+1)
	
	return chan_info_df


def padtrim(buf, num):
	num -= len(buf)
	if num>=0:
		# pad the input to the specified length
		buffer = (str(buf) + ' ' * num)
	else:
		# trim the input to the specified length
		buffer = (buf[0:num])
	
	return buffer

#%%


data_dir=r'/media/greydon/lhsc_data/datasets/emory_seeg/derivatives'
isub='montages'


for isub in [x for x in os.listdir(data_dir) if os.path.isdir(os.path.join(data_dir,x))]:
	
	files=glob.glob(os.path.join(data_dir,isub,'*.mtg'))
	
	for ifile in files:
		out_info=get_montage(ifile)
		out_info.to_csv(os.path.splitext(ifile)[0]+'.tsv',sep='\t',float_format='%.3f',index=False)



#%%


subs=np.unique([x.split('_')[0] for x in os.listdir(data_dir) if os.path.isfile(os.path.join(data_dir,x))])

for isub in subs:
	if not os.path.exists(os.path.join(data_dir,isub)):
		os.makedirs(os.path.join(data_dir,isub))
	
	for ifile in glob.glob(data_dir+os.path.sep+f'{isub}*.mtg'):
		out_info=get_montage(ifile)
		new_out=os.path.join(data_dir,isub,os.path.basename(ifile))
		if not os.path.exists(new_out):
			shutil.move(ifile,new_out)
			out_info.to_csv(os.path.splitext(new_out)[0]+'.tsv',sep='\t',float_format='%.3f',index=False)
