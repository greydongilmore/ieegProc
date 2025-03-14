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


def get_montage(fname):
	
	ignore_keys={'IChannelId','IInputId','ISiteId','ITypeId','OChannelId','OTypeId','GroupId'}
	
	mtg_file = np.fromfile(fname, dtype='uint8')
	mtg_file_tmp = "".join([struct.unpack('s', x)[0].decode('ISO-8859-1') for x in mtg_file])
	mtg_file_tmp = re.findall(r'\(.\(..*?\)\)', mtg_file_tmp)
	
	#chan info
	chans_info = [x.replace('(.','(') for x in mtg_file_tmp if x.startswith('(.(."ChanIndex"')]
	
	chan_info=[]
	for ichan in range(len(chans_info)):
		chan_info_tmp=[eval(re.findall(r'\(.*?\)',chans_info[ichan])[0].replace('((','('))]+[eval(x) for x in re.findall(r'\(.*?\)',chans_info[ichan])[2:] if not any( y in x for y in ignore_keys)]
		chan_info.append({key: value for (key, value) in chan_info_tmp})
		
	chan_info=pd.DataFrame(chan_info)
	
	return chan_info

import xml.etree.ElementTree as ET
import pandas as pd

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

new_labels=pd.read_csv(r'/home/greydon/Downloads/sub-P019.txt',sep='\t',header=None)
new_labels=new_labels.set_index(0).T.to_dict('records')[0]

tree = ET.parse('/home/greydon/Downloads/my_montage.mtg')
root = tree.getroot()
for label in root.findall('./signalcomposition/signal/label'):
	if label.text.strip() in list(new_labels):
		new_lbl=padtrim(new_labels[label.text.strip()], len(label.text))
		del new_labels[label.text.strip()]
		
		label.text=new_lbl

for alias in root.findall('./signalcomposition/alias'):
	if label.text.strip() in list(new_labels):
		new_lbl=padtrim(new_labels[label.text.strip()], len(label.text))
		del new_labels[label.text.strip()]
		
		alias.text=''

ET.ElementTree(root).write('/home/greydon/Downloads/my_montage_update.mtg', encoding="UTF-8", xml_declaration=False)

#%%

data_dir=r'/home/greydon/Documents/data/emory_seeg/derivatives'
isub='montage'


for isub in [x for x in os.listdir(data_dir) if os.path.isdir(os.path.join(data_dir,x))]:
	
	files=glob.glob(os.path.join(data_dir,isub,'*.mtg'))
	#files=glob.glob(r'/media/stereotaxy/3E7CE0407CDFF11F/data/iEEG/resources/montage/'+f'{isub}/*.mtg')
	
	for ifile in files:
		out_info=get_montage(ifile)
		out_info.to_csv(os.path.splitext(ifile)[0]+'.tsv',sep='\t',float_format='%.3f',index=False)

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
