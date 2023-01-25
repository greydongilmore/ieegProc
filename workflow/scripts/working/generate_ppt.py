#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 15:26:06 2022

@author: greydon
"""

import os
import glob
from pptx import Presentation
from pptx.util import Inches,Pt
from bs4 import BeautifulSoup as Soup
from pptx.dml.color import RGBColor
from pptx.oxml.xmlchemy import OxmlElement
from pptx.enum.text import PP_ALIGN
import pandas as pd
from pptx.enum.text import MSO_AUTO_SIZE
import datetime


def add_slide(presentation, layout, title_dict):
	slide = presentation.slides.add_slide(layout) # adding a slide
	for ititle in list(title_dict):
		title = slide.shapes.add_textbox(*title_dict[ititle]['position'])
		
		tf=title.text_frame
		tf.auto_size = MSO_AUTO_SIZE.NONE
		tf.word_wrap=False
		
		p = tf.paragraphs[0]
		p.alignment= PP_ALIGN.CENTER
		p.text = ititle
		p.font.size = Pt(title_dict[ititle]['font_size'])
		p.font.color.rgb = title_dict[ititle]['color']
	
	return slide

color_map={
	"gray": (102,102,102),
	"grey": (102,102,102),
	"black": (0,0,0),
	"red": (255,0,0),
	"blue": (42,96,153),
	"purple": (128,0,128),
	"orange": (255,128,0),
	"yellow": (255,255,0),
	"brown": (43,34,2),
	"green": (0,169,51),
	"white": (255,255,255),
}

#%%


debug = True
if debug:
	class dotdict(dict):
		"""dot.notation access to dictionary attributes"""
		__getattr__ = dict.get
		__setattr__ = dict.__setitem__
		__delattr__ = dict.__delitem__

	class Namespace:
		def __init__(self, **kwargs):
			self.__dict__.update(kwargs)

	isub = 'sub-P103'
	data_dir = r'/media/data/data/SEEG/derivatives'

	input = dotdict({
			'shopping_list': f'{data_dir}/seega_scenes/{isub}/*shopping_list.xlsx',
	 })
	
	snakemake = Namespace(input=input)


df_elec_raw = pd.read_excel(glob.glob(snakemake.input.shopping_list)[0],header=None)
df_elec=df_elec_raw.iloc[4:,:].reset_index(drop=True)
df_elec.columns=df_elec_raw.iloc[3]
df_elec=df_elec.iloc[0:df_elec.iloc[:,1].isnull().idxmax()]


pt_pin='PIN'
pin_idx=[i for i,x in enumerate(df_elec_raw.iloc[1].values) if x =='PIN']
if pin_idx:
	pt_pin = df_elec_raw.iloc[1,pin_idx[0]+1]


sx_date = 'yyyy-mm-dd'
sx_idx=[i for i,x in enumerate(df_elec_raw.iloc[2].values) if x =='Date']
if sx_idx:
	sx_date = df_elec_raw.iloc[2,sx_idx[0]+1].strftime('%Y-%m-%d')


lastname="lastname"
firstname="firstname"
name_idx=[i for i,x in enumerate(df_elec_raw.iloc[0].values) if x =='Name']
if name_idx:
	if ',' in df_elec_raw.iloc[0,name_idx[0]+1]:
		lastname,firstname=df_elec_raw.iloc[0,name_idx[0]+1].split(',')
	else:
		firstname,lastname=df_elec_raw.iloc[0,name_idx[0]+1].split(' ')
	
	firstname=firstname.strip()
	lastname=lastname.strip()


prs=Presentation()
prs.slide_width = Inches(16)
prs.slide_height = Inches(9)

blank_slide_layout = prs.slide_layouts[6]
fill=blank_slide_layout.background.fill
fill.solid()
fill.fore_color.rgb=RGBColor(0,0,0)

# Title slide
title_dict={
	f"{lastname}, {firstname}":{
		"font_size":52,
		"color": RGBColor(255,255,255),
		"position": (Inches(3), Inches(2), Inches(10), Inches(1))
		},
	f"{pt_pin}":{
		"font_size":36,
		"color": RGBColor(255,255,255),
		"position": (Inches(5), Inches(3), Inches(6), Inches(.8))
		},
	f"Implantation Date:\n{sx_date}":{
		"font_size":36,
		"color": RGBColor(255,255,255),
		"position": (Inches(5), Inches(4.5), Inches(6), Inches(1.5))
		}
	}

title_slide=add_slide(prs, blank_slide_layout, title_dict)
title_slide.name="title slide"

# Shopping list
shopping_list_slide=add_slide(prs, blank_slide_layout, {})
shopping_list_slide.name="shopping list"


# Errors
title_dict={
	"Errors":{
		"font_size":48,
		"color": RGBColor(255,255,255),
		"position": (Inches(3), Inches(.5), Inches(10), Inches(1))
		}
	}


error_slide=add_slide(prs, blank_slide_layout, title_dict)
error_slide.name="errors"


for _, row in df_elec.iterrows():
	if not row['Electrode label'] =='aborted':
		title_dict={
			row['Target']:{
				"font_size": 48,
				"color": RGBColor(255,255,255),
				"position": (Inches(3), Inches(.5), Inches(10), Inches(1))
				}
			}
		
		elec_slide=add_slide(prs, blank_slide_layout, title_dict)
		elec_slide.name=row['Target']
		
		if isinstance(row['Electrode label'],int):
			elec_color = 'black'
			elec_text = f"{row['Electrode label']}".zfill(3)
		else:
			elec_color = ''.join([x for x in row['Electrode label'] if x.isalpha()]).lower()
			elec_text = row['Electrode label']
		
		textbox = elec_slide.shapes.add_textbox(Inches(13.5), Inches(4.5), Inches(2), Inches(.5))
		textbox.fill.solid()
		if elec_color in ("yellow","green","white"):
			textbox.fill.fore_color.rgb = RGBColor(0, 0, 0)
		else:
			textbox.fill.fore_color.rgb = RGBColor(255, 255, 255)
		tf=textbox.text_frame
		tf.auto_size = MSO_AUTO_SIZE.NONE
		tf.word_wrap=False
		
		p = tf.paragraphs[0]
		p.alignment= PP_ALIGN.CENTER
		p.text = elec_text
		p.font.size = Pt(24)
		p.font.bold=True
		p.font.color.rgb = RGBColor(color_map[elec_color][0],color_map[elec_color][1],color_map[elec_color][2])
		
		line = textbox.line
		line.color.rgb = RGBColor(255, 0, 0)
		line.width = Inches(0.04)

out_fname = f"{lastname.replace(' ','')}_{firstname}_{sx_date}_maps.pptx"
prs.save(f'{data_dir}/seega_scenes/{isub}/{out_fname}')




