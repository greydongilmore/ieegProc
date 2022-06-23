#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 21:18:35 2022

@author: greydon
"""

from nilearn import plotting,surface
import nibabel as nb
import numpy as np
import matplotlib.pyplot as plt

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
	
	isub='sub-P092'
	data_dir=r'/home/greydon/Documents/data/SEEG/derivatives/atlasreg'
	
	input=dotdict({'t1':f'/home/greydon/Documents/data/SEEG/derivatives/fastsurfer/sub-P092/fastsurfer/mri/orig.mgz',
				'pet':f'{data_dir}/{isub}/{isub}_space-T1w_desc-rigid_pet.nii.gz',
				'gm':f'{data_dir}/{isub}/{isub}_label-GM_desc-atropos3seg_probseg.nii.gz',
				'wm':f'{data_dir}/{isub}/{isub}_label-GM_desc-atropos3seg_probseg.nii.gz',
				'mask':f'{data_dir}/{isub}/{isub}_desc-brain_from-MNI152NLin2009cSym_reg-affine_mask.nii.gz',
				'seg':f'{data_dir}/{isub}/{isub}_desc-atroposKseg_dseg.nii.gz',
				})
	
	output=dotdict({
		't1_grad':f'{data_dir}/{isub}/{isub}_desc-magnitude_T1w.nii.gz',
		't1_intensity':f'{data_dir}/{isub}/{isub}_desc-intensity_T1w.nii.gz',
		'png_mask':f'{data_dir}/{isub}/qc/{isub}_desc-brain_maskqc.png',
		'png_seg':f'{data_dir}/{isub}/qc/{isub}_desc-segmentation_segqc.png'
	})
	
	snakemake = Namespace(output=output, input=input)


lh_pial_file=r'/home/greydon/Documents/data/SEEG/derivatives/fastsurfer/sub-P092/fastsurfer/surf/lh.pial'
rh_pial_file=r'/home/greydon/Documents/data/SEEG/derivatives/fastsurfer/sub-P092/fastsurfer/surf/rh.pial'
lh_sulc_file=r'/home/greydon/Documents/data/SEEG/derivatives/fastsurfer/sub-P092/fastsurfer/surf/lh.sulc'
rh_sulcl_file=r'/home/greydon/Documents/data/SEEG/derivatives/fastsurfer/sub-P092/fastsurfer/surf/rh.sulc'
lh_infl_file=r'/home/greydon/Documents/data/SEEG/derivatives/fastsurfer/sub-P092/fastsurfer/surf/lh.inflated'
rh_infl_file=r'/home/greydon/Documents/data/SEEG/derivatives/fastsurfer/sub-P092/fastsurfer/surf/rh.inflated'


lh_data = nb.freesurfer.read_geometry(lh_pial_file)
rh_data = nb.freesurfer.read_geometry(rh_pial_file)
lh_sulc_data=surface.load_surf_data(lh_sulc_file)
rh_sulc_data=surface.load_surf_data(rh_sulcl_file)
lh_infl_data = nb.freesurfer.read_geometry(lh_infl_file)
rh_infl_data = nb.freesurfer.read_geometry(rh_infl_file)


all_ver=np.concatenate([lh_data[0],rh_data[0]],axis=0)
tmp_facer=rh_data[1]+lh_data[0].shape[0]
all_face=np.concatenate([lh_data[1],tmp_facer],axis=0)
surf_mesh=[all_ver, all_face]
bg_map = np.concatenate((lh_sulc_data, rh_sulc_data))


fig=plotting.plot_surf(surf_mesh, surf_map=bg_map,bg_map=None,
					   view='lateral', cmap = plt.cm.Greys_r, engine='plotly', 
					   ouput_file="/home/greydon/Downloads/text.html",alpha=1,
					   bg_on_data=True,darkness=0.8,title='Cortical Surface',symmetric_cmap=True)

fig.show()
fig.figure.write_html("/home/greydon/Downloads/text.html")



surf_mesh_dict={}
surf_mesh_dict['pial_left']=lh_data
surf_mesh_dict['pial_right']=rh_data
surf_mesh_dict['sulc_left']=lh_sulc_data
surf_mesh_dict['sulc_right']=rh_sulc_data
surf_mesh_dict['infl_left']=lh_infl_data
surf_mesh_dict['infl_right']=rh_infl_data

fig=plotting.plot_img_on_surf(nb.load(snakemake.input.pet),surf_mesh_dict)








