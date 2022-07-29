#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 08:58:27 2022

@author: greydon
"""

import os
from nipype.interfaces import freesurfer
from nipype.interfaces.freesurfer import FSCommand,petsurfer,Info
import subprocess
import os
import numpy as np

def readRegMatrix(trsfPath):
	with open(trsfPath) as (f):
		return np.loadtxt(f.readlines())

isub = 'sub-P093'
data_dir = '/home/greydon/Documents/data/SEEG'

SUBJECTS_DIR=f'{data_dir}/derivatives/fastsurfer'
os.environ['FREESURFER_HOME']='/usr/local/freesurfer/7.2.0'
os.environ['SUBJECTS_DIR']=SUBJECTS_DIR


out_dir = f'{data_dir}/derivatives/petsurfer/{isub}'
if not os.path.exists(out_dir):
	os.makedirs(out_dir)

pet_file = f'{data_dir}/derivatives/atlasreg/{isub}/{isub}_space-T1w_desc-rigid_pet.nii.gz'



fs_cmd=f"cd {SUBJECTS_DIR}/{isub}/mri&&/usr/local/freesurfer/7.2.0/mri_ca_register -align-after -nobigventricles -mask brainmask.mgz -T transforms/talairach.lta norm.mgz /usr/local/freesurfer/7.2.0/average/RB_all_2016-05-10.vc700.gca transforms/talairach.m3z"
subprocess.run(fs_cmd, shell=True)

fs_cmd=f'SUBJECTS_DIR={SUBJECTS_DIR} &&gtmseg --s {isub} --no-seg-stats --xcerseg'
subprocess.run(fs_cmd, shell=True)


#brainmask.mgz to nii.gz
fs_cmd=f'/usr/local/freesurfer/7.2.0/bin/mri_convert --in_type mgz --out_type nii {SUBJECTS_DIR}/{isub}/mri/brainmask.mgz {SUBJECTS_DIR}/{isub}/mri/brainmask.nii.gz'
subprocess.run(fs_cmd, shell=True)

fs_cmd=f'/usr/local/freesurfer/7.2.0/bin/mri_convert --in_type mgz --out_type nii {SUBJECTS_DIR}/{isub}/mri/aparc+aseg.mgz {SUBJECTS_DIR}/{isub}/mri/aparc+aseg.nii.gz'
subprocess.run(fs_cmd, shell=True)

# register PET to T1w
fs_cmd=f'SUBJECTS_DIR={SUBJECTS_DIR} &&mri_coreg --s {isub} '+\
	f'--mov {pet_file} '+\
	f'--reg {out_dir}/pet_coreg.reg.lta'
subprocess.run(fs_cmd, shell=True)

#convert transform matrix from lta to xfm to tfm
fs_cmd=f'/usr/local/freesurfer/7.2.0/bin/lta_convert --invert --inlta {out_dir}/pet_coreg.reg.lta --outmni {out_dir}/pet_coreg.xfm --src {SUBJECTS_DIR}/{isub}/mri/brainmask.mgz --trg {pet_file}'
subprocess.run(fs_cmd, shell=True)

with open(f'{out_dir}/pet_coreg.xfm', 'r') as fid:
	L = fid.readlines()

transformMatrix  = np.r_[np.stack([[np.float(s) for s in ll.replace(';\n','').replace('\n','').split() if s] for ll in L[-3:]]),np.zeros((1,4))]
Parameters = " ".join([str(x) for x in np.concatenate((transformMatrix[0:3,0:3].reshape(9), transformMatrix[0:3,3]*np.array([-1,-1,1])))])
output_matrix_txt = f'{out_dir}/pet_coreg.tfm'

with open(output_matrix_txt, 'w') as fid:
	fid.write("#Insight Transform File V1.0\n")
	fid.write("#Transform 0\n")
	fid.write("Transform: AffineTransform_double_3_3\n")
	fid.write("Parameters: " + Parameters + "\n")
	fid.write("FixedParameters: 0 0 0\n")

#resample PET based on T1w brainmask
fs_cmd=f'SUBJECTS_DIR={SUBJECTS_DIR} &&mri_vol2vol --mov {pet_file} '+\
	f'--targ $SUBJECTS_DIR/{isub}/mri/brainmask.mgz '+\
	f'--reg {out_dir}/pet_coreg.reg.lta '+\
	f'--o {out_dir}/pet_coreg.nii.gz'
subprocess.run(fs_cmd, shell=True)


out_gtmseg=f'{SUBJECTS_DIR}/{isub}/mri/gtmseg.mgz'
out_gtmpvc=f'{out_dir}/gtmpvc_pet.output'

fs_cmd='mri_gtmpvc --auto-mask 1.000000 0.100000 --default-seg-merge '+\
	f'--i {pet_file} --km-hb 11 12 50 51 --km-ref 8 47 --no-rescale --psf 4.000000 --o pvc '+\
	f'--reg {out_dir}/pet_coreg.reg.lta --save-input --seg {out_gtmseg} --o {out_gtmpvc}'


#%%


fsl_cmd=f'fslmaths {out_dir}/pet_coreg.nii.gz -mas {SUBJECTS_DIR}/{isub}/mri/brainmask.nii.gz {out_dir}/pet_brain.nii.gz'
subprocess.run(fsl_cmd, shell=True)

#### perform normalization
# cerebellar white matter mask
fs_cmd=f'/usr/local/freesurfer/7.2.0/bin/mri_binarize --i {SUBJECTS_DIR}/{isub}/mri/aparc+aseg.mgz --o {out_dir}/cereb_gm_mask.nii --match 8 47'
subprocess.run(fs_cmd, shell=True)


# normalize PET
fsl_cmd=f'fslmeants -i {out_dir}/pet_brain.nii.gz -m {out_dir}/cereb_gm_mask.nii -o {out_dir}/cereb_gm_intensity.txt'
subprocess.run(fsl_cmd, shell=True)


fsl_cmd=f'cereb_intensity=$(cat {out_dir}/cereb_gm_intensity.txt)'+'&&cereb_intensity=${cereb_intensity%.*}&&'+\
	f"fslmaths {out_dir}/pet_brain.nii.gz -div $cereb_intensity {out_dir}/pet_brain_norm.nii.gz"
subprocess.run(fsl_cmd, shell=True)


#### perform PVC
# get binarized gray and white matter images from freesurfer
fs_cmd1=f'/usr/local/freesurfer/7.2.0/bin/mri_binarize --i {SUBJECTS_DIR}/{isub}/mri/aparc+aseg.mgz --o {out_dir}/wm_mask.nii --all-wm'
fs_cmd2=f'/usr/local/freesurfer/7.2.0/bin/mri_binarize --i {SUBJECTS_DIR}/{isub}/mri/aparc+aseg.mgz --o {out_dir}/gm_mask.nii --gm'
fsl_cmd=f'fslmerge -t {out_dir}/gm_wm_mask.nii.gz {out_dir}/gm_mask.nii {out_dir}/wm_mask.nii'

combined='&&'.join([fs_cmd1,fs_cmd2,fsl_cmd])
subprocess.run(combined, shell=True)


fsl_cmd=f'fslcpgeom {out_dir}/pet_brain_norm.nii.gz {out_dir}/gm_wm_mask.nii.gz -d'
subprocess.run(fsl_cmd, shell=True)

petpvc_cmd=f'/home/greydon/Documents/software/PETPVC-1.2.9/bin/petpvc -i {out_dir}/pet_brain_norm.nii.gz -o {out_dir}/pet_brain_norm_PVC.nii.gz -m {out_dir}/gm_wm_mask.nii.gz -x 6.0 -y 6.0 -z 6.0 --pvc MG'
subprocess.run(petpvc_cmd, shell=True)


for Hemisphere in ('lh', 'rh'):
	for isurf in ["pial", "white", "thickness"]:
		in_file = f"{SUBJECTS_DIR}/{isub}/surf/{Hemisphere}.{isurf}"
		Structure="CORTEX_LEFT" if Hemisphere =='lh' else "CORTEX_RIGHT"
		if 'thickness' in isurf:
			out_file = f'{SUBJECTS_DIR}/{isub}/surf/{Hemisphere}.{isurf}.shape.gii'
			if not os.path.exists(out_file):
				white_file = f'{SUBJECTS_DIR}/{isub}/surf/{Hemisphere}.white'
				wb_cmd = f"/usr/local/freesurfer/7.2.0/bin/mris_convert -c {in_file} {white_file} {out_file}"
				subprocess.run(wb_cmd, shell=True)
				
				wb_cmd = f"wb_command -set-structure {out_file} {Structure}"
				subprocess.run(wb_cmd, shell=True)
				
				wb_cmd = f"wb_command -metric-palette {out_file} MODE_AUTO_SCALE_PERCENTAGE -pos-percent 2 98 -palette-name Gray_Interp -disp-pos true -disp-neg true -disp-zero true"
				subprocess.run(wb_cmd, shell=True)
		else:
			out_file = f'{SUBJECTS_DIR}/{isub}/surf/{Hemisphere}.{isurf}.surf.gii'
			if not os.path.exists(out_file):
				wb_cmd = f"/usr/local/freesurfer/7.2.0/bin/mris_convert {in_file} {out_file}"
				subprocess.run(wb_cmd, shell=True)
				
				wb_cmd = f"wb_command -set-structure {out_file} {Structure}"
				subprocess.run(wb_cmd, shell=True)
	
	midthickness_file = f"{out_dir}/{Hemisphere}.midthickness.surf.gii"
	white_file = f'{SUBJECTS_DIR}/{isub}/surf/{Hemisphere}.white.surf.gii'
	pial_file = f'{SUBJECTS_DIR}/{isub}/surf/{Hemisphere}.pial.surf.gii'
	wb_cmd = f"wb_command -surface-average {midthickness_file} -surf {white_file} -surf {pial_file}"
	subprocess.run(wb_cmd, shell=True)
	
	# initialize surfaces
	PET_volume=f'{out_dir}/pet_brain_norm_PVC.nii.gz'
	PET_volume_bin=f'{out_dir}/pet_brain_norm_PVC_bin.nii.gz'
	
	pial_file = f'{SUBJECTS_DIR}/{isub}/surf/{Hemisphere}.pial.surf.gii'
	white_file = f'{SUBJECTS_DIR}/{isub}/surf/{Hemisphere}.white.surf.gii'
	thickness_file = f'{out_dir}/{Hemisphere}.midthickness.surf.gii'
	
	PET_surf_out=f'{out_dir}/{Hemisphere}.shape.gii'
	
	# create volume-roi by binarizing PET_volume
	# otherwise 0 voxels which lie between the surfaces will used for the average too
	fs_cmd=f'/usr/local/freesurfer/7.2.0/bin/mri_binarize --i {PET_volume} --o {PET_volume_bin} --min 0.00001'
	subprocess.run(fs_cmd, shell=True)

	# volume to surface mapping
	wb_cmd=f'wb_command -volume-to-surface-mapping {PET_volume} {thickness_file} {PET_surf_out} '+\
			f'-ribbon-constrained {white_file} {pial_file} -volume-roi {PET_volume_bin}'
	subprocess.run(wb_cmd, shell=True)
	

fs_cmd=f'/usr/local/freesurfer/7.2.0/bin/mri_binarize --i {SUBJECTS_DIR}/{isub}/mri/aparc+aseg.nii.gz --o {out_dir}/subcort_mask.nii.gz --subcort-gm'
subprocess.run(fs_cmd, shell=True)

fsl_cmd=f'fslmaths {SUBJECTS_DIR}/{isub}/mri/aparc+aseg.nii.gz -mas {out_dir}/subcort_mask.nii.gz {out_dir}/subcort_gm.nii.gz'
subprocess.run(fsl_cmd, shell=True)


fsl_cmd=f'/usr/local/freesurfer/7.2.0/bin/mri_binarize --i {out_dir}/subcort_gm.nii.gz'+ \
	f'--o {out_dir}/subcort_gm.nii.gz'+ \
	'--replace 8 361 --replace 10 362 --replace 11 363 --replace 12 364 '+ \
	'--replace 13 365 --replace 17 366 --replace 18 367 --replace 26 368 '+ \
	'--replace 28 369 --replace 47 370 --replace 49 371 --replace 50 372 '+ \
	'--replace 51 373 --replace 52 374 --replace 53 375 --replace 54 376 '+ \
	'--replace 58 377 --replace 60 378 --replace 16 379'
subprocess.run(fsl_cmd, shell=True)

fsl_cmd=f'fslmeants -i {out_dir}/pet_brain_norm_PVC.nii.gz --label={out_dir}/subcort_gm.nii.gz -o {out_dir}/pet_load.subcortical.txt'
subprocess.run(fsl_cmd, shell=True)


# there is a issue with itk in the toolbox, giving error:
#           "Description: itk::ERROR: MullerGartnerImageFilter(0x7fe7ebc89ea0): Inputs do not occupy the same physical space!
#           InputImage Origin: [-9.0000000e+01, 1.2600000e+02, -7.2000000e+01],
#           InputImage_1 Origin: [-9.0000000e+01, 1.2599999e+02, -7.2000000e+01]
#	        Tolerance: 6.9999999e-07"
# "Inputs do not occupy the same physical space!" --> but this in reality is just rounding error of the floats, see th numbers above.
# quick and dirty solution: https://github.com/UCL/PETPVC/issues/31 use fslcpgeom <source> <destination> and copy the header from one nii file on the other

fslcpgeom $PET_results/PET2T1w_brain_norm.nii.gz $PET_results/gm_wm_mask.nii.gz -d #-d ... don't copy image dimensions
petpvc -i $PET_results/PET2T1w_brain_norm.nii.gz -o $PET_results/PET2T1w_brain_norm_PVC.nii.gz -m $PET_results/gm_wm_mask.nii.gz -x 6.0 -y 6.0 -z 6.0 --pvc MG


