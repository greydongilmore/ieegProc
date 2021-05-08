#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import subprocess
import numpy as np
import pandas as pd
import argparse
import glob
import shutil
import nibabel
import yaml
from sklearn.model_selection import train_test_split

c3d_path='/opt/c3d/bin/c3d'

debug = False

if debug:
    class Namespace:
        def __init__(self, **kwargs):
            self.__dict__.update(kwargs)
    
    input_dir = '/media/veracrypt6/projects/iEEG/working_dir/out/deriv'
    output_dir = '/media/veracrypt6/projects/iEEG/imaging/clinical/deriv/model_pred'
    
    args = Namespace(input_dir=input_dir, output_dir=output_dir)

script_dir = os.path.dirname(os.path.realpath(sys.argv[0]))
config_path = os.path.join(os.path.dirname(os.path.dirname(script_dir)), 'config/config.yml')

with open(config_path) as file:
    config = yaml.load(file, Loader=yaml.FullLoader)

c3d_path=os.path.join(config['c3d_path'],'c4d')

def run_command(cmdLineArguments):
    process = subprocess.Popen(cmdLineArguments, stdout=subprocess.PIPE, shell=True, stderr=subprocess.PIPE)
    stdout = process.communicate()[0]
    p_status = process.wait()
    
def main():
    
    nifti_files=glob.glob(os.path.join(config['input_dir'])+'/*.nii.gz')
    x_train ,x_test = train_test_split(nifti_files,test_size=0.2)
    
    dir_setup={'test_data','test_labels','train_data','train_labels'}
    for idir in dir_setup:
        if not os.path.exists(os.path.join(config['output_dir'],idir)):
            os.makedirs(os.path.join(config['output_dir'],idir))
    
    for ifile in x_train:
        
        sub = os.path.basename(ifile).split('_')[0]
        
        input_landmarks=glob.glob(os.path.join(config['input_dir'],'input_labels')+f'/{sub}*_contacts.nii.gz')
        
        if not input_landmarks:
        
        else:
            
            if not os.path.exists(os.path.join(args.output_dir,'data','niftis')):
                os.makedirs(os.path.join(args.output_dir,'data','niftis'))
            
            if not os.path.exists(os.path.join(args.output_dir,'data','labels')):
                os.makedirs(os.path.join(args.output_dir,'data','labels'))
            
            output_landmarks=os.path.join(args.output_dir,'data','labels', os.path.basename(input_landmarks[0]))
            
            seg = nibabel.load(input_landmarks[0])
            data = seg.get_fdata()
            data[data==0]=2
            
            ni_img = nibabel.Nifti1Image(data, seg.affine)
            nibabel.save(ni_img, output_landmarks)

            shutil.copyfile(ifile, os.path.join(args.output_dir,'data','niftis',os.path.basename(ifile)))
