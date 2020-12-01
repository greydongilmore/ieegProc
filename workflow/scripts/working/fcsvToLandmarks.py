#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import subprocess
import numpy as np
import pandas as pd
import argparse
import glob
import shutil
import nibabel

c3d_path='/opt/c3d/bin/c3d'

debug = False

if debug:
    class Namespace:
        def __init__(self, **kwargs):
            self.__dict__.update(kwargs)
    
    input_dir = '/media/veracrypt6/projects/iEEG/working_dir/out/deriv'
    output_dir = '/media/veracrypt6/projects/iEEG/imaging/clinical/deriv/model_pred'
    
    args = Namespace(input_dir=input_dir, output_dir=output_dir)
    
def run_command(cmdLineArguments):
    process = subprocess.Popen(cmdLineArguments, stdout=subprocess.PIPE, shell=True, stderr=subprocess.PIPE)
    stdout = process.communicate()[0]
    p_status = process.wait()
    
def main(args):
    
    for ifile in glob.glob(args.input_dir+'/*/*/*_space-T1w_desc-affine_ct.nii.gz'):
        
        sub = os.path.basename(ifile).split('_')[0]
        input_landmarks=glob.glob(args.input_dir+f'/*/*/{sub}*_contacts.nii.gz')
        
        if input_landmarks:
            
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
    
###############################################################################
if __name__ == "__main__":
    
    ###############################################################################
    # command-line arguments with flags
    ###############################################################################
    parser = argparse.ArgumentParser(description="Run c3d to convert 3D Slicer fcsv files into binary landmark masks.")
    
    parser.add_argument("-i", dest="input_dir", help="Input data directory with nifti files and fcsv files.")
    parser.add_argument("-o", dest="output_dir", help="Output directory to store label volumes)")
    
    args = parser.parse_args()
    
    main(args)