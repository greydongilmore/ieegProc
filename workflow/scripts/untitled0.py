#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 22:29:53 2022

@author: greydon
"""
import SimpleITK as sitk

grad_t1 = sitk.ReadImage(r'/home/greydon/Documents/data/SEEG/derivatives/atlasreg/sub-P096/sub-P096_desc-magnitude_T1w.nii.gz')
grad_t1 = sitk.ScalarToRGBColormap(grad_t1, sitk.ScalarToRGBColormapImageFilter.Hot)

sitk.WriteImage(grad_t1,r'/home/greydon/Documents/data/SEEG/derivatives/atlasreg/sub-P096/sub-P096_desc-magnitude_T1w_test.nii.gz')
	
fig, ax = plt.subplots(1, 1)
tracker = IndexTracker(ax, sitk.GetArrayFromImage(grad_t1)[25], points=None,rotate_img=True, rotate_points=False)
fig.canvas.mpl_connect('scroll_event', tracker.on_scroll)
plt.show()