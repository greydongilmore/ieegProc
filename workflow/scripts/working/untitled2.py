#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import nibabel as nib
import subprocess
import glob
import skimage.morphology as morph 
import subprocess
from scipy.spatial import ConvexHull, Delaunay
import matplotlib.pyplot as plt

class IndexTracker:
	def __init__(self, ax, img_data,points=None, rotate_img=False, rotate_points=False, title=None):
		self.ax = ax
		self.scatter=None
		self.points=None
		self.title=title
		no_index= True
		if rotate_img:
			self.img_data = np.fliplr(np.rot90(img_data.copy(),3))
			#if self.points is not None:
			#	self.points = np.rot90(points, k=1)
		else:
			self.img_data = img_data.copy()
		
		rows, cols, self.slices = self.img_data.shape
		self.ind = self.slices//2
		if points is not None:
			self.points=points.copy()
			if rotate_points:
				self.points[:,[0,1]]=self.points[:,[1,0]]
			while no_index==True:
				if any(self.points[:,2]==self.ind):
					no_index=False
					point_plot=np.vstack([np.mean(self.points[(self.points[:,2]==self.ind)*(self.points[:,3]==x),:2],0) for x in np.unique(self.points[self.points[:,2]==self.ind,3])])
					self.scatter,=ax.plot(point_plot[:,1],point_plot[:,0], marker="o", markersize=12, c = 'yellow', fillstyle='none', markeredgewidth=1, linestyle = 'None')
				else:
					self.ind+=1
				
		self.im = ax.imshow(self.img_data[:, :, self.ind], origin='lower')
		self.update()

	def on_scroll(self, event):
		print("%s %s" % (event.button, event.step))
		if event.button == 'up':
			self.ind = (self.ind + 1) % self.slices
		else:
			self.ind = (self.ind - 1) % self.slices
		self.update()

	def update(self):
		
		if self.points is not None:
			if any(self.points[:,2]==self.ind):
				point_plot=np.vstack([np.mean(self.points[(self.points[:,2]==self.ind)*(self.points[:,3]==x),:2],0) for x in np.unique(self.points[self.points[:,2]==self.ind,3])])
				self.scatter.set_xdata(point_plot[:,1])
				self.scatter.set_ydata(point_plot[:,0])
		self.im.set_data(self.img_data[:, :, self.ind])
		
		if self.title is not None:
			plot_title=self.title
		else:
			plot_title='slice %s' % self.ind
		
		self.ax.set_title(plot_title, fontdict={'fontsize': 18, 'fontweight': 'bold'})
		self.ax.tick_params(axis='x', labelsize=14)
		self.ax.tick_params(axis='y', labelsize=14)
		self.im.axes.figure.canvas.draw()


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

	input = dotdict(
		{
		"t1w_file": "/home/greydon/Documents/data/SEEG accuracy/bids/sub-P103/sub-P103_T1w.nii.gz",
		"ct_file": "/home/greydon/Documents/data/SEEG accuracy/bids/sub-P103/sub-P103_ct.nii.gz",
		}
	)
	
	output = dotdict(
		{
		"ct_brain_file": "/home/greydon/Documents/data/SEEG accuracy/derivatives/contact_localization/sub-P103/sub-P103_desc-masked_ct.nii.gz",
		"ct_file_wo_bone": "/home/greydon/Documents/data/SEEG accuracy/derivatives/contact_localization/sub-P103/sub-P103_desc-thresh_ct.nii.gz",
		"brain_file": "/home/greydon/Documents/data/SEEG accuracy/derivatives/contact_localization/sub-P103/sub-P103_desc-masked_T1w.nii.gz",
		}
	)
	
	snakemake = Namespace(input=input,output=output)


ct_img_obj = nib.load(snakemake.input.ct_file)
ct_img_data = ct_img_obj.get_fdata()
ct_img_data[ct_img_data>200] = 200
ct_img_data[ct_img_data<=0] = 0
ct_img_wo_bone = nib.Nifti1Image(ct_img_data, affine=ct_img_obj.affine)
nib.save(ct_img_wo_bone, snakemake.output.ct_file_wo_bone)

t1w_img_obj = nib.load(snakemake.input.t1w_file)
t1w_img_data = t1w_img_obj.get_fdata()


bet_command = '/usr/local/fsl/bin/bet "{}" "{}" -f 0.5 -R -m -n'.format(snakemake.input.t1w_file, snakemake.output.brain_file)

bet_cmd = ' '.join(['/usr/local/fsl/bin/bet2', 
	f'"{snakemake.input.t1w_file}"',
	f'"{snakemake.output.brain_file}"',
	'-f 0.6'
])

bet_process = subprocess.Popen(bet_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
output, errors = bet_process.communicate()


ct_img_masked_obj = nib.load(snakemake.output.ct_brain_file)
ct_img_masked_data = ct_img_masked_obj.get_fdata()


fig, ax = plt.subplots(1, 1)
tracker2 = IndexTracker(ax, t1w_img_data, points=None,rotate_img=True, rotate_points=False)
fig.canvas.mpl_connect('scroll_event', tracker2.on_scroll)
plt.show()