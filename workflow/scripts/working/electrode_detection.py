# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 06:54:17 2020

@author: Greydon
"""

import nibabel as nib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D
import numpy as np
import pandas as pd
from skimage import morphology
from skimage import measure
from sklearn.cluster import KMeans
from sklearn.cluster import AffinityPropagation
from sklearn.cluster import DBSCAN


def find_markers(img_data, threshold):
	final_location = pd.DataFrame([])
	for slice_idx in range(img_data.shape[2]):
		slice_img = img_data[:,:,slice_idx]
		
		# threshold the image
		thresh_img = np.where(slice_img<threshold,1.0,0.0)
		eroded = morphology.erosion(thresh_img, np.ones([3,3]))
		contours = measure.find_contours(eroded, 0.3)
		
		final_loc = []
		for contour in contours:
			# find xy coordinates in contour blob
			left = tuple([contour[contour[:,0].argmin()][0], contour[contour[:,0].argmin()][1]])
			right = tuple([contour[contour[:,0].argmax()][0], contour[contour[:,0].argmax()][1]])
			
			centx = int(np.sqrt(((right[0] + left[0])**2)/4))
			centy = int(np.sqrt(((right[1] + left[1])**2)/4 ))
			if np.all([slice_img[centx,centy] > threshold, contour.shape[0] < 40]):
				final_loc.append([centx, centy, slice_idx,slice_img[centx,centy]])
		
		if final_loc:
			final_location = pd.concat([final_location, pd.DataFrame(final_loc)], axis = 0, ignore_index=True)
	
	final_location = np.c_[final_location.iloc[:,0].values, final_location.iloc[:,1].values, final_location.iloc[:,2].values]
	
	return final_location


class cluster():
	
	def __init__(self, pcd_list=None, method=None):
		self.pcd_list = pcd_list
		
		method_dic={
			"kmeans":{
				'method': KMeans,
				'options':{
					'init':'random',
					'n_clusters':3,
					'n_jobs':3,
					'n_init':10
					}
				},
			"affinity":{
				'method': AffinityPropagation,
				'options':{
					'preference':-1
					}
				},
			"dbscan":{
				'method': DBSCAN,
				'options':{
					'eps':0.0001
					}
				}
			}
		
		self.method = method_dic[method]
		
	def clusterData(self):
		"""
		Take an input array of 3D points (x,y,z) and cluster based on the desired method.
		pcd_list is a list of points
		method_str is a string describing the method to use to cluster
		options is the options dictionary used by some scikit clustering
		functions
		"""
		
		# Keep two dimensions only to compute cluster (remove elevation)
		pcd_list = np.array(self.pcd_list)[:,:2]
	
		# Build the estimator with the given options
		estimator =self.method['method'](**self.method['options'])
	
		# Fit the estimator
		estimator.fit(self.pcd_list)
	
		# Get the labels and return the labeled points
		labels = estimator.labels_
		self.clusters = np.append(self.pcd_list, np.array([labels]).T, 1)
	
		return self.clusters

	def getCluster(self, i):
		"""
		Return the points belonging to a given cluster.
		clusters: list of labeled points where labels are the cluster ids
		i: id of the cluster to get
		label_field: field where the cluster id is stored in the list of labeled
		points
		"""
		return [c.tolist() for c in self.clusters if c[2] == i]

	def clustersSize(self):
		"""
		Return the size of the clusters given a list of labeled points.
		"""
		from collections import Counter
		
		labels = self.clusters[:,2]
		counter = Counter(labels)
		
		return counter.most_common()

#%%

filen = r'/home/greydon/Documents/GitHub/seeg2bids-pipeline/resources/sub-P221_ses-post_acq-DBS_ct.nii.gz'
filen = r'/media/veracrypt6/projects/iEEG/working_dir/out/results/sub-049/sub-049_desc-masked_from-atropos3seg_ct.nii.gz'

# threshold for CT
threshold=2000

img = nib.load(filen)
img_data = img.get_fdata()

# extract cluster blobs above threshold
final_location = find_markers(img_data, threshold)


# cluster the points (mthod can be: kmeans, affinity, dbscan)
clust = cluster(final_location, method="affinity")
clusters = clust.clusterData()


# plot the clusters
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(*clusters[:,0:3].T, s=2, c=clusters[:,3])








