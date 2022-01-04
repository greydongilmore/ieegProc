from nilearn import plotting, image
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Qt5Agg')
import nibabel as nib
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
import json
import time
import numpy as np


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
	
	isub="070"
	datap=r'/media/veracrypt6/projects/SEEG/derivatives/atlasreg/'
	
	input=dotdict({
		'img':datap+f'sub-P{isub}/sub-P{isub}_desc-masked_from-atropos3seg_T1w.nii.gz',
		'seg4d':datap+f'sub-P{isub}/sub-P{isub}_desc-atropos3seg_probseg.nii.gz',
		'mapping':datap+f'sub-P{isub}/sub-P{isub}_desc-atropos3seg_mapping.json',
	})
	
	output=dotdict({
		#'html':'/home/greydon/Downloads/' + f'sub-P{isub}_from-contrast_to-noncontrast_regqc.html',
		'html':'/home/greydon/Downloads/' + f'sub-P{isub}_desc-affine_from-subject_to-MNI152NLin2009cSym_regqc.html',
		#'png':'/home/greydon/Downloads/' + f'sub-P{isub}_from-contrast_to-noncontrast_regqc.png'
		'png':'/home/greydon/Downloads/' + f'sub-P{isub}_desc-affine_from-subject_to-MNI152NLin2009cSym_regqc.png'
	})
	
	snakemake = Namespace(output=output, input=input)

	title = 'sub-P001'

#html_view = plotting.view_img(stat_map_img=snakemake.input.seg,bg_img=snakemake.input.img,
#                              opacity=0.5,cmap='viridis',dim=-1,threshold=0.5,
#                              symmetric_cmap=False,title='sub-{subject}'.format(**snakemake.wildcards))
#
#html_view.save_as_html(snakemake.output.html)


ref_img=nib.load(snakemake.input.img)
ref_resamp = nib.nifti1.Nifti1Image(ref_img.get_fdata(), affine=ref_img.affine,header=ref_img.header)
ref_resamp = image.resample_img(ref_img, target_affine=np.eye(3), interpolation='continuous')

flo_img=nib.load(snakemake.input.seg4d)
flo_resamp = nib.nifti1.Nifti1Image(flo_img.get_fdata(), affine=flo_img.affine,header=flo_img.header)
flo_resamp = image.resample_img(flo_img, target_affine=np.eye(3), interpolation='continuous')


with open(snakemake.input.mapping, "r+") as fid:
	mapping_data = json.load(fid)

mapping_data = {str(y):x for x,y in mapping_data.items()}
mapping_data['3']=mapping_data['0']
del mapping_data['0']

coords = plotting.find_xyz_cut_coords(ref_resamp)

colors_dict=[(102,204,238),(34,136,51),(238,102,119),(170,51,119),(204,51,17),(222,143,5),(213,94,0)]
colors_map=[]
for imap in range(len(mapping_data)):
	colors_map.append(colors_dict[imap])

colors_map=[list(np.array(x)/255) for x in colors_map]


display = plotting.plot_prob_atlas(bg_img=ref_resamp,maps_img=flo_resamp, view_type='continuous',display_mode='ortho', draw_cross=False, alpha=.2,cmap=LinearSegmentedColormap.from_list("",colors_map),vmin=0,vmax=4,cut_coords=coords,colorbar=True)

new_yticks=[]
for j, lab in enumerate(mapping_data.values()):
	new_yticks.append(matplotlib.text.Text(0, (j + 1.5), text=lab))

display._cbar.ax.set_yticklabels(new_yticks)
display.savefig(snakemake.output.png,dpi=300)
display.close()

