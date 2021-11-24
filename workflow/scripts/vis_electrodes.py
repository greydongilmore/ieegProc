
import pandas as pd
import numpy as np
import matplotlib
import re
import matplotlib.pyplot as plt


matplotlib.use('Agg')

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
	
	input=dotdict({'fcsv':'/media/veracrypt6/projects/iEEG/imaging/clinical/deriv/seega_coordinates/sub-P061/sub-P061_space-native_SEEGA.fcsv',
				'xfm_ras':'/media/veracrypt6/projects/iEEG/imaging/clinical/deriv/atlasreg/sub-P061/sub-P061_desc-affine_from-subject_to-MNI152NLin2009cSym_type-ras_xfm.txt'
				})
	output=dotdict({'html':'/home/greydon/Downloads/sub-P061_space-MNI152NLin2009cSym_desc-affine_electrodes.html',
				'png':'/home/greydon/Downloads/sub-P061_space-MNI152NLin2009cSym_desc-affine_electrodevis.png'
				})
	
	snakemake = Namespace(output=output, input=input)

def determine_groups(iterable):
	values = []
	for item in iterable:
		if re.findall(r"([a-zA-Z]+)([0-9]+)([a-zA-Z]+)", item):
			temp = "".join(list(re.findall(r"([a-zA-Z]+)([0-9]+)([a-zA-Z]+)", item)[0]))
		else:
			temp="".join(x for x in item if not x.isdigit())
			
		values.append(temp)
	
	vals,indexes,count = np.unique(values, return_index=True, return_counts=True)
	values = [values[index] for index in sorted(indexes)]
	
	return values,count

#read fcsv electrodes file
df = pd.read_table(snakemake.input.fcsv,sep=',',header=2)

groups,n_members=determine_groups(df['label'].tolist())
cmap = plt.get_cmap('rainbow')
colors=np.repeat(cmap(np.linspace(0, 1, len(groups))), n_members, axis=0)

labels=[str(x) for x in  range(colors.shape[0])]

#load transform from subj to template
sub2template= np.loadtxt(snakemake.input.xfm_ras)

#plot electrodes transformed (affine) to MNI space, with MNI glass brain
from nilearn import plotting

print(plotting.__file__)

coords = df[['x','y','z']].to_numpy()
#print(coords.shape)

#to plot in mni space, need to transform coords
tcoords = np.zeros(coords.shape)
for i in range(len(coords)):

    vec = np.hstack([coords[i,:],1])
    tvec = np.linalg.inv(sub2template) @ vec.T
    tcoords[i,:] = tvec[:3]

html_view = plotting.view_markers(tcoords, marker_size=6.0, marker_color=colors, marker_labels=df['label'].tolist())
#html_view.open_in_browser()
html_view.save_as_html(snakemake.output.html)

#plot subject native space electrodes with glass brain
adjacency_matrix = np.zeros([len(coords),len(coords)])

node_label=np.repeat(groups, n_members, axis=0)

group = np.array([1,3,2,1,3])
cdict = {1: 'red', 2: 'blue', 3: 'green'}

_, idx = np.unique(colors, return_index=True, axis=0)

label_dict=dict(zip(groups,colors[np.sort(idx)].tolist()))

display = plotting.plot_connectome(adjacency_matrix, tcoords, node_color=colors, node_size=3)
display.savefig(snakemake.output.png,dpi=300)
display.close()
