import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')

#read fcsv electrodes file
df = pd.read_table(snakemake.input.fcsv,sep=',',header=2)

#load transform from subj to template
sub2template= np.loadtxt(snakemake.input.xfm_ras)

#plot electrodes transformed (affine) to MNI space, with MNI glass brain
from nilearn import plotting

coords = df[['x','y','z']].to_numpy()
#print(coords.shape)

#to plot in mni space, need to transform coords
tcoords = np.zeros(coords.shape)
for i in range(len(coords)):

    vec = np.hstack([coords[i,:],1])
    tvec = np.linalg.inv(sub2template) @ vec.T
    tcoords[i,:] = tvec[:3]


html_view = plotting.view_markers(tcoords)
html_view.save_as_html(snakemake.output.html)

#plot subject native space electrodes with glass brain
adjacency_matrix = np.zeros([len(coords),len(coords)])


display = plotting.plot_connectome(adjacency_matrix, tcoords)
display.savefig(snakemake.output.png)
display.close()

