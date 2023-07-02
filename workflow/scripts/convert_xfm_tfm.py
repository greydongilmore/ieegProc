import numpy as np
import pandas as pd
from scipy.io import loadmat


# filen=r'/home/greydon/Documents/data/SEEG_peds/derivatives/atlasreg/sub-P010/sub-P010_acq-noncontrast_desc-rigid_from-noncontrast_to-contrast_type-ras_xfm.mat'
# filen_out=r'/home/greydon/Documents/data/SEEG_peds/derivatives/atlasreg/sub-P010/sub-P010_acq-noncontrast_desc-rigid_from-noncontrast_to-contrast_type-ras_xfm_1.txt'
# filen=r'/home/greydon/Documents/data/SEEG/derivatives/atlasreg/sub-P097/sub-P097_desc-rigid_from-ct_to-T1w_type-ras_xfm.txt'


transformMatrix = np.loadtxt(snakemake.input.xfm)
lps2ras=np.diag([-1, -1, 1, 1])
ras2lps=np.diag([-1, -1, 1, 1])
transform_lps=np.dot(ras2lps, np.dot(np.linalg.inv(transformMatrix),lps2ras))

Parameters = " ".join([str(x) for x in np.concatenate((transform_lps[0:3,0:3].reshape(9), transform_lps[0:3,3]))])
#output_matrix_txt = filen.split('.txt')[0] + '.tfm'

with open(snakemake.output.tfm, 'w') as fid:
	fid.write("#Insight Transform File V1.0\n")
	fid.write("#Transform 0\n")
	fid.write("Transform: AffineTransform_double_3_3\n")
	fid.write("Parameters: " + Parameters + "\n")
	fid.write("FixedParameters: 0 0 0\n")

# for irow in df1.iterrows():
# 	if irow[1].values[0].lower().startswith('parameters'):
# 		temp=irow[1].values[0].lower().split('parameters:')[-1].split(' ')
# 		temp=[float(x) for x in temp if x !='']
# 		
# 		data = np.eye(4)
# 		data[0:3, 3] = np.array(temp)[9:12]
# 		data[:3, :3] = np.array(temp)[0:9].reshape((3, 3))
# 		
# 		np.savetxt(snakemake.output.txt, np.linalg.inv(data), fmt = '%.3f')
