import numpy as np
import pandas as pd
from scipy.io import loadmat


# filen=r'/home/greydon/Documents/data/SEEG_peds/derivatives/atlasreg/sub-P010/sub-P010_acq-noncontrast_desc-rigid_from-noncontrast_to-contrast_type-ras_xfm.mat'
# filen_out=r'/home/greydon/Documents/data/SEEG_peds/derivatives/atlasreg/sub-P010/sub-P010_acq-noncontrast_desc-rigid_from-noncontrast_to-contrast_type-ras_xfm_1.txt'
# filen=r'/home/greydon/Documents/data/SEEG/derivatives/atlasreg/sub-P097/sub-P097_desc-rigid_from-ct_to-T1w_type-ras_xfm.txt'


sub2template= np.loadtxt(snakemake.input.xfm)
lps2ras=np.diag([-1, -1, 1, 1])
ras2lps=np.diag([-1, -1, 1, 1])
sub2template=np.dot(ras2lps, np.dot(sub2template,lps2ras))

Parameters = " ".join([str(x) for x in np.concatenate((sub2template[0:3,0:3].reshape(9), sub2template[0:3,3]))])

with open(snakemake.output.xfm_inv, 'w') as fid:
	fid.write("#Insight Transform File V1.0\n")
	fid.write("#Transform 0\n")
	fid.write("Transform: AffineTransform_double_3_3\n")
	fid.write("Parameters: " + Parameters + "\n")
	fid.write("FixedParameters: 0 0 0\n")

