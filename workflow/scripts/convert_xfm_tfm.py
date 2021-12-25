import numpy as np

#filen=r'/home/greydon/Documents/data/SEEG/derivatives/atlasreg/sub-P078/sub-P078_acq-noncontrast_desc-affine_from-T1w_to-T1w_type-ras_xfm.txt'
#transformMatrix = np.loadtxt(filen)

transformMatrix = np.loadtxt(snakemake.input.xfm)
lps2ras=np.diag([-1, -1, 1, 1])
ras2lps=np.diag([-1, -1, 1, 1])
transform_lps=np.dot(ras2lps, np.dot(transformMatrix,lps2ras))

Parameters = " ".join([str(x) for x in np.concatenate((transform_lps[0:3,0:3].reshape(9), transform_lps[0:3,3]))])
#output_matrix_txt = filen.split('.txt')[0] + '.tfm'

with open(snakemake.output.tfm, 'w') as fid:
	fid.write("#Insight Transform File V1.0\n")
	fid.write("#Transform 0\n")
	fid.write("Transform: AffineTransform_double_3_3\n")
	fid.write("Parameters: " + Parameters + "\n")
	fid.write("FixedParameters: 0 0 0\n")
