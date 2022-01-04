import matplotlib
matplotlib.use('Agg')
import ants
from nilearn import plotting,image
import nibabel as nib
import numpy as np

# snakemake.input.
# snakemake.output.
# html_view.open_in_browser()

template = ants.image_read(ants.get_ants_data('mni'))

ct_img=nib.load(snakemake.input.ct)
if (np.isnan(ct_img.get_fdata())).any():
	ct_img=nib.Nifti1Image(np.nan_to_num(ct_img.get_fdata()), header=ct_img.header, affine=ct_img.affine)
	nib.save(ct_img,snakemake.input.ct)

ct_ants = ants.image_read(snakemake.input.ct)
mask_ants = ants.image_read(snakemake.input.mask)

ct_ants_reg = ants.registration(template, ct_ants, type_of_transform='QuickRigid')
ct_ants_reg_applied=ants.apply_transforms(template, ct_ants, transformlist=ct_ants_reg['fwdtransforms'])
ct_resample = ants.to_nibabel(ct_ants_reg_applied)

mask_ants_reg_applied = ants.apply_transforms(ct_ants, mask_ants, transformlist=ct_ants_reg['fwdtransforms'])
mask_resample = ants.to_nibabel(mask_ants_reg_applied)

mask_params = {
			'symmetric_cmap': True,
			'cut_coords':[0,0,0],
			'dim': 1,
			'cmap':'viridis',
			'opacity':0.7
			}

html_view = plotting.view_img(stat_map_img=mask_resample,bg_img=ct_resample,**mask_params)
html_view.save_as_html(snakemake.output.html)

