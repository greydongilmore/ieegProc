from nilearn import plotting, image
from nilearn.datasets import load_mni152_template
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

# ref = r'/media/veracrypt6/projects/iEEG/working_dir/out/deriv/atlasreg/sub-P001/results/sub-P001_T1w.nii.gz'
# flo = r'/media/veracrypt6/projects/iEEG/working_dir/out/deriv/atlasreg/sub-P001/results/sub-P001_desc-masked_from-atropos3seg_ct.nii.gz'
# title = 'sub-P001'
# html=r'/media/veracrypt6/projects/iEEG/working_dir/out/deriv/atlasreg/sub-P001/qc/sub-P001_desc-masked_from-ct_to-T1w_regqc.html'
# png='/media/veracrypt6/projects/iEEG/working_dir/out/deriv/atlasreg/sub-P001/qc/sub-P001_desc-masked_from-ct_to-T1w_regqc.png'
# html_view.open_in_browser()

template = load_mni152_template()

ref_resamp = image.resample_img(snakemake.input.ref, target_affine=template.affine, interpolation='nearest')
flo_resamp = image.resample_to_img(snakemake.input.flo, ref_resamp)

title='sub-{subject}'.format(**snakemake.wildcards)

html_view = plotting.view_img(stat_map_img=ref_resamp,bg_img=flo_resamp,
                              opacity=0.4,cmap='viridis',dim=0,
                              symmetric_cmap=False,title=title,
                              cut_coords=[0,0,0])

html_view.save_as_html(snakemake.output.html)


display = plotting.plot_anat(flo_resamp,display_mode='ortho',cut_coords=[0,0,0])
display.add_contours(ref_resamp, colors='r')
display.savefig(snakemake.output.png,dpi=250)
display.close()

