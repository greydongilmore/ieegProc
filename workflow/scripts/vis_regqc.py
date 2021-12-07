from nilearn import plotting, image
from nilearn.datasets import load_mni152_template
import matplotlib.pyplot as plt
import matplotlib
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
	
	isub="078"
	input=dotdict({
		'flo':r'/media/veracrypt6/projects/iEEG/imaging/clinical/derivatives/atlasreg/'+f'sub-P{isub}/sub-P{isub}_desc-atropos3seg_probseg.nii.gz',
		'ref':r'/media/veracrypt6/projects/iEEG/imaging/clinical/derivatives/atlasreg/'+f'sub-P{isub}/sub-P{isub}_desc-n4_T1w.nii.gz'
	})
	
	output=dotdict({
		'html':'/home/greydon/Downloads/' + f'sub-P{isub}_desc-masked_from-ct_to-T1w_regqc.html',
		'png':'/home/greydon/Downloads/' + f'sub-P{isub}_desc-masked_from-ct_to-T1w_regqc.png'
	})
	
	snakemake = Namespace(output=output, input=input)

	title = 'sub-P001'

# html_view.open_in_browser()

template = load_mni152_template()

ref_resamp = image.resample_img(snakemake.input.ref, target_affine=template.affine, interpolation='nearest')
flo_resamp = image.resample_to_img(snakemake.input.flo, ref_resamp)

title='sub-{subject}'.format(**snakemake.wildcards)

html_view = plotting.view_img(stat_map_img=ref_resamp, bg_img=flo_resamp,
                              opacity=0.4,cmap='viridis',dim=0,
                              symmetric_cmap=False,title=title,
                              cut_coords=[0,0,0])

html_view.save_as_html(snakemake.output.html)


display = plotting.plot_anat(flo_resamp,display_mode='ortho',cut_coords=[0,0,0])
display.add_contours(ref_resamp, colors='r')
display.savefig(snakemake.output.png,dpi=250)
display.close()

#html_view.save_as_html(snakemake.output.html)