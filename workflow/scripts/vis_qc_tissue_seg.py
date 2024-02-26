from nilearn import plotting
import matplotlib.pyplot as plt
import matplotlib
from scipy import ndimage
from nilearn import plotting, image
import matplotlib.pyplot as plt
import matplotlib
import nibabel as nib
from nibabel.affines import apply_affine
import numpy as np
import base64
import os
from io import BytesIO
import base64
from svgutils.transform import SVGFigure, GroupElement,fromstring
from svgutils.compose import Unit
from tempfile import TemporaryDirectory
from pathlib import Path
from uuid import uuid4
import re
import numpy as np
from matplotlib import gridspec
from matplotlib.colors import LinearSegmentedColormap
np.set_printoptions(precision=6,suppress=True)

from nilearn.datasets import load_mni152_template


def svg2str(display_object, dpi):
	"""Serialize a nilearn display object to string."""
	from io import StringIO

	image_buf = StringIO()
	display_object.frame_axes.figure.savefig(
		image_buf, dpi=dpi, format="svg", facecolor="k", edgecolor="k"
	)
	return image_buf.getvalue()

def extract_svg(display_object, dpi=300):
	"""Remove the preamble of the svg files generated with nilearn."""
	image_svg = svg2str(display_object, dpi)
	
	image_svg = re.sub(' height="[0-9]+[a-z]*"', "", image_svg, count=1)
	image_svg = re.sub(' width="[0-9]+[a-z]*"', "", image_svg, count=1)
	image_svg = re.sub(
		" viewBox", ' preseveAspectRation="xMidYMid meet" viewBox', image_svg, count=1
	)
	start_tag = "<svg "
	start_idx = image_svg.find(start_tag)
	end_tag = "</svg>"
	end_idx = image_svg.rfind(end_tag)
	
	# rfind gives the start index of the substr. We want this substr
	# included in our return value so we add its length to the index.
	end_idx += len(end_tag)
	return image_svg[start_idx:end_idx]

def clean_svg(fg_svgs, bg_svgs, ref=0):
	# Find and replace the figure_1 id.
	svgs = bg_svgs+fg_svgs
	roots = [f.getroot() for f in svgs]
	
	sizes = []
	for f in svgs:
		viewbox = [float(v) for v in f.root.get("viewBox").split(" ")]
		width = int(viewbox[2])
		height = int(viewbox[3])
		sizes.append((width, height))
	nsvgs = len([bg_svgs])
	
	sizes = np.array(sizes)

	# Calculate the scale to fit all widths
	width = sizes[ref, 0]
	scales = width / sizes[:, 0]
	heights = sizes[:, 1] * scales

	# Compose the views panel: total size is the width of
	# any element (used the first here) and the sum of heights
	fig = SVGFigure(Unit(f"{width}px"), Unit(f"{heights[:nsvgs].sum()}px"))

	yoffset = 0
	for i, r in enumerate(roots):
		r.moveto(0, yoffset, scale_x=scales[i])
		if i == (nsvgs - 1):
			yoffset = 0
		else:
			yoffset += heights[i]
	
	# Group background and foreground panels in two groups
	if fg_svgs:
		newroots = [
			GroupElement(roots[:nsvgs], {"class": "background-svg"}),
			GroupElement(roots[nsvgs:], {"class": "foreground-svg"}),
		]
	else:
		newroots = roots
	
	fig.append(newroots)
	fig.root.attrib.pop("width", None)
	fig.root.attrib.pop("height", None)
	fig.root.set("preserveAspectRatio", "xMidYMid meet")

	with TemporaryDirectory() as tmpdirname:
		out_file = Path(tmpdirname) / "tmp.svg"
		fig.save(str(out_file))
		# Post processing
		svg = out_file.read_text().splitlines()

	# Remove <?xml... line
	if svg[0].startswith("<?xml"):
		svg = svg[1:]

	# Add styles for the flicker animation
	if fg_svgs:
		svg.insert(
			2,
			"""\
<style type="text/css">
@keyframes flickerAnimation%s { 0%% {opacity: 1;} 100%% { opacity:0; }}
.foreground-svg { animation: 1s ease-in-out 0s alternate none infinite running flickerAnimation%s;}
.foreground-svg:hover { animation-play-state: running;}
</style>"""
			% tuple([uuid4()] * 2),
		)
	
	return svg

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
		'wm':datap+f'sub-P{isub}/sub-P{isub}_label-WM_desc-atropos3seg_probseg.nii.gz',
		'gm':datap+f'sub-P{isub}/sub-P{isub}_label-GM_desc-atropos3seg_probseg.nii.gz',
		'csf':datap+f'sub-P{isub}/sub-P{isub}_label-CSF_desc-atropos3seg_probseg.nii.gz',
	})
	
	output=dotdict({
		#'html':'/home/greydon/Downloads/' + f'sub-P{isub}_from-contrast_to-noncontrast_regqc.html',
		'html':'/home/greydon/Downloads/' + f'sub-P{isub}_desc-affine_from-subject_to-MNI152NLin2009cSym_regqc.html',
		#'png':'/home/greydon/Downloads/' + f'sub-P{isub}_from-contrast_to-noncontrast_regqc.png'
		'png':'/home/greydon/Downloads/' + f'sub-P{isub}_desc-affine_from-subject_to-MNI152NLin2009cSym_regqc.png'
	})
	
	snakemake = Namespace(output=output, input=input)

	title = 'sub-P001'


#%%

ref_img=nib.load(snakemake.input.img)
ref_resamp = nib.nifti1.Nifti1Image(ref_img.get_fdata(), affine=ref_img.affine,header=ref_img.header)
ref_resamp = image.resample_img(ref_img, target_affine=np.eye(3), interpolation='continuous')

coords = plotting.find_xyz_cut_coords(ref_resamp)

wm_seg=nib.load(snakemake.input.wm)
wm_seg = nib.nifti1.Nifti1Image(wm_seg.get_fdata(), affine=wm_seg.affine,header=wm_seg.header)

gm_seg=nib.load(snakemake.input.gm)
gm_seg = nib.nifti1.Nifti1Image(gm_seg.get_fdata(), affine=gm_seg.affine,header=gm_seg.header)

csf_seg=nib.load(snakemake.input.csf)
csf_seg = nib.nifti1.Nifti1Image(csf_seg.get_fdata(), affine=csf_seg.affine,header=csf_seg.header)


fig, axes = plt.subplots(3, 1,figsize=(16,12))
fig.tight_layout(pad=2)

# make figure of thalamic contours
plotting.plot_roi(roi_img=csf_seg, bg_img=ref_resamp, display_mode='ortho', draw_cross=False, cut_coords=coords,cmap=LinearSegmentedColormap.from_list("",['black','red']),axes=axes[0])
plotting.plot_roi(roi_img=wm_seg, bg_img=ref_resamp, display_mode='ortho', draw_cross=False, cut_coords=coords,cmap=LinearSegmentedColormap.from_list("",['black','yellow']),axes=axes[1])
plotting.plot_roi(roi_img=gm_seg, bg_img=ref_resamp, display_mode='ortho', draw_cross=False, cut_coords=coords,cmap=LinearSegmentedColormap.from_list("",['black','green']),axes=axes[2])

axes[0].set_title('CSF', fontdict={'fontsize': 20, 'fontweight': 'bold'})
axes[1].set_title('White Matter', fontdict={'fontsize': 20, 'fontweight': 'bold'})
axes[2].set_title('Gray Matter', fontdict={'fontsize': 20, 'fontweight': 'bold'})

fig.savefig(snakemake.output.png,dpi=300)
plt.close(fig)


