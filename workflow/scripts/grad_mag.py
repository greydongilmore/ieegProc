import ants
import numpy as np
from collections import Counter
from nilearn import plotting
from svgutils.transform import SVGFigure, GroupElement,fromstring
from svgutils.compose import Unit
from tempfile import TemporaryDirectory
from pathlib import Path
from uuid import uuid4
import re
import matplotlib.pyplot as plt
from io import BytesIO
import base64
import os
import nibabel as nb
from nilearn.maskers import NiftiMasker
import nilearn


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
	
	return svg,sizes

def threshold_percentile(x, lower_q, upper_q):
	x = x.numpy()
	lq = np.percentile(x, lower_q)
	uq = np.percentile(x, upper_q)
	x = x[np.logical_and(x>lq, x<=uq)]
	return x.flatten()

def peakfinder(gm, wm, lower_q, upper_q):
	gm_peak = Counter(threshold_percentile(gm, lower_q, upper_q)).most_common(1)[0][0]
	wm_peak = Counter(threshold_percentile(wm, lower_q, upper_q)).most_common(1)[0][0]
	bg = 0.5 * (gm_peak + wm_peak)
	return bg

def compute_RI(image, bg, mask):
	ri = np.zeros_like(image)
	bgm = np.stack(np.where(np.logical_and(image < bg, mask == 1)), axis=1)
	bgm_ind = bgm[:,0], bgm[:,1], bgm[:,2]
	bgp = np.stack(np.where(np.logical_and(image > bg, mask == 1)), axis=1)
	bgp_ind = bgp[:,0], bgp[:,1], bgp[:,2]

	ri[bgm_ind] = 100 * (1 - (bg - image[bgm_ind]) / bg )
	ri[bgp_ind] = 100 * (1 + (bg - image[bgp_ind]) / bg )
	return ri

def qc_mosaic_plot(anat_img, overlay_img, out_filen, dim=0, threshold=0):
	
	plt.ioff()
	
	display = plotting.plot_anat(anat_img, 
								 display_mode='mosaic',
								 black_bg=True,
								 cmap=plt.cm.Greys_r,
								 cut_coords=None, 
								 threshold=0, 
								 dim=dim,
								 draw_cross=False)
	fg_svgs = [fromstring(extract_svg(display,400))]
	
	mask_dat = overlay_img.get_fdata()[:,:,:]
	mask_bin = mask_dat > 0
	regions_masked_dat = anat_img.get_fdata() * mask_bin
	regions_masked_img = nilearn.image.new_img_like(anat_img, regions_masked_dat, affine=anat_img.affine, copy_header=True)
	display.add_overlay(regions_masked_img,alpha=0.7,cmap=plt.cm.coolwarm_r)
	
	bg_svgs = [fromstring(extract_svg(display,400))]
	display.savefig(out_filen,dpi=400)
	display.close()
	
	display = plotting.plot_anat(anat_img, 
								 display_mode='mosaic',
								 black_bg=True,
								 cmap=plt.cm.Greys_r,
								 cut_coords=None, 
								 threshold=0, 
								 dim=dim,
								 draw_cross=False)
	display.add_overlay(overlay_img,threshold=threshold,alpha=0.7,cmap=plt.cm.coolwarm_r)
	
	cleaned_svg,sizes=clean_svg(fg_svgs,bg_svgs)
	width = sizes[0, 0]
	scales = width / sizes[:, 0]
	heights = sizes[:, 1] * scales
	
	fig = SVGFigure(Unit(f"{width}px"), Unit(f"{heights[:1].sum()}px"))
	fig.append(GroupElement([fromstring(extract_svg(display,400)).getroot()], {"class": "background-svg"}))
	fig.root.attrib.pop("width", None)
	fig.root.attrib.pop("height", None)
	fig.root.set("preserveAspectRatio", "xMidYMid meet")
	display.close()
	
	with TemporaryDirectory() as tmpdirname:
		out_file = Path(tmpdirname) / "tmp.svg"
		fig.save(str(out_file))
		svg_final = out_file.read_text().splitlines()
	
	htmlbase='<!DOCTYPE html> <html lang="en"> <head> <title>Slice viewer</title>  <meta charset="UTF-8" /> </head> <body>'
	htmlend='</body> </html>'
	htmlfull=htmlbase  + "\n".join(cleaned_svg + svg_final) + htmlend
	
	# Write HTML String to file.html
	with open(out_filen.replace('png','html'), "w") as file:
		file.write(htmlfull)
	
	plt.ion()

#%%

subs=[
	"sub-P080",
	"sub-P079",
	"sub-P078",
	"sub-P077",
	"sub-P076",
	"sub-P074",
	"sub-P073",
	"sub-P072",
	"sub-P070",
	"sub-P066",
	"sub-P065",
	"sub-P062",
	"sub-P053",
	"sub-P052",
	"sub-P050",
	"sub-P046",
	"sub-P036",
	"sub-P021",  
  ]


for isub in subs:
	debug = True
	if debug:
		class dotdict(dict):
			"""dot.notation access to dictionary attributes"""
			__getattr__ = dict.get
			__setattr__ = dict.__setitem__
			__delattr__ = dict.__delitem__
		
		class Namespace:
			def __init__(self, **kwargs):
				self.__dict__.update(kwargs)
		
		
		data_dir=r'/media/veracrypt6/projects/SEEG/derivatives/atlasreg'
		
		input=dotdict({'t1':f'{data_dir}/{isub}/{isub}_desc-n4_T1w.nii.gz',
					'gm':f'{data_dir}/{isub}/{isub}_label-GM_desc-atropos3seg_probseg.nii.gz',
					'wm':f'{data_dir}/{isub}/{isub}_label-GM_desc-atropos3seg_probseg.nii.gz',
					'mask':f'{data_dir}/{isub}/{isub}_desc-brain_from-MNI152NLin2009cSym_reg-affine_mask.nii.gz',
					'seg':f'{data_dir}/{isub}/{isub}_desc-atroposKseg_dseg.nii.gz',
					})
		
		output=dotdict({
			't1_grad':f'{data_dir}/{isub}/{isub}_desc-magnitude_T1w.nii.gz',
			't1_intensity':f'{data_dir}/{isub}/{isub}_desc-intensity_T1w.nii.gz',
			'png_mask':f'{data_dir}/{isub}/qc/{isub}_desc-brain_maskqc.png',
			'png_seg':f'{data_dir}/{isub}/qc/{isub}_desc-segmentation_segqc.png'
		})
		
		snakemake = Namespace(output=output, input=input)
	
	
	t1_file=ants.image_read(snakemake.input.t1)
	mask_file=ants.image_read(snakemake.input.mask)
	seg_file=ants.image_read(snakemake.input.seg)
	
	if not os.path.exists(snakemake.output.t1_grad):
		grad_t1 = ants.iMath(t1_file, "Grad", 1)
		ants.image_write(grad_t1, snakemake.output.t1_grad)
	else:
		grad_t1=ants.image_read(snakemake.output.t1_grad)
	
	if not os.path.exists(snakemake.output.t1_intensity):
		gm_file=ants.image_read(snakemake.input.gm)
		wm_file=ants.image_read(snakemake.input.wm)
		t1_n4_gm = t1_file * t1_file.new_image_like(gm_file.numpy())
		t1_n4_wm = t1_file * t1_file.new_image_like(wm_file.numpy())
		bg_t1 = peakfinder(t1_n4_gm, t1_n4_wm, 1, 99.5)
		t1_ri = compute_RI(t1_file.numpy(), bg_t1, mask_file.numpy())
		tmp = t1_file.new_image_like(t1_ri)
		ri = ants.smooth_image(tmp, sigma=3, FWHM=True)
		ants.image_write(ri, snakemake.output.t1_intensity)
	else:
		ri=ants.image_read(snakemake.output.t1_intensity)
	
	lmin = float(t1_file.numpy().min())
	lmax = float(t1_file.numpy().max())
	norm_img=np.floor((t1_file.numpy()-lmin)/(lmax-lmin) * 255.)
	
	t1_file_nb = ants.to_nibabel(t1_file)
	t1_file_norm = nb.Nifti1Image(norm_img, t1_file_nb.affine, t1_file_nb.header)
	
	print(np.percentile(t1_file_norm.get_fdata(), 99), np.mean(t1_file_norm.get_fdata()))
	
	anat_img=t1_file_norm
	overlay_img= ants.to_nibabel(mask_file)
	out_file=snakemake.output.png_mask
	
	qc_mosaic_plot(t1_file_norm, ants.to_nibabel(mask_file), snakemake.output.png_mask)
	qc_mosaic_plot(t1_file_norm, ants.to_nibabel(seg_file), snakemake.output.png_seg, threshold=1)


