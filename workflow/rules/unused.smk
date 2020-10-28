#-- affine directly
rule rigid_aladin:
    input: 
        flo = bids(root='results',subject='{subject}',suffix='T1w.nii.gz'),
        ref = config['template_t1w'],
    output: 
        warped_subj = bids(root='results',subject='{subject}',suffix='T1w.nii.gz',space='{template}',desc='rigid'),
        xfm_ras = bids(root='results',subject='{subject}',suffix='xfm.txt',from_='subject',to='{template}',desc='rigid',type_='ras'),
    container: config['singularity']['neuroglia']
    group: 'preproc'
    shell:
        'reg_aladin -rigOnly -flo {input.flo} -ref {input.ref} -res {output.warped_subj} -aff {output.xfm_ras}'


#=======  below not used, too computationally heavy..

#split atlas discrete seg into multiple images, one for each label, and smooth each one
#then normalize to sum to 1
rule create_atlas_probseg:
    input:
        dseg = config['template_atlas_dseg_nii'],
    params:
        smoothing_kernel = '3mm',
    output:
        probseg = bids(root='results',prefix='tpl-{template}/tpl-{template}',atlas='{atlas}',suffix='probseg.nii.gz')
    shadow: 'minimal' #runs in isolated dir so we can use intermediate files (e.g. combined.nii), and throw them out after
    container: config['singularity']['neuroglia']
    shell:
        'c3d -verbose {input.dseg} -split  -foreach -smooth {params.smoothing_kernel} -endfor -oo smoothed_%03d.nii && '
        'c3d -verbose smoothed_*.nii -accum -add -endaccum -o sum.nii && '
        'fslmerge -t combined.nii smoothed_*.nii && '
        'fslmaths combined.nii -div sum.nii {output.probseg}'

rule warp_probseg_from_template:
    input: 
        probseg = bids(root='results',prefix='tpl-{template}/tpl-{template}',atlas='{atlas}',suffix='probseg.nii.gz'),
        ref = bids(root='results',subject='{subject}',suffix='T1w.nii.gz'),
        inv_composite = bids(root='results',suffix='InverseComposite.h5',from_='subject',to='{template}',subject='{subject}'),
    output:
        probseg = bids(root='results',subject='{subject}',suffix='probseg.nii.gz',atlas='{atlas}',from_='{template}'),
    container: config['singularity']['neuroglia']
    group: 'preproc'
    threads: 8
    resources:
        mem_mb = 16000
    shadow: 'minimal'
    shell: 
        'FSLOUTPUTTYPE=NIFTI fslsplit {input.probseg} in_seg -t && '
        'for seg in `ls in_seg*.nii`; do '
        'ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} '
        'antsApplyTransforms -d 3 --interpolation Linear -i ${{seg}} -o out_seg${{seg##in_seg}} -r {input.ref} ' 
            ' -t {input.inv_composite}; done && ' #use inverse xfm (going from template to subject)
        'FSLOUTPUTTYPE=NIFTI_GZ fslmerge -t {output.probseg} out_seg*.nii'

#--- atropos unused:

rule tissue_segk4:
    input: 
        t1 = bids(root='results',subject='{subject}',desc='n4',from_='{template}', suffix='T1w.nii.gz'),
        mask = bids(root='results',subject='{subject}',suffix='mask.nii.gz',from_='{template}',reg='SyN',desc='brain'),
    params:
        prior_fmt = 'results/sub-{subject}/prior_%d.nii.gz',
        prior_f = '0.{prior_f}',
        mrf_f = '0.{mrf_f}'    
    output:
        seg = bids(root='results',subject='{subject}',from_='{template}',suffix='dseg.nii.gz',prior='{prior_f}',mrf='{mrf_f}',desc='atroposk4')
    log: bids(root='logs',subject='{subject}',from_='{template}',suffix='log.txt',prior='{prior_f}',mrf='{mrf_f}',desc='atropos') 
    container: config['singularity']['neuroglia']
    shell:
        'Atropos -d 3 -a {input.t1} -i PriorProbabilityImages[4,{params.prior_fmt},{params.prior_f}] -x {input.mask} -o [{output.seg}] &> {log}'

rule tissue_seg:
    input: 
        t1 = bids(root='results',subject='{subject}',desc='n4',from_='{template}', suffix='T1w.nii.gz'),
        mask = bids(root='results',subject='{subject}',suffix='mask.nii.gz',from_='{template}',reg='SyN',desc='brain'),
    params:
        prior_fmt = 'results/sub-{subject}/prior_%d.nii.gz',
        prior_f = '0.{prior_f}',
        mrf_f = '0.{mrf_f}'    
    output:
        seg = bids(root='results',subject='{subject}',from_='{template}',suffix='dseg.nii.gz',prior='{prior_f}',mrf='{mrf_f}',desc='atropos')
    log: bids(root='logs',subject='{subject}',from_='{template}',suffix='log.txt',prior='{prior_f}',mrf='{mrf_f}',desc='atropos') 
    container: config['singularity']['neuroglia']
    shell:
        'Atropos -v -d 3 -a {input.t1} -i PriorProbabilityImages[3,{params.prior_fmt},{params.prior_f}] -x {input.mask} -o [{output.seg}] &> {log}'

rule tissue_seg_priorinit:
    input: 
        t1 = bids(root='results',subject='{subject}',desc='n4',from_='{template}', suffix='T1w.nii.gz'),
        mask = bids(root='results',subject='{subject}',suffix='mask.nii.gz',from_='{template}',reg='SyN',desc='brain'),
    params:
        init = 'PriorProbabilityImages[{k},results/sub-{subject}/prior_%d.nii.gz,0]',
    output:
        seg = bids(root='results',subject='{subject}',from_='{template}',suffix='dseg.nii.gz',method='PriorProbabilityImages',k='{k}',desc='atropos')
    container: config['singularity']['neuroglia']
    shell:
        'Atropos -v -d 3 -a {input.t1} -i {params.init} -x {input.mask} -o [{output.seg}]'



