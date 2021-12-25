def get_k_tissue_classes(wildcards):
    if wildcards.subject in config['subject_k_tissue_classes']:
        return config['subject_k_tissue_classes'][wildcards.subject]
    else:
        return config['default_k_tissue_classes']

#this performs Atropos with k-means as initialization
rule tissue_seg_kmeans_init:
    input:
        t1 = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'), subject=subject_id,desc='n4', suffix='T1w.nii.gz'),
        mask = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='mask.nii.gz',from_='{template}'.format(template=config['template']),reg='affine',desc='brain'),
    params:
        k = get_k_tissue_classes,
        posterior_fmt = 'posteriors_%d.nii.gz',
        posterior_glob = 'posteriors_*.nii.gz',
    output:
        seg = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='dseg.nii.gz',desc='atroposKseg'),
        posteriors = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='probseg.nii.gz',desc='atroposKseg'),
    shadow: 'minimal'
    #container: config['singularity']['neuroglia']
    group: 'preproc'
    shell:
        'Atropos -d 3 -a {input.t1} -i KMeans[{params.k}] -x {input.mask} -o [{output.seg},{params.posterior_fmt}] && '
        'fslmerge -t {output.posteriors} {params.posterior_glob} ' #merge posteriors into a 4d file (intermediate files will be removed b/c shadow)

rule map_channels_to_tissue:
    input:
        tissue_priors = expand(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='probseg.nii.gz',label='{tissue}',from_='{template}'.format(template=config['template']),reg='affine'),
                            tissue=config['tissue_labels'],allow_missing=True),
        seg_channels_4d = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='probseg.nii.gz',desc='atroposKseg'),
    output:
        mapping_json = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='mapping.json',desc='atropos3seg'),
        tissue_segs = expand(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='probseg.nii.gz',label='{tissue}',desc='atropos3seg'),
                            tissue=config['tissue_labels'],allow_missing=True),
    group: 'preproc'
    script: '../scripts/map_channels_to_tissue.py'
       
rule tissue_seg_to_4d:
    input:
        tissue_segs = expand(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='probseg.nii.gz',label='{tissue}',desc='atropos3seg'),
                            tissue=config['tissue_labels'],allow_missing=True),
    output:
        tissue_seg = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='probseg.nii.gz',desc='atropos3seg')
    group: 'preproc'
    #container: config['singularity']['neuroglia']
    shell:
        'fslmerge -t {output} {input}'

rule brainmask_from_tissue:
    input:
        tissue_seg = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='probseg.nii.gz',desc='atropos3seg')
    params:
        threshold = 0.5
    output:
        mask = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='mask.nii.gz',from_='atropos3seg',desc='brain')
    #container: config['singularity']['neuroglia']
    group: 'preproc'
    shell:
        #max over tissue probs, threshold, binarize, fill holes
       'fslmaths {input} -Tmax -thr {params.threshold} -bin -fillh {output}'     
 

#TODO: make lesion mask from the holes in the brainmask (instead of just filling them..) -- could be a nice way to exclude contrast enhanced vessels
