def get_electrodes_filename(wildcards): 
    if wildcards.subject in config['subject_electrodes_custom']:
        return config['out_dir'] + config['subject_electrodes_custom'][wildcards.subject]
    else:
        return config['out_dir'] + config['subject_electrodes']

rule electrode_coords:
    input:
        seega_scene = config['subject_seega_scene']
    output:
        seega_fcsv = bids(root=join(config['out_dir'],'deriv','seega_coordinates'),subject='{subject}',space='native', suffix='SEEGA.fcsv'),
    group: 'preproc'
    script: '../scripts/working/elec_labels_coords.py'

rule vis_electrodes:
    input: 
        fcsv = get_electrodes_filename,
        t1w = bids(root=join(config['out_dir'],'results'),subject='{subject}',desc='n4', suffix='T1w.nii.gz'),
        xfm_ras = bids(root=join(config['out_dir'],'results'),subject='{subject}',suffix='xfm.txt',from_='subject',to='{template}',desc='{desc}',type_='ras'),
    output:
        html = report(bids(root=join(config['out_dir'],'qc'),subject='{subject}',suffix='electrodes.html',desc='{desc}',space='{template}'),
                caption='../reports/electrodes_vis.rst',
                category='Electrodes in template space',
                subcategory='{desc} reg to {template}'),
        png = report(bids(root=join(config['out_dir'],'qc'),subject='{subject}',suffix='electrodevis.png',desc='{desc}',space='{template}'),
                caption='../reports/electrodes_vis.rst',
                category='Electrodes in template space',
                subcategory='{desc} reg to {template}'),
    group: 'preproc'
    script: '../scripts/vis_electrodes.py'

rule label_electrodes_atlas:
    input: 
        fcsv = get_electrodes_filename,
        dseg_tsv = config['template_atlas_dseg_tsv'],
        dseg_nii = bids(root=join(config['out_dir'],'results'),subject='{subject}',suffix='dseg.nii.gz',desc='dilated',atlas='{atlas}',from_='{template}',reg='SyN'),
        tissue_seg = expand(bids(root=join(config['out_dir'],'results'),subject='{subject}',suffix='probseg.nii.gz',label='{tissue}',desc='atropos3seg'),
                            tissue=config['tissue_labels'],allow_missing=True),
    output:
        tsv = bids(root=join(config['out_dir'],'results'),subject='{subject}',suffix='electrodes.tsv',atlas='{atlas}',desc='dilated',from_='{template}'),
#        tsv = report(bids(root='results',subject='{subject}',suffix='electrodes.tsv',desc='{atlas}',from_='{template}'),
#                caption='../reports/electrodes_vis.rst',
#                category='Electrodes Labelled',
#                subcategory='Atlas: {atlas}, Template: {template}')           
    group: 'preproc'
    script: '../scripts/label_electrodes_atlas.py'


