def get_electrodes_filename(wildcards): 
    if wildcards.subject in config['subject_electrodes_custom']:
        return config['out_dir'] + config['subject_electrodes_custom'][wildcards.subject]
    else:
        return config['out_dir'] + config['subject_electrodes']

rule qc_reg:
    input:
        ref = config['template_t1w'],
        flo = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='T1w.nii.gz',space='{template}',desc='{desc}'),
    output:
        png = report(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='regqc.png',from_='subject', to='{template}',desc='{desc}',include_subject_dir=False),
                caption='../reports/regqc.rst',
                category='Registration QC',
                subcategory='{desc} {template}'),
        html = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='regqc.html',from_='subject', to='{template}', desc='{desc}',include_subject_dir=False),
#        html = report(bids(root='qc',subject=subject_id,suffix='regqc.html',from_='subject', to='{template}', desc='{desc}'),
#                caption='../reports/regqc.rst',
#                category='Registration QC',
#                subcategory='{desc} {template}'),
    group: 'preproc'
    script: '../scripts/vis_regqc.py'

if config['atlas_reg']['greedy']['active']:
    final_outputs.extend(expand(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='regqc.png',from_='subject', to='{template}',desc='{desc}',include_subject_dir=False),
                            subject=subjects, desc=['affine','greedydeform','SyN'],template=config['template']))
else:
    final_outputs.extend(expand(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='regqc.png',from_='subject', to='{template}',desc='{desc}',include_subject_dir=False),
                        subject=subjects, desc=['affine','SyN'],template=config['template']))

if config['noncontrast_t1']['present']:
    rule qc_reg_noncontrast:
        input:
            flo = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,acq='contrast',suffix='T1w.nii.gz'),
            ref = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,acq='noncontrast',space='T1w',desc='rigid',suffix='T1w.nii.gz'),
        output:
            png = report(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='regqc.png',from_='contrast', to='noncontrast',include_subject_dir=False),
                    caption='../reports/regqc.rst',
                    category='Registration QC',
                    subcategory='{desc} T1w'),
            html = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='regqc.html',from_='contrast', to='noncontrast', include_subject_dir=False),
    #        html = report(bids(root='qc',subject=subject_id,suffix='regqc.html',from_='subject', to='{template}', desc='{desc}'),
    #                caption='../reports/regqc.rst',
    #                category='Registration QC',
    #                subcategory='{desc} {template}'),
        group: 'preproc'
        script: '../scripts/vis_regqc.py'

    final_outputs.extend(expand(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='regqc.png',from_='contrast', to='noncontrast',include_subject_dir=False), 
                        subject=subjects))

if config['post_ct']['present']:
    rule qc_reg_ct:
        input:
            ref = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='T1w.nii.gz'),
            flo = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='ct.nii.gz',space='T1w',desc='rigid')
        output:
            png = report(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='regqc.png',from_='ct', to='T1w',desc='rigid',include_subject_dir=False),
                    caption='../reports/regqc.rst',
                    category='Registration QC',
                    subcategory='{desc} T1w'),
            html = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='regqc.html',from_='ct', to='T1w', desc='rigid',include_subject_dir=False),
    #        html = report(bids(root='qc',subject=subject_id,suffix='regqc.html',from_='subject', to='{template}', desc='{desc}'),
    #                caption='../reports/regqc.rst',
    #                category='Registration QC',
    #                subcategory='{desc} {template}'),
        group: 'preproc'
        script: '../scripts/vis_regqc.py'

    final_outputs.extend(expand(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='regqc.png',from_='ct', to='T1w',desc='rigid',include_subject_dir=False), 
                        subject=subjects))

if config['pet']['present']:
    rule qc_reg_pet:
        input:
            ref = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='T1w.nii.gz'),
            flo = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='pet.nii.gz',space='T1w', desc='rigid')
        output:
            png = report(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='regqc.png',from_='pet', to='T1w',desc='rigid',include_subject_dir=False),
                    caption='../reports/regqc.rst',
                    category='Registration QC',
                    subcategory='{desc} T1w'),
            html = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='regqc.html',from_='pet', to='T1w', desc='rigid',include_subject_dir=False),
    #        html = report(bids(root='qc',subject=subject_id,suffix='regqc.html',from_='subject', to='{template}', desc='{desc}'),
    #                caption='../reports/regqc.rst',
    #                category='Registration QC',
    #                subcategory='{desc} {template}'),
        group: 'preproc'
        script: '../scripts/vis_regqc.py'

    final_outputs.extend(expand(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='regqc.png',from_='pet', to='T1w',desc='rigid',include_subject_dir=False), 
                        subject=subjects))

#if config['pet']['present'] and config['fastsurfer']['run']:
#    rule qc_pet_surf:
#        input:
#            pet = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='pet.nii.gz',space='T1w', desc='rigid'),
#            lh_pial = join(config['out_dir'], 'derivatives', config['fastsurfer']['sid'], 'sub-' + subject_id,config['fastsurfer']['sid'],'surf','lh.pial'),
#            rh_pial = join(config['out_dir'], 'derivatives', config['fastsurfer']['sid'], 'sub-' + subject_id,config['fastsurfer']['sid'],'surf','rh.pial'),
#        output:
#            pet_png = report(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='brain4views.png',desc='PET',include_subject_dir=False),
#                    caption='../reports/regqc.rst',
#                    category='Registration QC',
#                    subcategory='{desc} T1w'),
#        group: 'preproc'
#        script: '../scripts/pet_surf_brain4views.py'
#
#    final_outputs.extend(expand(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='brain4views.png',desc='PET',include_subject_dir=False), 
#                        subject=subjects))
rule qc_probseg:
    input:
        img = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'), subject=subject_id, desc='masked', from_='atropos3seg', suffix='T1w.nii.gz'),
        seg4d = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='probseg.nii.gz',desc='atropos3seg'),
        mapping = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='mapping.json',desc='atropos3seg'),
    params:
        ct_png=bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='regqc.png',from_='ct', to='T1w',desc='masked',include_subject_dir=False)
    output:
        png = report(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='probseg.png',desc='brainmask',label='atropos3seg',include_subject_dir=False),
                caption='../reports/segqc.rst',
                category='Segmentation QC',
                subcategory='3-class Tissue Segmentation'),
    group: 'preproc'
    script: '../scripts/vis_qc_probseg.py'

final_outputs.extend(expand(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='probseg.png', desc='brainmask',label='atropos3seg',include_subject_dir=False),
                        subject=subjects ))

rule qc_dseg:
    input:
        img = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'), subject=subject_id, desc='masked', from_='atropos3seg', suffix='T1w.nii.gz'),
        seg = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'), subject=subject_id,suffix='dseg.nii.gz',atlas='{atlas}',from_='{template}',reg='SyN'),
    output:
        png = report(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='dseg.png',atlas='{atlas}', from_='{template}',desc='brainmask', include_subject_dir=False),
                caption='../reports/segqc.rst',
                category='Segmentation QC',
                subcategory='{atlas} Atlas from {template}'),
        html = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='dseg.html',atlas='{atlas}', from_='{template}',desc='brainmask',include_subject_dir=False),
#        html = report(bids(root='qc',subject=subject_id,suffix='dseg.html',atlas='{atlas}', from_='{template}'),
#                caption='../reports/segqc.rst',
#                category='Segmentation QC',
#                subcategory='{atlas} Atlas from {template}'),
    group: 'preproc'
    script: '../scripts/vis_qc_dseg.py'

final_outputs.extend(expand(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='dseg.png',atlas='{atlas}', from_='{template}',desc='brainmask',include_subject_dir=False),
                        subject=subjects, atlas=config['atlases'],template=config['template']))

rule qc_dseg_dilated:
    input:
        img = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'), subject=subject_id, desc='masked', from_='atropos3seg', suffix='T1w.nii.gz'),
        seg = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='dseg.nii.gz',atlas='{atlas}',from_='{template}',desc='dilated',reg='SyN'),
    output:
        png = report(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='dseg.png',atlas='{atlas}', from_='{template}',desc='brainmask',label='dilated',include_subject_dir=False),
                caption='../reports/segqc.rst',
                category='Segmentation QC',
                subcategory='{atlas} Atlas from {template} dilated'),
        html = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='dseg.html',atlas='{atlas}', from_='{template}',desc='brainmask',label='dilated',include_subject_dir=False),
#        html = report(bids(root='qc',subject=subject_id,suffix='dseg.html',atlas='{atlas}', from_='{template}'),
#                caption='../reports/segqc.rst',
#                category='Segmentation QC',
#                subcategory='{atlas} Atlas from {template}'),
    group: 'preproc'
    script: '../scripts/vis_qc_dseg.py'

final_outputs.extend(expand(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='dseg.png',atlas='{atlas}', from_='{template}',desc='brainmask',label='dilated',include_subject_dir=False),
                        subject=subjects, atlas=config['atlases'],template=config['template']))

rule qc_tissue_class:
    input:
        img = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'), subject=subject_id, desc='masked', from_='atropos3seg', suffix='T1w.nii.gz'),
        wm = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='probseg.nii.gz',label='WM',desc='atropos3seg'),
        gm = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='probseg.nii.gz',label='GM',desc='atropos3seg'),
        csf = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='probseg.nii.gz',label='CSF',desc='atropos3seg'),
    output:
        png = report(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id, suffix='probseg.png', desc='brainmask', label='classes', include_subject_dir=False),
                caption='../reports/segqc.rst',
                category='Segmentation QC',
                subcategory='Tissue classification'),
    group: 'preproc'
    script: '../scripts/vis_qc_tissue_seg.py'

final_outputs.extend(expand(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='probseg.png', desc='brainmask', label='classes',include_subject_dir=False),
                        subject=subjects))


if config['seeg_contacts']['present']:
    rule vis_contacts:
        input:
            ct = bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='ct.nii.gz',space='T1w',desc='rigid'),
            mask = bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='contacts.nii.gz',space='ct',desc='mask'),
        output:
            html = report(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='contacts.html',desc='mask',space='ct',include_subject_dir=False),
                    caption='../reports/contacts_vis.rst',
                    category='Contacts in CT space',
                    subcategory='landmarks mask'),
        group: 'preproc'
        script: '../scripts/vis_contacts.py'

    rule vis_electrodes:
        input: 
            fcsv = get_electrodes_filename,
            t1w = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,desc='n4', suffix='T1w.nii.gz'),
            xfm_ras = bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='xfm.txt',from_='subject',to='{template}',desc='affine',type_='ras'),
        params:
            contacts= bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='contacts.html',desc='mask',space='ct',include_subject_dir=False)
        output:
            html = report(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='electrodes.html',desc='affine',space='{template}'),
                    caption='../reports/electrodes_vis.rst',
                    category='Electrodes in template space',
                    subcategory='{desc} reg to {template}'),
            png = report(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='electrodevis.png',desc='affine',space='{template}',include_subject_dir=False),
                    caption='../reports/electrodes_vis.rst',
                    category='Electrodes in template space',
                    subcategory='{desc} reg to {template}'),
        group: 'preproc'
        script: '../scripts/vis_electrodes.py'

    final_outputs.extend(expand(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='electrodevis.png',desc='affine',space='{template}',include_subject_dir=False),
                        subject=subjects, desc=['rigid'],template=config['template']))
    
    final_outputs.extend(expand(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='electrodes.html',desc='affine',space='{template}',include_subject_dir=False),
                        subject=subjects, desc=['rigid'],template=config['template']))
    
    final_outputs.extend(expand(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='contacts.html',desc='mask',space='ct',include_subject_dir=False),
            subject=subjects))