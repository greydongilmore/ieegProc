
def get_electrodes_filename(wildcards): 
    if wildcards.subject in config['subject_electrodes_custom']:
        return config['out_dir'] + config['subject_electrodes_custom'][wildcards.subject]
    else:
        return expand(bids(root=join(config['out_dir'], 'derivatives',config['seeg_contacts']['dirname_coords'].split('/')[2]), subject=config['subject_prefix']+f'{wildcards.subject}', space='native', suffix='SEEGA.fcsv'))[0]

def get_reference_t1(wildcards):
    if config['contrast_t1']['present']:
        ref_file=expand(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'), subject=config['subject_prefix']+'{subject}', acq='contrast', suffix='T1w.nii.gz'),subject=wildcards.subject)
    elif not config['contrast_t1']['present'] and config['noncontrast_t1']['present']:
        ref_file=expand(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'), subject=config['subject_prefix']+'{subject}', acq='noncontrast', suffix='T1w.nii.gz'),subject=wildcards.subject)
    return ref_file[0]

rule qc_reg_t1:
    input:
        ref = get_age_appropriate_template_name(expand(subject_id,subject=subjects),'t1w'),
        flo = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='T1w.nii.gz', space=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'),desc='{desc}'),
    output:
        png = report(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='regqc.png',from_='subject', to=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'),desc='{desc}',include_subject_dir=False),
                caption='../reports/regqc.rst',
                category='Registration QC',
                subcategory=lambda wildcards, input: f"{get_age_appropriate_template_name(input,'space')}"),
        html = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='regqc.html',from_='subject', to=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'), desc='{desc}', include_subject_dir=False),
#        html = report(bids(root='qc',subject=subject_id,suffix='regqc.html',from_='subject', to=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'), desc='{desc}'),
#                caption='../reports/regqc.rst',
#                category='Registration QC',
#                subcategory='{desc} {template}'),
    group: 'preproc'
    script: '../scripts/vis_regqc.py'

if config['nonlin_reg']['algo']=='greedy':
    final_outputs.extend(expand(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='regqc.png',from_='subject', to=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'),desc='{desc}',include_subject_dir=False),
                            subject=subjects, desc=['affine','nonlin']))
else:
    final_outputs.extend(expand(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='regqc.png',from_='subject', to=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'),desc='{desc}',include_subject_dir=False),
                        subject=subjects, desc=['affine','nonlin']))

if config['contrast_t1']['present'] and config['noncontrast_t1']['present']:
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
    #        html = report(bids(root='qc',subject=subject_id,suffix='regqc.html',from_='subject', to=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'), desc='{desc}'),
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
            ref = get_reference_t1,
            flo = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='ct.nii.gz',space='T1w',desc='rigid'),
        output:
            png = report(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='regqc.png',from_='ct', to='T1w',desc='rigid',include_subject_dir=False),
                    caption='../reports/regqc.rst',
                    category='Registration QC',
                    subcategory='{desc} T1w'),
            html = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='regqc.html',from_='ct', to='T1w', desc='rigid',include_subject_dir=False),
    #        html = report(bids(root='qc',subject=subject_id,suffix='regqc.html',from_='subject', to=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'), desc='{desc}'),
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
            ref = get_reference_t1,
            flo = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='pet.nii.gz',space='T1w', desc='rigid'),
        output:
            png = report(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='regqc.png',from_='pet', to='T1w',desc='rigid',include_subject_dir=False),
                    caption='../reports/regqc.rst',
                    category='Registration QC',
                    subcategory='{desc} T1w'),
            html = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='regqc.html',from_='pet', to='T1w', desc='rigid',include_subject_dir=False),
    #        html = report(bids(root='qc',subject=subject_id,suffix='regqc.html',from_='subject', to=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'), desc='{desc}'),
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
        ct_png=bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='regqc.png',from_='ct', to='T1w',desc='masked',include_subject_dir=False),
        tissue_classes=config['default_k_tissue_classes'],
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
        seg = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'), subject=subject_id,suffix='dseg.nii.gz',atlas='{atlas}',from_=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'),desc='nonlin'),
    output:
        png = report(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='dseg.png',atlas='{atlas}', from_=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'),desc='brainmask', include_subject_dir=False),
                caption='../reports/segqc.rst',
                category='Segmentation QC',
                ),
        html = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='dseg.html',atlas='{atlas}', from_=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'),desc='brainmask',include_subject_dir=False),
#        html = report(bids(root='qc',subject=subject_id,suffix='dseg.html',atlas='{atlas}', from_=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space')),
#                caption='../reports/segqc.rst',
#                category='Segmentation QC',
#                subcategory='{atlas} Atlas from {template}'),
    group: 'preproc'
    script: '../scripts/vis_qc_dseg.py'

final_outputs.extend(expand(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='dseg.png',atlas='{atlas}', from_=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'),desc='brainmask',include_subject_dir=False),
                        subject=subjects, atlas=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'atlas')))

rule qc_dseg_dilated:
    input:
        img = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'), subject=subject_id, desc='masked', from_='atropos3seg', suffix='T1w.nii.gz'),
        seg = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='dseg.nii.gz',atlas='{atlas}',from_=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'),desc='nonlin',label='dilated'),
    output:
        png = report(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='dseg.png',atlas='{atlas}', from_=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'),desc='brainmask',label='dilated',include_subject_dir=False),
                caption='../reports/segqc.rst',
                category='Segmentation QC'),
        html = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='dseg.html',atlas='{atlas}', from_=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'),desc='brainmask',label='dilated',include_subject_dir=False),
#        html = report(bids(root='qc',subject=subject_id,suffix='dseg.html',atlas='{atlas}', from_=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space')),
#                caption='../reports/segqc.rst',
#                category='Segmentation QC',
#                subcategory='{atlas} Atlas from {template}'),
    group: 'preproc'
    script: '../scripts/vis_qc_dseg.py'

final_outputs.extend(expand(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='dseg.png',atlas='{atlas}', from_=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'),desc='brainmask',label='dilated',include_subject_dir=False),
                        subject=subjects, atlas=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'atlas')))

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
            xfm_ras = bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='xfm.txt',from_='subject',to=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'),desc='affine',type_='ras'),
        params:
            contacts= bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='contacts.html',desc='mask',space='ct',include_subject_dir=False)
        output:
            html = report(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='electrodes.html',desc='affine',space=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space')),
                    caption='../reports/electrodes_vis.rst',
                    category='Electrodes in template space',
                    subcategory="reg to {get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space')}"),
            png = report(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='electrodevis.png',desc='affine',space=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'),include_subject_dir=False),
                    caption='../reports/electrodes_vis.rst',
                    category='Electrodes in template space',
                    subcategory=lambda wildcards, input: f"reg to {get_age_appropriate_template_name(input,'space')}"),
        group: 'preproc'
        script: '../scripts/vis_electrodes.py'

    final_outputs.extend(expand(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='electrodevis.png',desc='affine',space=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'),include_subject_dir=False),
                        subject=subjects, desc=['rigid']))
    
    final_outputs.extend(expand(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='electrodes.html',desc='affine',space=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'),include_subject_dir=False),
                        subject=subjects, desc=['rigid']))
    
    final_outputs.extend(expand(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='contacts.html',desc='mask',space='ct',include_subject_dir=False),
            subject=subjects))