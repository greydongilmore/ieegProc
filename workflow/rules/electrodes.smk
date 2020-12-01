
def get_electrodes_filename(wildcards): 
    if wildcards.subject in config['subject_electrodes_custom']:
        return config['out_dir'] + config['subject_electrodes_custom'][wildcards.subject]
    else:
        return config['out_dir'] + config['subject_electrodes']

rule electrode_coords:
    input:
        seega_scene = config['subject_seega_scene']
    params:
        sub=subject_id
    output:
        seega_fcsv = bids(root=join(config['out_dir'],'deriv','seega_coordinates'),subject=subject_id,space='native', suffix='SEEGA.fcsv'),
    group: 'preproc'
    script: '../scripts/working/elec_labels_coords.py'

rule label_electrodes_atlas:
    input: 
        fcsv = get_electrodes_filename,
        dseg_tsv = config['template_atlas_dseg_tsv'],
        dseg_nii = bids(root=join(config['out_dir'], 'deriv', 'atlasreg'),subject=subject_id,suffix='dseg.nii.gz',desc='dilated',atlas='{atlas}',from_='{template}',reg='SyN'),
        tissue_seg = expand(bids(root=join(config['out_dir'], 'deriv', 'atlasreg'),subject=subject_id,suffix='probseg.nii.gz',label='{tissue}',desc='atropos3seg'),
                            tissue=config['tissue_labels'],allow_missing=True),
    output:
        tsv = bids(root=join(config['out_dir'], 'deriv', 'atlasreg'),subject=subject_id,suffix='electrodes.tsv',atlas='{atlas}',desc='dilated',from_='{template}'),
#        tsv = report(bids(root='results',subject=subject_id,suffix='electrodes.tsv',desc='{atlas}',from_='{template}'),
#                caption='../reports/electrodes_vis.rst',
#                category='Electrodes Labelled',
#                subcategory='Atlas: {atlas}, Template: {template}')           
    group: 'preproc'
    script: '../scripts/label_electrodes_atlas.py'

rule contact_landmarks:
    input: 
        fcsv = get_electrodes_filename,
    output:
        txt = bids(root=join(config['out_dir'], 'deriv', 'atlasreg'),subject=subject_id,suffix='landmarks.txt',space='ct'),    
    group: 'preproc'
    run:
        df = pd.read_table(input.fcsv,sep=',',header=2)
        coords = df[['x','y','z']].to_numpy()
        with open (output.txt, 'w') as fid:
            for i in range(len(coords)):
                fid.write(' '.join(str(i) for i in np.r_[np.round(coords[i,:],3),int(1)])+ "\n")

rule mask_contacts:
    input: 
        ct = bids(root=join(config['out_dir'],'deriv', 'atlasreg'),subject=subject_id,suffix='ct.nii.gz',space='T1w',desc='affine'),
        txt = bids(root=join(config['out_dir'], 'deriv', 'atlasreg'),subject=subject_id,suffix='landmarks.txt',space='ct'),
    output:
        mask = bids(root=join(config['out_dir'],'deriv', 'atlasreg'),subject=subject_id,suffix='contacts.nii.gz',space='ct',desc='mask'),
    group: 'preproc'
    shell:
        'c3d {input.ct} -scale 0 -landmarks-to-spheres {input.txt} 1 -o {output.mask}'

rule vis_contacts:
    input:
        ct = bids(root=join(config['out_dir'],'deriv', 'atlasreg'),subject=subject_id,suffix='ct.nii.gz',space='T1w',desc='affine'),
        mask = bids(root=join(config['out_dir'],'deriv', 'atlasreg'),subject=subject_id,suffix='contacts.nii.gz',space='ct',desc='mask'),
    output:
        html = report(bids(root=join(config['out_dir'], 'deriv', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='contacts.html',desc='mask',space='ct',include_subject_dir=False),
                caption='../reports/contacts_vis.rst',
                category='Contacts in CT space',
                subcategory='landmarks mask'),
    group: 'preproc'
    script: '../scripts/vis_contacts.py'

rule vis_electrodes:
    input: 
        fcsv = get_electrodes_filename,
        t1w = bids(root=join(config['out_dir'], 'deriv', 'atlasreg'),subject=subject_id,desc='n4', suffix='T1w.nii.gz'),
        xfm_ras = bids(root=join(config['out_dir'], 'deriv', 'atlasreg'),subject=subject_id,suffix='xfm.txt',from_='subject',to='{template}',desc='{desc}',type_='ras'),
    params:
        contacts= bids(root=join(config['out_dir'], 'deriv', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='contacts.html',desc='mask',space='ct',include_subject_dir=False)
    output:
        html = report(bids(root=join(config['out_dir'], 'deriv', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='electrodes.html',desc='{desc}',space='{template}'),
                caption='../reports/electrodes_vis.rst',
                category='Electrodes in template space',
                subcategory='{desc} reg to {template}'),
        png = report(bids(root=join(config['out_dir'], 'deriv', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='electrodevis.png',desc='{desc}',space='{template}',include_subject_dir=False),
                caption='../reports/electrodes_vis.rst',
                category='Electrodes in template space',
                subcategory='{desc} reg to {template}'),
    group: 'preproc'
    script: '../scripts/vis_electrodes.py'
