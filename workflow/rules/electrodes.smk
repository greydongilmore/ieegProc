
def get_electrodes_filename(wildcards): 
    if wildcards.subject in config['subject_electrodes_custom']:
        return config['out_dir'] + config['subject_electrodes_custom'][wildcards.subject]
    else:
        return config['out_dir'] + config['subject_electrodes']

def get_fcsv_files(wildcards):
    file=glob(bids(root=join(config['out_dir'], 'derivatives','seega_coordinates'), subject=config['subject_prefix']+f'{wildcards.subject}', space='native', suffix='*.fcsv'))
    print(f'fcsv file: {file}')
    return file

def get_fcsv_files_acpc(wildcards):
    file=glob(bids(root=join(config['out_dir'], 'derivatives','seega_scenes','sub-'+config['subject_prefix']+f'{wildcards.subject}'), suffix='acpc.fcsv'))
    print(f'fcsv file: {file}')
    return file

def get_noncontrast_T1w(wildcards):
    file=glob(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=config['subject_prefix']+f'{wildcards.subject}', acq='noncontrast', space='T1w',desc='rigid',suffix='T1w.nii.gz'))
    print(f'fcsv file: {file}')
    return file

def get_contrast_T1w(wildcards):
    file=glob(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'), subject=config['subject_prefix']+f'{wildcards.subject}', acq='contrast', suffix='T1w.nii.gz'))
    print(f'fcsv file: {file}')
    return file

def get_segs(wildcards):
    file=glob(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'), subject=config['subject_prefix']+f'{wildcards.subject}', label='*', desc='atropos3seg', suffix='probseg.nii.gz'))
    print(f'fcsv file: {file}')
    return file

def get_atlas_segs(wildcards):
    file=glob(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=config['subject_prefix']+f'{wildcards.subject}', atlas=config['atlases'][0], 
        from_=config['template'],reg='SyN',suffix='dseg.nii.gz'))
    print(f'fcsv file: {file}')
    return file

def get_ct_file(wildcards):
    file=glob(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=config['subject_prefix']+f'{wildcards.subject}', space='T1w', desc='rigid', suffix='ct.nii.gz'))
    print(f'ct file: {file}')
    return file

def get_pet_file(wildcards):
    file=glob(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=config['subject_prefix']+f'{wildcards.subject}', space='T1w', desc='rigid', suffix='pet.nii.gz'))
    print(f'pet file: {file}')
    return file

rule electrode_coords:
    input:
        seega_scene = config['subject_seega_scene']
    params:
        sub=subject_id
    output:
        seega_fcsv = bids(root=join(config['out_dir'],'derivatives','seega_coordinates'),subject=subject_id,space='native', suffix='SEEGA.fcsv'),
    group: 'preproc'
    script: '../scripts/working/elec_labels_coords.py'

rule label_electrodes_atlas:
    input: 
        fcsv = bids(root=join(config['out_dir'],'derivatives','seega_coordinates'),subject=subject_id,space='native', suffix='SEEGA.fcsv'),
        dseg_tsv = config['template_atlas_dseg_tsv'],
        dseg_nii = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='dseg.nii.gz',desc='dilated',atlas='{atlas}',from_='{template}',reg='SyN'),
        tissue_seg = expand(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='probseg.nii.gz',label='{tissue}',desc='atropos3seg'),
                            tissue=config['tissue_labels'],allow_missing=True),
    output:
        tsv = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='electrodes.tsv',atlas='{atlas}',desc='dilated',from_='{template}'),
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
        txt = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='landmarks.txt',space='ct'),    
    group: 'preproc'
    run:
        df = pd.read_table(input.fcsv,sep=',',header=2)
        coords = df[['x','y','z']].to_numpy()
        with open (output.txt, 'w') as fid:
            for i in range(len(coords)):
                fid.write(' '.join(str(i) for i in np.r_[np.round(coords[i,:],3),int(1)])+ "\n")

rule mask_contacts:
    input: 
        ct = bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='ct.nii.gz',space='T1w',desc='rigid'),
        txt = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='landmarks.txt',space='ct'),
    output:
        mask = bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='contacts.nii.gz',space='ct',desc='mask'),
    group: 'preproc'
    shell:
        'c3d {input.ct} -scale 0 -landmarks-to-spheres {input.txt} 1 -o {output.mask}'

rule generate_slicer_directory:
    input:
        fcsv_files=get_fcsv_files,
        fcsv_acpc=get_fcsv_files_acpc,
        ct = bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='ct.nii.gz',space='T1w',desc='rigid'),
        noncontrast_t1w=get_noncontrast_T1w,
        contrast_t1w=get_contrast_T1w,
        segs=get_segs,
        atlas_segs=bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,atlas=config['atlases'][0], from_=config['template'],suffix='dseg.nii.gz',reg='SyN'),
        ct_file=get_ct_file,
        pet_file=get_pet_file,
    output:
        touch_slicer=touch(bids(root=join(config['out_dir'], 'derivatives', 'slicer_scene'), subject=subject_id, suffix='slicer.done')),
    script:
        '../scripts/slicer_dir.py'
        



final_outputs.extend(
    expand(
        bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),
            subject=subject_id,
            suffix='electrodes.tsv',
            atlas='{atlas}',
            desc='dilated',
            from_='{template}'
        ),
        subject=subjects,
        atlas=config['atlases'],
        template=config['template']
    )
)



final_outputs.extend(
    expand(
        bids(root=join(config['out_dir'], 'derivatives', 'slicer_scene'),
            subject=subject_id,
            suffix='slicer.done'
        ),
        subject=subjects
    )
)
