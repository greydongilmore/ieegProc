def get_pre_t1_filename(wildcards):
    files=glob(bids(root=join(config['out_dir'], 'bids'), subject=config['subject_prefix']+f'{wildcards.subject}', datatype='anat', session='pre', acq=config['noncontrast_t1']['acq'], run='*', suffix='T1w.nii.gz'))
    if len(files) <=1:
        file=expand(bids(root=join(config['out_dir'], 'bids'), subject=config['subject_prefix']+'{subject}', datatype='anat', session='pre', acq=config['noncontrast_t1']['acq'], run='01', suffix='T1w.nii.gz'),subject=wildcards.subject)
        if not exists(file[0]):
            file=expand(bids(root=join(config['out_dir'], 'bids'), subject=config['subject_prefix']+'{subject}', datatype='anat', session='pre', acq=config['noncontrast_t1']['acq'], run='02', suffix='T1w.nii.gz'),subject=wildcards.subject)
            if not exists(file[0]):
                file=expand(bids(root=join(config['out_dir'], 'bids'), subject=config['subject_prefix']+'{subject}', datatype='anat', session='pre', run='01', suffix='T1w.nii.gz'),subject=wildcards.subject)
    else:
        files.sort(key=lambda f: int(re.sub('\D', '', f)))
        file=files[-1]
    if file:
        if not exists(file[0]):
            file=expand(bids(root=join(config['out_dir'], 'bids'), subject=config['subject_prefix']+'{subject}', datatype='anat', session='pre', run='01', suffix='T1w.nii.gz'),subject=wildcards.subject)
    if file:
        print(f'Pre T1w non-contrast file: {basename(file[0])}')
    return file

def get_electrodes_filename(wildcards): 
    if wildcards.subject in config['subject_electrodes_custom']:
        return config['out_dir'] + config['subject_electrodes_custom'][wildcards.subject]
    else:
        return config['out_dir'] + config['subject_electrodes']

def get_transform_filename(wildcards):
    file=[]
    if config['contrast_t1']['present'] and config['noncontrast_t1']['present']:
        file=expand(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=config['subject_prefix']+'{subject}',suffix='xfm.txt',from_='noncontrast',to='contrast',desc='rigid',type_='ras'),subject=wildcards.subject),
    if len(file) >0:
        file=file[0][0]
    print(file)
    return file

if config['fastsurfer']['version'] =='dev':
    rule fastsurfer_seg:
        input: 
            t1 = get_pre_t1_filename,
        params:
            fastsurfer_run = config['fastsurfer']['home'],
            sid = config['fastsurfer']['sid'],
            batch = config['fastsurfer']['batch'],
            threads = config['fastsurfer']['threads'],
            order = config['fastsurfer']['order'],
            py = config['fastsurfer']['py'],
            fastsurfer_out = directory(join(config['out_dir'], 'derivatives', 'fastsurfer')),
            subjid=expand('sub-' + subject_id,subject=subjects),
        output:
            touch_fastsurfer=touch(join(config['out_dir'], 'logs', 'sub-' + subject_id + "_fastsurfer.done")),
            t1_fname = join(config['out_dir'],'derivatives','fastsurfer','sub-' + subject_id, 'mri','orig.mgz'),
        #threads:config['fastsurfer']['threads']
        shell:
            "export FASTSURFER_HOME={params.fastsurfer_run} &&PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:4096 {params.fastsurfer_run}/run_fastsurfer.sh \
            --t1 {input.t1} --sd {params.fastsurfer_out} --sid {params.subjid} --py {params.py} --viewagg_device cpu --fsaparc --parallel"
else:
    rule fastsurfer_seg:
        input: 
            t1 = get_pre_t1_filename,
        params:
            fastsurfer_run = config['fastsurfer']['home'],
            sid = config['fastsurfer']['sid'],
            batch = config['fastsurfer']['batch'],
            threads = config['fastsurfer']['threads'],
            order = config['fastsurfer']['order'],
            py = config['fastsurfer']['py'],
            fastsurfer_out = directory(join(config['out_dir'], 'derivatives', 'fastsurfer')),
            subjid=expand('sub-' + subject_id,subject=subjects),
        output:
            touch_fastsurfer=touch(join(config['out_dir'], 'logs', 'sub-' + subject_id + "_fastsurfer.done")),
            t1_fname = join(config['out_dir'],'derivatives','fastsurfer','sub-' + subject_id, 'mri','orig.mgz'),
        #threads:config['fastsurfer']['threads']
        shell:
            "export FASTSURFER_HOME={params.fastsurfer_run} &&PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:4096 {params.fastsurfer_run}/run_fastsurfer.sh \
            --t1 {input.t1} --sd {params.fastsurfer_out} --sid {params.subjid} --py {params.py} --viewagg_device cpu --fsaparc --parallel"
    
final_outputs.extend(expand(rules.fastsurfer_seg.output.touch_fastsurfer, subject=subjects))

if config['seeg_contacts']['present']:
    rule vis_electrodes_native:
        input: 
            fcsv = get_electrodes_filename,
            t1_fname = join(config['out_dir'],'derivatives','fastsurfer','sub-' + subject_id, 'mri','orig.mgz'),
            xfm_noncontrast = get_transform_filename,
        params:
            lh_pial = join(config['out_dir'],'derivatives','fastsurfer','sub-' + subject_id, 'surf','lh.pial'),
            rh_pial = join(config['out_dir'],'derivatives','fastsurfer','sub-' + subject_id, 'surf','rh.pial.T1'),
            lh_sulc = join(config['out_dir'],'derivatives','fastsurfer','sub-' + subject_id, 'surf','lh.sulc'),
            rh_sulc = join(config['out_dir'],'derivatives','fastsurfer','sub-' + subject_id, 'surf','rh.sulc'),
        output:
            html = report(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='electrodes.html',space='native'),
                    caption='../reports/electrodes_vis.rst',
                    category='Electrodes in template space',
                    subcategory='{desc} reg to {template}'),
        group: 'preproc'
        script: '../scripts/vis_electrodes_native.py'

    final_outputs.extend(expand(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='electrodes.html',space='native',include_subject_dir=False),
                        subject=subjects))
