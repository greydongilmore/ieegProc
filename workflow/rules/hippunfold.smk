def get_noncontrast_filename(wildcards):
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

def get_pre_t1_filename(wildcards):
    files=glob(bids(root=join(config['out_dir'], 'bids'), subject=config['subject_prefix']+f'{wildcards.subject}', datatype='anat', session='pre', run='*', suffix='T1w.nii.gz'))
    if len(files)==0:
        file=expand(bids(root=join(config['out_dir'], 'bids'), subject=config['subject_prefix']+'{subject}', datatype='anat', session='pre', acq=config['noncontrast_t1']['acq'], run='01', suffix='T1w.nii.gz'), subject=wildcards.subject)
    else:
        files.sort(key=lambda f: int(re.sub('\D', '', f)))
        file=files[-1]
    
    if file:
        print(f'Pre T1w contrast file: {basename(file)}')
    return file

if config['noncontrast_t1']['present']:
    rule import_subj_t1_hippunfold:
        input: 
            in_t1w=get_noncontrast_filename,
        params:
            in_dir = directory(join(config['out_dir'], 'derivatives', 'hippunfold_in')),
            in_dir_sub = directory(join(config['out_dir'], 'derivatives', 'hippunfold_in','sub-'+subject_id)),
        output: 
            out_t1w=bids(root=join(config['out_dir'], 'derivatives','hippunfold_in'), subject=subject_id, datatype='anat', suffix='T1w.nii.gz')
        group: 'preproc'
        shell: "mkdir -p {params.in_dir} && mkdir -p {params.in_dir_sub} && mkdir -p {params.in_dir_sub}/anat && cp {input.in_t1w} {output.out_t1w}"

else:
    rule import_subj_t1_hippunfold:
        input: 
            in_t1w=get_pre_t1_filename,
        params:
            in_dir = directory(join(config['out_dir'], 'derivatives', 'hippunfold_in')),
            in_dir_sub = directory(join(config['out_dir'], 'derivatives', 'hippunfold_in','sub-'+subject_id)),
        output: 
            out_t1w=bids(root=join(config['out_dir'], 'derivatives','hippunfold_in'), subject=subject_id, datatype='anat', suffix='T1w.nii.gz')
        group: 'preproc'
        shell: 
            "mkdir -p {params.in_dir} && mkdir -p {params.in_dir_sub} && mkdir -p {params.in_dir_sub}/anat && cp {input.in_t1w} {output.out_t1w}"

rule hippunfold_seg:
    input: 
        out_t1w=bids(root=join(config['out_dir'], 'derivatives','hippunfold_in'), subject=subject_id, datatype='anat', suffix='T1w.nii.gz')
    params:
        in_dir = directory(join(config['out_dir'], 'derivatives', 'hippunfold_in')),
        modality = config['hippunfold']['modality'],
        path_T1w = bids(root=join(config['out_dir'], 'derivatives','hippunfold_in'), subject=subject_id, datatype='anat', suffix='T1w.nii.gz'),
        participant_label = subject_id,
        hippunfold_out = directory(join(config['out_dir'], 'derivatives', 'hippunfold')),
        hippunfold_container= config['singularity']['hippunfold'],
    output:
        touch_hippunfold=touch(join(config['out_dir'], 'logs', 'sub-' + subject_id + "_hippunfold.done")),
        t1_fname = join(config['out_dir'], 'derivatives','hippunfold','hippunfold','sub-' + subject_id, 'anat','sub-' + subject_id + "_desc-preproc_T1w.nii.gz"),
    #threads:config['fastsurfer']['threads']
    shell:
        "singularity run -e {params.hippunfold_container} {params.in_dir} {params.hippunfold_out} participant --force_output --participant_label {params.participant_label}\
        --modality {params.modality} --cores 4"
    
final_outputs.extend(expand(rules.hippunfold_seg.output.touch_hippunfold, subject=subjects))

