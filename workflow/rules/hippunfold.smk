
rule hippunfold_prep:
    input:
        in_t1w= bids(root=join(config['out_dir'],'derivatives', 'seega_scenes'),subject=subject_id,acq='noncontrast',suffix='T1w.nii.gz'),
    params:
        in_dir = directory(join(config['out_dir'], 'derivatives', 'hippunfold_in')),
        in_dir_sub = directory(join(config['out_dir'], 'derivatives', 'hippunfold_in','sub-'+subject_id)),
    output:
        out_t1w=bids(root=join(config['out_dir'], 'derivatives','hippunfold_in'), subject=subject_id, datatype='anat', suffix='T1w.nii.gz')
    shell:
        "mkdir -p {params.in_dir} && mkdir -p {params.in_dir_sub} && mkdir -p {params.in_dir_sub}/anat && cp {input.in_t1w} {output.out_t1w}"

rule hippunfold_seg:
    input: 
        in_dir = directory(join(config['out_dir'], 'derivatives', 'hippunfold_in')),
    params:
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
        "singularity run -e {params.hippunfold_container} {input.in_dir} {params.hippunfold_out} participant --force_output --participant_label {params.participant_label}\
        --modality {params.modality} --cores 4"
    
final_outputs.extend(expand(rules.hippunfold_seg.output.touch_hippunfold, subject=subjects))

