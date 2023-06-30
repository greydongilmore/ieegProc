
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

def get_postop_filename(wildcards):
    files=glob(bids(root=join(config['out_dir'], 'bids'), subject=config['subject_prefix']+f'{wildcards.subject}', datatype='ct', session='post', acq='Electrode', run='*', suffix='ct.nii.gz'))
    if len(files) <=1:
        file=expand(bids(root=join(config['out_dir'], 'bids'), subject=config['subject_prefix']+'{subject}', datatype='ct', session='post', acq='Electrode', run='01', suffix='ct.nii.gz'),subject=wildcards.subject)
    else:
        files.sort(key=lambda f: int(re.sub('\D', '', f)))
        file=files[config['post_ct']['position']]
    return file

def get_pet_filename(wildcards):
    files=glob(bids(root=join(config['out_dir'], 'bids'), subject=config['subject_prefix']+f'{wildcards.subject}', datatype='pet', session='pre', task=config['pet']['task'], run='*', suffix='pet.nii.gz'))
    if len(files) <=1:
        file=expand(bids(root=join(config['out_dir'], 'bids'), subject=config['subject_prefix']+'{subject}', datatype='pet', session='pre', task=config['pet']['task'], run='01', suffix='pet.nii.gz'),subject=wildcards.subject)
    else:
        files.sort(key=lambda f: int(re.sub('\D', '', f)))
        file=files[config['pet']['position']]
    return file

def get_reference_t1(wildcards):
    if config['contrast_t1']['present']:
        ref_file=expand(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'), subject=config['subject_prefix']+'{subject}', acq='contrast', suffix='T1w.nii.gz'),subject=wildcards.subject)
    elif not config['contrast_t1']['present'] and config['noncontrast_t1']['present']:
        ref_file=expand(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'), subject=config['subject_prefix']+'{subject}', acq='noncontrast', suffix='T1w.nii.gz'),subject=wildcards.subject)
    return ref_file

if config['contrast_t1']['present']:
    rule import_subj_t1:
        input: get_pre_t1_filename,
        output: bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,acq='contrast',suffix='T1w.nii.gz'),
        group: 'preproc'
        shell: "echo {input} &&cp {input} {output}"

if config['noncontrast_t1']['present']:
    rule import_subj_t1_noncontrast:
        input: get_noncontrast_filename,
        output: bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,acq='noncontrast',suffix='T1w.nii.gz'),
        group: 'preproc'
        shell: "echo {input} &&cp {input} {output}"

if config['contrast_t1']['present'] and not config['noncontrast_t1']['present']:
    rule cp_contrast_to_T1w:
        input: bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,acq='contrast',suffix='T1w.nii.gz'),
        output: bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='T1w.nii.gz'),
        shell:
            'cp {input} {output}'

elif not config['contrast_t1']['present'] and config['noncontrast_t1']['present']:
    rule cp_contrast_to_T1w:
        input: bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,acq='noncontrast',suffix='T1w.nii.gz'),
        output: bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='T1w.nii.gz'),
        shell:
            'cp {input} {output}'

elif config['contrast_t1']['present'] and config['noncontrast_t1']['present']:
    if config['noncontrast_t1']['algo'] =='reg_aladin':
        rule rigonly_aladin_contrast:
            input: 
                flo = bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,acq='noncontrast',suffix='T1w.nii.gz'),
                ref = bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,acq='contrast',suffix='T1w.nii.gz'),
            params:
                dof=config['linear_reg']['reg_aladin']['dof']
            output: 
                warped_subj = bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,acq='noncontrast',suffix='T1w.nii.gz',space='T1w',desc='rigidInterp'),
                xfm_ras = bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='xfm.txt',from_='noncontrast',to='contrast',desc='rigid',type_='ras'),
            #container: config['singularity']['neuroglia']
            group: 'preproc'
            shell:
                'reg_aladin -flo {input.flo} -ref {input.ref} {params.dof} -interp 0 -res {output.warped_subj} -aff {output.xfm_ras} -speeeeed'
                #'flirt -in {input.flo} -ref {input.ref} -out {output.warped_subj} -omat {output.xfm_ras} -dof 6'
    elif config['noncontrast_t1']['algo'] =='greedy':
        rule rigonly_greedy_contrast:
            input: 
                flo = bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,acq='noncontrast',suffix='T1w.nii.gz'),
                ref = bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,acq='contrast',suffix='T1w.nii.gz'),
            params:
                n_iterations_affine=config['linear_reg']['greedy']['n_iterations_affine'],
                dof=config['linear_reg']['greedy']['dof'],
            output: 
                warped_subj = bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,acq='noncontrast',suffix='T1w.nii.gz',space='T1w',desc='rigidInterp'),
                xfm_ras = bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='xfm.txt',from_='noncontrast',to='contrast',desc='rigid',type_='ras'),
            group: 'preproc'
            shell:
                'greedy -d 3 -threads 4 -a -ia-image-centers -m MI -dof {params.dof} -i {input.ref} {input.flo} -o {output.xfm_ras} -n {params.n_iterations_affine} &&'
                'greedy -d 3 -threads 4 -rf {input.ref} -rm {input.flo} {output.warped_subj} -r {output.xfm_ras}'

    rule apply_noninterp_transform_noncontrast:
        input:
            xfm=bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='xfm.txt',from_='noncontrast',to='contrast',desc='rigid',type_='ras'),
            flo = bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,acq='noncontrast',suffix='T1w.nii.gz'),
        output:
            warped_subj = bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,acq='noncontrast',suffix='T1w.nii.gz',space='T1w',desc='rigid'),
        group: 'preproc'
        script: 
            '../scripts/apply_transform_noninterp.py'

    rule convert_T1w_xfm_tfm:
        input:
            xfm=bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='xfm.txt',from_='noncontrast',to='contrast',desc='rigid',type_='ras'),
        output:
            tfm=bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='xfm.tfm',from_='noncontrast',to='contrast',desc='rigid',type_='ras'),
        group: 'preproc'
        script: 
            '../scripts/convert_xfm_tfm.py'

    rule cp_noncontrast_to_T1w:
        input: bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,acq='noncontrast',suffix='T1w.nii.gz',space='T1w',desc='rigid'),
        output: bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='T1w.nii.gz'),
        shell:
            'cp {input} {output}'

    final_outputs.extend(expand(bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='xfm.tfm',from_='noncontrast',to='contrast',desc='rigid',type_='ras'),
                        subject=subjects))

if config['post_ct']['present']:
    rule import_subj_ct:
        input: get_postop_filename,
        output: bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='ct.nii.gz')
        group: 'preproc'
        shell: 'cp {input} {output}'

    if config['post_ct']['algo'] =='reg_aladin':
        rule rigonly_aladin_ct:
            input: 
                flo = bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='ct.nii.gz'),
                ref = get_reference_t1,
            params:
                dof=config['linear_reg']['reg_aladin']['dof']
            output: 
                warped_subj = bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='ct.nii.gz',space='T1w',desc='rigidInterp'),
                xfm_ras = bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='xfm.txt',from_='ct',to='T1w',desc='rigid',type_='ras'),
            #container: config['singularity']['neuroglia']
            group: 'preproc'
            shell:
                'reg_aladin -flo {input.flo} -ref {input.ref} {params.dof} -interp 0 -res {output.warped_subj} -aff {output.xfm_ras} -speeeeed'
                #'flirt -in {input.flo} -ref {input.ref} -out {output.warped_subj} -omat {output.xfm_ras} -dof 6'
    elif config['post_ct']['algo'] =='greedy':
        rule rigonly_greedy_ct:
            input: 
                flo = bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='ct.nii.gz'),
                ref = get_reference_t1,
            params:
                n_iterations_affine=config['linear_reg']['greedy']['n_iterations_affine'],
                dof=config['linear_reg']['greedy']['dof'],
            output: 
                warped_subj = bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='ct.nii.gz',space='T1w',desc='rigidInterp'),
                xfm_ras = bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='xfm.txt',from_='ct',to='T1w',desc='rigid',type_='ras'),
            group: 'preproc'
            shell:
                'greedy -d 3 -threads 4 -a -ia-image-centers -m MI -dof {params.dof} -i {input.ref} {input.flo} -o {output.xfm_ras} -n {params.n_iterations_affine} &&'
                'greedy -d 3 -threads 4 -rf {input.ref} -rm {input.flo} {output.warped_subj} -r {output.xfm_ras}'

    rule apply_noninterp_transform_ct:
        input:
            xfm=bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='xfm.txt',from_='ct',to='T1w',desc='rigid',type_='ras'),
            flo = bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='ct.nii.gz'),
        output:
            warped_subj = bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='ct.nii.gz',space='T1w',desc='rigid'),
        group: 'preproc'
        script: 
            '../scripts/apply_transform_noninterp.py'

    rule convert_ct_xfm_tfm:
        input:
            xfm=bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='xfm.txt',from_='ct',to='T1w',desc='rigid',type_='ras'),
        output:
            tfm=bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='xfm.tfm',from_='ct',to='T1w',desc='rigid',type_='ras'),
        group: 'preproc'
        script: 
            '../scripts/convert_xfm_tfm.py'

    final_outputs.extend(expand(bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='xfm.tfm',from_='ct',to='T1w',desc='rigid',type_='ras'),
                        subject=subjects))

if config['pet']['present']:
    rule import_subj_pet:
        input: get_pet_filename,
        output: bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='pet.nii.gz')
        group: 'preproc'
        shell: 'cp {input} {output}'

    rule rigonly_aladin_pet:
        input: 
            flo = bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='pet.nii.gz'),
            ref = get_reference_t1,
        params:
            dof=config['linear_reg']['reg_aladin']['dof']
        output: 
            warped_subj = bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='pet.nii.gz',space='T1w',desc='rigidInterp'),
            xfm_ras = bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='xfm.txt',from_='pet',to='T1w',desc='rigid',type_='ras'),
        #container: config['singularity']['neuroglia']
        group: 'preproc'
        shell:
            'reg_aladin -flo {input.flo} -ref {input.ref} {params.dof} -res {output.warped_subj} -aff {output.xfm_ras} -speeeeed'
            #'flirt -in {input.flo} -ref {input.ref} -out {output.warped_subj} -omat {output.xfm_ras} -dof 6'

    rule apply_noninterp_transform_pet:
        input:
            xfm=bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='xfm.txt',from_='pet',to='T1w',desc='rigid',type_='ras'),
            flo = bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='pet.nii.gz'),
        output:
            warped_subj = bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='pet.nii.gz',space='T1w',desc='rigid'),
        group: 'preproc'
        script: 
            '../scripts/apply_transform_noninterp.py'

    rule convert_pet_xfm_tfm:
        input:
            xfm=bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='xfm.txt',from_='pet',to='T1w',desc='rigid',type_='ras'),
        output:
            tfm=bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='xfm.tfm',from_='pet',to='T1w',desc='rigid',type_='ras'),
        group: 'preproc'
        script: 
            '../scripts/convert_xfm_tfm.py'

    final_outputs.extend(expand(bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='xfm.tfm',from_='pet',to='T1w',desc='rigid',type_='ras'),
                        subject=subjects))

if config['affine_reg']['algo']=='reg_aladin':
    rule affine_xfm:
        input: 
            flo = bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='T1w.nii.gz'),
            ref = get_age_appropriate_template_name(expand(subject_id,subject=subjects),'t1w'),
        output: 
            warped_subj = bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='T1w.nii.gz',space=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'),desc='affine'),
            affine_xfm_ras = bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='xfm.txt',from_='subject',to=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'),desc='affine',type_='ras'),
        #container: config['singularity']['neuroglia']
        group: 'preproc'
        shell:
            'reg_aladin -flo {input.flo} -ref {input.ref} -res {output.warped_subj} -aff {output.affine_xfm_ras} -speeeeed'

elif config['affine_reg']['algo']=='greedy':
    rule affine_xfm:
        input: 
            flo = bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='T1w.nii.gz'),
            ref = get_age_appropriate_template_name(expand(subject_id,subject=subjects), 't1w'),
        output: 
            affine_xfm_ras = bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='xfm.txt',from_='subject',to=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'),desc='affine',type_='ras'),
            xfm_deform = bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='xfm.nii.gz',from_='subject',to=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'),desc='warp',type_='ras'),
            xfm_deform_inv = bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='xfm.nii.gz',from_='subject',to=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'),desc='warpInverse',type_='ras'),
            warped_subj_affine = bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='T1w.nii.gz',space=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'),desc='affine'),
            warped_subj_greedy = bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='T1w.nii.gz',space=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'),desc='greedydeform'),
        params:
            n_iterations_affine=config['affine_reg']['greedy']['n_iterations_affine'],
            n_iterations_deform=config['affine_reg']['greedy']['n_iterations_deform'],
            grad_sigma=config['affine_reg']['greedy']['grad_sigma'],
            warp_sigma=config['affine_reg']['greedy']['warp_sigma'],
        #container: config['singularity']['neuroglia']
        group: 'preproc'
        shell:
            'greedy -d 3 -threads 4 -a -ia-image-centers -m MI -i {input.ref} {input.flo} -o {output.affine_xfm_ras} -n {params.n_iterations_affine}&&'
            'greedy -d 3 -threads 4 -m MI -i {input.ref} {input.flo} -it {output.affine_xfm_ras}  -o {output.xfm_deform} -oinv {output.xfm_deform_inv} -n {params.n_iterations_deform} -s {params.grad_sigma} {params.warp_sigma} &&'
            'greedy -d 3 -threads 4 -rf {input.ref} -rm {input.flo} {output.warped_subj_affine} -r {output.affine_xfm_ras}&&'
            'greedy -d 3 -threads 4 -rf {input.ref} -rm {input.flo} {output.warped_subj_greedy} -r {output.xfm_deform} {output.affine_xfm_ras}'

rule convert_affine_xfm_tfm:
    input:
        xfm=bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='xfm.txt',from_='subject',to=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'),desc='affine',type_='ras'),
    output:
        tfm=bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='xfm.tfm',from_='subject',to=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'),desc='affine',type_='ras'),
    group: 'preproc'
    script: 
        '../scripts/convert_xfm_tfm.py'

final_outputs.extend(expand(bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='xfm.tfm',from_='subject',to=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'),desc='affine',type_='ras'),
                        subject=subjects))

rule convert_xfm_ras2itk:
    input:
        tfm=bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='xfm.txt',from_='subject',to=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'),desc='affine',type_='ras'),
        xfm=bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='xfm.txt',from_='subject',to=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'),desc='{desc}',type_='ras'),
    output:
        bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='xfm.txt',from_='subject',to=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'),desc='{desc}',type_='itk'),
    #container: config['singularity']['neuroglia']
    group: 'preproc'
    shell:
        'c3d_affine_tool {input.xfm}  -oitk {output}'

rule warp_brainmask_from_template_affine:
    input: 
        mask = get_age_appropriate_template_name(expand(subject_id,subject=subjects),'mask'),
        ref = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='T1w.nii.gz'),
        xfm = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='xfm.txt',from_='subject',to=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'),desc='{desc}',type_='itk'),
    output:
        mask = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='mask.nii.gz',from_=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'),reg='{desc}',desc='brain'),
    #container: config['singularity']['neuroglia']
    group: 'preproc'
    shell: 'antsApplyTransforms -d 3 --interpolation NearestNeighbor -i {input.mask} -o {output.mask} -r {input.ref} '
            ' -t [{input.xfm},1] ' #use inverse xfm (going from template to subject)

rule warp_tissue_probseg_from_template_affine:
    input: 
        probseg = get_age_appropriate_template_name(expand(subject_id,subject=subjects),'tissue_probseg'),
        ref = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='T1w.nii.gz'),
        xfm = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='xfm.txt',from_='subject',to=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'),desc='{desc}',type_='itk'),
    output:
        probseg = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='probseg.nii.gz',label='{tissue}',from_=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'),reg='{desc}'),
    #container: config['singularity']['neuroglia']
    group: 'preproc'
    threads: 1
    resources:
        mem_mb = 16000
    shell: 
        'ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} '
        'antsApplyTransforms -d 3 --interpolation Linear -i {input.probseg} -o {output.probseg} -r {input.ref} '
            ' -t [{input.xfm},1]' #use inverse xfm (going from template to subject)

rule n4biasfield:
    input: 
        t1 = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='T1w.nii.gz'),
    output:
        t1 = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,desc='n4', suffix='T1w.nii.gz'),
    threads: 8
    #container: config['singularity']['neuroglia']
    group: 'preproc'
    shell:
        'ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} '
        'N4BiasFieldCorrection -d 3 -i {input.t1} -o {output}'
        #'N4BiasFieldCorrection -d 3 -i {input.t1} -x {input.mask} -o {output}'

rule mask_template_t1w:
    input:
        t1 = get_age_appropriate_template_name(expand(subject_id,subject=subjects),'t1w'),
        mask = get_age_appropriate_template_name(expand(subject_id,subject=subjects),'mask'),
    output:
        t1 = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix=lambda wildcards, input: f"tpl-{get_age_appropriate_template_name(input,'space')}/tpl-{get_age_appropriate_template_name(input,'space')}",desc='masked',suffix='T1w.nii.gz')
    #container: config['singularity']['neuroglia']
    group: 'preproc'
    shell:
        'fslmaths {input.t1} -mas {input.mask} {output}'

if config['deface']:
    rule deface_subject:
        input: 
            t1=bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,desc='n4', suffix='T1w.nii.gz'),
            ct=bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='ct.nii.gz',space='T1w',desc='rigid'),
            mask = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='mask.nii.gz',from_='atropos3seg',desc='brain'),
        params:
            template=config['mean_reg2mean'],
            facemask=config['facemask']
        output:
            t1_masked = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='T1w.nii.gz',from_='atropos3seg',desc='masked'),
            template_to_t1=bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='template.nii.gz',space='subject',desc='affine'),
            template_to_t1_matrix = bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='xfm.txt',from_='template',desc='affine',type_='ras'),
            template_to_t1_itk = bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='xfm.txt',from_='template',desc='affine',type_='itk'),
            warped_mask=bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='facemask.nii.gz',space='subject',desc='affine'),
            #warped_mask_matrix = bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='xfm.txt',from_='facemask',to='subject',desc='affine',type_='ras'),
            #warped_mask_itk = bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='xfm.txt',from_='facemask',to='subject',desc='affine',type_='itk'),
            defaced_t1=bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='T1w.nii.gz',desc='defaced'),
            defaced_ct=bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='ct.nii.gz',desc='defaced',space='T1w'),
        group: 'preproc'
        #shell: 'cp {input} {output.t1w} &&\
        #        flirt -in {params.template} -ref {input} -out {output.template_reg} -omat {output.out_matrix} -cost mutualinfo &&\
        #        flirt -in {params.facemask} -ref {input} -init {output.template_reg} -applyxfm -out {output.warped_mask} -omat {output.out_mask_matrix} &&\
        #        fslmaths {input} -mas {output.warped_mask} {output.deface}'

        shell: 'fslmaths {input.t1} -mas {input.mask} {output.t1_masked} &&\
                reg_aladin -flo {params.template} -ref {input.t1} -res {output.template_to_t1} -aff {output.template_to_t1_matrix} -interp 0 -speeeeed &&\
                c3d_affine_tool {output.template_to_t1_matrix}  -oitk {output.template_to_t1_itk} &&\
                antsApplyTransforms -d 3 -i {params.facemask} -o {output.warped_mask} -r {input.t1} -t [{output.template_to_t1_itk},0] &&\
                fslmaths {input.t1} -mas {output.warped_mask} {output.defaced_t1} &&\
                fslmaths {input.ct} -mas {output.warped_mask} {output.defaced_ct}'
else:
    rule mask_subject_t1w:
        input:
            t1 = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,desc='n4', suffix='T1w.nii.gz'),
            mask = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='mask.nii.gz',from_='atropos3seg',desc='brain')
        output:
            t1 = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='T1w.nii.gz',from_='{atropos3seg}',desc='masked'),
        #container: config['singularity']['neuroglia']
        group: 'preproc'
        shell:
            'fslmaths {input.t1} -mas {input.mask} {output}'

if config['post_ct']['present']:
    rule mask_subject_ct:
        input:
            ct = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,desc='rigid',space='T1w', suffix='ct.nii.gz'),
            mask = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='mask.nii.gz',from_='atropos3seg',desc='brain')
        output:
            ct = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='ct.nii.gz',from_='atropos3seg',desc='masked'),
        #container: config['singularity']['neuroglia']
        group: 'preproc'
        shell:
            'fslmaths {input.ct} -mas {input.mask} {output.ct}'

if  config['warp_reg']['algo']=='ants':
    rule ants_syn_affine_init:
        input: 
            flo = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='T1w.nii.gz',from_='atropos3seg',desc='masked'),
            ref = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix=lambda wildcards, input:  f"tpl-{get_age_appropriate_template_name(input,'space')}/tpl-{get_age_appropriate_template_name(input,'space')}",desc='masked',suffix='T1w.nii.gz'),
            init_xfm = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='xfm.txt',from_='subject',to=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'),desc='affine',type_='itk'),
        params:
            out_prefix = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='',from_='subject',to=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space')),
            base_opts = '--write-composite-transform -d {dim} --float 1 '.format(dim=config['warp_reg']['ants']['dim']),
            intensity_opts = config['warp_reg']['ants']['intensity_opts'],
            init_transform = lambda wildcards, input: '-r {xfm}'.format(xfm=input.init_xfm),
            linear_multires = '-c [{reg_iterations},1e-6,10] -f {shrink_factors} -s {smoothing_factors}'.format(
                                    reg_iterations = config['warp_reg']['ants']['linear']['reg_iterations'],
                                    shrink_factors = config['warp_reg']['ants']['linear']['shrink_factors'],
                                    smoothing_factors = config['warp_reg']['ants']['linear']['smoothing_factors']),
            linear_metric = lambda wildcards, input: '-m MI[{template},{target},1,32,Regular,0.25]'.format( template=input.ref,target=input.flo),
            deform_model = '-t {deform_model}'.format(deform_model = config['warp_reg']['ants']['deform']['transform_model']),
            deform_multires = '-c [{reg_iterations},1e-9,10] -f {shrink_factors} -s {smoothing_factors}'.format(
                                    reg_iterations = config['warp_reg']['ants']['deform']['reg_iterations'],
                                    shrink_factors = config['warp_reg']['ants']['deform']['shrink_factors'],
                                    smoothing_factors = config['warp_reg']['ants']['deform']['smoothing_factors']),
            deform_metric = lambda wildcards, input: '-m {metric}[{template},{target},1,4]'.format(
                                    metric=config['warp_reg']['ants']['deform']['sim_metric'],
                                    template=input.ref, target=input.flo)
        output:
            out_composite = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='Composite.h5',from_='subject',to=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space')),
            out_inv_composite = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='InverseComposite.h5',from_='subject',to=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space')),
            warped_flo = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='T1w.nii.gz',space=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'),desc='SyN',label='brain'),
        threads: 8
        resources:
            mem_mb = 16000, # right now these are on the high-end -- could implement benchmark rules to do this at some point..
            time = 60 # 1 hrs
        #container: config['singularity']['neuroglia']
        group: 'preproc'
        shell: 
            'ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} '
            'antsRegistration {params.base_opts} {params.intensity_opts} '
            '{params.init_transform} ' #initial xfm  -- rely on this for affine
        #    '-t Rigid[0.1] {params.linear_metric} {params.linear_multires} ' # rigid registration
        #    '-t Affine[0.1] {params.linear_metric} {params.linear_multires} ' # affine registration
            '{params.deform_model} {params.deform_metric} {params.deform_multires} '  # deformable registration
            '-o [{params.out_prefix},{output.warped_flo}]'
elif  config['warp_reg']['algo']=='greedy':
    rule greedy_syn_affine_init:
        input: 
            flo = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='T1w.nii.gz',from_='atropos3seg',desc='masked'),
            ref = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix=lambda wildcards, input:  f"tpl-{get_age_appropriate_template_name(input,'space')}/tpl-{get_age_appropriate_template_name(input,'space')}",desc='masked',suffix='T1w.nii.gz'),
            init_xfm = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='xfm.txt',from_='subject',to=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'),desc='affine',type_='itk'),

rule warp_dseg_from_template:
    input: 
        dseg = get_age_appropriate_template_name(expand(subject_id,subject=subjects),'atlas_dseg_nii'),
        ref = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'), subject=subject_id,suffix='T1w.nii.gz'),
        inv_composite = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='InverseComposite.h5',from_='subject',to=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space')),
    output:
        dseg = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='dseg.nii.gz',atlas='{atlas}',from_=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'),reg='SyN'),
    #container: config['singularity']['neuroglia']
    group: 'preproc'
    threads: 1
    resources:
        mem_mb = 16000
    shell: 
        'ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} '
        'antsApplyTransforms -d 3 --interpolation NearestNeighbor -i {input.dseg} -o {output.dseg} -r {input.ref} '
            ' -t {input.inv_composite} ' #use inverse xfm (going from template to subject)

rule warp_t1w_to_template:
    input: 
        t1 = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='T1w.nii.gz'),
        ref = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),prefix=lambda wildcards, input:  f"tpl-{get_age_appropriate_template_name(input,'space')}/tpl-{get_age_appropriate_template_name(input,'space')}",desc='masked',suffix='T1w.nii.gz'),
        composite = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='Composite.h5',from_='subject',to=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space')),
    output:
        warped_t1 = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='T1w.nii.gz',space=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'),desc='SyN'),
    #container: config['singularity']['neuroglia']
    group: 'preproc'
    threads: 1
    resources:
        mem_mb = 16000
    shell: 
        'ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} '
        'antsApplyTransforms -d 3 --interpolation Linear -i {input.t1} -o {output.warped_t1} -r {input.ref} '
            ' -t {input.composite} ' #use inverse xfm (going from template to subject)

rule warp_tissue_probseg_from_template:
    input: 
        probseg = get_age_appropriate_template_name(expand(subject_id,subject=subjects),'tissue_probseg'),
        ref = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='T1w.nii.gz'),
        inv_composite = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='InverseComposite.h5',from_='subject',to=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space')),
    output:
        probseg = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='probseg.nii.gz',label='{tissue}',from_=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'),reg='SyN'),
    #container: config['singularity']['neuroglia']
    group: 'preproc'
    threads: 1
    resources:
        mem_mb = 16000
    shell: 
        'ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} '
        'antsApplyTransforms -d 3 --interpolation Linear -i {input.probseg} -o {output.probseg} -r {input.ref} '
            ' -t {input.inv_composite} ' #use inverse xfm (going from template to subject)

rule warp_brainmask_from_template:
    input: 
        mask = get_age_appropriate_template_name(expand(subject_id,subject=subjects),'mask'),
        ref = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='T1w.nii.gz'),
        inv_composite = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='InverseComposite.h5',from_='subject',to=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space')),
    output:
        mask = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='mask.nii.gz',from_=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'),reg='SyN',label='brain'),
    #container: config['singularity']['neuroglia']
    group: 'preproc'
    threads: 1
    resources:
        mem_mb = 16000
    shell: 
        'ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} '
        'antsApplyTransforms -d 3 --interpolation NearestNeighbor -i {input.mask} -o {output.mask} -r {input.ref} '
            ' -t {input.inv_composite} ' #use inverse xfm (going from template to subject)

rule dilate_brainmask:
    input:
        mask = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='mask.nii.gz',from_=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'),reg='{desc}',desc='brain'),
    output:
        mask = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='mask.nii.gz',from_=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'),reg='{desc}',desc='braindilated'),
    #container: config['singularity']['neuroglia']
    group: 'preproc'
    shell:
        'fslmaths {input} -dilD {output}'

#dilate labels N times to provide more of a fudge factor when assigning GM labels
rule dilate_atlas_labels:
    input:
        dseg = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='dseg.nii.gz',atlas='{atlas}',from_=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'),reg='SyN'),
    params:
        dil_opt =  ' '.join([ '-dilD' for i in range(config['n_atlas_dilate'])])
    output:
        dseg = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='dseg.nii.gz',atlas='{atlas}',from_=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'),desc='dilated',reg='SyN'),
    #container: config['singularity']['neuroglia']
    group: 'preproc'
    shell:
        'fslmaths {input} {params.dil_opt} {output}'

final_outputs.extend(expand(bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'),subject=subject_id,suffix='T1w.nii.gz',space=get_age_appropriate_template_name(expand(subject_id,subject=subjects),'space'),desc='SyN',label='brain'),subject=subjects))