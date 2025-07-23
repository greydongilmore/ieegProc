
if config['meld']['run']:
	rule import_meld:
		input: get_meld_filename,
		output: bids(root=join(config['out_dir'],'derivatives', 'meld_in','input'),subject=subject_id, datatype=config['meld_vol']['datatype'], suffix=config['meld_vol']['suffix']+config['meld_vol']['ext'], include_session_dir=False)
		group: 'preproc'
		shell: 'cp {input} {output}'
	rule import_meld_params:
		input: 
			subjid=subject_id,
			demographics_file = config['meld_config']['demographics_file'],
			sid = config['fastsurfer_config']['sid'],
			batch = config['fastsurfer_config']['batch'],
		output:
			demographics_out: bids(root=join(config['out_dir'],'derivatives', 'meld_in','input'),subject=subject_id, datatype=config['meld_vol']['datatype'], suffix=config['meld_vol']['suffix']+config['meld_vol']['ext'], include_session_dir=False)
	

	rule fastsurfer_seg:
		input: 
			t1 = rules.import_meld.output,
		params:
			fastsurfer_run = config['fastsurfer_config']['home'],
			sid = config['fastsurfer_config']['sid'],
			batch = config['fastsurfer_config']['batch'],
			threads = config['fastsurfer_config']['threads'],
			vox_size = config['fastsurfer_config']['vox_size'],
			py = config['fastsurfer_config']['py'],
			fastsurfer_out = directory(join(config['out_dir'], 'derivatives', 'fastsurfer')),
			subjid=subject_id,
		output:
			touch_fastsurfer=touch(join(config['out_dir'], 'logs', 'sub-' + subject_id + "_fastsurfer.done")),
			t1_fname = join(config['out_dir'],'derivatives','fastsurfer','sub-' + subject_id, 'mri','orig.mgz'),
			segs = join(config['out_dir'],'derivatives','fastsurfer','sub-' + subject_id, 'mri','aparc+aseg.mgz'),
		group: 'preproc'
		threads: 8
		shell:
			"export FASTSURFER_HOME={params.fastsurfer_run} &&PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:4096 {params.fastsurfer_run}/run_fastsurfer.sh \
			--t1 {input.t1} --sd {params.fastsurfer_out} --threads {params.threads} --vox_size {params.vox_size} --sid sub-{params.subjid} --py {params.py} --viewagg_device cpu --fsaparc --parallel --allow_root"