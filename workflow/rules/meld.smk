
if config['meld']['run']:
	rule import_meld:
		input: 
			img_in=get_meld_filename,
			demographics_file = config['meld_config']['demographics_file'],
			meld_params = config['meld_config']['meld_params'],
			meld_models = config['meld_config']['meld_models'],
			freesurfer_lic = config['meld_config']['freesurfer_license'],
			dataset_description = config['meld_config']['dataset_description'],
			meld_bids_config = config['meld_config']['meld_bids_config'],
		params:
			subjid='sub-' +subject_id,
			demographics_out= join(config['out_dir'],'derivatives', 'meld_in','demographics_file.csv'),
			meld_params_out=join(config['out_dir'],'derivatives', 'meld_in', 'meld_params'),
			meld_models_out=join(config['out_dir'],'derivatives', 'meld_in', 'models'),
			freesurfer_lic_out=join(config['out_dir'],'derivatives', 'meld_in', 'license.txt'),
			dataset_description_out=join(config['out_dir'],'derivatives', 'meld_in', 'input','dataset_description.json'),
			meld_bids_config_out=join(config['out_dir'],'derivatives', 'meld_in', 'input', 'meld_bids_config.json'),
		output:
			img_out=bids(root=join(config['out_dir'],'derivatives', 'meld_in','input'),subject=subject_id, datatype=config['meld_vol']['datatype'], suffix=config['meld_vol']['suffix']+config['meld_vol']['ext'], include_session_dir=False),
		run:
			import pandas as pd
			import shutil,os
			df = pd.read_table(input.demographics_file,sep=',',header=0)
			df['ID']=params.subjid
			df['Harmo code (harmonisation code, put “noHarmo” if not using harmonisation)']="noHarmo"
			df['Group (“patient” or “control”)']="patient"
			df['Age at preoperative (in years)']=''
			df['Sex (“female” or “male”)']=''
			df['Scanner (“3T” for 3Tesla or “15T” for 1.5T)']=''
			df.to_csv(params.demographics_out, sep=',', index=False)
			shutil.copy2(input.img_in, output.img_out)
			if not os.path.exists(params.meld_params_out):
				shutil.copytree(input.meld_params, params.meld_params_out)
			if not os.path.exists(params.meld_models_out):
				shutil.copytree(input.meld_models, params.meld_models_out)
			if not os.path.exists(params.freesurfer_lic_out):
				shutil.copy2(input.freesurfer_lic, params.freesurfer_lic_out)
			if not os.path.exists(params.dataset_description_out):
				shutil.copy2(input.dataset_description, params.dataset_description_out)
			if not os.path.exists(params.meld_bids_config_out):
				shutil.copy2(input.meld_bids_config, params.meld_bids_config_out)
	rule meld_run:
		input:
			img_in=rules.import_meld.output.img_out,
		params:
			subjid='sub-' +subject_id,
			meld_in=directory(join(config['out_dir'],'derivatives', 'meld_in')),
			meld_container= config['singularity']['meld'],
			freesurfer_lic = rules.import_meld.params.freesurfer_lic_out,
			demographics_file = rules.import_meld.params.demographics_out,
		output:
			touch_meld=touch(join(config['out_dir'], 'logs', 'sub-' + subject_id + "_meld.done")),
		group: 'preproc'
		threads: 8
		shell:
			'export SINGULARITY_BINDPATH={params.meld_in}:/data,{params.freesurfer_lic}:/license.txt:ro&&export SINGULARITYENV_FS_LICENSE=/license.txt&&'
			'singularity exec {params.meld_container} /bin/bash -c "cd /app && source \$FREESURFER_HOME/FreeSurferEnv.sh && python scripts/new_patient_pipeline/new_pt_pipeline.py -id {params.subjid} -demos {params.demographics_file} --fastsurfer --parallelise"'

	final_outputs.extend(expand(bids(root=join(config['out_dir'],'derivatives', 'meld_in','input'),subject=subject_id, datatype=config['meld_vol']['datatype'], suffix=config['meld_vol']['suffix']+config['meld_vol']['ext'], include_session_dir=False),subject=subjects))
	final_outputs.extend(expand(rules.meld_run.output.touch_meld, subject=subjects))
