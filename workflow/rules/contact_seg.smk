if config['contact_seg']['run']:
	rule import_contact_seg:
		input: 
			fcsv_in= get_electrodes_coords(subject_id, coords_type='actual'),
			xfm_ras = bids(root=join(config['out_dir'],'derivatives', 'atlasreg'),subject=subject_id,suffix='ses-'+config['post_image']['session']+'_xfm.txt',from_=config['post_image']['datatype'],to=config['contrast_t1']['suffix'],desc='rigid',type_='ras'),
			post_img=rules.apply_noninterp_transform_post.output.warped_subj,
			pre_img=get_pre_t1_filename,
		params:
			contactseg_out=directory(join(config['out_dir'], 'derivatives', 'contactseg_out')),
		output:
			fcsv_out= bids(root=join(config['out_dir'],'derivatives', 'contactseg_in'),subject=subject_id, run=config['post_image']['run'], suffix='planned.fcsv', include_session_dir=False),
			xfm__out= bids(root=join(config['out_dir'],'derivatives', 'contactseg_in'),subject=subject_id, run=config['post_image']['run'], suffix='xfm.txt', include_session_dir=False),
			img_out=bids(root=join(config['out_dir'],'derivatives', 'contactseg_in'),subject=subject_id, datatype=config['post_image']['datatype'], session=config['post_image']['session'], acq=config['post_image']['acq'], run=config['post_image']['run'], suffix=config['post_image']['suffix']+config['post_image']['ext'], include_session_dir=True),
			pre_out=bids(root=join(config['out_dir'], 'derivatives', 'contactseg_in'), subject=subject_id, datatype=config['contrast_t1']['datatype'], session=config['contrast_t1']['session'], acq=config['contrast_t1']['acq'], run='02', suffix=config['contrast_t1']['suffix']+config['contrast_t1']['ext'], include_session_dir=True),
		run:
			import shutil,os
			shutil.copy2(input.post_img, output.img_out)
			shutil.copy2(input.pre_img, output.pre_out)
			if not os.path.exists(output.fcsv_out):
				shutil.copy2(input.fcsv_in, output.fcsv_out)
			if not os.path.exists(output.xfm__out):
				shutil.copy2(input.xfm_ras, output.xfm__out)
	
	rule contact_seg_run:
		input:
			dir_in=dirname(dirname(rules.import_contact_seg.output.fcsv_out)),
			reg_in=rules.import_contact_seg.output.xfm__out
		params:
			contact_seg_home=config['contact_seg_config']['home'],
			contact_seg_cmd=config['contact_seg_config']['cmd'],
			subjid=subject_id,
		output:
			touch_contact_seg=touch(join(config['out_dir'], 'logs', 'sub-' + subject_id + "_contact_seg.done")),
		group: 'preproc'
		threads: 4
		shell:
			'cd {params.contact_seg_home}&&{params.contact_seg_cmd}&&'
			'contactseg {input.dir_in} {rules.import_contact_seg.params.contactseg_out} participant --cores 4 --use_gpu --label --transform --skip-bids-validation  --participant-label {params.subjid}'

	final_outputs.extend(expand(bids(root=join(config['out_dir'],'derivatives', 'contactseg_in'),subject=subject_id, run=config['post_image']['run'], suffix='planned.fcsv', include_session_dir=False),subject=subjects))
	final_outputs.extend(expand(rules.contact_seg_run.output.touch_contact_seg, subject=subjects))
