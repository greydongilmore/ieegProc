
rule pet_asymmetry_copy:
    input:
        t1 = get_pre_t1_filename,
        pet=get_pet1_filename,
        flair=get_flair_filename,
    params:
        in_json= config['pet_asymmetry_config']['in_json'],
        in_dir= config['pet_asymmetry_config']['home'],
        in_cmd= config['pet_asymmetry_config']['cmd'],
    output:
        out_t1w_nii=join(config['out_dir'], 'derivatives','pet_asymmetry','sub-' + subject_id, 'in',"T1.nii"),
        out_pet_nii=join(config['out_dir'], 'derivatives','pet_asymmetry','sub-' + subject_id, 'in',"PET.nii"),
        out_flair_nii=join(config['out_dir'], 'derivatives','pet_asymmetry','sub-' + subject_id, 'in',"FLAIR.nii"),
        out_json=join(config['out_dir'], 'derivatives','pet_asymmetry','sub-' + subject_id, 'in','sub-' + subject_id+".json"),
        #out_report=join(config['out_dir'], 'derivatives','pet_asymmetry','sub-' + subject_id, 'out','report',"report.html"),
    log: join(config['out_dir'], 'derivatives','pet_asymmetry','sub-' + subject_id, 'out',"pet_ai_log.txt") 
    run:
        import shutil,os,gzip,json,re

        with gzip.open(input.t1, 'rb') as f_in:
            with open(output.out_t1w_nii, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        
        with gzip.open(input.flair, 'rb') as f_in:
            with open(output.out_flair_nii, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

        with gzip.open(input.pet, 'rb') as f_in:
            with open(output.out_pet_nii, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        
        with open(params.in_json) as json_file:
            json_info = json.load(json_file)
        
        json_info['T1']=output.out_t1w_nii
        json_info['PET']=output.out_pet_nii
        json_info['FLAIR']=output.out_flair_nii
        json_info['output_dir']=os.path.join(os.path.dirname(os.path.dirname(output.out_json)),'out')
        json_output = json.dumps(json_info, indent=4)
        
        with open(output.out_json, 'w') as (fid):
            fid.write(json_output)
            fid.write('\n')

        import matlab.engine
        eng = matlab.engine.start_matlab()
        eng.cd(params.in_dir, nargout=0)
        try:
            eng.run_PET_AI(output.out_json)
        except RuntimeError as e:
            print(e)

final_outputs.extend(expand(rules.pet_asymmetry_copy.output.out_t1w_nii, subject=subjects))
#final_outputs.extend(expand(rules.pet_asymmetry_copy.output.out_report, subject=subjects))
