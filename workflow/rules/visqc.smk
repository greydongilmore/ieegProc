
rule qc_reg:
    input:
        ref = config['template_t1w'],
        flo = bids(root=join(config['out_dir'], 'deriv', 'atlasreg'),subject=subject_id,suffix='T1w.nii.gz',space='{template}',desc='{desc}'),
    output:
        png = report(bids(root=join(config['out_dir'], 'deriv', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='regqc.png',from_='subject', to='{template}',desc='{desc}',include_subject_dir=False),
                caption='../reports/regqc.rst',
                category='Registration QC',
                subcategory='{desc} {template}'),
        html = bids(root=join(config['out_dir'], 'deriv', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='regqc.html',from_='subject', to='{template}', desc='{desc}',include_subject_dir=False),
#        html = report(bids(root='qc',subject=subject_id,suffix='regqc.html',from_='subject', to='{template}', desc='{desc}'),
#                caption='../reports/regqc.rst',
#                category='Registration QC',
#                subcategory='{desc} {template}'),
    group: 'preproc'
    script: '../scripts/vis_regqc.py'

rule qc_reg_ct:
    input:
        ref = bids(root=join(config['out_dir'], 'deriv', 'atlasreg'),subject=subject_id,suffix='T1w.nii.gz'),
        flo = bids(root=join(config['out_dir'], 'deriv', 'atlasreg'),subject=subject_id,suffix='ct.nii.gz',from_='atropos3seg',desc='masked')
    output:
        png = report(bids(root=join(config['out_dir'], 'deriv', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='regqc.png',from_='ct', to='T1w',desc='masked',include_subject_dir=False),
                caption='../reports/regqc.rst',
                category='Registration QC',
                subcategory='{desc} T1w'),
        html = bids(root=join(config['out_dir'], 'deriv', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='regqc.html',from_='ct', to='T1w', desc='masked',include_subject_dir=False),
#        html = report(bids(root='qc',subject=subject_id,suffix='regqc.html',from_='subject', to='{template}', desc='{desc}'),
#                caption='../reports/regqc.rst',
#                category='Registration QC',
#                subcategory='{desc} {template}'),
    group: 'preproc'
    script: '../scripts/vis_regqc.py'

rule qc_probseg:
    input:
        img = bids(root=join(config['out_dir'], 'deriv', 'atlasreg'),subject=subject_id,desc='n4', suffix='T1w.nii.gz'),
        seg4d = bids(root=join(config['out_dir'], 'deriv', 'atlasreg'),subject=subject_id,suffix='probseg.nii.gz',desc='atropos3seg'),
    params:
        ct_png=bids(root=join(config['out_dir'], 'deriv', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='regqc.png',from_='ct', to='T1w',desc='masked',include_subject_dir=False)
    output:
        png = report(bids(root=join(config['out_dir'], 'deriv', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='probseg.png',desc='atropos3seg',include_subject_dir=False),
                caption='../reports/segqc.rst',
                category='Segmentation QC',
                subcategory='3-class Tissue Segmentation'),
    group: 'preproc'
    script: '../scripts/vis_qc_probseg.py'

rule qc_dseg:
    input:
        img = bids(root=join(config['out_dir'], 'deriv', 'atlasreg'),subject=subject_id,desc='n4', suffix='T1w.nii.gz'),
        seg = bids(root=join(config['out_dir'], 'deriv', 'atlasreg'),subject=subject_id,suffix='dseg.nii.gz',atlas='{atlas}',from_='{template}',reg='SyN'),
    output:
        png = report(bids(root=join(config['out_dir'], 'deriv', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='dseg.png',atlas='{atlas}', from_='{template}',include_subject_dir=False),
                caption='../reports/segqc.rst',
                category='Segmentation QC',
                subcategory='{atlas} Atlas from {template}'),
        html = bids(root=join(config['out_dir'], 'deriv', 'atlasreg'),prefix='sub-'+subject_id+'/qc/sub-'+subject_id,suffix='dseg.html',atlas='{atlas}', from_='{template}',include_subject_dir=False),
#        html = report(bids(root='qc',subject=subject_id,suffix='dseg.html',atlas='{atlas}', from_='{template}'),
#                caption='../reports/segqc.rst',
#                category='Segmentation QC',
#                subcategory='{atlas} Atlas from {template}'),
    group: 'preproc'
    script: '../scripts/vis_qc_dseg.py'



