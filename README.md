# electrodeProc

A Snakemake workflow for processing SEEG electrode imaging data in BIDS-style datasets.

## Overview

`electrodeProc` is a neuroimaging workflow for processing intracranial electrode datasets, with emphasis on SEEG imaging, anatomical registration, tissue segmentation, electrode localization, atlas labeling, and visual quality control.

The workflow is built around Snakemake and is configured through `config/config.yml`. It supports modular execution, allowing users to turn major processing stages on or off depending on the dataset and analysis goal.

## Main Features

* BIDS-style subject discovery
* Subject-level anatomical registration
* Tissue segmentation
* SEEG electrode coordinate processing
* Atlas-based electrode labeling
* Native and template-space visualization
* Registration quality control
* Segmentation quality control
* Optional FastSurfer integration
* Optional fMRIPrep integration
* Optional HippUnfold integration
* Optional MELD integration
* Optional PET asymmetry workflow
* Optional contact segmentation workflow


## Installation

Clone the repository:

```bash
git clone https://github.com/greydongilmore/electrodeProc.git
cd electrodeProc
git checkout development
```

Create and activate a Python environment:

```bash
python -m venv .venv
source .venv/bin/activate
```

Install Python dependencies:

```bash
pip install -r requirements.txt
```

For larger datasets or cluster use, a Conda or Mamba environment is recommended.

## External Dependencies

This workflow may require external neuroimaging tools depending on which modules are enabled.

Common dependencies include:

* Snakemake
* ANTs
* FSL
* Convert3D / c3d
* Singularity or Apptainer
* FreeSurfer, if using FreeSurfer-based modules
* FastSurfer, if `fastsurfer.run` is enabled
* fMRIPrep, if `fmriprep.run` is enabled
* HippUnfold, if `hippunfold.run` is enabled
* MELD, if `meld.run` is enabled

Several container paths and external tool paths are configured in `config/config.yml`.

## Input Data

The workflow expects a BIDS-style dataset.

At minimum, the dataset should contain a BIDS folder with subject directories:

```text
bids/
├── sub-001/
├── sub-002/
└── sub-003/
```

The workflow can identify subjects in one of two ways:

1. From a `participants_run.tsv` file in the configured BIDS directory
2. From subject folders in the BIDS dataset

The configured subject prefix is defined in `config/config.yml`.

## Configuration

Set the main BIDS directory:

```yaml
bids_dir: /path/to/project
```

Enable or disable workflow modules:

```yaml
registration:
  run: True

segmentation:
  run: True

visqc:
  run: True

fastsurfer:
  run: False

hippunfold:
  run: False

meld:
  run: False

pet_asymmetry:
  run: False

contact_seg:
  run: False
```

Configure input image types in the imaging volume section. Example:

```yaml
contrast_t1:
  present: True
  session: 'pre'
  datatype: 'anat'
  run: '02'
  suffix: 'T1w'
  ext: '.nii.gz'
  algo: greedy

post_image:
  present: True
  session: 'post'
  datatype: 'ct'
  acq: 'Electrode'
  run: '01'
  suffix: 'ct'
  ext: '.nii.gz'
  algo: greedy
```

## Running the Workflow

Perform a dry run first:

```bash
snakemake -nr
```

Run locally using available CPU cores:

```bash
snakemake -j8
```

## Main Workflow Modules

### Registration

Handles subject-level and template-level registration.

Supported registration backends include:

* `greedy`
* `reg_aladin`
* `ants`

Registration parameters are configured in `config/config.yml`.

### Segmentation

Runs tissue segmentation workflows and generates tissue-class outputs.

Default tissue labels include:

```text
GM
WM
CSF
```

### Electrode Processing

Processes SEEG coordinate files and related electrode localization outputs when electrode inputs are enabled.

Expected coordinate patterns are configured in `config/config.yml`.

### Visual QC

Generates visual outputs for checking:

* Registration quality
* Tissue segmentation
* Atlas segmentation
* Electrode localization
* Native-space and template-space outputs

### Optional Connected Pipelines

The workflow includes optional hooks for:

* FastSurfer
* fMRIPrep
* HippUnfold
* MELD
* PET asymmetry analysis
* Contact segmentation

Enable these only after confirming the required containers, licenses, and external paths are available.

## Outputs

Outputs are written under the configured project directory, typically inside derivative folders.

Expected outputs may include:

* Registered anatomical images
* Transform files
* Tissue segmentations
* Electrode coordinate derivatives
* Atlas labels
* Visual QC images
* HTML or image-based reports
* Pipeline-specific derivative folders

Exact outputs depend on which modules are enabled.

## Citation

If you use this workflow in a publication or project, cite the repository:

```text
Gilmore G. electrodeProc. GitHub repository.
https://github.com/greydongilmore/electrodeProc
```

## License

This project is distributed under the MIT License. See `LICENSE` for details.

## Author

* Greydon Gilmore @greydongilmore
