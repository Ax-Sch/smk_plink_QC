# Data Cleaning Pipeline

This is a simple data cleaning pipeline built with Snakemake for the binary PLINK file set: `.bed`, `.bim`, and `.fam`. The pipeline automates the data cleaning process to ensure reproducibility. It assumes that the files are in binary PLINK format with phenotypes already added.

## Table of Contents

- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Configuration](#configuration)
- [Acknowledgements](#acknowledgements)

## Requirements

- [Git](https://git-scm.com/)
- [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
- [Gunzip](https://www.gnu.org/software/gzip/)

## Installation

### 1. Clone the Repository

First, clone the repository to your local machine:

```sh
git clone https://github.com/tanavader0/smk_plink_QC.git
cd smk_plink_QC
```

### 2.  Set Up the Conda Environment

Create and activate the conda environment:
```sh
conda env create -f workflow/envs/snakemake8.yaml
conda activate snakemake8
```

### 3. Usage
1. Configure Input Files and Parameters

Navigate to the config.yaml file in the config directory. Here, you can adjust parameters according to your requirements. For example, you can set the path to your input files:
```sh
input_plink: "path/to/your/input_files/prepared.fam"
```
The input files for this pipeline are genotype data in .bim, .bed, and .fam formats, with sex information and case/control status already included as phenotypes.

2. Run the Pipeline

Perform a dry run to ensure everything is set up correctly:
```sh
snakemake -np
```
If the dry run is successful, run the pipeline using the following command (replace 1 with the number of cores you want to use):
```sh
snakemake --cores 1 --use-conda
```

### 4. Configuration
All configuration options can be found in the config.yaml file. This file allows you to customize paths, parameters, and settings to fit your specific needs.

### 5. Acknowledgements
I would like to thank my supervisors, Dr. Axel Schmidt ([@Ax-Sch](https://github.com/Ax-Sch)) and Dr. Kerstin Ludwig, for their guidance and support. Additionally, I would like to thank the developers of the software used within this repository.
