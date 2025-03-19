# THIS IS CURRENTLY UNDER ACTIVE DEVELOPMENT




# Data Cleaning Pipeline

This is a simple data cleaning pipeline built with Snakemake for the binary PLINK file set: `.bed`, `.bim`, and `.fam`. The pipeline automates the data cleaning process to ensure reproducibility. It assumes that the files are in binary PLINK format with phenotypes such as sex and case/control status already added.

## Table of Contents

- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Output](#output)
- [Acknowledgements](#acknowledgements)

## Requirements

- [Git](https://git-scm.com/)
- [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
- [Gunzip](https://www.gnu.org/software/gzip/)

## Installation

### 1. Clone the Repository

First, clone the repository to your local machine:

```sh
git clone https://github.com/Ax-Sch/smk_plink_QC 
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

All configuration options can be found in the config.yaml file within the directory config. This file allows you to customize paths, parameters, and settings to fit your specific needs. For example, you can set the path to your input files, and Human genome reference version:
```sh
input_plink: "path/to/your/input_files/files.fam"
genome_ref:
    version: b37 or b38
```
An toy data set is located within the folder example_data to test if the pipeline is working in principle.

As mentioned earlier, the input files for this pipeline are genotype data in .bim, .bed, and .fam formats, with sex information and case/control status already included as phenotypes.

As part of the pipeline, Principal Component Analysis (PCA) is performed to filter out individuals who do not fall within the defined ranges for the ancestry you are interested in. The boundaries for the first two principal components (PC1 and PC2) are specified in the configuration file under the section "pca_ancestry_filters:". These ranges can be adjusted as necessary to meet the requirements of specific analyses.


2. Run the Pipeline

Perform a dry run to ensure everything is set up correctly:
```sh
snakemake -np
```
If the dry run is successful, run the pipeline using the following command (replace 1 with the number of cores you want to use):
```sh
snakemake --cores 1 --use-conda --conda-frontend conda
```

### 4. Output
All output directories and their corresponding snakemake rules have a capital letter as prefix (A-Z) to separate the corresponding steps:
- A_Prepare_correct_x
- B_VariantCallrate1 
- C_SampleCallrate 
- D_Filter_sex_checked 
- D_Graph_Sex_check
- D_Sex_check
- E_Check_heterozygosity 
- E_Filter_het_samples 
- E_Get_heterozygosity 
- F_VariantCallrate2 
- G_Check_MissDiff_HWE 
- G_Filter_MissDiff_HWE
- G_Get_MissDiff_HWE 
- H_Change_ID_for_1000G_PCA
- H_Download_1000G_chromosomes
- H_Download_1000G_sample_info 
- H_Download_fasta_files 
- H_Filter_plink_for_ancestry
- H_Make_PCA_plots 
- H_Merge_data_w_1000G_run_PCA_step3 
- H_Prepare_1000G_for_ancestry_PCA_step1
- H_Prepare_1000G_for_ancestry_PCA_step2
- H_Run_pca_filter 
- H_Unzip_fasta
- I_Kinship_analysis 
- I_Kinship_check2 
- I_Kinship_check2_R 
- I_kinship_analysis_R 
- I_kinship_scatter_plot 
- I_remove_relateds


### 5. Acknowledgements
I would like to thank my supervisors, Dr. Axel Schmidt ([@Ax-Sch](https://github.com/Ax-Sch)) and Dr. Kerstin Ludwig, for their guidance and support. Additionally, I would like to thank the developers of the software used within this repository.
