#!/bin/bash
#SBATCH -J Smk_plink_QC                      # Job name, rename accordingly(optional)
#SBATCH --time=02-00:00:00               # Walltime, Probably ADJUST THIS to the constraints of your slurm setup
#SBATCH --partition=medium               # short, long or medium on our cluster
#SBATCH --threads-per-core=1
#SBATCH --mem=4G
#SBATCH --ntasks=1                       # 1 task
#SBATCH --cpus-per-task=1
#SBATCH --account=ag_ukihg_ludwig       # account - modify / delete; as needed on your cluster

snakemake --profile config/snakemake_slurm_profile/ --rerun-incomplete --use-conda
