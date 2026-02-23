#!/bin/bash
#SBATCH --job-name=AUS1_bismark_extractor
#SBATCH --output=AUS1_bismark_extractor.log
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mail-user=cfcollins1@sheffield.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL 

# source the correct programs
source /users/bop22cfc/miniconda3/bin/activate
source activate /users/bop22cfc/miniconda3/envs/bismark

sample_folder="/mnt/parscratch/users/bop22cfc/methylation-work-flow/5_run_bismark"

bismark_methylation_extractor \
-p \
--no_overlap \
--report \
--multicore 10 \
--comprehensive \
--bedGraph \
$sample_folder/AUS1_1_val_1_bismark_bt2_pe.bam
