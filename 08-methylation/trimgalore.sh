#!/bin/bash
#SBATCH --job-name=trim_galore
#SBATCH --output=trim_galore_%a.log  # %a will be replaced with the array index
#SBATCH --mem=16G
#SBATCH --array=0-3
#SBATCH --time=24:00:00
#SBATCH --mail-user=cfcollins1@sheffield.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

source /users/bop22cfc/miniconda3/bin/activate
source activate /users/bop22cfc/miniconda3/envs/trim-galore

sample_folder="/mnt/parscratch/users/bop22cfc/methylation-work-flow/2_access_data/catenated"

# get a list of all the forward reads
SAMPLE1=($(ls $sample_folder/*_1.fq))
# get a list of all the reverse reads
SAMPLE2=($(ls $sample_folder/*_2.fq))

# set up an index for the reads/tasks
INDEX=${SLURM_ARRAY_TASK_ID}

# run trim galore on the 4 sets
trim_galore --paired --fastqc \
${SAMPLE1[$INDEX]} \
${SAMPLE2[$INDEX]}
