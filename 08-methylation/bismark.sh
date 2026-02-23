#!/bin/bash
#SBATCH --job-name=bismark_run
#SBATCH --output=bismark_run_%a.log
#SBATCH --mem=64G
#SBATCH --array=0-2

#SBATCH --time=24:00:00
#SBATCH --mail-user=cfcollins1@sheffield.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL 

# source the correct programs
source /users/bop22cfc/miniconda3/bin/activate
source activate /users/bop22cfc/miniconda3/envs/bismark

# set the paths
sample_folder="/mnt/parscratch/users/bop22cfc/methylation-work-flow/4_trim_galore"
accession="AUS1"

# get a list of all the forward reads
SAMPLE1=($(ls $sample_folder/$accession/*_1.fq))
# get a list of all the reverse reads
SAMPLE2=($(ls $sample_folder/$accession/*_2.fq))

# set up an index for the reads/tasks
INDEX=${SLURM_ARRAY_TASK_ID}

bismark \
--path_to_bowtie /users/bop22cfc/miniconda3/envs/bowtie2/bin \
--samtools_path /users/bop22cfc/miniconda3/envs/samtools/bin \
--sam \
--genome_folder /mnt/parscratch/users/bop22cfc/methylation-work-flow/1_bismark_genome_preparation/$accession \
-1 ${SAMPLE1[$INDEX]} \
-2 ${SAMPLE2[$INDEX]}
