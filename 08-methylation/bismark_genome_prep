#!/bin/bash
#SBATCH --job-name=KWT-bismark-genome-prep
#SBATCH --output=KWT-bismark-genome-prep.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16GB
#SBATCH --time=24:00:00
#SBATCH --mail-user=bop22cfc@sheffield.ac.uk
#SBATCH --mail-type=FAIL

source /users/bop22cfc/miniconda3/bin/activate
source activate /users/bop22cfc/miniconda3/envs/bismark

bismark_genome_preparation --path_to_aligner /users/bop22cfc/miniconda3/envs/bowtie2 --verbose /mnt/parscratch/users/bop22cfc/bismark/KWT
