#!/bin/bash
#SBATCH --job-name=trimmomatic
#SBATCH --output=ZAMtrimmomatic_%a.log 
#SBATCH --mem=4G
#SBATCH --array=0-2
#SBATCH --time=24:00:00
#SBATCH --mail-user=cfcollins1@sheffield.ac.uk
#SBATCH --mail-type=FAIL

source /users/bop22cfc/miniconda3/bin/activate
source activate /users/bop22cfc/miniconda3/envs/trimmomatic

sample_folder="/mnt/parscratch/users/bop22cfc/methylation-work-flow/2_access_data/D33_ZAM"
accession_folder="D33_ZAM"
accession="ZAM"

# set the path to trimmomatic
ProgramPath="/mnt/parscratch/users/bop22cfc/methylation-work-flow/4_trimmomatic/Trimmomatic-0.39"

# get a list of all the forward reads
SAMPLE1=($(ls $sample_folder/*_1.fq))
# get a list of all the reverse reads
SAMPLE2=($(ls $sample_folder/*_2.fq))

# set up an index for the reads/tasks
INDEX=${SLURM_ARRAY_TASK_ID}

# run trimmomatic on the 3 sets
java -jar $ProgramPath/trimmomatic-0.39.jar PE -phred33 \
${SAMPLE1[$INDEX]} ${SAMPLE2[$INDEX]} \
${SAMPLE1[$INDEX]}.out_paired.fq ${SAMPLE1[$INDEX]}.out_unpaired.fq \
${SAMPLE2[$INDEX]}.out_paired.fq ${SAMPLE2[$INDEX]}.out_unpaired.fq \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
