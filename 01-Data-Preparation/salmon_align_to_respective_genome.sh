#!/bin/bash
#$ -N salmon-KWT
#$ -l h_rt=96:00:00
#$ -l rmem=12G
#$ -pe openmp 1 # threads                                                           
#$ -v OP_NUM_THREADS=4
#$ -e salmon_KWT.error
#$ -o salmon_KWT.log

source /usr/local/extras/Genomics/.bashrc

# Define paths
genome_index="/mnt/fastdata/bop22cfc/expression_project/Salmon/run_Salmon/index_ref_cds/KWT_index"  # Provide the path to the Salmon index for your reference genome/
output_dir="/mnt/fastdata/bop22cfc/expression_project/Salmon/run_Salmon/KWT"     # Choose a directory to store the output files
threads=8 


# Loop through samples
samples=("L1" "L55" "L59" "L6" "L84" "L85")

for sample in "${samples[@]}"; do
    # Run Salmon for each sample
    salmon quant -i $genome_index -l A -1 /shared/dunning_lab/User/boa18bta-backup/RNA-seq/new_data_alignment/02-trimmed/${sample}_1.paired.fastq -2 /shared/dunning_lab/User/boa18bta-backup/RNA-seq/new_data_alignment/02-trimmed/${sample}_2.paired.fastq \
                 -p $threads --validateMappings -o $output_dir/${sample}_salmon

done
