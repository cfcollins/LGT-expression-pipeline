#!/bin/bash
#$ -N salmon-AUS1
#$ -l h_rt=96:00:00
#$ -l rmem=12G
#$ -pe openmp 1 # threads                                                           
#$ -v OP_NUM_THREADS=4
#$ -e salmon_AUS1.error
#$ -o salmon_AUS1.log

source /usr/local/extras/Genomics/.bashrc

# Define paths
genome_index="/mnt/fastdata/bop22cfc/expression_project/Salmon/run_Salmon/index_ref_cds/AUS1_index"  
output_dir="/mnt/fastdata/bop22cfc/expression_project/Salmon/align-to-AUS1-ref/AUS1"    
threads=8 


# Loop through samples
samples=("L13" "L14" "L46" "L48" "L54" "L86")

for sample in "${samples[@]}"; do
    # Run Salmon for each sample
    salmon quant -i $genome_index -l A -1 /shared/dunning_lab/User/boa18bta-backup/RNA-seq/new_data_alignment/02-trimmed/${sample}_1.paired.fastq -2 /shared/dunning_lab/User/boa18bta-backup/RNA-seq/new_data_alignment/02-trimmed/${sample}_2.paired.fastq \
                 -p $threads --validateMappings -o $output_dir/${sample}_salmon

done

