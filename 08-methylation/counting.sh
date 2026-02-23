#!/bin/bash
#SBATCH --job-name=counting                             
#SBATCH --mem=32G
#SBATCH --time=72:00:00
#SBATCH --mail-user=cfcollins1@sheffield.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

accession="AUS1"
working_dir=/mnt/parscratch/users/bop22cfc/methylation-work-flow/7_genes_of_interest/AUS1
type="upstream"

cat ${working_dir}/${accession}-genes.bed | cut -f 5 | while read line ; do cat ${type}-methyl/"$line".bed | wc -l | sed 's/^/'$line'\t/g' >> ${type}_total_sites.txt ; done
cat ${working_dir}/${accession}-genes.bed | cut -f 5 | while read line ; do cat ${type}-methyl/"$line".bed | awk '$4>=50' | wc -l | sed 's/^/'$line'\t/g' >> ${type}_methylated_sites.txt ; done
paste ${type}_total_sites.txt ${type}_methylated_sites.txt | cut -f 1,2,4 > ${accession}-${type}-results.txt
