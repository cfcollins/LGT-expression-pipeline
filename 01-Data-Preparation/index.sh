#!/bin/bash
#$ -N Salmon-index
#$ -l h_rt=96:00:00
#$ -l rmem=12G
#$ -pe openmp 1 # threads                                                           
#$ -v OP_NUM_THREADS=4

source /usr/local/extras/Genomics/.bashrc

salmon index -t AANG_MRL48.all.maker.CDS.fasta -i /mnt/fastdata/bop22cfc/expression_project/Salmon/run_Salmon/index_ref_cds/MRL_index -k 21
salmon index -t ASEM_C4_v1.0.CDS.fasta -i /mnt/fastdata/bop22cfc/expression_project/Salmon/run_Salmon/index_ref_cds/AUS1_index -k 21
salmon index -t ASEM_L04B_v1.0.all.maker.cds.fasta -i /mnt/fastdata/bop22cfc/expression_project/Salmon/run_Salmon/index_ref_cds/L04B_index -k 21
salmon index -t ASEM_ZAM15-05-10.v1.0.all.maker.cds.fasta -i /mnt/fastdata/bop22cfc/expression_project/Salmon/run_Salmon/index_ref_cds/ZAM_index -k 21
salmon index -t KWT_v1.0.all.maker.CDS.fasta -i /mnt/fastdata/bop22cfc/expression_project/Salmon/run_Salmon/index_ref_cds/KWT_index -k 21
salmon index -t Setaria_italica.Setaria_italica_v2.0.cds.all.fa -i /mnt/fastdata/bop22cfc/expression_project/Salmon/run_Salmon/index_ref_cds/SET_index -k 21
salmon index -t Sorghum_bicolor.Sorghum_bicolor_NCBIv3.cds.all.fa -i /mnt/fastdata/bop22cfc/expression_project/Salmon/run_Salmon/index_ref_cds/THE_index -k 21


