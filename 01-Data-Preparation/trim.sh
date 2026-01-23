#!/bin/bash
#$ -N RNA-seq_boa18bta
#$ -l h_rt=24:00:00
#$ -l rmem=12G
#$ -pe smp 1
#$ -t 1-46 # max number of tasks
#$ -tc 6 # number running concurrently
#$ -o /fastdata/boa18bta/RNA-seq/logs/trim-logs
#$ -e /fastdata/boa18bta/RNA-seq/errors/trim-errors
source /usr/local/extras/Genomics/.bashrc

#trimmomatic \
#PE -phred33 \
#"/shared/dunning_lab/Shared/SequencingProjects/Alloteropsis/RNA-seq-2022/X204SC22110619-Z01-F001/01.RawData/L1/L1_2.fq.gz" "/shared/dunning_lab/Shared/SequencingProjects/Alloteropsis/RNA-seq-2022/X204SC22110619-Z01-F001/01.RawData/L1/L1_1.fq.gz" \
#02-trimmed/L1_2_paired.fastq 02-trimmed/L1_2_unpaired.fastq \
#02-trimmed/L1_1_paired.fastq 02-trimmed/L1_1_unpaired.fastq \
#ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
#LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:50

#### Trimming and Filering ## currently set up for dry run
sample_directory=/shared/dunning_lab/Shared/SequencingProjects/Alloteropsis/RNA-seq-2022/X204SC22110619-Z01-F001/01.RawData
wd=/fastdata/boa18bta/RNA-seq

cd ${wd}

for directory in  ${sample_directory}/L*
do
for file1 in $directory/*1.fq.gz
do
cd ${wd}/$(basename $directory)
mkdir 02-trimmed
file2=${file1%1.fq.gz}2.fq.gz
trimmomatic PE -phred33 \
$file1 $file2 \
02-trimmed/$(basename -s .fq.gz $file1).paired.fastq 02-trimmed/$(basename -s .fq.gz $file1).unpaired.fastq.gz \
02-trimmed/$(basename -s .fq.gz $file2).paired.fastq 02-trimmed/$(basename -s .fq.gz $file2).unpaired.fastq.gz \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:50
done
done
