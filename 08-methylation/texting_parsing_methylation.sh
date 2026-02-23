############ SCRIPT FOR PARSING INFORMATION OUT OF METHYLATION OUTPUTS #################
######################## Cat Collins + Luke Dunning ####################################
############################### March 2024 #############################################

### SETTING UP

# Define variables for paths and accessions
# gff=/mnt/autofs/dunning_lab/Shared/reference_genomes/Alloteropsis_semialata_GCA_004135705.1/ASEM_AUS1_v1.0.gff
# gff=/mnt/autofs/christin_lab1/genome_data/Assemberlies/REFERENCE/KWT/KWT_v1.0.all.gff
# gff=/mnt/autofs/christin_lab1/genome_data/Assemberlies/REFERENCE/L04B/ASEM_L04B_v1.0.gff
 gff=/mnt/autofs/christin_lab1/genome_data/Assemberlies/REFERENCE/ZAM-05-10/ASEM_ZAM15-05-10.v1.0.gff
accession="AUS1"
base_dir="/mnt/parscratch/users/bop22cfc/methylation-work-flow"
gene_dir="${base_dir}/7_genes_of_interest/${accession}"
bismark_dir="${base_dir}/6_bismark_extractor/${accession}"
output_dir="${gene_dir}/1_whole_gene_region"

# Move to working directory
cd "$gene_dir"

## Activate Bedtools
source /users/bop22cfc/miniconda3/bin/activate
source activate /users/bop22cfc/miniconda3/envs/bedtools/

## Create some folders
mkdir -p 1_whole_gene_region
mkdir -p 2_exons_only
mkdir -p 3_upstream
mkdir -p 4_downstream

### SETTING UP THE GENES LOCATIONS AND REDUCING THE METHYLATION BED FILE 

# Create a BED file from the GFF for every gene 
### NOTE: THIS NEEDS TO BE CHANGED DEPENDING ON ACCESSION!!!
# AUS1: grep "\sgene\s" ${gff} | sed 's/Name=/\t/g' | sed 's/-gene//g' | cut -f 1,4,5,7,10 > ${accession}-genes.bed
# KWT, L04 & ZAM: grep "\sgene\s" ${gff} | cut -f 1,4,5,7,9 | sed 's/Name=/\t/g' | sed 's/|size/\t/g' | cut -f 1,3,4,5,8 > ${accession}-genes.bed

# Edit the bedGraph file so it removes the first line (a title line)
tail -n +2 "${bismark_dir}/${accession}_bedgraph" > ${accession}-mod.1

# Combine the bedGraph file and the coverage file and reduce 
# (because each file has some useful information that it would be good to keep)
paste ${accession}-mod.1 ${bismark_dir}/${accession}.cov | cut -f 1,2,3,8,9,10 > ${accession}-mod.2
# for AUS1 
${accession}-mod.2 > ${accession}-mod.3
# for non-AUS1 accessions remove the size bit 
# cat ${accession}-mod.2 | sed 's/|size/\t/g' | cut -f 1,3,4,5,6,7 > ${accession}-mod.3

# Reduce this file using Bedtools to make it more manageable
# This finds the intersection between the two files and directs the output into a new file 
# Keeps the basepairs that are within these genes
bedtools intersect -a ${accession}-mod.3 -b ${accession}-genes.bed -nonamecheck > ${accession}-trim-methyl.bed

# This filters the file so that only cytosines with 10 or above reads aligned are included
# It does this by adding together columns 5 and 6 and seeing if it is greater or equal to 10
cat ${accession}-trim-methyl.bed | awk '$5+$6>9' > ${accession}-trim-methyl-depth10.bed 

# This filters the file so that only cytosines with 50 or below reads aligned are included (as above)
cat ${accession}-trim-methyl-depth10.bed | awk '$5+$6<51' > ${accession}-trim-methyl-depth10-50.bed 

### SECTION 1: METHYLATION COORDINATES FOR THE WHOLE GENE REGION

# Make a new directory for gene files
mkdir -p ${output_dir}/${accession}-gene-bed

# This creates a loop to generate a file for every gene in the genome containing the BED file information for that gene
cat ${accession}-genes.bed | cut -f 5 | while read line ; do grep "$line" ${accession}-genes.bed > ${output_dir}/${accession}-gene-bed/"$line".bed ; done

# Make a new directory 
mkdir -p ${output_dir}/${accession}-gene-methyl

# This loops through each gene file made above and uses Bedtools intersect to find the cytosines for that gene 
# It inputs the coverage information for these cytosines into the relevant gene file 
# Warning: This takes a very very long time 
# Splitting the gene-bed folder into 10 subfolders and running it separately on each speeds it up a lot 
cat ${accession}-genes.bed \
| cut -f 5 \
| while read line ; \
	do bedtools intersect \
	-a ${accession}-trim-methyl-depth10-50.bed \
	-b 1_whole_gene_region/${accession}-gene-bed/"$line".bed \
	-nonamecheck >> 1_whole_gene_region/${accession}-gene-methyl/"$line" ; 
 done 

## some code for splitting into subfolders
cd 1_whole_gene_region/${accession}-gene-bed
for i in {01..10}; do mkdir "$i"; done 
# change this for each subfolder
subfolder="01"
ls -1 | head -n 6710 | tail -n +11 | xargs -I {} mv "{}" "${subfolder}/"

# This counts the total number of cytosines in this gene with reads aligned
cat ${accession}-genes.bed | cut -f 5 | while read line ; do cat 1_whole_gene_region/${accession}-gene-methyl/"$line".bed | wc -l | sed 's/^/'$line'\t/g' >> 1_whole_gene_region/total_sites.txt ; done
# This counts the totwl number of cytosines where 50% and over reads are methylated 
cat ${accession}-genes.bed | cut -f 5 | while read line ; do cat gene-methyl/"$line" | awk '$4>=50'  | wc -l | sed 's/^/'$line'\t/g' >> 1_whole_gene_region/methylated_sites.txt ; done
# This compiles these two files into a results file
paste 1_whole_gene_region/total_sites.txt 1_whole_gene_region/methylated_sites.txt | cut -f 1,2,4 > 1_whole_gene_region/AUS1-results.txt

### SECTION 2: METHYLATION COORDINATES FOR EXONS ONLY 

exon_dir=${gene_dir}/2_exons_only

## The code below is for the GFF in the christin_lab1 folder 
## It involves extracting the exon without the gene IDs (as for some reason they don't appear for the exons in this GFF
# gff=/mnt/autofs/christin_lab1/genome_data/Assemberlies/REFERENCE/AUS1/alloteropsis_semialata_final_assembly/ASEM_C4_v1.0.all.gff
# grep "\sgene\s" ${gff} | cut -f 9 | sed 's/;Name=/\t/g' | sed 's/ID=//g' | sed 's/-RA-gene//g' > CDS_to_gene.txt
# grep "\sexon\s" ${gff} >> Exon.txt
# mkdir CDS-bed
# cat CDS_to_gene.txt | while read line ; do var1=$(echo "$line" | cut -f 1 ) ; var2=$(echo "$line" | cut -f 2 ) ; grep "$var1-mRNA" Exon.txt | cut -f 1,4,5,7 | sed 's/$/\t'$var2'/g'  >>  CDS-bed/"$var2".bed; done

# however the gff in dunning_lab does have the gene IDs so I've just used this code

# extract BED file for exons
cat ${gff} | grep "\sexon\s" | sed 's/Name=/\t/g' | sed 's/-gene:/\t/g' | cut -f 1,4,5,7,10 > ${accession}_exons.bed # adapt based on accession 
# cat ${gff} | grep "\sexon\s" | sed 's/ID=/\t/g' | sed 's/-RA/\t/g' | sed 's/|size/\t/g' | cut -f 1,5,6,8,11 > ${accession}_exons.bed

# extract gene IDs
cat ${gff} | grep "\sgene\s" | sed 's/Name=/\t/g' | cut -f 10 | sed 's/-gene//g' > ${accession}-geneIDs.txt # adapt based on accession 
# cat ${gff} | grep "\sgene\s" | sed 's/ID=/\t/g' | sed 's/|size/\t/g' | cut -f 11 > ${accession}-geneIDs.txt

mkdir ${exon_dir}/${accession}-CDS-bed

# Make a loop that puts the exon lines for each gene into the same file 
cat ${accession}-geneIDs.txt | while read line ; do grep "$line" ${accession}_exons.bed >> 2_exons_only/${accession}-CDS-bed/"$line".bed; done

mkdir ${exon_dir}/${accession}-CDS-methyl

# again this takes a while
# split up (see above)
cat AUS1-GENEs.bed | cut -f 4 | while read line ; do bedtools intersect -a Trim-methyl-depth10-50.bed -b CDS-bed/"$line".bed >> CDS-methyl/"$line" ; done
# script saved for this

cat AUS1-GENEs.bed | cut -f 4 | while read line ; do cat CDS-methyl/"$line" | wc -l | sed 's/^/'$line'\t/g' >> CDStotal_sites.txt ; done
cat AUS1-GENEs.bed | cut -f 4 | while read line ; do cat CDS-methyl/"$line" | awk '$4>=50'  | wc -l | sed 's/^/'$line'\t/g' >> CDSmethylated_sites.txt ; done

paste CDStotal_sites.txt CDSmethylated_sites.txt > CDSAUS1-results.txt


### SECTION 3: UPSTREAM AND DOWNSTREAM 

## SET UP 

genome_file=/mnt/autofs/dunning_lab/Shared/reference_genomes/Alloteropsis_semialata_GCA_004135705.1/ASEM_C4_v1.0.fasta
#need genome file for bedtools

## Activate Samtools
source /users/bop22cfc/miniconda3/bin/activate
source activate /users/bop22cfc/miniconda3/envs/samtools/

# index the genome 
samtools faidx  ${genome_file} -o ${accession}.fai.fa

## Activate Bedtools
source /users/bop22cfc/miniconda3/bin/activate
source activate /users/bop22cfc/miniconda3/envs/bedtools/

# Identify flanking regions 
bedtools flank -i ${accession}-genes.bed -g ${accession}.fai.fa -b 1000 > flanks

grep "\s+\s" AUS1-genes.bed >> AUS1-GENEs_forward.bed  
grep "\s-\s" AUS1-genes.bed >> AUS1-GENEs_reverse.bed 
mkdir upstream-bed
mkdir downstream-bed

cat AUS1-GENEs_forward.bed \ 
| while read line ;  \ 
	do GENE=$(echo "$line" | cut -f 5 ) ; \
	start=$(echo "$line" | cut -f 2 ) ; \
	stop=$(echo "$line" | cut -f 3 ) ; \
	grep "$GENE" flanks | grep "$start" > upstream-bed/"$GENE".bed ; \
	grep "$GENE" flanks | grep "$stop" > downstream-bed/"$GENE".bed;  \
done
# cat AUS1-GENEs_forward.bed | while read line ; do GENE=$(echo "$line" | cut -f 5 ) ; start=$(echo "$line" | cut -f 2 ) ; stop=$(echo "$line" | cut -f 3 ) ; grep "$GENE" flanks | grep "$start" > upstream-bed/"$GENE".bed ; grep "$GENE" flanks | grep "$stop" > downstream-bed/"$GENE".bed; done

cat AUS1-GENEs_reverse.bed | while read line ;  do GENE=$(echo "$line" | cut -f 5 ) ; start=$(echo "$line" | cut -f 3 ) ; stop=$(echo "$line" | cut -f 2 ) ; grep "$GENE" flanks | grep "$start" > upstream-bed/"$GENE".bed ; grep "$GENE" flanks | grep "$stop" > downstream-bed/"$GENE".bed;  done

## then do the bedtools intersect and the counting like you did above for exons of whole gene (do upstream and downstream seperately)


cat flanks | awk '{ $5 = $3 - $2 } 1' | grep -v "\s1000$" | sed 's/ /\t/g' > flanks_less_1000bp.txt

#there are 101 flanking regions <1000 bp - do you want to remove these? can help out if we do


accession="AUS1"
working_dir=/mnt/parscratch/users/bop22cfc/methylation-work-flow/7_genes_of_interest/AUS1
type="upstream"

cat ${working_dir}/${accession}-genes.bed | cut -f 5 | while read line ; do cat ${type}-methyl/"$line".bed | wc -l | sed 's/^/'$line'\t/g' >> ${type}_total_sites.txt ; done
cat ${working_dir}/${accession}-genes.bed | cut -f 5 | while read line ; do cat ${type}-methyl/"$line".bed | awk '$4>=50' | wc -l | sed 's/^/'$line'\t/g' >> ${type}_methylated_sites.txt ; done
paste ${type}_total_sites.txt ${type}_methylated_sites.txt | cut -f 1,2,4 > ${accession}-${type}-results.txt
