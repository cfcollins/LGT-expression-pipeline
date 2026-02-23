# 08 Methylation

### Generating methylation data

DNA was extracted from the four A. semialata accessions and sent to Novogene for whole genome bisulfite library preparation and sequencing. 

The bisulfite reads were trimmed using TrimGalore v0.6.10

Script: trimgalore.sh

### Bismark

Bisulfate sequencing data was analysed using Bismark 

Step 1: Index the reference genomes for each accession

Script: bismark_genome_prep.sh

Step 2: Run Bismark

Script: bismark.sh

Step 3: Bismark Methylation Extractor

Script: bismark_extractor.sh

### Text parsing to match methylation data to gene information

The file below includes coding notes for parsing information out of methylation outputs - please note it should not be run as a script and should be considered line-by-line 

text_parsing_methylation.sh

The pipeline within these notes is as follows:
- Step 1: editing the bedGraph and coverage file to be useable for analysis
- Step 2: calculating methylated cytosines for whole gene region
- Step 3: calculating methylated cytosines for exons only
- Step 4: calculating methylated cytosines for upstream and downstream the gene

This is an example counting script

Script: counting.sh

### Comparing methylation between LGT and vertically inherited genes

Script: methylation_LGT_vs_Native.R

### Comparing methylation to expression levels in laterally acquired genes

**Paralogs summed:** methylation_lower_silenced_vs_everything_else_SUM_ALL.R

**Paralogs summed - leaf/sheath only:** methylation_lower_silenced_vs_everything_else_SUM_leaf_sheath.R

**Most expressed paralog used:** methylation_lower_silenced_vs_everything_else_TOP_GENE_ALL.R

**Most expressed paralog used - leaf/sheath only:** methylation_lower_silenced_vs_everything_else_TOP_GENE_leaf_sheath.R

**Single genes only:** methylation_lower_silenced_vs_everything_else_singletons.R

### Creating Figure 4

Script: methylation_figure.R
