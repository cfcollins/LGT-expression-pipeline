# 01 Data Preparation

## RNA-Seq

### Trimming
Script: trim.sh

### Alignment and counts with Salmon
1. Created a Salmon index

Script: index.sh

2. Aligned to AUS1 genome (for clustering analysis)

Script: salmon_align_to_AUS1.sh

3. Aligned to respective genome

Script: salmon_align_to_respective_genome.sh

### Creating count matrices 
Script: count_matrix_collation.R

Script for per accession: creating_full_count_per_accession.R

### PCA to check clustering of expression data
Script: pca.R

### Analysing paralogs
Script: paralogs.R
