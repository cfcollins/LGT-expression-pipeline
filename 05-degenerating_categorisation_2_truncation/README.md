# 05 Degenerating Categorisation step 2: truncation 

To categorise the laterally acquired genes, I conducted the following steps:
- are they expressed significantly less than the vertically inherited and the donor orthologs? (detailed in previous section 04-degenerating_categorisation_1) if yes -> degenerating. if no -> next step
- are they truncated (<70% of other orthologous genes, lacking start/stop codon)? (detailed in this section) if yes -> degenerating. if no -> next step
- these genes were classified as 'putatively stable' and went through further categorisation (detailed in 06-categorising_putatively_stable)

### Generating truncation data
I used BLAST to find orthologous genes in good quality reference genomes (e.g. crop species) and then used Geneious to check each gene sets to see if the laterally acquired / vertically inherited genes were less than 70% gene length in comparison to the other orthologs. 

My coding notes are in the file below. Note: this script will not run on its own, it is command line lines of code for BLAST and text parsing within a HPC environment

sequence_blast.sh

I then tested whether truncated vertically inherited genes have full length LGTs

Script: truncation_w_full_length.R

### Comparing truncation between laterally acquired and vertically inherited genes
Script: truncation_LGT_vs_Native.R

### Do truncated genes have lower expression levels?
Script: truncation_effect_on_logCount.R

### Comparing truncation between laterally acquired genes previous classified as 'degenerating' in the previous section, to those who were not

**Paralogs summed:** SUM_ALL_truncation_vs_expression_category.R

**Paralogs summed - leaf/sheath only:** SUM_leaf_sheath_truncation_vs_expression_category.R

**Most expressed paralog used:** Top_gene_all_truncation_vs_expression_category.R

**Most expressed paralog used - leaf/sheath only:** Top_gene_leaf_sheath_truncation_vs_expression_category.R

**Single genes only:** Singletons_truncation_vs_expression_category.R

### Combining with the previous section to expand the degenerating category 

**Paralogs summed:** SUM_ALL_combine_truncation_w_expression.R

**Paralogs summed - leaf/sheath only:** SUM_LEAF_SHEATH_combine_truncation_w_expression.R

**Most expressed paralog used:** TOP_GENE_ALL_combine_truncation_w_expression.R

**Most expressed paralog used - leaf/sheath only:** TOP_GENE_LEAF_SHEATH_combine_truncation_w_expression.R

**Single genes only:** SINGLETONS_combine_truncation_w_expression.R


