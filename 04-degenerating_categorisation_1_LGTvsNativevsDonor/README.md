# 04 Degenerating Categorisation step 1: LGT vs vertically inherited vs donor
To categorise the laterally acquired genes, I conducted the following steps:
- are they expressed significantly less than the vertically inherited and the donor orthologs? (detailed in this section) if yes -> degenerating. if no -> next step
- are they truncated (<70% of other orthologous genes, lacking start/stop codon)? (detailed in next section 05-degenerating_categorisation_2) if yes -> degenerating. if no -> next step
- these genes were classified as 'putatively stable' and went through further categorisation (detailed in 06-categorising_putatively_stable)


### Overall comparison for LGT, vertically inherited and donor genes
I used a linear mixed effects model to compare all LGT gene expression with vertically inherited and donor ortholog expression, accounting for hierarchy in paralog groups and accessions. Script also generates Figure 2b (although edited in Inkscape). 
 
Script: LMM_LGT_vs_nat_vs_don.R

### Remove inactive homologs
Laterally acquired genes where both the vertically inherited and donor orthologs have negligible gene expression is not biologically meaningful. Therefore, these genes were removed from analysis. 

Script: remove_inactive_homologs.R

### Is the laterally acquired gene expressed less than both its vertically inherited and donor homolog?
This script uses Wilcoxon signed-rank tests for each laterally acquired gene to identify whether it is significantly expressed less than both its vertically inherited and donor homologs. This analysis was repeated five times, for each way paralogs were handled in analysis and for when only leaf/sheath tissue was used.

**Paralogs summed:** LGT_less_than_natives_donors_SUM_ALL.R

**Paralogs summed - leaf/sheath only:** LGT_less_than_natives_donors_SUM_LEAF_SHEATH.R

**Most expressed paralog used:** LGT_less_than_natives_donors_TOP_GENE_ALL.R

**Most expressed paralog used - leaf/sheath only:** LGT_less_than_natives_donors_TOP_GENE_LEAF_SHEATH.R

**Single genes only:** LGT_less_than_natives_donors_SINGLETONS.R
