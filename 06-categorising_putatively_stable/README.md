# 06 Categorising Putatively Stable

To categorise the laterally acquired genes, I conducted the following steps:
- are they expressed significantly less than the vertically inherited and the donor orthologs? (detailed in previous section 04-degenerating_categorisation_1) if yes -> degenerating. if no -> next step
- are they truncated (<70% of other orthologous genes, lacking start/stop codon)? (detailed in previous section 05-degenerating_categorisation_2) if yes -> degenerating. if no -> next step
- these genes were classified as 'putatively stable' and went through further categorisation (detailed in this section)

We used Kruskal-Wallis and Dunn's test to categorise the putatively stable laterally acquired genes based on expression levels. We then use binomial tests to identify well-represented categories. 

**Paralogs summed:** SUM_ALL_categorising_putatively_stable.R

**Paralogs summed - leaf/sheath only:** SUM_LEAF_SHEATH_categorising_putatively_stable.R

**Most expressed paralog used:** TOP_GENE_ALL_categorising_putatively_stable.R

**Most expressed paralog used - leaf/sheath only:** TOP_GENE_LEAF_SHEATH_categorising_putatively_stable.R

**Single genes only:** SINGLETONS_categorising_putatively_stable.R
