# 09 Selective Constraight Analysis 

### Generating dN/dS data 

Full coding notes for setting up alignments and running PAML are in the file below:

full_paml_coding_notes.sh

Please note that this is to run line-by-line rather than all in one go. 

I used alignment files from Raimondeau et al (2023). I reduced, cleaned and make trees for these files for PAML on Geneious.

### M0 model (one ratio)

I ran created a control file for each gene

This is template control file: M0_control_file_template.txt

I used this script to make all the control files: make_control_files.sh

I ran PAML with this script: run_paml.sh

### Branch model

Template control file: branch_model_control_file_template.txt

Then I reused these scripts: make_control_files.sh and run_paml.sh

### Branch site model

Template control file: branch_site_model_control_file_tempalate.txt

You need to compare the branch site model results to a null

Template null control file: branch_site_model_null_control_file_template.txt

### Stats

The coding notes within full_paml_coding_notes.sh include details for how to extract the data you need for analysis (e.g. omega and lnL)

I then analysed this data in R

Comparing LGT with vertically inherited genes: LGT_vs_native_dn_ds.R

Identifying statistical significance between LGT and vertically inherted: fishers_exact_dn_ds.R

Comparing degenerating with putatitively stable LGT genes: dn_ds_deg_vs_stab.R

Including truncation data: dn_ds_truncation.R
