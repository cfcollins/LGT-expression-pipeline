## line by line code for running PAML for LGTs / vertically inhertited genes
## note: this is intended to be run line-by-line, NOT as a script

# export alignments from Geneious as fasta files 
mkdir 01-alignments
# move these fasta files into 01-alignments
mkdir 02-fixed_titles

ls -1 01-alignments/ | while read line ; do cat 01-alignments/${line} | sed 's/_/-/g'| sed '/LGT/c\>LGT' | sed '/Alloteropsis/c\>Native' | sed '/Aristida-congesta/c\>ACON' | sed '/Dactyloctenium-aegyptium/c\>DAEG' | sed '/Oropetium-thomaeum/c\>OTHO' | sed '/Zoysia-japonica/c\>ZJAP' | sed '/Danthonia-californica/c\>DCAL' | sed '/Eriachne-aristidea/c\>EARI' | sed '/Brachypodium-distachyon/c\>BDIS' | sed '/Poa-sp/c\>Poa' | sed '/Hordeum-vulgare/c\>HVUL' | sed '/Triticum-aestivum/c\>TAES' | sed '/Stipagrostis-hirtigluma/c\>SHIR' | sed '/Oryza-sativa/c\>OSAT' | sed '/Eragrostis-tef/c\>ETEF' | sed '/Aristida-purpurea/c\>APUR' | sed '/Steinchisma-decipiens/c\>SDEC' | sed '/Pasplum-vaginatum/c\>PVAG' | sed '/Pasplum-fimbriatum/c\>PFIM' | sed '/Hymenachne-amplexicaulis/c\>HAMP' | sed '/Zea-mays/c\>ZMAY' | sed '/Sorghum-bicolor/c\>SBIC' | sed '/Dichanthium-sericeum/c\>DSER' | sed '/Arundinella-hirta/c\>AHIR' | sed '/Danthoniopsis-dinteri/c\>DDIN' | sed '/Otachyriinae-sp/c\>O_sp' | sed '/Chasmanthium-latifolium/c\>CLAT' | sed '/Arundinella-hookeri/c\>AHOO' | sed '/Leersia-perrieri/c\>LPER' | sed '/Digitaria-pentzii/c\>DPEN' | sed '/Chionochloa-macra/c\>CMAC' | sed '/Lasiacis-sorghoidea/c\>LSOR' | sed '/Dichanthelium-oligosanthes/c\>DOLI' | sed '/Cyrtococcum-patens/c\>CPAT' | sed '/Homopholis-proluta/c\>HPRO' | sed '/Acroceras-zizanioides/c\>AZIZ' | sed '/Steinchisma-sp/c\>S_sp' | sed '/Echinochloa-stagnina/c\>ESTA' | sed '/Panicum-vigratum/c\>PVIG' | sed '/Echinochloa-frumentacea/c\>EFRU' | sed '/Sacciolepis-indica/c\>SIND' | sed '/Panicum-queenslandicum/c\>PQUE' | sed '/Panicum-hallii/c\>PHAL' | sed '/Neurachne-alopecuroidea/c\>NALO' | sed '/Megathyrsus-maximus/c\>MMAX' | sed '/Echinochloa-crus/c\>ECRU' | sed '/Setaria-italica/c\>SITA' | sed '/Cenchrus-americanus/c\>CAME' | sed '/Panicum-pygmaeum/c\>PPYG' | sed '/Panicum-miliaceum/c\>PMIL' | sed '/Digitaria-ciliaris/c\>DCIL' | sed '/Dichanthelium-clandestinum/c\>DCLA' | sed '/Paraneurachne-muelleri/c\>PMUE' | sed '/Setaria-barbata/c\>SBAR' | sed '/Panicum-repens/c\>PREP' | sed '/Panicum-capillare/c\>PCAP' | sed '/Urochloa-fusca/c\>UFUS' | sed '/Echinochloa-esculenta/c\>EESC' | sed '/Anthephora-pubescens/c\>APUB' | sed '/Stenotaphrum-secundatum/c\>SSEC' | sed '/Panicum-trichanthum/c\>PTRI' | sed '/Panicum-coloratum/c\>PCOL' | sed '/Oplismenus-burmannii/c\>OBUR' | sed '/Urochloa-brizantha/c\>UBRI' | sed '/Sacciolepis-striata/c\>SSTR' | sed '/Panicum-virgatum/c\>PVIR' | sed '/Dichanthelium-scoparium/c\>DSCO' | sed '/Acroceras-calcicola/c\>ACAL' > 02-fixed_titles/${line}  ; done

# move the edited fasta files back into Geneious
# generate tree files for each 
# export alignments from Geneious as phylip files 
mkdir 03-phy_alignments
# move phylip alignments to 03-alignments
mkdir 04-phy_alignments_w_spaces
# Geneious only puts one space in phylip files for some reason - add another space
ls -1 03-phy_alignments/ | while read line ; do cat 03-phy_alignments/${line} | sed 's/ /  /g' > 04-phy_alignments_w_spaces/${line} ; done

mkdir 05-trees
# Geneious adds extra bits to the tree files, fix this
cd 05-trees
for file in *".phy PhyML Tree"*; do mv "$file" "${file//.phy PhyML Tree/}" ; done
cd ..

# make a list of files for reference for next steps
ls -1 04-phy_alignments_w_spaces/ | sed 's/.phy//g' > list_of_files.txt
# set up PAML
source /users/bop22cfc/miniconda3/bin/activate
source activate /users/bop22cfc/miniconda3/envs/paml

#################
### M0 (one ratio)
# model=0 and NSsites=0
mkdir 06-M0_model
mkdir 06-M0_model/01-control_files
mkdir 06-M0_model/02-output

# template M0 control file
seqfile = 04-phy_alignments_w_spaces/xxx.phy     * Path to alignment file
treefile = 05-trees/xxx.newick       * Path to tree file
outfile = 06-M0_model/02-output/xxx.txt        * Output file

noisy = 3                    * 0,1,2,3,9: how much rubbish on the screen
verbose = 1                  * 1: detailed output, 0: concise output
runmode = 0                  * 0: user tree; 1: semi-automatic; 2: automatic

seqtype = 1                  * 1:codons; 2:AAs; 3:codon-translated

CodonFreq = 2                * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
model = 0                    * 0: one omega; 1: nearly-neutral; 2: branch model

NSsites = 0                  * 0: one omega for all sites

fix_kappa = 0                * 1: kappa fixed; 0: kappa to be estimated
kappa = 2                    * initial or fixed kappa

fix_omega = 0                * 1: omega fixed; 0: omega to be estimated
omega = 1                    * initial omega

fix_alpha = 1                * 0: estimate gamma shape parameter; 1: no gamma
alpha = 0.0                  * initial alpha, 0:infinity (constant rate)

Malpha = 0                   * different alphas for genes
ncatG = 10                   * # of categories in dG of NSsites models

getSE = 0                    * 0: dont want them, 1: want S.E.s of estimates
RateAncestor = 0             * (1/0): rates (alpha>0) or ancestral states (alpha=0)
Small_Diff = 5e-7            * Default value.
cleandata = 1                * Remove sites with ambiguity data (1) or not (0)
method = 0                   * 0: simultaneous; 1: one branch at a time

# create this as "template.txt"
nano 06-M0_model/M0_template.txt

# here's a script for making the control files

##
#!/bin/bash

template="06-M0_model/M0_template.txt"

# Read each string in list.txt and create a new file
while read -r line; do
    # Define the new file name using the string
    newfile="${line}.ctl"
    
    # Replace 'xxx' in the template and save to the new file
    sed "s/xxx/$line/g" "$template" > 06-M0_model/01-control_files/"$newfile"
done < list_of_files.txt
##

# create this as "make_control_files.sh"
nano 06-M0_model/make_control_files.sh
sbatch 06-M0_model/make_control_files.sh

# and a script for running paml
##
#!/bin/bash

source /users/bop22cfc/miniconda3/bin/activate
source activate /users/bop22cfc/miniconda3/envs/paml

# Run codeml on all files in the 07-control_files directory
for file in 06-M0_model/01-control_files/*; do
  codeml "$file"
done
##
# make this into a script called 'run_paml.sh'
nano 06-M0_model/run_paml.sh
sbatch 06-M0_model/run_paml.sh

#######################
#######################

### Branch model
# model=2 and NSsites=0
mkdir 07-branch_model
# need to add a #1 to LGT and #2 to Native for each tree
mkdir 07-branch_model/01-tree_LGT_labelled
ls -1 05-trees/ | while read line ; do cat 05-trees/${line} | sed 's/LGT/LGT#1/g' | sed 's/Native/Native#2/g' > 07-branch_model/01-tree_LGT_labelled/${line} ; done
mkdir 07-branch_model/02-control_files
mkdir 07-branch_model/03-output

# template branch model file:
seqfile = 04-phy_alignments_w_spaces/xxx.phy     * Path to alignment file
treefile = 07-branch_model/01-tree_LGT_labelled/xxx.newick       * Path to tree file
outfile = 07-branch_model/03-output/xxx.txt        * Output file

noisy = 3                    * 0,1,2,3,9: how much rubbish on the screen
verbose = 1                  * 1: detailed output, 0: concise output
runmode = 0                  * 0: user tree; 1: semi-automatic; 2: automatic

seqtype = 1                  * 1:codons; 2:AAs; 3:codon-translated

CodonFreq = 2                * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
model = 2                    * 0: one omega; 1: nearly-neutral; 2: branch model

NSsites = 0                  * 0: one omega for all sites

fix_kappa = 0                * 1: kappa fixed; 0: kappa to be estimated
kappa = 2                    * initial or fixed kappa

fix_omega = 0                * 1: omega fixed; 0: omega to be estimated
omega = 1                    * initial omega

fix_alpha = 1                * 0: estimate gamma shape parameter; 1: no gamma
alpha = 0.0                  * initial alpha, 0:infinity (constant rate)

Malpha = 0                   * different alphas for genes
ncatG = 10                   * # of categories in dG of NSsites models

getSE = 0                    * 0: dont want them, 1: want S.E.s of estimates
RateAncestor = 0             * (1/0): rates (alpha>0) or ancestral states (alpha=0)
Small_Diff = 5e-7            * Default value.
cleandata = 1                * Remove sites with ambiguity data (1) or not (0)
method = 0                   * 0: simultaneous; 1: one branch at a time

# create this as "template.txt"
nano 07-branch_model/branch_model_template.txt

# here's a script for making the control files

##
#!/bin/bash

template="07-branch_model/branch_model_template.txt"

# Read each string in list.txt and create a new file
while read -r line; do
    # Define the new file name using the string
    newfile="${line}.ctl"
    
    # Replace 'xxx' in the template and save to the new file
    sed "s/xxx/$line/g" "$template" > 07-branch_model/02-control_files/"$newfile"
done < list_of_files.txt

###################################
# create this as "make_control_files.sh"
nano 07-branch_model/make_control_files.sh
sbatch 07-branch_model/make_control_files.sh

# and a script for running paml
##
#!/bin/bash

source /users/bop22cfc/miniconda3/bin/activate
source activate /users/bop22cfc/miniconda3/envs/paml

# Run codeml on all files in the 07-control_files directory
for file in 07-branch_model/02-control_files/*; do
  codeml "$file"
done
##
# make this into a script called 'run_paml.sh'
nano 07-branch_model/run_paml.sh
sbatch 07-branch_model/run_paml.sh


##################################
##################################
### Branch-site model
# model=2 and NSsites=2
mkdir 08-branch_site_model

## first for LGT genes
mkdir 08-branch_site_model/01-tree_LGT_labelled
ls -1 05-trees/ | while read line ; do cat 05-trees/${line} | sed 's/LGT/LGT#1/g' > 08-branch_site_model/01-tree_LGT_labelled/${line} ; done
mkdir 08-branch_site_model/02-LGT-control_files
mkdir 08-branch_site_model/03-LGT-output

# now make control files
seqfile = 04-phy_alignments_w_spaces/xxx.phy     * Path to alignment file
treefile = 08-branch_site_model/01-tree_LGT_labelled/xxx.newick       * Path to tree file
outfile = 08-branch_site_model/03-LGT-output/xxx.txt        * Output file

noisy = 3                    * 0,1,2,3,9: how much rubbish on the screen
verbose = 1                  * 1: detailed output, 0: concise output
runmode = 0                  * 0: user tree; 1: semi-automatic; 2: automatic

seqtype = 1                  * 1:codons; 2:AAs; 3:codon-translated

CodonFreq = 2                * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
model = 2                    * 0: one omega; 1: nearly-neutral; 2: branch model

NSsites = 2                  * 0: one omega for all sites

fix_kappa = 0                * 1: kappa fixed; 0: kappa to be estimated
kappa = 2                    * initial or fixed kappa

fix_omega = 0                * 1: omega fixed; 0: omega to be estimated
omega = 1                    * initial omega

fix_alpha = 1                * 0: estimate gamma shape parameter; 1: no gamma
alpha = 0.0                  * initial alpha, 0:infinity (constant rate)

Malpha = 0                   * different alphas for genes
ncatG = 10                   * # of categories in dG of NSsites models

getSE = 0                    * 0: dont want them, 1: want S.E.s of estimates
RateAncestor = 0             * (1/0): rates (alpha>0) or ancestral states (alpha=0)
Small_Diff = 5e-7            * Default value.
cleandata = 1                * Remove sites with ambiguity data (1) or not (0)
method = 0                   * 0: simultaneous; 1: one branch at a time

# create this as "template.txt"
nano 08-branch_site_model/LGT_branch_side_model_template.txt

# here's a script for making the control files

##
#!/bin/bash

template="08-branch_site_model/LGT_branch_side_model_template.txt"

# Read each string in list.txt and create a new file
while read -r line; do
    # Define the new file name using the string
    newfile="${line}.ctl"
    
    # Replace 'xxx' in the template and save to the new file
    sed "s/xxx/$line/g" "$template" > 08-branch_site_model/02-LGT-control_files/"$newfile"
done < list_of_files.txt

###################################
# create this as "make_control_files.sh"
nano 08-branch_site_model/LGT_make_control_files.sh
sbatch 08-branch_site_model/LGT_make_control_files.sh

# and a script for running paml
##
#!/bin/bash

source /users/bop22cfc/miniconda3/bin/activate
source activate /users/bop22cfc/miniconda3/envs/paml

# Run codeml on all files in the 07-control_files directory
for file in 08-branch_site_model/02-LGT-control_files/*; do
  codeml "$file"
done
##
# make this into a script called 'run_paml.sh'
nano 08-branch_site_model/LGT_run_paml.sh
sbatch 08-branch_site_model/LGT_run_paml.sh

#################
## and now native
mkdir 08-branch_site_model/04-tree_native_labelled
ls -1 05-trees/ | while read line ; do cat 05-trees/${line} | sed 's/Native/Native#1/g' > 08-branch_site_model/04-tree_native_labelled/${line} ; done
mkdir 08-branch_site_model/05-native-control_files
mkdir 08-branch_site_model/06-native-output

# now make control files
seqfile = 04-phy_alignments_w_spaces/xxx.phy     * Path to alignment file
treefile = 08-branch_site_model/04-tree_native_labelled/xxx.newick       * Path to tree file
outfile = 08-branch_site_model/06-native-output/xxx.txt        * Output file

noisy = 3                    * 0,1,2,3,9: how much rubbish on the screen
verbose = 1                  * 1: detailed output, 0: concise output
runmode = 0                  * 0: user tree; 1: semi-automatic; 2: automatic

seqtype = 1                  * 1:codons; 2:AAs; 3:codon-translated

CodonFreq = 2                * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
model = 2                    * 0: one omega; 1: nearly-neutral; 2: branch model

NSsites = 2                  * 0: one omega for all sites

fix_kappa = 0                * 1: kappa fixed; 0: kappa to be estimated
kappa = 2                    * initial or fixed kappa

fix_omega = 0                * 1: omega fixed; 0: omega to be estimated
omega = 1                    * initial omega

fix_alpha = 1                * 0: estimate gamma shape parameter; 1: no gamma
alpha = 0.0                  * initial alpha, 0:infinity (constant rate)

Malpha = 0                   * different alphas for genes
ncatG = 10                   * # of categories in dG of NSsites models

getSE = 0                    * 0: dont want them, 1: want S.E.s of estimates
RateAncestor = 0             * (1/0): rates (alpha>0) or ancestral states (alpha=0)
Small_Diff = 5e-7            * Default value.
cleandata = 1                * Remove sites with ambiguity data (1) or not (0)
method = 0                   * 0: simultaneous; 1: one branch at a time

# create this as "template.txt"
nano 08-branch_site_model/native_branch_side_model_template.txt

# here's a script for making the control files

##
#!/bin/bash

template="08-branch_site_model/native_branch_side_model_template.txt"

# Read each string in list.txt and create a new file
while read -r line; do
    # Define the new file name using the string
    newfile="${line}.ctl"
    
    # Replace 'xxx' in the template and save to the new file
    sed "s/xxx/$line/g" "$template" > 08-branch_site_model/05-native-control_files/"$newfile"
done < list_of_files.txt

###################################
# create this as "make_control_files.sh"
nano 08-branch_site_model/native_make_control_files.sh
sbatch 08-branch_site_model/native_make_control_files.sh

# and a script for running paml
##
#!/bin/bash

source /users/bop22cfc/miniconda3/bin/activate
source activate /users/bop22cfc/miniconda3/envs/paml

# Run codeml on all files in the 07-control_files directory
for file in 08-branch_site_model/05-native-control_files/*; do
  codeml "$file"
done
##
# make this into a script called 'run_paml.sh'
nano 08-branch_site_model/native_run_paml.sh
sbatch 08-branch_site_model/native_run_paml.sh
####################
# now running the null branch-site model
# model=2 and NSsites=2
# fix_omega = 1

## first for LGT genes
mkdir 08-branch_site_model/07-LGT-null-control_files
mkdir 08-branch_site_model/08-LGT-null-output

# now make control files
seqfile = 04-phy_alignments_w_spaces/xxx.phy     * Path to alignment file
treefile = 08-branch_site_model/01-tree_LGT_labelled/xxx.newick       * Path to tree file
outfile = 08-branch_site_model/08-LGT-null-output/xxx.txt        * Output file

noisy = 3                    * 0,1,2,3,9: how much rubbish on the screen
verbose = 1                  * 1: detailed output, 0: concise output
runmode = 0                  * 0: user tree; 1: semi-automatic; 2: automatic

seqtype = 1                  * 1:codons; 2:AAs; 3:codon-translated

CodonFreq = 2                * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
model = 2                    * 0: one omega; 1: nearly-neutral; 2: branch model

NSsites = 2                  * 0: one omega for all sites

fix_kappa = 0                * 1: kappa fixed; 0: kappa to be estimated
kappa = 2                    * initial or fixed kappa

fix_omega = 1                * 1: omega fixed; 0: omega to be estimated
omega = 1                    * initial omega

fix_alpha = 1                * 0: estimate gamma shape parameter; 1: no gamma
alpha = 0.0                  * initial alpha, 0:infinity (constant rate)

Malpha = 0                   * different alphas for genes
ncatG = 10                   * # of categories in dG of NSsites models

getSE = 0                    * 0: dont want them, 1: want S.E.s of estimates
RateAncestor = 0             * (1/0): rates (alpha>0) or ancestral states (alpha=0)
Small_Diff = 5e-7            * Default value.
cleandata = 1                * Remove sites with ambiguity data (1) or not (0)
method = 0                   * 0: simultaneous; 1: one branch at a time

# create this as "template.txt"
nano 08-branch_site_model/LGT_branch_site_model_null_template.txt

# here's a script for making the control files

##
#!/bin/bash

template="08-branch_site_model/LGT_branch_site_model_null_template.txt"

# Read each string in list.txt and create a new file
while read -r line; do
    # Define the new file name using the string
    newfile="${line}.ctl"
    
    # Replace 'xxx' in the template and save to the new file
    sed "s/xxx/$line/g" "$template" > 08-branch_site_model/07-LGT-null-control_files/"$newfile"
done < list_of_files.txt

###################################
# create this as "make_control_files.sh"
nano 08-branch_site_model/LGT_make_null_control_files.sh
sbatch 08-branch_site_model/LGT_make_null_control_files.sh

# and a script for running paml
##
#!/bin/bash

#SBATCH --job-name=codeml_array      # Job name
#SBATCH --output=logs/codeml_%A_%a.out  # Standard output and error log
#SBATCH --error=logs/codeml_%A_%a.err   # Error log
#SBATCH --array=0-115%10               # Array range; replace <N-1> with the number of files minus one
#SBATCH --ntasks=1                    # Number of tasks (1 per job)
#SBATCH --time=01:00:00               # Time limit hrs:min:sec
#SBATCH --mem=1G                      # Memory per CPU core

# Activate environment
source /users/bop22cfc/miniconda3/bin/activate
source activate /users/bop22cfc/miniconda3/envs/paml

# List all files in the directory
files=(08-branch_site_model/07-LGT-null-control_files/*)

# Select the file based on the SLURM_ARRAY_TASK_ID
file=${files[$SLURM_ARRAY_TASK_ID]}

# Run codeml on the selected file
codeml "$file"

##
# make this into a script called 'run_paml.sh'
nano 08-branch_site_model/LGT_null_run_paml.sh
sbatch 08-branch_site_model/LGT_null_run_paml.sh

#################
## and now native null
mkdir 08-branch_site_model/09-native-null-control_files
mkdir 08-branch_site_model/10-native-null-output

# now make control files
seqfile = 04-phy_alignments_w_spaces/xxx.phy     * Path to alignment file
treefile = 08-branch_site_model/04-tree_native_labelled/xxx.newick       * Path to tree file
outfile = 08-branch_site_model/10-native-null-output/xxx.txt        * Output file

noisy = 3                    * 0,1,2,3,9: how much rubbish on the screen
verbose = 1                  * 1: detailed output, 0: concise output
runmode = 0                  * 0: user tree; 1: semi-automatic; 2: automatic

seqtype = 1                  * 1:codons; 2:AAs; 3:codon-translated

CodonFreq = 2                * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
model = 2                    * 0: one omega; 1: nearly-neutral; 2: branch model

NSsites = 2                  * 0: one omega for all sites

fix_kappa = 0                * 1: kappa fixed; 0: kappa to be estimated
kappa = 2                    * initial or fixed kappa

fix_omega = 1                * 1: omega fixed; 0: omega to be estimated
omega = 1                    * initial omega

fix_alpha = 1                * 0: estimate gamma shape parameter; 1: no gamma
alpha = 0.0                  * initial alpha, 0:infinity (constant rate)

Malpha = 0                   * different alphas for genes
ncatG = 10                   * # of categories in dG of NSsites models

getSE = 0                    * 0: dont want them, 1: want S.E.s of estimates
RateAncestor = 0             * (1/0): rates (alpha>0) or ancestral states (alpha=0)
Small_Diff = 5e-7            * Default value.
cleandata = 1                * Remove sites with ambiguity data (1) or not (0)
method = 0                   * 0: simultaneous; 1: one branch at a time

# create this as "template.txt"
nano 08-branch_site_model/native_branch_site_null_model_template.txt

# here's a script for making the control files

##
#!/bin/bash

template="08-branch_site_model/native_branch_site_null_model_template.txt"

# Read each string in list.txt and create a new file
while read -r line; do
    # Define the new file name using the string
    newfile="${line}.ctl"
    
    # Replace 'xxx' in the template and save to the new file
    sed "s/xxx/$line/g" "$template" > 08-branch_site_model/09-native-null-control_files/"$newfile"
done < list_of_files.txt

###################################
# create this as "make_control_files.sh"
nano 08-branch_site_model/native_make_null_control_files.sh
sbatch 08-branch_site_model/native_make_null_control_files.sh

# and a script for running paml
##
#!/bin/bash

#SBATCH --job-name=codeml_array      # Job name
#SBATCH --output=logs3/codeml_%A_%a.out  # Standard output and error log
#SBATCH --error=logs3/codeml_%A_%a.err   # Error log
#SBATCH --array=0-46%10               # Array range; replace <N-1> with the number of files minus one
#SBATCH --ntasks=1                    # Number of tasks (1 per job)
#SBATCH --time=01:00:00               # Time limit hrs:min:sec
#SBATCH --mem=1G                      # Memory per CPU core

# Activate environment
source /users/bop22cfc/miniconda3/bin/activate
source activate /users/bop22cfc/miniconda3/envs/paml

# List all files in the directory
files=(08-branch_site_model/09-native-null-control_files/*)

# Select the file based on the SLURM_ARRAY_TASK_ID
file=${files[$SLURM_ARRAY_TASK_ID]}

# Run codeml on the selected file
codeml "$file"
##
# make this into a script called 'run_paml.sh'
nano 08-branch_site_model/native_null_run_paml.sh
sbatch 08-branch_site_model/native_null_run_paml.sh



###############

# fixing broken alignments / runs 
mkdir 09-fixing_errors
## fixing things
broken="LGT-020_ASEM_L04B_06933"
broken="LGT-021_ASEM_ZAM15-05-10_29632"
broken="LGT-059_ASEM_ZAM15-05-10_02275"
broken="LGT-114_ASEM_ZAM15-05-10_05026"

# for M0
cat 03-phy_alignments/${broken}.phy | sed 's/ /  /g' > 04-phy_alignments_w_spaces/${broken}.phy
# for branch model
# need to add a #1 to LGT and #2 to Native for each tree
cat 05-trees/${broken}.newick | sed 's/LGT/LGT#1/g' | sed 's/Native/Native#2/g' > 07-branch_model/01-tree_LGT_labelled/${broken}.newick 

# for branch site model - LGT
cat 05-trees/${broken}.newick | sed 's/LGT/LGT#1/g' > 08-branch_site_model/01-tree_LGT_labelled/${broken}.newick

# for branch side model - native
cat 05-trees/${broken}.newick | sed 's/Native/Native#1/g' > 08-branch_site_model/04-tree_native_labelled/${broken}.newick

#####################
## compiling 

mkdir 10-summary_stats
# M0 omegas
grep "omega (dN/dS) = " 06-M0_model/02-output/* | sed 's/06-M0_model\/02-output\///g' | sed 's/.txt:omega (dN\/dS) = //g' | sed 's/ /\t/g' | sed 's/_ASEM/ ASEM/g' | sed 's/_AANG/ AANG/g' > 10-summary_stats/M0_omegas.txt

# M0 lnL
grep "lnL" 06-M0_model/02-output/* | sed 's/06-M0_model\/02-output\///g' | sed 's/.txt:/ /g' | sed 's/      +0.000000//g' | sed 's/_ASEM/ ASEM/g' | sed 's/_AANG/ AANG/g' | sed 's/lnL(ntime: //g' | sed 's/ np: //g' | sed 's/):   \|):  / /g' > 10-summary_stats/M0_lnL.txt

# branch model omegas
grep "w (dN/dS) for branches:" 07-branch_model/03-output/* | sed 's/07-branch_model\/03-output\///g' | sed 's/.txt:w (dN\/dS) for branches: //g' | sed 's/ /\t/g' | sed 's/_ASEM/ ASEM/g' | sed 's/_AANG/ AANG/g' > 10-summary_stats/branch_model_omegas.txt

# branch model lnL
grep "lnL" 07-branch_model/03-output/* | sed 's/07-branch_model\/03-output\///g' | sed 's/.txt:/ /g' | sed 's/      +0.000000//g' | sed 's/_ASEM/ ASEM/g' | sed 's/_AANG/ AANG/g' | sed 's/lnL(ntime: //g' | sed 's/ np: //g' | sed 's/):   \|):  / /g' > 10-summary_stats/branch_model_lnL.txt

# branch model LGT only omegas
grep "w (dN/dS) for branches:" 07-branch_model/06-LGT-output/* | sed 's/07-branch_model\/06-LGT-output\///g' | sed 's/.txt:w (dN\/dS) for branches: //g' | sed 's/ /\t/g' | sed 's/_ASEM/ ASEM/g' | sed 's/_AANG/ AANG/g' > 10-summary_stats/branch_model_LGT_only_omegas.txt

# branch model LGT only lnL
grep "lnL" 07-branch_model/06-LGT-output/* | sed 's/07-branch_model\/06-LGT-output\///g' | sed 's/.txt:/ /g' | sed 's/      +0.000000//g' | sed 's/_ASEM/ ASEM/g' | sed 's/_AANG/ AANG/g' | sed 's/lnL(ntime: //g' | sed 's/ np: //g' | sed 's/):   \|):  / /g' > 10-summary_stats/branch_model_LGT_only_lnL.txt


# branch site model LGT results
grep -A5 "MLEs of dN/dS (w) for site classes (K=4)" 08-branch_site_model/03-LGT-output/* | sed 's/08-branch_site_model\/03-LGT-output\///g' | sed 's/.txt-/ /g' | sed 's/.txt:MLEs of dN\/dS (w) for site classes (K=4)//g' | sed 's/site class       /site_class /g' | sed 's/proportion       /proportion /g' | sed 's/background w     /background_w /g' | sed 's/foreground w     /foreground_w /g' | sed 's/      0/class_0/g' | sed 's/      1/class_1/g' | sed 's/      2a/class_2a/g' | sed 's/     2b/class_2b/g' | sed 's/  / /g' | sed 's/ /\t/g' | awk 'NR % 7 >= 4 && NR % 7 <= 6' | sed 's/_ASEM/ ASEM/g' | sed 's/_AANG/ AANG/g' > 10-summary_stats/branch_site_LGT_results.txt 

# branch site model native results
grep -A5 "MLEs of dN/dS (w) for site classes (K=4)" 08-branch_site_model/06-native-output/* | sed 's/08-branch_site_model\/06-native-output\///g' | sed 's/.txt-/ /g' | sed 's/.txt:MLEs of dN\/dS (w) for site classes (K=4)//g' | sed 's/site class       /site_class /g' | sed 's/proportion       /proportion /g' | sed 's/background w     /background_w /g' | sed 's/foreground w     /foreground_w /g' | sed 's/      0/class_0/g' | sed 's/      1/class_1/g' | sed 's/      2a/class_2a/g' | sed 's/     2b/class_2b/g' | sed 's/  / /g' | sed 's/ /\t/g' | awk 'NR % 7 >= 4 && NR % 7 <= 6' | sed 's/_ASEM/ ASEM/g' | sed 's/_AANG/ AANG/g' > 10-summary_stats/branch_site_native_results.txt

# branch site model LGT lnL
grep "lnL" 08-branch_site_model/03-LGT-output/* | sed 's/08-branch_site_model\/03-LGT-output\///g' | sed 's/.txt:/ /g' | sed 's/_ASEM/ ASEM/g' | sed 's/_AANG/ AANG/g' | sed 's/ lnL(ntime://g' | sed 's/ np: //g' | sed 's/): //g' | sed 's/      +0.000000//g' | sed 's/  / /g' > 10-summary_stats/branch_site_LGT_lnL.txt

# branch site model native lnL
grep "lnL" 08-branch_site_model/06-native-output/* | sed 's/08-branch_site_model\/06-native-output\///g' | sed 's/.txt:/ /g' | sed 's/_ASEM/ ASEM/g' | sed 's/_AANG/ AANG/g' | sed 's/ lnL(ntime://g' | sed 's/ np: //g' | sed 's/): //g' | sed 's/      +0.000000//g' | sed 's/  / /g' > 10-summary_stats/branch_site_native_lnL.txt

# branch site model LGT null results
grep -A5 "MLEs of dN/dS (w) for site classes (K=4)" 08-branch_site_model/08-LGT-null-output/* | sed 's/08-branch_site_model\/08-LGT-null-output\///g' | sed 's/.txt-/ /g' | sed 's/.txt:MLEs of dN\/dS (w) for site classes (K=4)//g' | sed 's/site class       /site_class /g' | sed 's/proportion       /proportion /g' | sed 's/background w     /background_w /g' | sed 's/foreground w     /foreground_w /g' | sed 's/      0/class_0/g' | sed 's/      1/class_1/g' | sed 's/      2a/class_2a/g' | sed 's/     2b/class_2b/g' | sed 's/  / /g' | sed 's/ /\t/g' | awk 'NR % 7 >= 4 && NR % 7 <= 6' | sed 's/_ASEM/ ASEM/g' | sed 's/_AANG/ AANG/g' > 10-summary_stats/branch_site_LGT_null_results.txt

# branch site model native null results
grep -A5 "MLEs of dN/dS (w) for site classes (K=4)" 08-branch_site_model/10-native-null-output/* | sed 's/08-branch_site_model\/10-native-null-output\///g' | sed 's/.txt-/ /g' | sed 's/.txt:MLEs of dN\/dS (w) for site classes (K=4)//g' | sed 's/site class       /site_class /g' | sed 's/proportion       /proportion /g' | sed 's/background w     /background_w /g' | sed 's/foreground w     /foreground_w /g' | sed 's/      0/class_0/g' | sed 's/      1/class_1/g' | sed 's/      2a/class_2a/g' | sed 's/     2b/class_2b/g' | sed 's/  / /g' | sed 's/ /\t/g' | awk 'NR % 7 >= 4 && NR % 7 <= 6' | sed 's/_ASEM/ ASEM/g' | sed 's/_AANG/ AANG/g' > 10-summary_stats/branch_site_native_null_results.txt

# branch site model LGT null lnL
grep "lnL" 08-branch_site_model/08-LGT-null-output/* | sed 's/08-branch_site_model\/08-LGT-null-output\///g' | sed 's/.txt:/ /g' | sed 's/_ASEM/ ASEM/g' | sed 's/_AANG/ AANG/g' | sed 's/ lnL(ntime://g' | sed 's/ np: //g' | sed 's/): //g' | sed 's/      +0.000000//g' | sed 's/  / /g' > 10-summary_stats/branch_site_LGT_null_lnL.txt

# branch site model native lnL
grep "lnL" 08-branch_site_model/10-native-null-output/* | sed 's/08-branch_site_model\/10-native-null-output\///g' | sed 's/.txt:/ /g' | sed 's/_ASEM/ ASEM/g' | sed 's/_AANG/ AANG/g' | sed 's/ lnL(ntime://g' | sed 's/ np: //g' | sed 's/): //g' | sed 's/      +0.000000//g' | sed 's/  / /g' > 10-summary_stats/branch_site_native_null_lnL.txt






