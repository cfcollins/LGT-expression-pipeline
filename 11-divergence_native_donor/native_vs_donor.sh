#### comparing dN for donor vs native
#!/bin/bash

# Input files
GENOME=Sorghum_bicolor.Sorghum_bicolor_NCBIv3.cds.all.fa
GENE_LIST=SOR_donors

# Make an output directory to keep things tidy
mkdir -p extracted_genes

# Read gene list line by line
while IFS=$'\t' read -r GROUP GENE; do
    # Output file name
    OUTFILE="extracted_genes/${GROUP}_${GENE}.fa"

    # Use awk to extract sequence from FASTA
    awk -v gene="$GENE" -v outfile="$OUTFILE" '
        BEGIN { found=0 }
        /^>/ {
            if (found) exit;
            found=index($0, gene) > 0;
        }
        found { print >> outfile }
    ' "$GENOME"

done < "$GENE_LIST"



for prefix in $(ls LGT-*.fa | cut -d'_' -f1 | sort -u); do
  cat ${prefix}_*.fa > new_files/${prefix}_AUS1.fa
done


###
# remove " alignment by translation"
for f in *" alignment by translation"*; do
  mv -- "$f" "${f// alignment by translation/}"
done

# add spaces
ls -1 finished_alignments/ | while read line ; do cat finished_alignments/${line} | sed 's/ASEM_ZAM/ZAM  /g' | sed 's/Sorghum_bi/Sorghum_bi  /g' | sed 's/Alloterops/AUS1  /g' | sed 's/Setaria_it/Setaria_it  /g' | sed 's/ASEM_KWT3_/KWT  /g' | sed 's/ASEM_L04B_/L04B  /g' | sed 's/AANG_MRL48/MRL  /g' > 02-fixed_alignments/${line}  ; done


seqfile = lgt_nat_aligned1.phy
outfile = yn00_output.txt
verbose = 1

##### make control files
#!/bin/bash

# Folder containing the .phy files
INPUT_DIR="02-fixed_alignments"

# Output folder for control files
CONTROL_DIR="03-control_files"
mkdir -p "$CONTROL_DIR"

# Loop through each .phy file
for filepath in "$INPUT_DIR"/*.phy; do
    # Get the base filename without the path and extension
    filename=$(basename "$filepath" .phy)

    # Create control file content
    cat <<EOF > "$CONTROL_DIR/${filename}.ctl"
seqfile = /mnt/parscratch/users/bop22cfc/native_vs_donor/02-fixed_alignments/${filename}.phy
outfile = /mnt/parscratch/users/bop22cfc/native_vs_donor/04-output/${filename}.txt
verbose = 1
EOF

done


ls -1 03-control_files/ | while read line ; do yn00 03-control_files/${line} ; done

