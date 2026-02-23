### coding notes for blast

# Access CDS files 
AUS1_cds_file=/mnt/autofs/dunning_lab/Shared/reference_genomes/Alloteropsis_semialata_GCA_004135705.1/ASEM_AUS1_v1.0.CDS.fa
KWT_cds_file=/mnt/autofs/christin_lab1/genome_data/Assemberlies/REFERENCE/KWT/KWT_v1.0.all.maker.CDS.fasta
L04B_cds_file=/mnt/autofs/christin_lab1/genome_data/Assemberlies/REFERENCE/L04B/ASEM_L04B_v1.0.all.maker.cds.fasta
ZAM_cds_file=/mnt/autofs/christin_lab1/genome_data/Assemberlies/REFERENCE/ZAM-05-10/ASEM_ZAM15-05-10.v1.0.all.maker.cds.fasta
MRL_cds_file=/mnt/autofs/christin_lab1/genome_data/Assemberlies/REFERENCE/MRL/AANG_MRL48.all.maker.CDS.fasta

# combine CDS files 
cat $AUS1_cds_file $KWT_cds_file $L04B_cds_file $ZAM_cds_file $MRL_cds_file > alloteropsis_genomes.cds.fa

# unwrap all the files 
cat alloteropsis_genomes.cds.fa | awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' > corrected_alloteropsis_genomes.cds.fa

# access the lists of LGTs and natives and move to this folder 
# may want to check the lists have no secret characters
cat -A list_of_LGTs.txt
dos2unix list_of_natives.txt

# gain the sequence files 
cat list_of_LGTs.txt | while read line ; do grep "$line" -A1 corrected_alloteropsis_genomes.cds.fa >> LGT_native_sequences.fasta; done
cat list_of_natives.txt | while read line ; do grep "$line" -A1 corrected_alloteropsis_genomes.cds.fa >> LGT_native_sequences.fasta; done

# blast!

# load environments
source /users/bop22cfc/miniconda3/bin/activate
source activate /users/bop22cfc/miniconda3/envs/blast

# download ref genome CDS
mkdir ref_genome_cds
cd ref_genome_cds
wget http://ftp.ensemblgenomes.org/pub/release-58/plants/fasta/sorghum_bicolor/cds/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.cds.all.fa
wget http://ftp.ensemblgenomes.org/pub/release-58/plants/fasta/brachypodium_distachyon/cds/Brachypodium_distachyon.Brachypodium_distachyon_v3.0.cds.all.fa.gz
wget http://ftp.ensemblgenomes.org/pub/release-58/plants/fasta/setaria_italica/cds/Setaria_italica.Setaria_italica_v2.0.cds.all.fa.gz
wget http://ftp.ensemblgenomes.org/pub/release-58/plants/fasta/zea_mays/cds/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.cds.all.fa
wget http://ftp.ensemblgenomes.org/pub/release-58/plants/fasta/oryza_sativa/cds/Oryza_sativa.IRGSP-1.0.cds.all.fa.gz

# edit these so the gene names are neater
ref_location=/mnt/parscratch/users/bop22cfc/sequence_fun/previous/ref_genomes
cat ${ref_location}/Brachypodium_distachyon.v3.0.cds.all.fa | sed 's/>/>Brach-/g' > Brachypodium_edit.cds.fa
cat ${ref_location}/Oryza_sativa.IRGSP-1.0.cds.all.fa | sed 's/>/>Oryza-/g' > Oryza_edit.cds.fa
cat ${ref_location}/Setaria_italica_v2.0.cds.all.fa | sed 's/>/>Setaria-/g' > Setaria_edit.cds.fa
cat ${ref_location}/Sorghum_bicolor_NCBIv3.cds.all.fa | sed 's/>/>Sorghum-/g' > Sorghum_edit.cds.fa
cat ${ref_location}/Zea_mays.NAM-5.0.cds.all.fa | sed 's/>/>Zea-/g' > Zea_edit.cds.fa

# it's combining weirdly so :
cat Brachypodium_edit.cds.fa Zea_edit.cds.fa > half1.fa
cat Setaria_edit.cds.fa Oryza_edit.cds.fa Sorghum_edit.cds.fa > half2.fa
cat half1.fa half2.fa > ref_sequences.fa
# unwrap
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' ref_sequences.fa > ref_corrected_sequences.fa

# make blast database
makeblastdb -in ref_corrected_sequences.fa -out ref_db -parse_seqids -dbtype nucl
# move back up a folder
cd ../

# BLAST 
blastn -query LGT_native_sequences.fasta -db ref_genome_cds/ref_db -out results.out -outfmt 6

# filter results so only matches with over 300bp are included
cat results.out | awk '$4>=300' > above_300bp.out
cat results.out | awk '$4<300' > below_300bp.out

# split up 
mkdir split_up
mkdir split_up_top

# move the sequences from each species into separate folders
grep "Brach-" above_300bp.out > split_up/Brach.out
grep "Oryza-" above_300bp.out > split_up/Oryza.out
grep "Setaria-" above_300bp.out > split_up/Setaria.out
grep "Sorghum-" above_300bp.out > split_up/Sorghum.out
grep "Zea-" above_300bp.out > split_up/Zea.out

# now find the top results for each file 
ls -1 split_up/ | while read line ; do awk '!seen[$1]++' split_up/$line > split_up_top/top_$line ; done

# for each move
mkdir sequence_files
cat split_up_top/top_Brach.out | cut -f 1 | sort | uniq  | while read line ; do grep "$line" split_up_top/top_Brach.out | cut -f 2 | sort | uniq >> sequence_files/"$line" ; done
cat split_up_top/top_Oryza.out | cut -f 1 | sort | uniq  | while read line ; do grep "$line" split_up_top/top_Oryza.out | cut -f 2 | sort | uniq >> sequence_files/"$line" ; done
cat split_up_top/top_Setaria.out | cut -f 1 | sort | uniq  | while read line ; do grep "$line" split_up_top/top_Setaria.out | cut -f 2 | sort | uniq >> sequence_files/"$line" ; done
cat split_up_top/top_Sorghum.out | cut -f 1 | sort | uniq  | while read line ; do grep "$line" split_up_top/top_Sorghum.out | cut -f 2 | sort | uniq >> sequence_files/"$line" ; done
cat split_up_top/top_Zea.out | cut -f 1 | sort | uniq  | while read line ; do grep "$line" split_up_top/top_Zea.out | cut -f 2 | sort | uniq >> sequence_files/"$line" ; done

# now include the headers in these files 
for file in sequence_files/*; do if [ -f "$file" ]; then filename=$(basename "$file"); echo "$filename" > "$file"_with_header; cat "$file" >> "$file"_with_header; fi; done
# move into its own folder
mkdir sequence_files/with_header
mv sequence_files/*_with_header sequence_files/with_header/

# combine together the original sequence files so we can look at Alloteropsis and the refs at the same time
cat LGT_native_sequences.fasta ref_genome_cds/ref_corrected_sequences.fa > all_sequences.fa

# now add sequences to all of these
mkdir finished_sequences

for filename in sequence_files/with_header/*; do while read -r line; do grep "$line" -A 1 all_sequences.fa >> "finished_sequences/$(basename "$filename")" ; done < "$filename"; done

# the AUS1 files are poorly named so rename these
for file in Alloteropsis_semialata_GCA_004135705.1_*; do
    # Extract the part of the filename after the prefix
    new_name="${file#Alloteropsis_semialata_GCA_004135705.1_}"
    # Rename the file
    mv "$file" "$new_name"
    echo "Renamed '$file' to '$new_name'"
done





#### kept old code below in case useful


# it's combining weirdly so :
cat Brachypodium_distachyon.Brachypodium_distachyon_v3.0.cds.all.fa Zea_mays.Zm-B73-REFERENCE-NAM-5.0.cds.all.fa > half1.fa
cat Setaria_italica.Setaria_italica_v2.0.cds.all.fa Oryza_sativa.IRGSP-1.0.cds.all.fa Sorghum_bicolor.Sorghum_bicolor_NCBIv3.cds.all.fa > half2.fa
cat half1.fa half2.fa > ref_sequences.fa

awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' ref_sequences.fa > ref_corrected_sequences.fa

# make blast database
makeblastdb -in ref_corrected_sequences.fa -out ref_db -parse_seqids -dbtype nucl

# BLAST 
blastn -query LGT_native_sequences.fasta -db ref_genome_cds/ref_db -out results.out -outfmt 6

# filter results so only matches with over 300bp are included
cat results.out | awk '$4>=300' > above_300bp.out
cat results.out | awk '$4<300' > below_300bp.out

# work out which gene names correspond to which species 
# e.g. for this collection:
# EER, EES, KXG, OQU = Sorghum 
# KQJ, PNT = Brachypodium 
# Os = Oryza
# KQL = Setaria 
# Zm = Zea
# unfortunately KQK is both Brachypodium and Setaria so will have to be considered separately 

#mkdir split_up
#mkdir split_up_top
#grep -E "EER|EES|KXG|OQU" results.out > split_up/sorghum.out
#grep -E "KQJ|PNT" results.out > split_up/most_of_brachypodium.out
#grep "Os" results.out > split_up/oryza.out
#grep "KQL" results.out > split_up/most_of_setaria.out
#grep "Zm" results.out > split_up/zea.out
#grep "KQK" results.out > split_up/KQK.out

# now find the top results for each file 
#ls -1 split_up/ | while read line ; do awk '!seen[$1]++' split_up/$line > split_up_top/top_$line ; done

#mkdir sequence_files
#cat split_up_top/*.out | cut -f 2 | while read line ; do grep -A1 "$line" ref_corrected_sequences.fa > sequence_files/"$line".fa ; done
# this will likely end up with less files than the genes in the list in split_up_top as some genes are repeated
# you can check this with this:
# cat split_up_top/* | cut -f 2 | awk '!seen[$1]++' | wc -l

# take top results from blast
cat above_300bp.out | awk '!seen[$1]++' > top_results.out 

# convert into files 
mkdir top_results_sequence_files

cat top_results.out | cut -f 2 | while read line ; do grep -A1 "$line" ref_corrected_sequences.fa > top_results_sequence_files/"$line".fa ; done


cat above_300bp.out | cut -f 1 | sort | uniq  | while read line ; do grep "$line" above_300bp.out | cut -f 2 | sort | uniq >> crap/"$line" ; done
for filename in * ; do cat "$filename" | while read line ; grep "$line" -A 1 yourcdsfile >> "$filenname".fa 
