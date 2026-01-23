### Cleaning Count Matrices from Salmon ###
### Cat Collins - October 2023 ############

# load libraries
library(readr)
library(BiocManager)
library(tximport)

# generate file paths to the quant files
dirs <- list.files("./clean_salmon_matrices/old_data/align_to_own_ref/SET/") # change for the accession
quant_files <- list.files(
  "./clean_salmon_matrices/old_data/align_to_own_ref/SET/", # change for the accession
  pattern="quant.sf", 
  recursive=TRUE, 
  full.names=TRUE)
dirs
quant_files

# inspect the first of the files
#quants <- read_tsv(quant_files[1])
#head(quants)

# we can demonstrate how TPM works
#rpk <- quants$NumReads / quants$EffectiveLength
#scale_factor <- sum(rpk) / 1e6
#tpm <- rpk / scale_factor

# extract abundance data
txi <- tximport(
  quant_files,
  type="salmon", 
  txIn=FALSE, 
  txOut = FALSE, 
  geneIdCol = "Name", 
  abundanceCol = "TPM")
names(txi)
head(txi$abundance)       

# change the column headers
quant_files # names should be same order as quant_files, so list first
new_headers <- c("X32", "X41", "X44", "X45", "55") # change column headers
colnames(txi$abundance) <- new_headers
TPM_data <- txi$abundance

# convert dataframe into a csv
write.csv(TPM_data, file = "clean_salmon_matrices/old_data/align_to_own_ref/all_count_matrices/SET_old_own_ref_TPM_count_matrix.csv", quote=FALSE) # change name


