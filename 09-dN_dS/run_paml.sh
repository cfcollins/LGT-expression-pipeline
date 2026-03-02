#!/bin/bash

source /users/bop22cfc/miniconda3/bin/activate
source activate /users/bop22cfc/miniconda3/envs/paml

# Run codeml on all files in the 07-control_files directory
for file in 06-M0_model/01-control_files/*; do
  codeml "$file"
done
