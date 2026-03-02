#!/bin/bash

template="06-M0_model/M0_template.txt"

# Read each string in list.txt and create a new file
while read -r line; do
    # Define the new file name using the string
    newfile="${line}.ctl"
    
    # Replace 'xxx' in the template and save to the new file
    sed "s/xxx/$line/g" "$template" > 06-M0_model/01-control_files/"$newfile"
done < list_of_files.txt
