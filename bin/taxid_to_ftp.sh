#!/bin/bash

# Input file containing taxon numbers
input_file=$1

# Default output file
output_file="ftp_table.tsv"

# Output file if name is declared
if [ $# -ge 2 ]; then
    output_file=$2
fi

binDir=$3

# Loop through each taxon number and process
while IFS= read -r taxid; do
    #echo $taxid
    #python getRef.py -t "$taxid" -v viral_assembly_summary.txt -b bacteria_assembly_summary.txt
    output=$(python "$binDir"/getRef.py -t "$taxid" -v viral_assembly_summary.txt -b bacteria_assembly_summary.txt)
    
    # Extract the first line as the organism name
    org_name=$(echo "$output" | awk 'NR==1{print $0}')

    # Extract the second line as the URL
    url=$(echo "$output" | awk 'NR==2{print $0}')

    #echo "taxid: $taxid"
    #echo "output: $output"
    #echo "orgname: $org_name"
    #echo "url: $url"    
    #echo -e "$taxid\t$url"
    echo -e "$taxid\t$org_name\t$url"
    echo -e "$taxid\t$org_name\t$url" >> "$output_file"
done < "$input_file"

# Edit the middle column of the new ftp table and save to a temp file
awk -F'\t' -v OFS='\t' '{ gsub(",","",$2); gsub("\\(|\\)","",$2); gsub("/","_",$2); gsub(" ","_",$2); gsub(",","_",$2); gsub(",","",$2); print }' "$output_file" > "$output_file.tmp"

# Replace the original file with the temporary file to overwrite it
mv "$output_file.tmp" "$output_file"