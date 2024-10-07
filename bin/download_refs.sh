#!/bin/bash

ftp_table=$1

# Assuming the file name is ftp_table.tsv
awk -F'\t' 'BEGIN {OFS=FS} {gsub(/ /, "_", $2); print}' $1 |
while IFS=$'\t' read -r taxid org_name url; do
    rsync -avP "$url" ref.fa.gz
    gzip -d ref.fa.gz && mv ref.fa "${org_name}.fasta"

    # Edit fna download link to gff
    modified_url=$(echo "$url" | sed 's/\.fna\.gz/\.gff\.gz/')
    rsync -avP "$modified_url" ref.gff.gz
    gzip -d ref.gff.gz && mv ref.gff "${org_name}.gff"
done