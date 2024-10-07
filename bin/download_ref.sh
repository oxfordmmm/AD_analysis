#/bin/bash

taxid_row=$1
#org_name=$1
#url=$2

    #IFS=$'\t' read -r taxid org_name url <<< "$taxid_row"
    org_name=$(cat $taxid_row | awk '{print $2}')
    url=$(cat $taxid_row | awk '{print $3}')

    echo $url

    # Download fasta
    rsync -avP "$url" ref.fa.gz
    gzip -d ref.fa.gz && mv ref.fa $org_name.fasta

    # Edit fna download link to gff
    modified_url=$(echo "$url" | sed 's/\.fna\.gz/\.gff\.gz/')

    # Download gff
    rsync -avP "$modified_url" ref.gff.gz
    gzip -d ref.gff.gz && mv ref.gff $org_name.gff