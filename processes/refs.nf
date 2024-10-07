// Download the virus and bacteria RefSeq metadata tables
process GET_ASSEMBLY_META{
    label 'online'

    output:
    path('viral_assembly_summary.txt'), emit: viral
    path('bacteria_assembly_summary.txt'), emit: bacteria

    script:
    """
    wget ftp://ftp.ncbi.nih.gov/genomes/refseq/viral/assembly_summary.txt
    mv assembly_summary.txt viral_assembly_summary.txt
    wget ftp://ftp.ncbi.nih.gov/genomes/refseq/bacteria/assembly_summary.txt
    mv assembly_summary.txt bacteria_assembly_summary.txt
    """
    stub:
    """
    touch viral_assembly_summary.txt bacteria_assembly_summary.txt
    """
}

process GET_ASSEMBLY_META_GENBANK{
    label 'online'

    output:
    path('viral_assembly_summary.txt'), emit: viral
    path('bacteria_assembly_summary.txt'), emit: bacteria

    script:
    """
    wget ftp://ftp.ncbi.nih.gov/genomes/genbank/viral/assembly_summary.txt
    mv assembly_summary.txt viral_assembly_summary.txt
    wget ftp://ftp.ncbi.nih.gov/genomes/genbank/bacteria/assembly_summary.txt
    mv assembly_summary.txt bacteria_assembly_summary.txt
    """
    stub:
    """
    touch viral_assembly_summary.txt bacteria_assembly_summary.txt
    """
}

// Use the list of taxids to create a table with columns: taxid, organism name, reference ftp link
// This process creates a row for each taxid; the pre_main.nf merges and publishes the resulting table
// process CREATE_FTP_TABLE{
//     conda "$params.envs/biopython"

//     cpus 1

//     input:
//     //tuple val(taxid), path('viral_assembly_summary.txt'), path('bacteria_assembly_summary.txt')
//     tuple val(taxid), path('viral_assembly_summary.txt'), path('bacteria_assembly_summary.txt'), \
//         path('viral_assembly_summary_genbank.txt'), path('bacteria_assembly_summary_genbank.txt')

//     output:
//     tuple val(taxid), path("${taxid}_row.tsv"), emit: tsv
//     tuple val(taxid), path("${taxid}_row.csv"), emit: csv

//     script:
//     """
//         # Extract metadata for taxid
//         results=\$(getRef.py -t $taxid -v viral_assembly_summary.txt -b bacteria_assembly_summary.txt)
//         org_name=\$(echo "\$results" | head -n 1)
//         url=\$(echo "\$results" | tail -n 1)

//         # Check if either url or org_name is empty
//         if [ -z "\$url" ] || [ -z "\$org_name" ]; then
//             echo "No data pulled from RefSeq - trying GenBank...
//             results=\$(getRef.py -t $taxid -v viral_assembly_summary_genbank.txt -b bacteria_assembly_summary_genbank.txt)
//             org_name=\$(echo "\$results" | head -n 1)
//             url=\$(echo "\$results" | tail -n 1)
//         else
//             echo "Data found in RefSeq Database."
//         fi

//         echo "Organism name: \$org_name"
//         echo "URL: \$url"

//         # Edit the organism names - we will rename the files using this later
//         org_name=\$(edit_org_names.sh "\$org_name")

//         # Print row to a tsv file
//         echo "$taxid \$org_name \$url"
//         echo "$taxid \$org_name \$url" > ${taxid}_row.tsv

//         tr ' ' ',' < ${taxid}_row.tsv > ${taxid}_row.csv
//     """

//     stub:
//     """
//     echo "$taxid orgname ftp-path"
//     """
// }

process CREATE_FTP_TABLE{
    conda "$params.envs/biopython"

    cpus 1

    input:
    tuple val(taxid), path('viral_assembly_summary.txt'), path('bacteria_assembly_summary.txt')

    output:
    tuple val(taxid), path("${taxid}_row.tsv"), emit: tsv
    tuple val(taxid), path("${taxid}_row.csv"), emit: csv

    script:
    """
        # Extract metadata for taxid
        results=\$(getRef.py -t $taxid -v viral_assembly_summary.txt -b bacteria_assembly_summary.txt)
        org_name=\$(echo "\$results" | head -n 1)
        url=\$(echo "\$results" | tail -n 1)

        echo $taxid
        echo \$results
        echo \$url
        echo \$org_name

        # Edit the organism names - we will rename the files using this later
        org_name=\$(edit_org_names.sh "\$org_name")

        # Print row to a tsv file
        echo "$taxid \$org_name \$url" > ${taxid}_row.tsv

        #Create csv
        tr ' ' ',' < ${taxid}_row.tsv > ${taxid}_row.csv
    """

    stub:
    """
    echo "$taxid orgname ftp-path"
    """
}

// Use ftp table to download references in parallel
process DOWNLOAD_REFS{
    conda "$params.envs/biopython"
    cpus 4

    label 'online'

    input:
    tuple val(taxid), path('taxid_row')

    output:
    tuple val(taxid), file('*.fasta'), emit: fasta
    tuple val(taxid), file('*.gff'), emit: gff

    script:
    """
    #IFS=\$'\t' read -r taxid org_name url <<< "$taxid_row"
    #echo "\$org_name \$taxid_row"
    download_ref.sh $taxid_row 
    """
    stub:
    """
    touch "$params.refs/refs/virus1.fasta"
    touch "$params.refs/refs/virus2.fasta"
    touch "$params.refs/refs/virus1.gff"
    touch "$params.refs/refs/virus2.gff"
    """
}

process EDIT_FASTA_HEADERS{
    conda "$params.envs/biopython"
    cpus 4

    publishDir "${params.refs}/refs", mode: 'copy'

    input:
    //tuple val(taxid), file('references.fasta')
    //tuple val(taxid), file('ftp_table_row')
    tuple val(taxid), file('references.fasta'), file('ftp_table_row')

    output:
    tuple val(taxid), file('*.fasta'), emit: fasta

    script:
    """
    cat ftp_table_row
    org_name=\$(awk '{print \$2}' ftp_table_row)
    org_name=\$org_name.fasta
    echo \$org_name
    edit_fasta_headers.py references.fasta \$org_name
    """

    stub:
    """
    touch "$params.refs/refs/virus1.fasta"
    touch "$params.refs/refs/virus2.fasta"
    """
}

process MERGE_REFS{
    conda "$params.envs/biopython"
    cpus 4

    publishDir "${params.refs}", mode: 'copy'

    input:
    path('*')

    output:
    file('meta_ref.fasta')

    script:
    """
    cat *.fasta > meta_ref.fasta
    #merge_refs.sh *.fasta
    cat $params.refs/spikes/*.fasta >> meta_ref.fasta
    cat $params.refs/other_refs/*.fasta >> meta_ref.fasta
    """

    stub:
    """
    touch "$params.refs/meta_ref.fasta
    """
}