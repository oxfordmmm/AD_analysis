process FASTQ_TO_FASTA {
    conda "$params.envs/biopython"

    input:
    tuple val(barcode), val(split), path(fastq)

    output:
    tuple val(barcode), val(split), path("${barcode}.fasta")

    script:
    """
    seqtk seq -a $fastq > ${barcode}.fasta
    """
}

/*
Uses blast to find all locations of the 22bp primer within the reads.
Returns a tsv where each row is a primer hit
*/
process BLASTN {
    tag {barcode  + ' ' + split}
    //cpus 4

    conda "$params.envs/blast"
    publishDir "$params.output/primer_blast", mode: 'copy', saveAs: { filename -> "${barcode}_${split}_primer_blast.tsv" }

    input:
    tuple val(barcode), val(split), path(reads)
    path(query) // This is the sequence you want to search for

    output:
    tuple val(barcode), val(split), path("${barcode}_blast.tsv") 

    script:
    """
    makeblastdb -in ${reads} -title "read db" -dbtype nucl -out read_db -max_file_sz '3.5GB'
    blastn -task blastn-short -query $query -db read_db -word_size 7 -gapopen 2 -out _blastn.tsv -max_target_seqs 10000000 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen"
    echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tslen" > ${barcode}_blast.tsv && cat _blastn.tsv >> ${barcode}_blast.tsv
    """
}

process GRAPH_PRIMER_HITS {
    conda "$params.envs/r"
    publishDir "$outdir", mode: 'copy'

    input:
    tuple val(barcode), path(blast_hits), path(kaiju_out), path(k2_out)
    val(outdir)

    output:
    tuple val(barcode), path("*.png"), emit: 'graphs'

    script:
    """
    Rscript $params.bin/primer_graphs.R $blast_hits $kaiju_out $k2_out
    for g in *.png
    do
        mv \$g ${barcode}_\$g
    done
    """
}


/*
Process takes reads and splits into "pseudoreads" based on primer locations.

The script split_by_primer.py, will produce several outputs:
    - non_primer_reads.fastq: all reads with no primer
    - primer_reads.fastq: all reads with primer found in them
    - pseudoreads.fastq: the pseudoreads made from splitting the reads with primer
    - read_stats.csv: a table with info about every read
    - report.txt: a short summary of the splitting
*/
process PHI_SPLIT_READS_BY_PRIMER {
    tag {barcode  + ' ' + split}

    conda "$params.envs/biopython"
    cpus = 4

    input:
    tuple val(barcode), val(split), path(reads), path(primer_blast)

    output:
    tuple val(barcode), val(split), path("${split}_split_reads.fastq"), emit: 'fastq' // This is non-primer reads and pseudoreads together
    tuple val(barcode), val(split), path("${split}_read_stats.csv"), emit: 'read_stats'

    script:
    """
    echo barcode $barcode 
    echo split $split
    #gzip -dc $reads > reads.fastq
    #python3 $params.bin/split_by_primer.py -r reads.fastq -b $primer_blast -o ${barcode} -m 30
    split_by_primer.py -r $reads -b $primer_blast -o ${split} -m 30
    
    #cat ${barcode}_non_primer_reads.fastq ${barcode}_pseudoreads.fastq > ${barcode}_split_reads.fastq
    cat ${split}_non_primer_reads.fastq ${split}_pseudoreads.fastq > ${split}_split_reads.fastq
    """
}

/*
Process will cut the primer off the end of all reads to produce trimmed reads.
If a read has no primer, it leaves it as is.
If a read has a primer in the middle then something is wrong and the read is failed
Output all fastq files so that we can view them
*/
process SISPA_TRIM_PRIMER {
    tag {barcode  + ' ' + split}
    conda "$params.envs/biopython"

    input:
    tuple val(barcode), val(split), path(reads), path(primer_blast)

    output:
    tuple val(barcode), val(split), path("${split}_trimmed.fastq"), emit: 'fastq'
    tuple val(barcode), val(split), path("${split}_removed.fastq"), emit: 'failed'
    tuple val(barcode), val(split), path("${split}_read_stats.csv"), emit: 'read_stats'

    script:
    """
    echo barcode $barcode
    trim_primer.py -r $reads -b $primer_blast -o ${split} -m 20
    """
}

process PHI_MERGE_SPLITTING_OUTPUTS{
    label 'online'

    publishDir "$params.output/splitting", mode: 'copy', pattern: "*.fastq.gz"
    publishDir "$params.output/alignment_info", mode: 'copy', pattern: '*.csv'
    
    input:
    tuple val(barcode), val(split), path('*')
    tuple val(barcode), val(split), path('*')

    output:
    tuple val(barcode), path("${barcode}_split_reads.fastq.gz"), emit: 'fastq' // This is non-primer reads and pseudoreads together
    tuple val(barcode), path("${barcode}_read_stats.csv"), emit: 'read_stats'

    script:
    """
    cat *.fastq > "${barcode}_split_reads.fastq"
    gzip "${barcode}_split_reads.fastq" 

    # cat *.txt > "${barcode}_report.txt"

    awk 'NR==1{print; next} FNR>1' *.csv > "${barcode}_read_stats.csv"
    """

}

process SISPA_MERGE_TRIMMING_OUTPUTS {
    tag {barcode}
    label 'online'

    publishDir "$params.output/trimmed_reads", mode: 'copy', pattern: "*.fastq.gz"
    publishDir "$params.output/alignment_info", mode: 'copy', pattern: '*.csv'
    
    input:
    tuple val(barcode), val(split), path('*')
    tuple val(barcode), val(split), path('*')
    tuple val(barcode), val(split), path('*')

    output:
    tuple val(barcode), path("${barcode}_trimmed.fastq.gz"), emit: 'fastq' // This is non-primer reads and pseudoreads together
    tuple val(barcode), path("${barcode}_removed.fastq.gz"), emit: 'removed'
    tuple val(barcode), path("${barcode}_read_stats.csv"), emit: 'read_stats'

    script:
    """
    cat *trimmed.fastq > "${barcode}_trimmed.fastq"
    cat *removed.fastq > "${barcode}_removed.fastq"
    gzip "${barcode}_trimmed.fastq"
    gzip "${barcode}_removed.fastq"
    awk 'NR==1{print; next} FNR>1' *.csv > "${barcode}_read_stats.csv"
    """

}

process PHI_SPLITTING_REPORT {
    conda "$params.envs/biopython"

    publishDir "$outdir", mode: 'copy'

    input:
    tuple val(barcode), path(read_stats_split)
    val(outdir)

    output:
    path("${barcode}_split_report.txt")

    script:
    """
    python3 $params.bin/count_stats.py $read_stats_split "${barcode}_split_report.txt"
    """
}

// This script needs writing
process SISPA_TRIMMING_REPORT {
    conda "$params.envs/biopython"
    tag {barcode}

    publishDir "$outdir", mode: 'copy'

    input:
    tuple val(barcode), path(read_stats_split)
    val(outdir)

    output:
    path("${barcode}_trim_report.txt")

    script:
    """
    python3 $params.bin/count_stats.py $read_stats_split "${barcode}_split_report.txt"
    """
}




