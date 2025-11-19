

params.fastqs = ''
params.meta_pathogens="$projectDir/../meta.csv"
params.method = 'sispa'

process FASTQ_TO_FASTA {
    //conda "$params.envs/biopython"
    tag {run + ' ' + barcode}

    input:
    tuple val(run), val(barcode), path(fastq), val(method)

    output:
    tuple val(run), val(barcode), path("${barcode}.fasta")

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
    tag {run + ' ' + barcode}
    //cpus 4

    //conda "$params.envs/blast"
    //publishDir "$params.output/primer_blast", mode: 'copy', saveAs: { filename -> "${barcode}_${split}_primer_blast.tsv" }

    input:
    tuple val(run), val(barcode), path(reads)
    path(query) // This is the sequence you want to search for

    output:
    tuple val(run), val(barcode), path("${barcode}_blast.tsv") 

    script:
    """
    makeblastdb -in ${reads} -title "read db" -dbtype nucl -out read_db -max_file_sz '3.5GB'
    blastn -task blastn-short -query $query -db read_db -word_size 7 -gapopen 2 -out _blastn.tsv -max_target_seqs 10000000 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen"
    echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tslen" > ${barcode}_blast.tsv && cat _blastn.tsv >> ${barcode}_blast.tsv
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
    tag {run + ' ' + barcode}

//    conda "$params.envs/biopython"
    cpus = 1

    input:
    tuple val(run), val(barcode), path(reads), path(primer_blast)

    output:
    tuple val(run), val(barcode), path("split_split_reads.fastq"), emit: 'fastq' // This is non-primer reads and pseudoreads together
    tuple val(run), val(barcode), path("split_read_stats.csv"), emit: 'read_stats'

    script:
    """
    #gzip -dc $reads > reads.fastq
    split_by_primer.py -r $reads -b $primer_blast -o split -m 30
    
    #cat ${barcode}_non_primer_reads.fastq ${barcode}_pseudoreads.fastq > ${barcode}_split_reads.fastq
    cat split_non_primer_reads.fastq split_pseudoreads.fastq > split_split_reads.fastq
    """
}

process CAT_FASTQS {
    tag {run + ' ' + barcode}
    //publishDir "fastqs/${run}/", mode: 'copy'

    //label 'cat'

    input:
    tuple val(run), val(barcode), path('???.fastq.gz')

    output:
    tuple val(run), val(barcode), path("${run}_${barcode}.fastq.gz")

    script:
    """
    cat *.fastq.gz > ${run}_${barcode}.fastq.gz
    """
}

process JOIN_FQS {
    tag {barcode}

    input:
    tuple val(barcode), file('R1.fq.gz'), file('R2.fq.gz')

    output:
    tuple val(barcode),val(barcode), path("${barcode}.fq.gz")

    script:
    """

    join_reads.py R1.fq.gz R2.fq.gz ${barcode}.fq.gz

    """
}

process KRAKEN2SERVER {
    tag {run + ' ' + barcode}
    publishDir "classifications/${run}/${barcode}/", mode: 'copy'
    maxForks 15


    input:
    tuple val(run), val(barcode), file('in.fastq.gz')

    output:
    tuple val(run), val(barcode), path("${run}_${barcode}_read_tax.txt"), emit: read_tax
    tuple val(run), val(barcode), path("${run}_${barcode}_report.txt"), emit: reports
    tuple val(run), val(barcode), path("${run}_${barcode}.txt"), emit: txt
    
    script:
    port=params.kraken2server_port
    host=params.kraken2server_host
    """
    kraken2_client -s in.fastq.gz \
        --port $port --host-ip $host \
        -r ${run}_${barcode}_report.txt \
        > ${run}_${barcode}.txt

    echo "readID taxID" > read_tax.txt
    awk '{print \$2,\$3}' ${run}_${barcode}.txt >> read_tax.txt
    mv read_tax.txt ${run}_${barcode}_read_tax.txt
    """
    stub:
    """
    touch kraken.txt.gz read_tax.txt "${barcode}_${batchID}_report.txt" "${barcode}_${batchID}.txt"
    """

}

process KRAKEN2{
    tag {run + ' ' + barcode}
    publishDir "classifications/${run}/${barcode}/", mode: 'copy', pattern:"*_report.txt"
    publishDir "reads/${run}/${barcode}/", mode: 'copy', pattern:"*_reads.txt"

    maxForks 6
    cpus 4


    input:
    tuple val(run), val(barcode), file('in.fastq.gz')

    output:
    //tuple val(run), val(barcode), path("${run}_${barcode}_read_tax.txt"), emit: read_tax
    tuple val(run), val(barcode), path("${run}_${barcode}_report.txt"), emit: reports
    tuple val(run), val(barcode), path("${run}_${barcode}_reads.txt"), emit: txt
    
    script:
    port=params.kraken2server_port
    host=params.kraken2server_host
    """
    kraken2 --db /mnt/nanostore/dbs/k2_eupathdb48_20230407/ \
        --report ${run}_${barcode}_report.txt \
        --threads $task.cpus \
        --report-minimizer-data \
        --gzip-compressed in.fastq.gz \
        > ${run}_${barcode}_reads.txt

    """
    stub:
    """
    touch kraken.txt.gz read_tax.txt "${barcode}_${batchID}_report.txt" "${barcode}_${batchID}.txt"
    """

}

process EXTRACT_RESULTS {
    tag {run + ' ' + barcode}
    publishDir "results/${run}/${barcode}/", mode: 'copy'

    label 'python'

    input:
    tuple val(run), val(barcode), path("${run}_${barcode}_report.txt"), path("seqkit_stats.tsv")

    output:
    path("${run}_${barcode}_report.tsv"), emit: tsv

    script:
    """
    extract_results.py -i ${run}_${barcode}_report.txt \
        -o ${run}_${barcode}_report.tsv \
        -t 10508 10509 \
            11118 11137 290028 277944 31631 2697049 \
            11308 11320 641809 211044 335341 11520 518987 \
            12058 12059 147711 147712 463676 138948 138949 138950 138951\
            11158 12730 2560525 11216 11224 162145 \
            2842319 12022 2946187 10886 11158 3052731 11050 64320 \
            351073 3049954 11250 \
            687329 \
            519 520 83558 2104 0\
        -r $run -b $barcode -s seqkit_stats.tsv
    """
}

process SEQKIT_STATS {
    tag {run + ' ' + barcode}
    publishDir "stats/${run}/", mode: 'copy'

    input:
    tuple val(run), val(barcode), path("${run}_${barcode}.gz")

    output:
    tuple val(run), val(barcode), path("${barcode}_seqkit_stats.txt")

    script:
    """
    seqkit stats -a -b -T ${run}_${barcode}.gz > ${barcode}_seqkit_stats.txt
    """
}

process COMPILE_RESULTS {
    publishDir "compiled_results/", mode: 'copy'

    label 'python'

    input:
    path "*.tsv"
    path('meta.csv')

    output:
    path "compiled_results.csv"

    script:
    """
    compile_results.py -i *.tsv -o compiled_results.csv -m meta.csv
    """
}


workflow {
    //fqs=Channel.fromPath("$params.fastqs/**/*.fastq.gz")
    //        .map {it ->
    //            tuple(it.getParent(), it.simpleName, it)
    //        }
    //        .view()

    fastqs=Channel.fromPath("$params.fastqs/**/*.fastq.gz")
            .map {it ->
                tuple(it.getParent().getParent().getParent().getParent().getParent().getName(),it.getParent().getName(), it)
            }
            .groupTuple(by: [0,1])

    fqs=Channel.fromPath("$params.fastqs/**/*.fastq.gz")
            .map {it ->
                tuple(it.getParent().getName(),it.baseName, it)
            }
            .groupTuple(by: [0,1])

    illumina_fastqs=Channel.fromFilePairs("$params.illuminafastqs/*{1,2}.fq.gz")
        .map{it -> tuple(it[0], it[1][0], it[1][1])}
        .view()


    meta_pathogens = Channel.fromPath(params.meta_pathogens)

    // split csv file in tuple of just run,method columns. Remove duplicates
    run_method = meta_pathogens
        .splitCsv(header:true)
        .map{ tuple(it['run'], it['method']) }
        .distinct()
        .view()
    //fqs=JOIN_FQS(illumina_fastqs)

    //fqs=CAT_FASTQS(fastqs)

    fqs.combine(run_method, by:0)
        .branch{ it -> 
                sispa_fqs: it[3] == 'sispa'
                phi_fqs:   it[3] == 'phi' 
                }.set{fq_by_method}

    
    read_fastas = FASTQ_TO_FASTA(fq_by_method.phi_fqs)
    BLASTN(read_fastas, "$params.primer")
    phi_fqs = fq_by_method
                .phi_fqs.map{ tuple(it[0], it[1], it[2]) }
                .combine(BLASTN.out, by:[0,1])
                //.view()

    PHI_SPLIT_READS_BY_PRIMER(phi_fqs)

    // remove method from sispa fqs
    sispa_fqs = fq_by_method.sispa_fqs.map{ tuple(it[0], it[1], it[2]) }

    // combine sispa and SPLIT_READS outputs
    fqs=sispa_fqs.concat(PHI_SPLIT_READS_BY_PRIMER.out.fastq)
    
    KRAKEN2SERVER(fqs)

    //KRAKEN2(fqs)
    SEQKIT_STATS(fqs)

    EXTRACT_RESULTS(KRAKEN2SERVER.out.reports.combine(SEQKIT_STATS.out, by:[0,1]))

    results=EXTRACT_RESULTS.out.tsv // filter our duplicates if any
        .unique()
        .collect()

    COMPILE_RESULTS(results, meta_pathogens)

    // collectS

}