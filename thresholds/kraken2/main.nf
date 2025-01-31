

params.fastqs = ''
params.meta_pathogens="$projectDir/../meta.csv"

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

process EXTRACT_RESULTS {
    tag {run + ' ' + barcode}
    publishDir "results/${run}/${barcode}/", mode: 'copy'

    label 'python'

    input:
    tuple val(run), val(barcode), path("${run}_${barcode}_report.txt")

    output:
    path("${run}_${barcode}_report.tsv"), emit: tsv

    script:
    """
    extract_results.py -i ${run}_${barcode}_report.txt \
        -o ${run}_${barcode}_report.tsv \
        -t 463676 147711 44130 518987 641809 211044 335341 47681 2760809 138948 42789 138949 \
            138951 138950 33757 11224 12730 11216 2560525 162145 11137 277944 186938 2697049 12814 31631 \
            290028 11520 \
            11320 12059 11118 10509 \
            519 520 83558 2104 \
            351073 64320 3052731 12022\
        -r $run -b $barcode
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
    //fastqs=Channel.fromPath("$params.fastqs/**/*.fastq.gz")
    //        .map {it ->
    //            tuple(it.getParent().getName() , it.simpleName, it)
    //        }

    //fastqs=Channel.fromPath("$params.fastqs/**/*.fastq.gz")
    //        .map {it ->
    //            tuple(it.getParent().getParent().getParent().getName(),it.getParent().getName(), it)
    //        }
    //        .groupTuple(by: [0,1])

    illumina_fastqs=Channel.fromFilePairs("$params.illuminafastqs/*{1,2}.fq.gz")
        .map{it -> tuple(it[0], it[1][0], it[1][1])}
        .view()


    //meta_pathogens = Channel.fromPath(params.meta_pathogens)
 
    fqs=JOIN_FQS(illumina_fastqs)

    //fqs=CAT_FASTQS(fastqs)

    KRAKEN2SERVER(fqs)

    //EXTRACT_RESULTS(KRAKEN2SERVER.out.reports)

    //results=EXTRACT_RESULTS.out.tsv
    //    .collect()

    //COMPILE_RESULTS(results, meta_pathogens)

}