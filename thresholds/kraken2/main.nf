

params.fastqs = ''


process KRAKEN2SERVER {
    tag {run + ' ' + barcode}
    publishDir "classifications/${run}/${barcode}/", mode: 'copy'
    maxForks 15

    label 'kraken2server'

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


workflow {
    fastqs=Channel.fromPath("$params.fastqs/**/*.fastq.gz")
            .map {it ->
                tuple(it.getParent().getName() , it.simpleName, it)
            }


    KRAKEN2SERVER(fastqs)

}