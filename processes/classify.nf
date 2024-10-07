// processes taken from Nick's crumpit2

process DOWNLOAD_TAXONOMY {
    output:
    tuple path('names.dmp'), path('nodes.dmp')

    script:
    """
    centrifuge-download -o tax_dump taxonomy
    mv tax_dump/* .
    """
    stub:
    """
    touch names.dmp nodes.dmp
    """
}

process MAKE_KTAXONOMY {

    input:
    tuple val(centname), file('cent_db'), path('names.dmp'), path('nodes.dmp')

    output:
    tuple val(centname), path('ktaxonomy.txt')

    script:
    if (params.classifier=='centrifuge') {
    """
    centrifuge-inspect --conversion-table cent_db/$centname > seqid2taxid.txt
    make_ktaxonomy.py --nodes nodes.dmp --names names.dmp --seqid2taxid seqid2taxid.txt -o ktaxonomy.txt
    """
    }
    else if (params.classifier=='kraken2') {
    """
    make_ktaxonomy.py --nodes nodes.dmp --names names.dmp --seqid2taxid cent_db/seqid2taxid.map -o ktaxonomy.txt
    """

    }
    stub:
    """
    touch ktaxonomy.txt
    """
}

process KRAKEN2 {
    conda "$params.envs/bracken"
    tag {barcode}

    //publishDir "$params.output/meta_classification/${outdir}/kraken2/${barcode}", mode: 'copy'
    publishDir "$outdir/kraken2", mode: 'copy'

    cpus 4
    //memory 9

    input:
    tuple val(barcode), path(fastq)
    val(outdir)

    output:
    tuple val(barcode), path('read_tax.txt'), emit: read_tax
    tuple val(barcode), path("${barcode}_kraken2_report.txt"), emit: reports
    tuple val(barcode), path("${barcode}_classification.txt"), emit: txt

    script:
    """
    kraken2 --gzip-compressed $fastq \
        --db ${params.kraken2db} --threads ${task.cpus} \
        --report ${barcode}_kraken2_report.txt \
        --output ${barcode}_classification.txt \
        # --report ${barcode}_kraken2_report.tsv
        

    echo "readID taxID" > read_tax.txt
    awk '{print \$2,\$3}' ${barcode}_kraken2_report.txt >> read_tax.txt
    """
    stub:
    """
    mkdir kraken2
    """
}

process KREPORT_BASES {
    tag {barcode}

    publishDir "$params.output/classifications/${barcode}/", mode: 'copy'

    input:
    tuple val(barcode), val(batchID), file('cent.txt'), val(centname), path('ktaxonomy.txt')

    output:
    tuple val(barcode), val(batchID), file("${barcode}_${batchID}_bases.txt"), emit: reports

    script:
    """
    create_kreport.py -i cent.txt -t ktaxonomy.txt \
        -o ${barcode}_bases.txt -c ${params.classifier}
    """

    stub:
    """
    touch ${barcode}_bases.txt
    """
}

process COMBINE_MULTI_KREPORTS {
    tag {barcode}
    errorStrategy 'ignore'

    publishDir "$params.output/${params.outdir}/kreports/${task.process.replaceAll(":","_").toLowerCase()}"

    input:
    tuple val(barcode), val(batchID), file("${barcode}_?????.txt")

    output:
    tuple val(barcode), file("${barcode}.txt")

    script:
    """
    combine_kreports.py -r  ${barcode}_*.txt \
        -o ${barcode}.txt --only-combined
    """
    stub:
    """
    touch ${barcode}.txt
    """
}

// this will mix barcodes in current formation
process PROPOGATE_CLASSIFICATIONS {
    tag {barcode}

    publishDir "$params.output/classifications/${barcode}/", mode: 'copy'

    input: 
    val data

    output:
    path "*_read_tax.txt"
    
    script:
    new_input = data.getAt(-1)
    if (data.size() > 3)
        { old_input = data.getAt(-2) }
            else { old_input='empty.txt' }
    barcode = data.getAt(0)
    batchID = data.getAt(1)
    """
    touch empty.txt
    cat $old_input $new_input > ${barcode}_${task.index}_read_tax.txt
    """
    stub:
    """
    echo "Task $task.index inputs: $data" > out_${task.index}.txt
    """
}