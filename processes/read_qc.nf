
process NANOSTAT {
    conda "$params.envs/qc"

    publishDir "$outdir/", mode: 'copy'

    input:
    tuple val(barcode), path(nano_reads)
    val(outdir)

    output:
    tuple val(barcode), path("${barcode}_nanostat.tsv")

    script:
    """
    NanoStat --tsv --fastq $nano_reads > ${barcode}_nanostat.tsv
    """
}

/*
Process needs internet!

Creates some plots with general read qc stats
*/
process NANOPLOT {
    conda "$params.envs/qc"

    publishDir "$outdir/", mode: 'copy'
    publishDir "$params.output/summary/nanoplot", mode: 'copy', pattern: "*.png"

    input:
    tuple val(barcode), path(nano_reads)
    val(outdir)

    output:
    tuple val(barcode), path("${barcode}_nanostats.tsv"), emit: 'stats'
    tuple val(barcode), path("${barcode}_graphs"), emit: 'graphs'
    tuple val(barcode), path("*.png"), emit: 'main_graphs'

    script:
    """
    NanoPlot --tsv_stats --fastq $nano_reads --plots kde --N50 -p ${barcode}_
    mv ${barcode}_NanoStats.txt ${barcode}_nanostats.tsv
    mkdir ${barcode}_graphs
    mv ${barcode}_*.html ${barcode}_graphs
    mv ${barcode}_*.png ${barcode}_graphs
    cp ${barcode}_graphs/*Non_weightedLog*.png ${barcode}_read_lengths.png 
    """
}

