process KRAKEN2 {
    conda "$params.envs/bracken"
    cpus=1
    publishDir "$outdir", mode: 'copy'

    input:
    tuple val(barcode), path(reads)
    val(outdir)
    val(confidence_theshold)

    output:
    tuple val(barcode), path("${barcode}_classification.tsv"), emit: classification
    tuple val(barcode), path("${barcode}_kraken2_report.tsv"), emit: report
    tuple val(barcode), path("${barcode}_bracken_report.tsv"), emit: bracken_report

    script:
    """
    kraken2 --threads $task.cpus --db $params.kraken_db $reads --confidence $confidence_theshold \
        --output ${barcode}_classification.tsv --report ${barcode}_kraken2_report.tsv
    bracken -d $params.kraken_db -i ${barcode}_kraken2_report.tsv -o ${barcode}_bracken_report.tsv
    """
}

process KAIJU {
    conda "$params.envs/kaiju"
    cpus=1
    publishDir "$outdir", mode: 'copy'

    input:
    tuple val(barcode), path(reads)
    val(outdir)

    output:
    tuple val(barcode), path("${barcode}_kaiju_summary.tsv"), emit: 'summary'
    tuple val(barcode), path("${barcode}_kaiju_out.tsv"), emit: 'output'

    script:
    """
    kaiju -v -z $task.cpus -t $params.kaiju_db/nodes.dmp -f $params.kaiju_db/viruses/kaiju_db_viruses.fmi -i $reads > ${barcode}_kaiju_out.tsv
    kaiju2table -v -e -l 'family,genus,species' -r species -t $params.kaiju_db/nodes.dmp -n $params.kaiju_db/names.dmp -o ${barcode}_kaiju_summary.tsv ${barcode}_kaiju_out.tsv
    """
}

process MERGE_META {
    conda "$params.envs/biopython"
    publishDir "$outdir", mode: 'copy'

    input:
    tuple val(barcode), path('kaiju_out.tsv'), path('k2_classification.tsv'), path('k2_report.tsv')
    val(outdir)

    output:
    tuple val(barcode), path("${barcode}_merged_meta.tsv"), emit: 'meta'

    script:
    """
    python3 $params.bin/combine_meta_classification.py kaiju_out.tsv k2_classification.tsv k2_report.tsv ${barcode}_merged_meta.tsv
    """
}