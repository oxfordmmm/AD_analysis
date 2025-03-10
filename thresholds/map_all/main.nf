

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

process MINIMAP2 {
    tag {run + ' ' + barcode + ' ' + ref}
    publishDir "bams/${run}/"

    //label 'cat'
    cpus 4

    input:
    tuple val(run), val(barcode), path("${run}_${barcode}.fastq.gz"),  val(ref), path("ref.fa")

    output:
    tuple val(run), val(barcode), path("${run}_${barcode}.bam")

    script:
    """
    minimap2 -t $task.cpus -a -x sr ref.fa ${run}_${barcode}.fastq.gz | samtools view -bS -q 40 -F 4 - | samtools sort -o ${run}_${barcode}.bam
    """
}

process SAMTOOLS_DEPTH {
    tag {run + ' ' + barcode }
    errorStrategy 'ignore'
    //publishDir "fastqs/${run}/", mode: 'copy'

    //label 'cat'

    input:
    tuple val(run), val(barcode), path("${run}_${barcode}.bam")

    output:
    tuple val(run), val(barcode), path("${run}_${barcode}.depth")

    script:
    """
    samtools depth -aa ${run}_${barcode}.bam > ${run}_${barcode}.depth
    """
}

process SAMTOOLS_STATS {
    tag {run + ' ' + barcode }
    publishDir "stats/${run}/", mode: 'copy'
    errorStrategy 'ignore'


    //label 'cat'

    input:
    tuple val(run), val(barcode), path("${run}_${barcode}.bam")

    output:
    tuple val(run), val(barcode), path("${run}_${barcode}.stats")

    script:
    """
    samtools stats  ${run}_${barcode}.bam > ${run}_${barcode}.stats
    """
}

process PLOT_BACTERIAL {
    tag {sampleName + ' ' + barcode}
    errorStrategy 'ignore'


    publishDir "bacteria_pdfs", mode: 'copy', pattern: '*.pdf'

    input:
    tuple val(run), val(barcode), path("${run}_${barcode}.depth")
    
    output:
    tuple val(run), val(barcode), path("${run}_${bacteria}.pdf"), emit: 'pdf', optional: true
    tuple val(run), val(barcode), path("${run}_${bacteria}.csv"), emit: 'csv', optional: true

    script:
    """
    plot_bacteria.py ${run}_${barcode}.depth ${run}_${barcode}
    """
}

workflow {
    //fastqs=Channel.fromPath("$params.fastqs/**/*.fastq.gz")
    //        .map {it ->
    //            tuple(it.getParent().getName() , it.simpleName, it)
    //        }

    fastqs=Channel.fromPath("$params.fastqs/**/*.fastq.gz")
            .map {it ->
                tuple(it.getParent().getName(), it.simpleName, it)
            }

    refs = Channel.fromPath("$params.refs")
            .map {it -> tuple(it.simpleName, it)
            }
 

    //CAT_FASTQS(fastqs)

    MINIMAP2(fastqs.combine(refs))

    SAMTOOLS_DEPTH(MINIMAP2.out)

    SAMTOOLS_STATS(MINIMAP2.out)

    PLOT_BACTERIAL(SAMTOOLS_DEPTH.out)

    //EXTRACT_RESULTS(KRAKEN2SERVER.out.reports)

    //results=EXTRACT_RESULTS.out.tsv
    //    .collect()

    //COMPILE_RESULTS(results, meta_pathogens)

}