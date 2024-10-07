
process COMBINE_CSVS {
    conda "$params.envs/pandas"

    publishDir "$outdir/", mode: 'copy'

    input:
    val(barcode_tsv_list)
    val(outfile)
    val(outdir)
    val(format_in)
    val(format_out)
    val(skip_rows)

    output:
    path("$outfile")

    script:
    """
    python3 $params.bin/combine_csvs.py '$barcode_tsv_list' $outfile \
        -i $format_in -o $format_out -s $skip_rows
    """
}
/* example use with barcode, "path"
    blasts = BLASTN.out
        .map({it -> 
            tuple(it[0], "\"${it[1]}\"")
        })
        .toSortedList()
    COMBINE_CSVS_BLAST(blasts, "promoters.csv", "$params.amr_genes_dir", \
        "tsv", "csv", 0)
*/


/*
Process takes fastq files (which can be gzipped), and converts to fastas
*/

process SUMMARISE_READ_STATS {
    conda "$params.envs/biopython"

    publishDir "$outdir", mode: 'copy'

    input:
    path(read_stats)
    val(outfile)
    val(outdir)

    output:
    path("$outfile")

    script:
    """
    python3 $params.bin/count_stats.py $read_stats $outfile
    """
}

process SUMMARISE_N_COUNT {
    conda "$params.envs/biopython"

    publishDir "$outdir", mode: 'copy'

    input:
    val(dir)
    val(outfile)
    val(outdir)

    output:
    path("$outfile")

    script:
    """
    python3 $params.bin/combine_error_counts.py $dir $outfile
    """
}

process SUMMARISE_ALIGNMENT_INFO {
    conda "$params.envs/biopython"

    publishDir "$outdir", mode: 'copy'

    input:
    val(dir)
    val(outdir)

    output:
    path("overall_summary.csv")

    script:
    """
    python3 $params.bin/combine_alignment_summaries.py $dir
    """
}