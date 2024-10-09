

params.virus = "none"
params.genome = "$projectDir/references/${params.virus}.fasta"
params.output = "$projectDir/sim_out"
params.num_reads = 1000

workflow {
    NANOSIM(params.genome, params.num_reads, "sim_${params.virus}", params.output)

    CIRCULARISE(NANOSIM.out.reads, "${params.virus}_circularised.fasta", params.output)
}

process NANOSIM {
    conda "$params.envs/nanosim"
    publishDir "$outdir", mode: 'copy'

    input:
    path(ref)
    val(num_reads)
    val(out_prefix)
    val(outdir)

    output:
    path("${out_prefix}.fa"), emit: reads
    path("${out_prefix}*"), emit: all

    script:
    """
    fasta_len=\$(awk '!/^>/' $ref | wc -m )
    len=\$((\$fasta_len - \$fasta_len/10))
    nanosim-h -o $out_prefix -n $num_reads --unalign-rate 0 --max-len \$len $ref
    """
}

process CIRCULARISE {
    conda "$params.envs/biopython"
    publishDir "$outdir", mode: 'copy'

    input:
    path(reads)
    val(outfile)
    val(outdir)

    output:
    path("$outfile")

    script:
    """
    python3 $projectDir/bin/loop_reads.py $reads $outfile
    """
}