/*
Process maps reads to reference to produce a sam file of alignments,
Samtools is used to sort and filter,
Tries to keep as many alignments as possible
*/
process MASK_REF {
    conda "$params.envs/blast"
    cpus 3

    input:
    tuple path('ref.fasta'), path('mask.csv')

    output:
    path('masked.fasta') , emit: fasta

    script:
    """
    mask_ref.py -r ref.fasta -m mask.csv -o masked.fasta
    """
}

process MINIMAP2 {
    tag {barcode}
    conda "$params.envs/map"
    cpus 3
    //publishDir "bams/minimap2"
    publishDir "bams/minimap2_unfiltered"

    input:
    tuple val(barcode), path('sequences.fastq'), path('ref.fasta')

    output:
    tuple val(barcode), path('ref.fasta'), path("${barcode}_minimap.bam"), path("${barcode}_minimap.bam.bai"), emit: bam

    script:
    """
    minimap2 -t $task.cpus -ax map-ont -N 1000 ref.fasta sequences.fastq | samtools view -bS - | \
        samtools sort -o ${barcode}_minimap.bam
    samtools index ${barcode}_minimap.bam
    """
}

process FILTER_BAM {
    conda "$params.envs/map"

    publishDir "bams/filtered"

    input:
    tuple val(barcode), path('ref.fasta'), path("${barcode}.bam"), path("${barcode}.bam.bai")
    val(min_pc_of_ref)
    val(min_read_length)
    val(outdir)

    output:
    tuple val(barcode), path('ref.fasta'), path("${barcode}.filt.bam"), path("${barcode}.bam.bai"), emit: bam
    tuple val(barcode), path("${barcode}.filt.bam"), path("${barcode}.filt.bam.bai"), emit: bam_no_ref

    script:
    """
    python3 $params.bin/filterBAM.py -i ${barcode}.bam -o ${barcode}.filt.bam -p $min_pc_of_ref -m $min_read_length
    samtools index ${barcode}.filt.bam
    """
}

/*
Takes all the alignments from minimap2 and filters them.
For phi a single read can have multiple alignments, but can't overlap,
    as we want them to be distinct alignments
It produces a filtered bam file, as well as some csvs with alignment info.
*/
process FILTER_PHI {
    conda "$params.envs/map"
    errorStrategy 'ignore'

    publishDir "$params.output/bams", pattern: "*.bam"
    publishDir "$params.output/alignment_info", pattern: "*.csv", mode: 'copy'
    publishDir "$params.output/reports", pattern: "*.txt", mode: 'copy'

    input:
    tuple val(barcode), path('ref.fasta'), path("${barcode}.bam"), path("${barcode}.bam.bai")
    val(min_align_length)
    val(max_read_overlap)

    output:
    tuple val(barcode), path('ref.fasta'), path("${barcode}.filt.bam"), path("${barcode}.filt.bam.bai"), emit: bam
    tuple val(barcode), path("${barcode}.filt.bam"), path("${barcode}.filt.bam.bai"), emit: bam_no_ref
    tuple val(barcode), path("${barcode}_all_alignments.csv"), emit: all_alignments
    tuple val(barcode), path("${barcode}_non_overlapping.csv"), emit: non_overlapping
    tuple val(barcode), path("*_report.txt"), emit: report

    script:
    """
    echo $barcode
    python3 ${params.bin}filter_phi_alignments.py -i ${barcode}.bam -o ${barcode}.filt.bam \
        -p ${barcode} -m $min_align_length -x $max_read_overlap
    samtools index ${barcode}.filt.bam
    """
}


/*
Process filters bam file from minimap2.
For SISPA we are just removing secondary alignments, are ones which are too short.
*/
process FILTER_SISPA {
    conda "$params.envs/map"
    errorStrategy 'ignore'
    
    //publishDir "$params.output/bams", pattern: "*.bam"
    publishDir "$params.output/alignment_info", pattern: "*_alignments.csv", mode: 'copy'
    publishDir "$params.output/reports", pattern: "*.txt", mode: 'copy'

    input:
    tuple val(barcode), path('ref.fasta'), path("${barcode}.bam"), path("${barcode}.bam.bai")
    val(min_align_length)

    output:
    tuple val(barcode), path('ref.fasta'), path("${barcode}.filt.bam"), path("${barcode}.filt.bam.bai"), emit: bam, optional: true
    tuple val(barcode), path("${barcode}.filt.bam"), path("${barcode}.filt.bam.bai"), emit: bam_no_ref, optional: true
    tuple val(barcode), path("${barcode}*.csv"), emit: stats, optional: true
    tuple val(barcode), path("${barcode}_alignments.csv"), emit: alignments, optional: true
    tuple val(barcode), path("*_report.txt"), emit: report, optional: true
    tuple val(barcode), path("${barcode}_alignment_counts.csv"), emit: alignment_counts

    script:
    """
    echo barcode $barcode
    filter_sispa_alignments.py -i ${barcode}.bam -o ${barcode}.filt.bam \
        -p ${barcode} -m $min_align_length -s $barcode
    samtools index ${barcode}.filt.bam
    """
}

process CALC_OVERLAPS {
    tag {barcode}
    conda "$params.envs/bedtools"
    errorStrategy 'ignore'

    //publishDir "$params.output/alignment_info/overlaps", pattern: "*.csv", mode: 'copy'

    input:
    tuple val(barcode), path('ref.fasta'), path("${barcode}.bam"), path("${barcode}.bam.bai")

    //output:
    //tuple val(barcode), path("${barcode}_overlaps.csv"), emit: overlaps

    script:
    """
    bedtools bamtobed -i ${barcode}.bam > reads.bed
    bedtools merge -i reads.bed > merged_reads.bed
    bedtools subtract -a reads.bed -b merged_reads.bed > non_overlapping_reads.bed
    """
}

process BCFTOOLS_CSQ {
    conda "$params.envs/map"
    publishDir "$params.output/vcfs", mode: 'copy'

    input:
    tuple val(sampleName), path('vcf'), path('ref.fasta'),path('ref.fasta.fai'), path('GFF3')

    output:
    tuple val(sampleName), path("${sampleName}.vcf"), emit: vcf

    script:
    """
    python3 $params.bin/replaceIUPAC.py -i ref.fasta -o ref_N.fasta
    bcftools csq \
	    -f ref_N.fasta \
	    -g GFF3 \
	    vcf \
	    -Ot -o ${sampleName}.vcf \
        -p a \
        -s - \
	    --force
    """
    stub:
    """
    touch ${sampleName}.vcf
    """
}

process GRAPH_PHI_ALIGNMENTS {
    conda "$params.envs/r"
    publishDir "$outdir", mode: 'copy'

    input:
    tuple val(barcode), path("read_stats.csv"), path("alignment_info.csv")
    val(outdir)

    output:
    tuple val(barcode), path("${barcode}_read_align_stats.csv"), emit: 'read_stats'
    tuple val(barcode), path("${barcode}_summary_stats.csv"), emit: 'summary_stats'
    tuple val(barcode), path("${barcode}_phi_plots"), emit: 'plots'

    script:
    """
    Rscript ${params.bin}graph_phi_alignments.R read_stats.csv alignment_info.csv
    mv read_align_stats.csv ${barcode}_read_align_stats.csv
    mv summary_stats.csv ${barcode}_summary_stats.csv

    mkdir ${barcode}_phi_plots
    mv *_plot.png ${barcode}_phi_plots
    """
}

process GRAPH_SISPA_ALIGNMENTS {
    conda "$params.envs/r"
    publishDir "$outdir", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(barcode), path("read_stats.csv"), path("alignment_info.csv")
    val(outdir)

    output:
    tuple val(barcode), path("${barcode}_read_align_stats.csv"), emit: 'read_stats'
    tuple val(barcode), path("${barcode}_summary_stats.csv"), emit: 'summary_stats'
    tuple val(barcode), path("${barcode}_sispa_plots"), emit: 'plots'

    script:
    """
    Rscript $params.bin/graph_sispa_alignments.R read_stats.csv alignment_info.csv $barcode
    mv read_align_stats.csv ${barcode}_read_align_stats.csv
    mv summary_stats.csv ${barcode}_summary_stats.csv

    mkdir ${barcode}_sispa_plots
    mv *_plot.png ${barcode}_sispa_plots
    """
}

process INDEX_REF {
    conda "$params.envs/map"

    input:
    path('ref.fa.gz')

    output:
    tuple path('ref_N.fa'), path('ref_N.fa.fai'), emit: index

    script:
    """
    zless ref.fa.gz > ref.fa
    python3 $params.bin/replaceIUPAC.py -i ref.fa -o ref_N.fa
    samtools faidx ref_N.fa
    """
}

/*
Clair3 is a tool for variant calling with ont data.
Produces a vcf with all the differences from reference
*/
process CLAIR3 {
    conda "$params.envs/clair3"
    cpus=6

    publishDir "$params.output/vcfs", mode: 'copy', saveAs: { filename -> "${sample}.vcf.gz"} 

    input:
    tuple val(sample), path('mapped.bam'), path('mapped.bam.bai'), path('ref.fa'), path('ref.fa.fai')

    output:
    tuple val(sample), path("clair_out/merge_output.vcf.gz"), emit: vcf

    script:
    MODEL_NAME='r941_prom_sup_g5014'
    """
    clair_path=\$(which run_clair3.sh)
    PREFIX=\$(dirname "\$clair_path")
    echo "dirname" \$PREFIX

    run_clair3.sh \
        --bam_fn=mapped.bam \
        --ref_fn=ref.fa \
        --threads=${task.cpus} \
        --platform="ont" \
        --model_path=\${PREFIX}/models/${MODEL_NAME} \
        --include_all_ctgs \
        --no_phasing_for_fa \
        --output=clair_out

    """
}

/*
Finds depth of read coverage across reference
*/
process GENOME_DEPTH {     
    tag {barcode}                                                      
    conda "$params.envs/map"
    publishDir "$params.output/depths", overwrite: true, mode: 'copy'  
                                                                                
    input:                                                                      
    tuple val(barcode), path('ref.fasta'), path("${barcode}.sorted.bam"), path("${barcode}.sorted.bam.bai"), path("alignment_counts.csv"), path("alignment_info.csv")
    each path("multi_segment_pathogens.csv")

    output:                                                                     
    tuple val(barcode), path("${barcode}_depth.csv"), emit: genomeDepths
    tuple val(barcode), path("${barcode}_depth.tsv"), emit: tsv 
    path("${barcode}_depth.csv"), emit: csv        
                                                                                
    script:                                                                     
    """                                                                         
    samtools depth -aa ${barcode}.sorted.bam > ${barcode}_depth.tsv                
    coverageStats.py ${barcode}_depth.tsv  ${barcode} alignment_counts.csv  multi_segment_pathogens.csv alignment_info.csv                    
                                                                                
    mv coverage_stats.csv ${barcode}_depth.csv                                 
    """                                                  
}

process DEPTH_SUMMARY {
    conda "$params.envs/map"
    publishDir "$params.output/depth_summary", mode: 'copy'

    input:
    path("depths???.csv")
    path("meta_pathogens.csv")
    path("pathogens_reduced.csv")

    output:
    path("*depth_summary.csv"), emit: summary
    //path("*depth_summary_pass.csv"), emit: pass_summary

    script:
    batch = params.batch
    """
    depth_summary.py -i *.csv -o ${params.batch}_depth_summary.csv -b $batch -mp meta_pathogens.csv -pr pathogens_reduced.csv
    """
}

process BCFTOOLS_CONSENSUS {
    conda "$params.envs/map"
    publishDir "$params.output/consensus", mode: 'copy'

    input:
    tuple val(sampleName), path('vcf'), path('ref.fasta'),path('ref.fasta.fai'),path('depth.txt')

    output:
    tuple val(sampleName), path("${sampleName}.fasta"), emit: fasta

    script:
    """
    python3 $params.bin/replaceIUPAC.py -i ref.fasta -o ref_N.fasta
    awk '{if (\$3<$params.min_depth) print \$1"\t"\$2}' depth.txt > mask.txt
    tabix -p vcf vcf
    bcftools consensus --mask mask.txt -f ref_N.fasta vcf > ${sampleName}.fasta
    """
    stub:
    """
    touch ${sampleName}.fasta
    """
}
process COUNT_NS {
    conda "$params.envs/biopython"
    publishDir "$params.output/n_counts", mode: 'copy'

    input:
    tuple val(sampleName), path(consensus_file), path('vcf.gz')

    output:
    tuple val(sampleName), path("${sampleName}_N_count.csv")

    script:
    """
    echo $sampleName
    gzip -dc vcf.gz > vcf 
    python3 $params.bin/count_errors.py $consensus_file vcf ${sampleName}_N_count.csv
    """
}

/*
Produces graph with coverage across the references
*/
process GRAPH_COVERAGE {
    conda "$params.envs/r"

    publishDir "$params.output/coverage_graphs", mode: 'copy'

    input:
    tuple val(sampleName), path(depth_file), val(outfile)
    val(group_size)
    val(low_depth_cutoff)
    val(high_cutoff)

    output:
    tuple val(sampleName), path("*.png")

    script:
    """
    Rscript $params.bin/graph_coverage.R $depth_file $group_size $low_depth_cutoff $high_cutoff $outfile $sampleName
    """
}

process PLOT_BACTERIAL {
    conda "$params.envs/pyplots"

    publishDir "$params.output/coverage_graphs/bacteria", mode: 'copy'

    input:
    tuple val(sampleName), path(depth_file) 
    
    output:
    tuple val(sampleName), path("${sampleName}.pdf"), emit: 'pdf'
    tuple val(sampleName), path("${sampleName}.csv"), emit: 'csv'

    script:
    """
    plot_bacteria.py $depth_file ${sampleName}
    """
}

process GRAPH_COVERAGE_SEPARATE {
    tag {sampleName}
    conda "$params.envs/r"
    errorStrategy 'ignore'

    publishDir "$params.output/coverage_graphs", mode: 'copy'

    input:
    tuple val(sampleName), path(depth_file)
    val(group_size)
    val(low_depth_cutoff)
    val(high_cutoff)

    output:
    tuple val(sampleName), path("*.png"), emit: 'png'
    tuple val(sampleName), path("*.pdf"), emit: 'pdf'


    script:
    """
    Rscript $params.bin/graph_coverage_separate.R $depth_file $group_size $low_depth_cutoff $high_cutoff $sampleName
    """
}

process GRAPH_COVERAGE_ALL_SAMPLES {
    conda "$params.envs/r"
    publishDir "$params.output/coverage_graphs_by_segment", mode: 'copy'
    cpus=3

    input:
    path(depth_files)
    val(group_size)
    val(low_depth_cutoff)
    val(high_cutoff)

    output:
    path("*")

    script:
    """
    graph_coverage_all_samples.R $depth_files $group_size $low_depth_cutoff $high_cutoff
    """
}

process GRAPH_COVERAGE_REPORT {
    conda "$params.envs/r2"
    //publishDir "$params.output/", mode: 'copy'
    cpus=3

    input:
    path(rmd)
    
    output:
    path("notebook.html")
    //ßpath("out_test.txt")

    script:
    """
    cp -L ${rmd} notebook.Rmd
    #head ${rmd} > out_test.txt
    Rscript -e "library(rmarkdown); \
                xfun::session_info(); \
                rmarkdown::render('notebook.Rmd', 'html_document')"
    
     #Rscript -e "library(rmarkdown); \
     #            library(rmarkdown); \
       #          library(tinytex); \
      #           rmarkdown::render('test_env.R','html_document')" 
    # Rscript -e 'rmarkdown::render("${rmd}", output_file="script.pdf", output_dir = getwd())'
    """
}


process GRAPH_COVERAGE_ALL_SAMPLES_OPTIMISED {
    conda "$params.envs/r"
    publishDir "$params.output/coverage_graphs_by_segment2", mode: 'copy'
    cpus=5

    input:
    path(depth_files)
    val(group_size)
    val(low_depth_cutoff)
    val(high_cutoff)

    output:
    path("*")

    script:
    """
    graph_coverage_all_samples_optimised.R $depth_files $group_size $low_depth_cutoff $high_cutoff
    """
}


// process TEST {
//     conda "$params.envs/r"
//     publishDir "$params.output/coverage_graphs_by_segment", mode: 'copy'

//     input:
//     path(depth_files)
//     val(group_size)
//     val(low_depth_cutoff)
//     val(high_cutoff)

//     output:
//     tuple path("*")

//     script:
//     """
//     test2.R $depth_files $group_size $low_depth_cutoff $high_cutoff
//     """
// }