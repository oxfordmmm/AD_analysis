// primer
include {FASTQ_TO_FASTA} from './processes/primer.nf'
include {BLASTN; PHI_SPLIT_READS_BY_PRIMER} from './processes/primer.nf'
include {SISPA_TRIM_PRIMER} from './processes/primer.nf'
include {GRAPH_PRIMER_HITS} from './processes/primer.nf'
include {NANOSTAT as NANOSTAT_RAW} from './processes/read_qc.nf'
include {NANOSTAT as NANOSTAT_FILTERED} from './processes/read_qc.nf'
include {NANOPLOT} from './processes/read_qc.nf'
include {PHI_MERGE_SPLITTING_OUTPUTS} from './processes/primer.nf'
include {SISPA_MERGE_TRIMMING_OUTPUTS} from './processes/primer.nf'
include {PHI_SPLITTING_REPORT} from './processes/primer.nf'
include {SISPA_TRIMMING_REPORT} from './processes/primer.nf'  // may need to change this
// include {COMBINE_CSVS} from './processes/general.nf'

// mapping
include {MINIMAP2} from './processes/map.nf'
include {FILTER_BAM} from './processes/map.nf'
include {FILTER_PHI} from './processes/map.nf'
include {FILTER_SISPA} from './processes/map.nf'
include {INDEX_REF} from './processes/map.nf'
include {CLAIR3} from './processes/map.nf'
include {GENOME_DEPTH} from './processes/map.nf' 
include {DEPTH_SUMMARY} from './processes/map.nf'
include {BCFTOOLS_CSQ} from './processes/map.nf' 
include {BCFTOOLS_CONSENSUS} from './processes/map.nf' 
include {GRAPH_COVERAGE} from './processes/map.nf' 
include {GRAPH_COVERAGE_ALL_SAMPLES} from './processes/map.nf' 
include {GRAPH_COVERAGE_SEPARATE} from './processes/map.nf'
include {COUNT_NS} from './processes/map.nf' 
include {GRAPH_PHI_ALIGNMENTS; GRAPH_SISPA_ALIGNMENTS} from './processes/map.nf' 
include {GRAPH_COVERAGE_ALL_SAMPLES_OPTIMISED} from './processes/map.nf'
include {GRAPH_COVERAGE_REPORT} from './processes/map.nf'

// summary
include {SUMMARISE_N_COUNT} from './processes/summaries.nf' 
include {SUMMARISE_ALIGNMENT_INFO} from './processes/summaries.nf' 



/*
For both:
plot raw read statss

For SISPA:
trim primers from ends. If primer hits are found in the middle of the read remove these.
Then run clean reads through kaiju/bracken
map clean reads to refs
make coverage graphs and report errors per 1000 bases

For PHI:
split by primers. Divide into no primers, 1/2 primer, 3-10 primers and many primers
run through kaiju/bracken but report number of reads with/without multiplicity
map split reads to refs
Check for reads which map to same spot multiple times (so probs missed primer)
Check to see if split reads map to the same spot as confirmation.

Could try for a read with 3+ primers. Blast/minimap2 gaps onto themselves to check for similarity
*/

// params.meta_ref=["${params.refs}/spikes.fasta"]


// This workflow simply runs NanoPlot. It requires internet so cannot be run on rescomp slurm cluster (needs -profile local)
workflow nanoplot {
    reads = Channel.fromPath("$params.input/*.fq.gz")
                .map {it ->
                    tuple(it.simpleName, it)
                }
    
    NANOPLOT(reads, "$params.output/nanoplots")
}


// New more generalised entries
workflow phi_workflow {
    references = Channel.fromList(params.meta_ref)

    phi_reads = Channel.fromPath("$params.input/*.fq.gz")
                .map {it ->
                    tuple(it.simpleName, it)
                }
    phi(phi_reads, references)
}

workflow sispa_workflow {
    references = Channel.fromList(params.meta_ref)

    sispa_reads = Channel.fromPath("$params.input/*.fq.gz")
                .map {it ->
                    tuple(it.simpleName, it)
                }
    sispa(sispa_reads, references)
}

workflow RV_workflow {
    references = Channel.fromList(params.RV_ref)

    RV_reads = Channel.fromPath("$params.RV_input/*.fq.gz")
                .map {it ->
                    tuple(it.simpleName, it)
                }

    //RV_reads.view()
    //references.view()

    RV(RV_reads, references)
}
  
workflow RV {
    take:
        reads
        references

    main:

    MINIMAP2( reads.combine(references) )  // map reads to reference
    
    FILTER_SISPA(MINIMAP2.out.bam, 20) //

    // pass to consensus workflow to produce consensus fasta and coverage graphs
    consensus(references, FILTER_SISPA.out.bam, FILTER_SISPA.out.bam_no_ref)
}

workflow corona_workflow {
    references = Channel.fromList(params.corona_ref)

    reads = Channel.fromPath("$params.corona_input/*.fq.gz")
                .map {it ->
                    tuple(it.simpleName, it)
                }

    // Map reads to reference
    MINIMAP2( reads.combine(references) )  
    
    // Filter the bam file
    FILTER_SISPA(MINIMAP2.out.bam, 20)

    // Pass to consensus workflow to produce consensus fasta and coverage graphs
    consensus(references, FILTER_SISPA.out.bam, FILTER_SISPA.out.bam_no_ref)
}

/*
This is the central workflow for dealing with SISPA reads.
1. Will convert reads to fastas in order to use blast to find primer in them.
2. Then it trims the reads to remove the primer parts
3. Then uses minimap2 to align trimmed reads to reference producing a bam file
4. This bam file needs to be filtered to get rid of secondary alignments,
    During the filtering step we also make a csv with alignment info.
5. We pass the filtered bam file to the consensus workflow
*/
workflow sispa {
    take:
        reads
        references

    main:
    read_fastas = FASTQ_TO_FASTA(reads)

    BLASTN(read_fastas, "$params.refs/primers.fasta")

    // Split reads into batches of n=2500 as SISPA_TRIM_PRIMER is slow when n is large
    reads.splitFastq(by: 2500, file: true, elem: [1])
            .map{ fileTuple -> tuple(fileTuple[0], file(fileTuple[1]).baseName, fileTuple[1])}
            .set { batch_fq }

    //SPLIT_READS_BY_PRIMER(batch_fq.combine(BLASTN.out, by: 0))
    //SISPA_TRIM_PRIMER(reads.join(BLASTN.out))

    SISPA_TRIM_PRIMER(batch_fq.combine(BLASTN.out, by: 0))

    //SISPA_TRIM_PRIMER.out.fastq.view()

    // Group the SISPA_TRIM_PRIMER outputs into their original samples
    SISPA_TRIM_PRIMER.out.fastq  // This is non-primer reads and pseudoreads together
        .groupTuple()
        .set { grouped_trimming_fastq }

    SISPA_TRIM_PRIMER.out.read_stats
        .groupTuple()
        .set { grouped_trimming_read_stats }

    SISPA_TRIM_PRIMER.out.failed
        .groupTuple()
        .set { grouped_trimming_failed }

    // Merge SISPA_TRIM_PRIMER outputs back into orgiginal samples
    SISPA_MERGE_TRIMMING_OUTPUTS( grouped_trimming_fastq, grouped_trimming_read_stats, grouped_trimming_failed )

    //     //MERGE_SPLITTING_OUTPUTS.out.fastq.view()
    //     //MERGE_SPLITTING_OUTPUTS.out.report.view()
    //     //MERGE_SPLITTING_OUTPUTS.out.read_stats.view()

    //SPLITTING_REPORT(MERGE_SPLITTING_OUTPUTS.out.read_stats,"$params.output/reports")
    // This needs editing before adding in


    //MINIMAP2(SISPA_TRIM_PRIMER.out.trimmed.combine(references))
    MINIMAP2( SISPA_MERGE_TRIMMING_OUTPUTS.out.fastq.combine(references) )
    
    FILTER_SISPA(MINIMAP2.out.bam, 20)

    sispa_summary_ch = SISPA_MERGE_TRIMMING_OUTPUTS.out.read_stats
            .join(FILTER_SISPA.out.alignments)

    GRAPH_SISPA_ALIGNMENTS(sispa_summary_ch, "$params.output/alignment_summary")

    // pass to consensus workflow to produce consensus fasta and coverage graphs
    consensus(references, FILTER_SISPA.out.bam, FILTER_SISPA.out.bam_no_ref, FILTER_SISPA.out.alignment_counts)
}

/*
This is the central workflow for dealing with PHI reads.
1. Will convert reads to fastas in order to use blast to find primer in them.
2. Then it splits the reads by primer to produce "pseudoreads"
3. Then uses minimap2 to align reads to reference producing a bam file
4. This bam file needs to be filtered. With Phi we expect to sometimes get
    multiple alignments for a single read, we just need to make sure these are 
    not overlapping. 
5. We pass the filtered bam file to the consensus workflow
*/

workflow phi {
    take:
        reads
        references
    
    main:
        read_fastas = FASTQ_TO_FASTA(reads)
        BLASTN(read_fastas, "$params.refs/primers.fasta")

        //reads.view(
        reads.splitFastq(by: 2500, file: true, elem: [1])
            .map{ fileTuple -> tuple(fileTuple[0], file(fileTuple[1]).baseName, fileTuple[1])}
            .set { batch_fq }
            
        //batch_fq.view()
        // [4, /gpfs3/well/bag/users/vbr851/projects/agnostic_diagnostic_v2/work/db/986cca81e162a07ae2854cda543fd6/4.1.fq]

        PHI_SPLIT_READS_BY_PRIMER(batch_fq.combine(BLASTN.out, by: 0))

        PHI_SPLIT_READS_BY_PRIMER.out.fastq  // This is non-primer reads and pseudoreads together
            .groupTuple()
            .set { grouped_splitting_fastq }

        PHI_SPLIT_READS_BY_PRIMER.out.read_stats
            .groupTuple()
            .set { grouped_splitting_read_stats }

        PHI_MERGE_SPLITTING_OUTPUTS( grouped_splitting_fastq, grouped_splitting_read_stats )

        //MERGE_SPLITTING_OUTPUTS.out.fastq.view()
        //MERGE_SPLITTING_OUTPUTS.out.read_stats.view()

        PHI_SPLITTING_REPORT(PHI_MERGE_SPLITTING_OUTPUTS.out.read_stats,"$params.output/reports")

        MINIMAP2( PHI_MERGE_SPLITTING_OUTPUTS.out.fastq.combine(references) )
        
        FILTER_PHI(MINIMAP2.out.bam, 20, 30)

        phi_summary_ch = PHI_MERGE_SPLITTING_OUTPUTS.out.read_stats
                .join(FILTER_PHI.out.non_overlapping)

        GRAPH_PHI_ALIGNMENTS(phi_summary_ch, "$params.output/alignment_summary")

        // // pass to consensus workflow to produce consensus fasta and coverage graphs
        consensus(references, FILTER_PHI.out.bam, FILTER_PHI.out.bam_no_ref)
}



/*
This section is modified from Nick's flu pipeline.
Uses the bam alignment file to do some basecalling with clair3, 
then goes on to produce a consensus fasta and depth coverage graph.
*/
workflow consensus {
    take:
        references
        bam
        bam_no_ref 
        alignment_counts

    main:

    // sam tools index
    INDEX_REF(references)

    // produces vcf file giving points of variance against reference. 
    // clair3 needed because ont output has more complex error profile
    CLAIR3(bam_no_ref.combine(INDEX_REF.out.index))

    // get depth of reads along ref
    GENOME_DEPTH(bam.combine(alignment_counts, by:0))

    // sumarise the depth for all samples
    depths=GENOME_DEPTH.out.csv
        .collect()

    DEPTH_SUMMARY(depths)

    // create consensus file by applying vcf to ref. Depth is used to mask low depth regions
    BCFTOOLS_CONSENSUS(CLAIR3.out.vcf.combine(INDEX_REF.out.index).combine(GENOME_DEPTH.out.tsv, by:0))
    
    // Calculate some basic stats about the resulting consensus fasta (count snps and N's etc)
    to_count_ns = BCFTOOLS_CONSENSUS.out.fasta.join(CLAIR3.out.vcf)
    COUNT_NS(to_count_ns)

    // Produce coverage graphs for the spiked references
    to_graph = GENOME_DEPTH.out.tsv.map {it -> 
        tuple(it[0], it[1], "${it[0]}_depth.png")
    }
    //GRAPH_COVERAGE(to_graph, 50, 5, 100)

    GRAPH_COVERAGE_SEPARATE(to_graph, 50, 5, 100)

    // all_depth_files = GENOME_DEPTH.out.tsv
    //                     .map{it -> it[1]}
    //                     .collect()

    //GRAPH_COVERAGE_ALL_SAMPLES(all_depth_files, 50, 5, 100)

    Channel.fromPath("${params.test_rmd}")
        .set{ ch_rmd }

    GRAPH_COVERAGE_REPORT(ch_rmd)

    //GRAPH_COVERAGE_ALL_SAMPLES_OPTIMISED(all_depth_files, 50, 5, 100)
}

/*
This workflow merges some outputs so that all the barcodes can be compared in single files.
Outputs from here will go into output/$batch/summary
*/
workflow summaries {
    SUMMARISE_N_COUNT("$params.output/n_counts", 'error_rates.csv', "$params.output/summary")
    SUMMARISE_ALIGNMENT_INFO("$params.output/alignment_summary", "$params.output/summary")
}