include {KRAKEN2} from './processes/classify.nf'
//include {KRAKEN2} from './processes/classification.nf'
//include {KAIJU} from './processes/classification.nf'
//include {MERGE_META} from './processes/classification.nf'

// if you dont want to do any trimming or splitting - hopefully for sispa trimming wont be significant
workflow rawreads {
    reads = Channel.fromPath("$params.input/*.fq.gz")
                .map {it ->
                    tuple(it.simpleName, it)
                }
    read_outdir = Channel.fromPath("$params.output/meta_classifications/rawreads").first()

    meta_classification(reads, read_outdir)
}


// dependent on output of main workflow
workflow pseudoreads {
    split_reads = Channel.fromPath("$params.output/splitting/*split_reads.fastq.gz")
                .map {it ->
                    tuple(it.simpleName.split('_')[0], it)
                }
    trimmed_reads = Channel.fromPath("$params.output/trimming/*trimmed.fastq")
                .map {it ->
                    tuple(it.simpleName.split('_')[0], it)
                }

    read_outdir = Channel.fromPath("$params.output/meta_classifications/pseudoreads").first()

    reads = split_reads.concat(trimmed_reads)
    
    meta_classification(reads, read_outdir)
}

workflow meta_classification {
    take:
        reads
        read_outdir

    main:
    //KRAKEN2(reads, outpath, 0)
    //KAIJU(reads, outpath)

    // combined = KAIJU.out.output
    //             .join(KRAKEN2.out.classification)
    //             .join(KRAKEN2.out.report)
    // MERGE_META(combined, outpath)

    KRAKEN2(reads, read_outdir)
}