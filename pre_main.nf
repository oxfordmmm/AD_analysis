// Most of this is taken from crumpit2

// include modules           
//include {DOWNLOAD_TAXONOMY} from '../modules/classify.nf'
//include {MAKE_KTAXONOMY} from '../modules/classify.nf'
//include {KRAKEN2} from '../modules/classify.nf'
//include {KREPORT_BASES} from '../modules/classify.nf'
//include {COMBINE_MULTI_KREPORTS} from '../modules/classify.nf'
//include {COMBINE_MULTI_KREPORTS as COMBINE_MULTI_KREPORTS_BASES} from '../modules/classify.nf'

//include {BATCHREADS} from './processes/refs.nf'

include {CREATE_FTP_TABLE} from './processes/refs.nf'
include {GET_ASSEMBLY_META} from './processes/refs.nf'
include {GET_ASSEMBLY_META_GENBANK} from './processes/refs.nf'
include {DOWNLOAD_REFS} from './processes/refs.nf'
include {EDIT_FASTA_HEADERS} from './processes/refs.nf'
include {MERGE_REFS} from './processes/refs.nf'

workflow get_refs {

    GET_ASSEMBLY_META()

    GET_ASSEMBLY_META_GENBANK()

    // Make a channel containing all the chosen taxids for references we want to download
    params.taxid_list = "$params.refs/taxid_list.txt"
    taxids = Channel.fromPath(params.taxid_list) \
        | splitCsv(header:false)
        //| map { row -> tuple(row, file("viral_assembly_summary.txt"), file("bacteria_assembly_summary.txt")) }

    //GETREFS(taxids.combine(GET_ASSEMBLY_META.out.viral)
    //                    .combine(GET_ASSEMBLY_META.out.bacteria))
    
    // get_refs_tuple = taxids.combine(GET_ASSEMBLY_META.out.viral)
    //                     .combine(GET_ASSEMBLY_META.out.bacteria)
    //                     .combine(GET_ASSEMBLY_META_GENBANK.out.bacteria)
    //                     .combine(GET_ASSEMBLY_META_GENBANK.out.viral)

    

    // CREATE_FTP_TABLE(get_refs_tuple)

    get_refs_tuple = taxids.combine(GET_ASSEMBLY_META.out.viral)
                        .combine(GET_ASSEMBLY_META.out.bacteria)

    

    CREATE_FTP_TABLE(get_refs_tuple)

    //CREATE_FTP_TABLE.out[1].view()

    // ftp_rows = CREATE_FTP_TABLE.out       
    //             .map{ it -> it[1]}

    // ftp_rows = CREATE_FTP_TABLE.out.tsv        
    //             .map{ it -> it[1]}

    ftp_rows = CREATE_FTP_TABLE.out.csv       
                .map{ it -> it[1]}
    
    // Create a table for reference
    ftp_rows.collectFile(name: 'ftp_table.csv', storeDir: "$params.refs")
    // It would be handy having a csv version too - change this


    DOWNLOAD_REFS(CREATE_FTP_TABLE.out.tsv)

    //EDIT_FASTA_HEADERS(DOWNLOAD_REFS.out.fasta, CREATE_FTP_TABLE.out)

    EDIT_FASTA_HEADERS( DOWNLOAD_REFS.out.fasta.combine(CREATE_FTP_TABLE.out.tsv, by: 0))

    all_fastas = EDIT_FASTA_HEADERS.out.fasta
                    .map{ it -> it[1]}
                    .collect()

    MERGE_REFS(all_fastas)
    
}

/*
// run nextflow on raw reads (untrimmed), so use this for now just on sispa
workflow raw_reads {
    reads = Channel.fromPath("$params.input/*.fq.gz")
                .map {it ->
                    tuple(it.simpleName, it)
                }
    read_outdir = Channel.fromPath("$params.output/raw_read_meta_classification").first()
    
    meta_classification(reads, read_outdir)
}

// ignore this for now, but we may tweak and use it later
// this is dependent on the output of main and uses trimmed (sispa) OR split (phi) reads
workflow pseudoreads {
    split_reads = Channel.fromPath("$params.output/splitting/*split_reads.fastq")
                .map {it ->
                    tuple(it.simpleName.split('_')[0], it)
                }
    trimmed_reads = Channel.fromPath("$params.output/trimming/*trimmed.fastq")
                .map {it ->
                    tuple(it.simpleName.split('_')[0], it)
                }
    read_outdir = Channel.fromPath("$params.output/pseudo_read_meta_classification").first()
    
    reads = split_reads.concat(trimmed_reads)
    get_ref(reads, read_outdir)
}

workflow get_ref {
    take:
        reads
        outpath

    main:

    DOWNLOAD_TAXONOMY()

    MAKE_KTAXONOMY(classifier_db.combine(DOWNLOAD_TAXONOMY.out))

    //KRAKEN2(reads)

    classifier=KRAKEN2(reads)
    cr=classifier

    crb=KREPORT_BASES(classifier.txt.combine(MAKE_KTAXONOMY.out))

    kreports=cr.reports
        .groupTuple(by:[0])

    kreports_bases=crb.reports
        .groupTuple(by:[0])
    
    COMBINE_MULTI_KREPORTS(kreports)

    COMBINE_MULTI_KREPORTS_BASES(kreports_bases)

    emit:
    txt = classifier.read_tax
    read_reports = COMBINE_MULTI_KREPORTS.out
    base_reports = COMBINE_MULTI_KREPORTS_BASES.out


    // This is from the mapping subworkflow
    BATCHREADS(kraken2.output.read_tax)
}

// now run classification
*/
