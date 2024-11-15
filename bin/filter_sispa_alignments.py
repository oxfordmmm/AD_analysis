#!/usr/bin/env python3
import sys
import pysam
import pandas as pd
from argparse import ArgumentParser


def get_good_alignments(in_file, out_file, min_length):
    input_bam = pysam.AlignmentFile(in_file, "rb")
    out_bam = pysam.AlignmentFile(out_file, "wb", template=input_bam)

    lens=input_bam.lengths
    names=input_bam.references
    ref_lenDict={}
    for n,l in zip(names, lens):
        ref_lenDict[n]=l

    alignments = []
    unmapped_count = 0
    secondary_count = 0
    supplementary_count = 0

    for read in input_bam.fetch(until_eof=True):  
        if read.is_unmapped:
            unmapped_count += 1
            continue
        if read.is_secondary:
            secondary_count += 1
        if read.is_supplementary:
            supplementary_count += 1

        read_len=read.infer_read_length()
        ref_len=ref_lenDict[read.reference_name]
        if read.query_alignment_length >= min_length:
            out_bam.write(read)
            read_stats = {
                'query' : read.query_name,
                'ref' : read.reference_name,
                'read_len' : read_len,
                'align_len' : read.query_alignment_length,
                'query_start' : read.query_alignment_start,
                'query_end' : read.query_alignment_end,
                'ref_start' : read.reference_start,
                'ref_end' : read.reference_end,
                'ref_len' : ref_len,
                'supp' : read.is_supplementary,
                'secondary' : read.is_secondary,
                'mapping_quality' : read.mapping_quality,
                'mismatches_gaps' : read.get_tag('NM'),
                'dp_alignment_score' : read.get_tag('AS'),
            }
            alignments.append(read_stats)

    input_bam.close()
    out_bam.close()
    return pd.DataFrame(alignments), unmapped_count, secondary_count, supplementary_count


if __name__ == '__main__':
    parser = ArgumentParser(description='Filter BAM file by read length proportion of reference')
    parser.add_argument('-i', '--inBam', required=True,
                        help='input bam file')
    parser.add_argument('-o', '--outBam', required=True,
                        help='output bam file')
    parser.add_argument('-p', '--out_prefix', required=True,
                        help='output prefix for other files')
    parser.add_argument('-m', '--min_length', required=True,
                        help='min alignment length required to pass filtering')    
    parser.add_argument('-s', '--sample_name', required=True,
                        help='sample name')
    args = parser.parse_args()
    input_bam = args.inBam
    output = args.outBam
    prefix = args.out_prefix
    min_length = int(args.min_length)

    # Simply filters reads with too short an alignment

    with open(f"{prefix}_filter_report.txt", "w") as file:
        alignments, unmapped_count, secondary_count, supplementary_count = get_good_alignments(input_bam, output, min_length)
        if len(alignments) > 0:
            alignments = alignments.sort_values(by=['query', 'query_start'])
            alignments.to_csv(f"{prefix}_all_alignments.csv", index=False)
            alignments_all  = alignments.copy()

            alignments = alignments.query('supp == False and secondary == False')
    
        # Get good alignments
        alignments.to_csv(f"{prefix}_alignments.csv", index=False)
        file.write(f"{unmapped_count} reads were not aligned. (May have been removed earlier)\n")
        file.write(f"{secondary_count} reads were secondary. \n")
        file.write(f"{supplementary_count} reads were supplementary. \n")
        file.write(f"{alignments.shape[0]} good alignments.\n")

    # count alignments per chrom
    #alignment_counts = alignments.groupby(['ref','query']).size().reset_index(name='counts')
    if len(alignments) > 0:
        alignment_counts = alignments_all.groupby(['ref'])[['query']].count()#.reset_index(name='counts')
        alignment_counts.rename(columns={'query':'non-unique counts'}, inplace=True)
        alignment_counts_unique = alignments_all.groupby(['ref'])[['query']].nunique()#.reset_index(name='unique counts')
        alignment_counts_unique.rename(columns={'query':'counts'}, inplace=True)
        alignment_counts = alignment_counts.merge(alignment_counts_unique, on='ref', how='left')
    else:
        # create empty dataframe with ref as index
        alignment_counts = pd.DataFrame(columns=['ref','non-unique counts','counts'])
        alignment_counts.set_index('ref', inplace=True)
    # add unmapped counts
    alignment_counts.loc['unmapped'] = [unmapped_count, unmapped_count]
    alignment_counts['Sample name'] = args.sample_name
    alignment_counts.reset_index(inplace=True)
    alignment_counts.rename(columns={'ref':'chrom'}, inplace=True)
    alignment_counts.to_csv(f"{prefix}_alignment_counts.csv", index=False)