#!/usr/bin/env python3
import sys
import pysam
import pandas as pd
from argparse import ArgumentParser


def get_good_alignments(in_file, min_length):
    input_bam = pysam.AlignmentFile(in_file, "rb")

    ref_lenDict={}
    for n,l in zip(input_bam.references, input_bam.lengths):
        ref_lenDict[n]=l

    alignments = []
    unmapped_count = 0

    for read in input_bam.fetch():        
        if read.is_unmapped:
            unmapped_count += 1
            continue

        read_len=read.infer_read_length()
        ref_len=ref_lenDict[read.reference_name]
        if read.query_alignment_length >= min_length:
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
    return pd.DataFrame(alignments), unmapped_count

def filter_dups(alignments, max_overlap):
    alignments['first_aln'] = alignments.groupby('query')['ref_start'].transform('min')
    alignments['num_alns'] = alignments.groupby('query').transform('size')
    alignments = alignments.sort_values(by=['first_aln', 'query', 'supp', 'dp_alignment_score'],
                                        ascending=[True, True, True, False])
    # df is ordered by query and then by score

    single_alns = alignments.query("num_alns == 1")
    multiple_alns = alignments.query("num_alns > 1")

    keep_rows = []
    current_range = []
    prev_query = ''
    for row in multiple_alns.itertuples(index=False):
        if prev_query != row.query:
            prev_query = row.query
            keep_rows.append(row)
            current_range = [(row.query_start, row.query_end)]
        else:
            s = row.query_start
            e = row.query_end
            overlap = False
            for x, y in current_range:
                # Check for overlap
                overlap_start = max(s, x)
                overlap_end = min(e, y)
                if overlap_end > overlap_start + max_overlap:
                    overlap = True
                    break
            
            if not overlap:
                keep_rows.append(row)
                current_range.append((s, e)) 

    non_overlap = pd.concat([single_alns, pd.DataFrame(keep_rows)])
    # Recalculate summary now overlaps removed
    non_overlap['first_aln'] = non_overlap.groupby('query')['ref_start'].transform('min')
    non_overlap['num_alns'] = non_overlap.groupby('query').transform('size')
    non_overlap = non_overlap.sort_values(by=['first_aln', 'query', 'supp', 'dp_alignment_score'],
                                        ascending=[True, True, True, False])
    return non_overlap

def output_clean_bam(in_file, out_file, nonoverlap):
    input_bam = pysam.AlignmentFile(in_file, "rb")
    output_bam = pysam.AlignmentFile(out_file, "wb", template=input_bam)

    non_overlap_ids = list(zip(
        nonoverlap['query'],
        nonoverlap['ref'],
        nonoverlap['query_start'],
        nonoverlap['query_end'],
        nonoverlap['ref_start']
    ))

    for read in input_bam.fetch():
        id = (read.query_name, read.reference_name,
              read.query_alignment_start, read.query_alignment_end,
              read.reference_start)
    
        if id in non_overlap_ids:
            output_bam.write(read)


    input_bam.close()
    output_bam.close()


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
    parser.add_argument('-x', '--max_overlap', required=True,
                        help='max_overlap allowed for secondary alignment')    
    args = parser.parse_args()
    input_bam = args.inBam
    output = args.outBam
    prefix = args.out_prefix
    min_length = int(args.min_length)
    max_overlap = int(args.max_overlap)

    with open(f"{prefix}_filter_report.txt", "w") as file:
    
        # Get alignments for good alignments
        alignments, unmapped_count = get_good_alignments(input_bam, min_length)
        if alignments.empty:
            file.write("alignments is empty, no reads map to reference!!")
        else:
            alignments.sort_values(by=['query', 'query_start']).to_csv(f"{prefix}_all_alignments.csv", index=False)
            file.write(f"{unmapped_count} reads were not aligned. (May have been removed earlier)\n")
            file.write(f"{alignments.shape[0]} alignments.\n")

            # Now need to prune the duplicate/overlapping alignments
            non_overlap = filter_dups(alignments, max_overlap)
            non_overlap.to_csv(f"{prefix}_non_overlapping.csv", index=False)
            file.write(f"{non_overlap.shape[0]} non overlapping alignments.\n")

            # Now go back to the bam and make a filtered copy bam
            output_clean_bam(input_bam, output, non_overlap)