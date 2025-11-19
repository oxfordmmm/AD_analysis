#!/usr/bin/env python3
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from statistics import mean,median
import gzip
import time

def split_by_primer(reads_file, blast_df, min_contig_length):
    reads = SeqIO.parse(gzip.open(reads_file, "rt"), "fastq")
    primer_read_ids = set(blast_df['sseqid'].drop_duplicates().tolist())

    t = time.time()
    blast_df['low'] = blast_df[['sstart', 'send']].min(axis=1)
    blast_df['high'] = blast_df[['sstart', 'send']].max(axis=1)
    # Want to remove random nonomers as well. These appear before the primer
    blast_df['forward'] = blast_df['sstart'] < blast_df['send']
    blast_df['low'] = blast_df['low'] + -9 * blast_df['forward']
    blast_df['high'] = blast_df['high'] + 9 * (1 - blast_df['forward'])

    blast_df = blast_df[['sseqid', 'low', 'high']]
    blast_df = blast_df.sort_values(by=['sseqid', 'low'])
    print(f"blastdf manipulation in {time.time() - t} seconds")


    primers_seqs = []
    non_primer_seqs = []
    split_records = []

    read_stats = []

    for record in reads:
        id = record.id
        has_primer = id in primer_read_ids

        if not has_primer:
            stats = {
                'id' : id,
                'read_length' : len(record.seq),
                'has_primer' : has_primer,
                'primer_count' : 0,
                'pseudoread_count' : 0,
                'pseudoread_mean_length' : 0,
                'pseudoread_median_length' : 0,
            }
            read_stats.append(stats)
            non_primer_seqs.append(record)
            continue


        hits = blast_df.query('sseqid == @id')

        current=0
        i=0
        pseudoread_lengths = []
        for low, high in zip(hits['low'], hits['high']):
            if low - current > min_contig_length:
                new_record = record[current:low]
                i += 1
                pseudoread_lengths.append(low - current + 1)
                new_record.id = f"{id}_split{i}"
                split_records.append(new_record)
            current = high + 1

        # Final split
        if len(record.seq) - current > min_contig_length:
            new_record = record[current:]
            new_record.id = f"{id}_split{i}"
            split_records.append(new_record)
            i += 1
            pseudoread_lengths.append(len(record.seq) - current + 1)
        
        stats = {
            'id' : id,
            'read_length' : len(record.seq),
            'has_primer' : has_primer,
            'primer_count' : hits.shape[0],
            'pseudoread_count' : i,
            'pseudoread_mean_length' : mean(pseudoread_lengths),
            'pseudoread_median_length' : median(pseudoread_lengths),
        }
        read_stats.append(stats)
        primers_seqs.append(record)

    return primers_seqs, non_primer_seqs, split_records, pd.DataFrame(read_stats)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', "--reads", help="fastq file of reads")
    parser.add_argument('-b', "--blast_out", help="tsv file of primer hits")
    parser.add_argument('-o', "--outfile_prefix", help="reads with primer")
    parser.add_argument('-m', "--min_contig_length", help="contigs shorter than this are removed after splitting", default='50')
    args = parser.parse_args()
    min_contig_length = int(args.min_contig_length)

    prefix = args.outfile_prefix
    blast_df = pd.read_csv(args.blast_out, sep='\t')

    t = time.time()
    primers_seqs, non_primer_seqs, split_records, read_stats = split_by_primer(args.reads, blast_df, min_contig_length)
    print(f"Splitting complete in {time.time() - t} seconds")
    SeqIO.write(primers_seqs, f"{prefix}_primer_reads.fastq", "fastq")
    SeqIO.write(non_primer_seqs, f"{prefix}_non_primer_reads.fastq", "fastq")
    SeqIO.write(split_records, f"{prefix}_pseudoreads.fastq", "fastq")
    read_stats.to_csv(f"{prefix}_read_stats.csv", index=False)

    pseudoreads = read_stats.query('has_primer == True')

    with open(f"{prefix}_split_report.txt", "w") as file:
        file.write(f"Non-primer reads: {len(non_primer_seqs)}, Primer reads: {len(primers_seqs)}, Split to {len(split_records)}.\n")
        if pseudoreads.shape[0] > 0:
            file.write(f"Splitting leads to an mean of {mean(pseudoreads['pseudoread_count'])}, and median of {median(pseudoreads['pseudoread_count'])} pseudoreads.\n")
        else:
            file.write(f"No primer reads found\n")