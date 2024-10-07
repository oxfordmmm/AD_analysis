#!/usr/bin/env pypy3
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def trim(reads, blast_df, margin_pc, min_read_length):
    sequences = SeqIO.parse(reads, "fastq")
    primer_read_ids = blast_df['sseqid'].drop_duplicates().tolist()
    trimmed_seqs = []
    removed_seqs = []
    read_stats = []

    all_hits = blast_df.copy()
    all_hits['low'] = all_hits[['sstart', 'send']].min(axis=1)
    all_hits['high'] = all_hits[['zsstart', 'send']].max(axis=1)
    # Want to remove random nonomers as well
    all_hits['forward'] = all_hits['sstart'] < all_hits['send']
    all_hits['low'] = all_hits['low'] + -9 * all_hits['forward']
    all_hits['high'] = all_hits['high'] + 9 * (1 - all_hits['forward'])
    print(all_hits)

    for seq_record in sequences:
        if seq_record.id not in primer_read_ids:
            trimmed_seqs.append(seq_record)
            stats = {
                'id' : seq_record.id,
                'read_length' : len(seq_record.seq),
                'has_primer' : False,
                'primer_count' : 0,
                'failed_trimming' : False,
                'fail_reason' : '',
            }
            read_stats.append(stats)
            continue

        seq_length = len(seq_record.seq)
        hits = all_hits.query('sseqid == @seq_record.id').copy()
        hits = hits.sort_values('low')

        if hits.shape[0] > 10:
            # Should have at most 2 primers, but can get strings of primers stuck together. However remove if too many
            removed_seqs.append(seq_record)
            stats = {
                'id' : seq_record.id,
                'read_length' : seq_length,
                'has_primer' : True,
                'primer_count' : hits.shape[0],
                'failed_trimming' : True,
                'fail_reason' : 'too_many_primers',
            }
            read_stats.append(stats)
            continue


        new_start = 0
        new_end = seq_length
        trim_fail = False
        for _, row in hits.iterrows():
            # Check if primer is at the front (up to margin_pc)
            if (100 * row['high'] / seq_length) < margin_pc:
                new_start = max(new_start, row['high'])
            # Check if primer is at the front (up to margin_pc)
            elif (100 * row['low'] / seq_length) > 100 - margin_pc:
                new_end = min(new_end, row['low'])
            else:
                # Primer must have been in the middle so remove
                removed_seqs.append(seq_record)
                trim_fail = True
                stats = {
                    'id' : seq_record.id,
                    'read_length' : seq_length,
                    'has_primer' : True,
                    'primer_count' : hits.shape[0],
                    'failed_trimming' : True,
                    'fail_reason' : 'primer_in_middle',
                }
                read_stats.append(stats)
                break
        
        if new_end - new_start < min_read_length:
            removed_seqs.append(seq_record)
            trim_fail = True
            stats = {
                'id' : seq_record.id,
                'read_length' : seq_length,
                'has_primer' : True,
                'primer_count' : hits.shape[0],
                'failed_trimming' : True,
                'fail_reason' : 'resulting_read_too_short',
            }
            read_stats.append(stats)

        if not trim_fail:
            new_record = seq_record[new_start:new_end]
            trimmed_seqs.append(new_record)
            stats = {
                'id' : seq_record.id,
                'read_length' : seq_length,
                'has_primer' : True,
                'primer_count' : hits.shape[0],
                'failed_trimming' : False,
                'fail_reason' : '',
            }
            read_stats.append(stats)

    return trimmed_seqs, removed_seqs, primer_read_ids, pd.DataFrame(read_stats)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', "--reads", help="fastq file of reads")
    parser.add_argument('-b', "--blast_out", help="tsv file of primer hits")
    parser.add_argument('-o', "--outfile_prefix", help="reads with primer trimmed")
    parser.add_argument('-m', "--margin_pc", type=int, help="primer must be at the start or end with this margin", default=20)
    parser.add_argument('-x', "--min_read_length", type=int, help="resulting read must be longer than this", default=50)
    args = parser.parse_args()


    prefix = args.outfile_prefix
    blast_df = pd.read_csv(args.blast_out, sep='\t')
    trimmed, failed, primer_read_ids, read_stats = trim(args.reads, blast_df, margin, min_read_length)

    SeqIO.write(trimmed, f"{prefix}_trimmed.fastq", "fastq")
    SeqIO.write(failed, f"{prefix}_removed.fastq", "fastq")
    read_stats.to_csv(f"{prefix}_read_stats.csv", index=False)

    read_stats.groupby(['has_primer', 'fail_reason']).size().reset_index(name='count').to_csv(f"{prefix}_trim_report.txt", index=False)
    # with open(f"{prefix}_report.txt", "w") as file:
    #     file.write(f"Failed: {len(failed)}, Success: {len(trimmed)}, Primer found in {len(primer_read_ids)} reads")