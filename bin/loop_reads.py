#!/usr/bin/env python3
import argparse
from Bio import SeqIO
import random

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta", help="fasta file of reads")
    parser.add_argument("outfile", help="summary file")
    args = parser.parse_args()

    new_records = []
    for record in SeqIO.parse(args.fasta, "fasta"):
        print(record.id)
        record.seq = "AAAAAAAAAAAA" + record.seq
        record_len = len(record.seq)

        start = random.randint(0, record_len)
        end = start + random.randint(100, record_len * 5)

        record.seq = (record.seq*6)[start:end]
        new_records.append(record)
    
    SeqIO.write(new_records, args.outfile, "fasta")