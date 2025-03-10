#!/usr/bin/env python3
import sys
import gzip

def join_reads(fastq1, fastq2, output_fastq):
    with gzip.open(fastq1, 'rt') as f1, gzip.open(fastq2, 'rt') as f2, gzip.open(output_fastq, 'wt') as out:
        while True:
            header1 = f1.readline().strip()
            seq1 = f1.readline().strip()
            plus1 = f1.readline().strip()
            qual1 = f1.readline().strip()

            header2 = f2.readline().strip()
            seq2 = f2.readline().strip()
            plus2 = f2.readline().strip()
            qual2 = f2.readline().strip()

            if not header1 or not header2:
                break

            joined_seq = seq1 + seq2
            joined_qual = qual1 + qual2

            out.write(f"{header1}\n{joined_seq}\n{plus1}\n{joined_qual}\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: join_reads.py <fastq1> <fastq2> <output_fastq>")
        sys.exit(1)

    fastq1 = sys.argv[1]
    fastq2 = sys.argv[2]
    output_fastq = sys.argv[3]

    join_reads(fastq1, fastq2, output_fastq)