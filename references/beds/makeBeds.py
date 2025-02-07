#!/usr/bin/env python
from Bio import SeqIO

# Read in the fasta file
fasta = SeqIO.parse("../meta_ref.fasta", "fasta")

bacteria_chroms = ['CP085971.1', 'NZ_CP025371.1', 'NC_005043.1', 'NZ_LR214945.1']

# make two bed files, one for bacteria and one for viruses using fasta to get the lengths
with open("bacteria.bed", "w") as bacteria, open("virus.bed", "w") as virus:
    for record in fasta:
        if record.id in bacteria_chroms:
            bacteria.write(f"{record.id}\t0\t{len(record.seq)}\n")
        else:
            virus.write(f"{record.id}\t0\t{len(record.seq)}\n")
