#!/usr/bin/env python3
import sys
import pysam
from argparse import ArgumentParser


def run(opts):
    insamfile = pysam.AlignmentFile(opts.inBam, "rb")
    outsamfile = pysam.AlignmentFile(opts.outBam, "wb", template=insamfile)

    lens=insamfile.lengths
    names=insamfile.references
    ref_lenDict={}
    for n,l in zip(names, lens):
        ref_lenDict[n]=l

    for read in insamfile.fetch():
        read_len=read.infer_read_length()
        ref_len=ref_lenDict[read.reference_name]
        read_pc = (read_len/ref_len)*100 
        if read_pc > int(opts.perLength) and read_pc < 115 and read_len >= int(opts.minLength):
            #print(read_len, ref_len)
            outsamfile.write(read)

    insamfile.close()
    outsamfile.close()

if __name__ == '__main__':
    parser = ArgumentParser(description='Filter BAM file by read length proportion of reference')
    parser.add_argument('-i', '--inBam', required=True,
                        help='input bam file')
    parser.add_argument('-o', '--outBam', required=True,
                        help='output bam file')
    parser.add_argument('-p', '--perLength', required=True,
                        help='percentant of reference read length required to pass filtering')    
    parser.add_argument('-m', '--minLength', required=True,
                        help='min read length required to pass filtering')    
    opts, unknown_args = parser.parse_known_args()
   
    run(opts)