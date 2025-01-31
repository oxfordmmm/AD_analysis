#!/usr/bin/env python3
import argparse
import pandas as pd
from Bio import SeqIO

def getBlast(blast,df):
    df2=pd.read_csv(blast,sep='\t',names=['qseqid','sseqid','pident','length','qstart','qend','sstart','send'])

    # remove hits that are not on the same contig
    df2=df2[df2['qseqid']==df2['sseqid']]

    # remove hits are are the entire length of the contig
    df2['slength']=df2['sstart']-df2['send']
    df2['slength']=df2['slength'].abs()
    df2=df2[df2['length']!=df2['slength']]

    # remove hits that are less than 50 bp
    df2=df2[df2['length']>=50]

    # remove hits that are less than 90% identical
    df2=df2[df2['pident']>=95]

    print(df2)

    # reduce to only the columns we need for masking
    df2=df2[['qseqid','qstart','qend']]
    df2.columns=['chrom','start','end']


    # merge the two dataframes
    df=pd.concat([df,df2])
    return df



def run(ref, mask, blast, output):
    df = pd.read_csv(mask)
    maskDict = {}
    
    seqs=[]
    if blast:
        df=getBlast(blast,df)
    for i, row in df.iterrows():
        maskDict.setdefault(row['chrom'], []).append({'start': row['start'], 'end': row['end']})
    print(maskDict)
    
    for seq in SeqIO.parse(ref, 'fasta'):
        if seq.id in maskDict:
            for mask in maskDict[seq.id]:
                print(f'Masking {seq.id} from {mask["start"]} to {mask["end"]}')
                start=mask['start']
                end=mask['end']
                seq.seq = seq.seq[:start] + 'N'*(end-start) + seq.seq[end:]

                
        seqs.append(seq)

    SeqIO.write(seqs, output, 'fasta')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Mask reference sequence')
    parser.add_argument('-r', '--reference', required=True, 
                        help='Reference sequence fasta file')
    parser.add_argument('-m', '--mask', required=True,
                        help='Masked regions csv file')
    parser.add_argument('-b', '--blast', required=False,   
                        help='blast output file')
    parser.add_argument('-o', '--output', required=True,
                        help='Output masked reference fasta file')
    args = parser.parse_args()

    run(args.reference, args.mask, args.blast, args.output)

