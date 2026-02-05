#!/usr/bin/env python3
from Bio import SeqIO
import pandas as pd
import argparse
import os

def ref_data(ref_fasta):
    """Extract reference sequence information from a FASTA file and save to CSV."""
    ref_sequences = []
    for record in SeqIO.parse(ref_fasta, "fasta"):
        ref_sequences.append({
            'Sequence_ID': record.id,
            'Length': len(record.seq)
        })

    df_ref = pd.DataFrame(ref_sequences)
    return df_ref

def format_table_S5(df):
    """Format to table S5 structure"""
    cols=['subject_id','Sequence_ID', 'pathogen_reduced', 'Length']
    df=df[cols]
    df=df.rename(columns={
        'subject_id': 'Accession',
        'pathogen_reduced': 'Routine PCR target',
        'Sequence_ID': 'Reference name'
    })
    # remove controls 
    controls=['MS2', 'orthoreovirus', 'murine_respirovirus', 'zika']
    df = df[~df['Routine PCR target'].isin(controls)]

    # combine multi segment references into the same row
    multi_segment_refs = ['Influenza A/H1-2009', 'Influenza A/H3', 'Influenza B', 'Influenza A/H5', 'Influenza A/H1']
    df_combined = pd.DataFrame()
    for ref in multi_segment_refs:
        df_ref = df[df['Routine PCR target'] == ref]
        if not df_ref.empty:
            combined_accession = ';'.join(df_ref['Accession'].tolist())
            combined_length = df_ref['Length'].sum()
            num_segments = df_ref.shape[0]
            df_combined = pd.concat([df_combined, pd.DataFrame({
                'Accession': [combined_accession],
                'Routine PCR target': [ref],
                'Reference name': [ref],
                'Length': [combined_length],
                'Number of segments': [num_segments]
            })], ignore_index=True)
            df = df[df['Routine PCR target'] != ref]
    df['Number of segments'] = 1
    df = pd.concat([df, df_combined], ignore_index=True)
    df = df.sort_values('Routine PCR target')
    df.to_csv('table_S5.csv', index=False)
    return df



def main(ref_fasta, output_file, meta_file, biofire_name_file, blastn_file):
    df_ref = ref_data(ref_fasta)
    
    # Load meta file to map sequence IDs to reference names
    df_meta = pd.read_csv(meta_file)
    
    # Merge reference data with meta information
    df_merged = pd.merge(df_ref, df_meta, left_on='Sequence_ID', right_on='reference', how='left')

    # Load Biofire reference name mapping
    df_biofire = pd.read_csv(biofire_name_file)

    # Merge to get Biofire reference names
    df_merged = pd.merge(df_merged, df_biofire, on='pathogen', how='left')
    
    # drop empty rows
    df_merged = df_merged.dropna(subset=['pathogen_reduced'])

    # drop duplicates
    df_merged = df_merged.drop_duplicates(subset=['Sequence_ID'], keep='first')

    # Load FTP table
    names=['query_id', 'subject_id', 'percent_identity', 'alignment_length', 'mismatches', 'gap_opens',
           'q_start', 'q_end', 's_start', 's_end', 'evalue', 'bit_score']
    dfb = pd.read_csv(blastn_file, sep='\t', names=names)
    
    # filter out hits with <100% identity
    dfb = dfb[dfb['percent_identity'] == 100.0]

    # take the longest alignment for each query_id
    dfb = dfb.sort_values('alignment_length', ascending=False).drop_duplicates(subset=['query_id'], keep='first')

    # Merge with main dataframe to get lengths
    df_merged = pd.merge(df_merged, dfb[['query_id', 'subject_id', 'alignment_length']], left_on='Sequence_ID', right_on='query_id', how='left')
    
    # Save to CSV
    df_merged.to_csv(output_file, index=False)

    # format table S5
    df_table_S5 = format_table_S5(df_merged)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get reference sequence information.")
    parser.add_argument('-r', '--ref_fasta', help='Reference FASTA file', required=True)
    parser.add_argument('-o', '--output_file', help='Output CSV file for reference info', required=True)
    parser.add_argument('-m', '--ref_meta', help='meta file look up table for seqid to reference name', required=True)
    parser.add_argument('-b', '--biofire_name', help='Biofire reference name mapping file', required=True)
    parser.add_argument('-n', '--blastn', help='TSV of blastn results in outfmt 6 format', required=True)
    args = parser.parse_args()

    main(args.ref_fasta, args.output_file, args.ref_meta, args.biofire_name, args.blastn)

    