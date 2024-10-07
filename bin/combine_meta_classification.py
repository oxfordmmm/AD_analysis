#!/usr/bin/env python3
import argparse
import pandas as pd

# Script combines the clasifications for each read from kraken2 and kaiju
# Will take the value with the more precise classification, selecting k2 in draws

# These are the taxon levels (unclassified, root, domain, kingdom etc), higher is more precise
taxon_mapping = {
    'U': 0,
    'R': 1,
    'D': 2,
    'K': 3,
    'P': 4,
    'C': 5,
    'O': 6,
    'F': 7,
    'G': 8,
    'S': 9
}

def combine_classification(kaiju_file, kraken2_file, k2_report, output_file):
    # Read Kaiju and Kraken2 classification files into pandas DataFrames
    kaiju_df = pd.read_csv(kaiju_file, sep='\t', header=None,
                           names=['kaiju_classified', 'read_id', 'kaiju_taxon', 'match_length', 'all_taxon', 'all_accessions', 'seq'])
    kaiju_df = kaiju_df[['kaiju_classified', 'read_id', 'kaiju_taxon', 'match_length']]
    kraken2_df = pd.read_csv(kraken2_file, sep='\t', header=None,
                             names=['k2_classified', 'read_id', 'k2_taxon', 'read_length', 'kmers'])
    
    # The kraken2 report is used to get the mappings from taxon id to taxon level and taxon name
    taxon_labels = pd.read_csv(k2_report, sep='\t', names=['pc', 'count', 'count_at_level', 'taxonomy_lvl', 'taxonomy_id', 'name'])[['name', 'taxonomy_id', 'taxonomy_lvl']]
    taxon_labels['name'] = taxon_labels['name'].str.replace(' ', '')
    # Sometimes the taxon level is G2 which means between G and S. For our purposes that is too precise so we remove the number with some regex
    taxon_labels['taxon_score'] = taxon_labels['taxonomy_lvl'].str.replace(r'\d+', '', regex=True).map(taxon_mapping)
    print(taxon_labels)
    
    # Need to separately add the taxon level/name to the k2 results and the kaiju results.
    k2_taxon_labels = taxon_labels.rename(columns = {'taxonomy_id' : 'k2_taxon', 'name' : 'k2_name', 'taxonomy_lvl' : 'k2_taxon_lvl', 'taxon_score' : 'k2_score'})
    kaiju_taxon_labels = taxon_labels.rename(columns = {'taxonomy_id' : 'kaiju_taxon', 'name' : 'kaiju_name', 'taxonomy_lvl' : 'kaiju_taxon_lvl', 'taxon_score' : 'kaiju_score'})

    df = pd.merge(kaiju_df, kraken2_df, on='read_id', how='outer') \
            .merge(k2_taxon_labels, on='k2_taxon', how='left') \
            .merge(kaiju_taxon_labels, on='kaiju_taxon', how='left')


    df['use_kaiju'] = df['kaiju_score'] > df['k2_score']
    df['classified'] = df['kaiju_classified'].where(df['use_kaiju'], df['k2_classified'])
    df['taxon'] = df['kaiju_taxon'].where(df['use_kaiju'], df['k2_taxon'])
    df['taxon_lvl'] = df['kaiju_taxon_lvl'].where(df['use_kaiju'], df['k2_taxon_lvl'])
    df['name'] = df['kaiju_name'].where(df['use_kaiju'], df['k2_name'])

    df = df[['read_id', 'read_length', 'classified', 'taxon', 'taxon_lvl', 'name', 'use_kaiju']]
    df.to_csv(output_file, sep='\t', index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Combine Kaiju and Kraken2 classification files.')
    parser.add_argument('kaiju_file', help='Path to the Kaiju classification file.')
    parser.add_argument('kraken2_file', help='Path to the Kraken2 classification file.')
    parser.add_argument('k2_report_file', help='Path to the braken file, used of taxon labels.')
    parser.add_argument('output_file', help='Path to the output combined file.')

    args = parser.parse_args()

    combine_classification(args.kaiju_file, args.kraken2_file, args.k2_report_file, args.output_file)
