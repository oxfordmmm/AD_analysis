#!/usr/bin/env python3
import argparse
import os
import pandas as pd

def merge_csv_files(input_dir):
    # Get all CSV files in the input directory
    csv_files = [file for file in os.listdir(input_dir) if file.endswith('summary_stats.csv')]

    # Initialize an empty DataFrame to store the merged results
    dfs = []

    # Iterate over each CSV file
    for file in csv_files:
        id = int(file.split('_')[0])
        file_path = os.path.join(input_dir, file)

        # Read the CSV file and extract the 'diff_rate' column
        df = pd.read_csv(file_path)
        df['barcode'] = id
        print(df)
        dfs.append(df)

    # Write the merged data to the output file
    df = pd.concat(dfs)
    barcodes = df.pop('barcode')
    df.insert(0, 'barcode', barcodes)
    df = df.sort_values('barcode')

    df = df.round({
        'read_length_Q1' : 0,
        'read_length_Q2' : 0,
        'read_length_Q3' : 0,
        'read_length_mean' : 0,
        'mapped_read_length_Q1' : 0,
        'mapped_read_length_Q2' : 0,
        'mapped_read_length_Q3' : 0,
        'mapped_read_length_mean' : 0,
        'pseudo_circle_size_Q1' : 0,
        'pseudo_circle_size_Q2' : 0,
        'pseudo_circle_size_Q3' : 0,
        'pseudo_circle_size_Q95' : 0,
        'pseudo_circle_size_mean' : 0,
        'mapped_pseudo_circle_size_Q1' : 0,
        'mapped_pseudo_circle_size_Q2' : 0,
        'mapped_pseudo_circle_size_Q3' : 0,
        'mapped_pseudo_circle_size_Q95' : 0,
        'mapped_pseudo_circle_size_mean' : 0,
        'align_circle_size_Q1' : 0,
        'align_circle_size_Q2' : 0,
        'align_circle_size_Q3' : 0,
        'align_circle_size_Q95' : 0,
        'align_circle_size_mean' : 0,
    })
    df.to_csv(f"overall_summary.csv", index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Merge CSV alignment summaries.')
    parser.add_argument('dir', help='Input directory containing phi CSV files')

    # Parse the arguments
    args = parser.parse_args()

    # check if dir exists
    if not os.path.isdir(args.dir):
        print("directory does not exist. Are you running a sispa batch??")
    else:
        # Call the merge_csv_files function with the provided arguments
        merge_csv_files(args.dir)
