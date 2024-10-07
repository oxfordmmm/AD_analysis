#!/usr/bin/env python3
import argparse
import os
import pandas as pd

def merge_csv_files(input_dir, output_file):
    # Get all CSV files in the input directory
    csv_files = [file for file in os.listdir(input_dir) if file.endswith('count.csv')]

    # Initialize an empty DataFrame to store the merged results
    dfs = []

    # Iterate over each CSV file
    for file in sorted(csv_files):
        id = int(file.split('_')[0])
        file_path = os.path.join(input_dir, file)

        # Read the CSV file and extract the 'diff_rate' column
        data = pd.read_csv(file_path)[['header', 'n_percent', 'diff_percent']]
        # data['n/diff'] = data.agg('{0[n_rate]} / {0[diff_rate]}'.format, axis=1)
        # data = data[['header', 'n/diff']]
        data = data.set_index('header').transpose().reset_index()
        data['barcode'] = id
        print(data)
        dfs.append(data)

    # Write the merged data to the output file
    df = pd.concat(dfs)
    barcodes = df.pop('barcode')
    df.insert(0, 'barcode', barcodes)
    df = df.rename(columns = {'index':'metric'}).sort_values(by = ['metric', 'barcode'] )
    df.to_csv(output_file, index=False)
    print(f"Merged data saved to '{output_file}'.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Merge CSV files based on the diff_rate column.')
    parser.add_argument('dir', help='Input directory containing CSV files')
    parser.add_argument('output', help='Output file for merged data')

    # Parse the arguments
    args = parser.parse_args()

    # Call the merge_csv_files function with the provided arguments
    merge_csv_files(args.dir, args.output)
