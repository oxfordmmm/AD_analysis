#!/usr/bin/env python3
import os
import sys
import pandas as pd
import argparse

def get_list(list_file, passed_only=True):
    df = pd.read_csv(list_file)
    if 'test normalized' not in df.columns:
        df['test normalized'] = df['test_type']
    df=df[df['test normalized'].isin(['BIOFIRE','ALINITY','CEPHEID'])].drop_duplicates(subset=['Run', 'barcode'])
    if passed_only:
        df=df[(df['run_pass']==True) & (df['PCs_passed']==1)]
    return df

def link_files(input_folder, output_file, df_list, set='validation'):
    """Create symbolic links for FASTQ files from a list to a directory."""
    
    
    for index, row in df_list.iterrows():
        run = row['Run']
        output_folder=os.path.join(output_file, f'{run}')
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        if set=='derivation':
            run=f'{run}.sup_barcodes_sup_fastq/'
        barcode = row['barcode']
        barcode=str(barcode).zfill(2)  # Ensure barcode is two digits
        dataFold = os.path.join(input_folder, run, 'data')

        if set=='derivation':
            dataFold = os.path.join(input_folder, run)
            fq_filename = f"SQK-RBK114-96_barcode{barcode}.fastq.gz"
        else:   
            fq_filename = f"barcode{barcode}.fastq.gz"
        fq_path = os.path.join(dataFold, fq_filename)
        
        if not os.path.isfile(fq_path):
            print(f"Warning: File {fq_path} does not exist. Skipping.", file=sys.stderr)
            continue
        
        link_name = os.path.join(output_folder, fq_filename)
        if not os.path.exists(link_name):
            os.symlink(fq_path, link_name)
        else:
            print(f"Link {link_name} already exists. Skipping.", file=sys.stderr)

def fq_cat(input_folder, output_file, df_list):
    """Concatenate multiple FASTQ files into a single FASTQ file."""
    with open(output_file, 'w') as outfile:
        for index, row in df_list.iterrows():
            run = row['Run']
            run=f'{run}.sup_barcodes_sup_fastq/'
            barcode = row['barcode']
            barcode=str(barcode).zfill(2)  # Ensure barcode is two digits
            #dataFold = os.path.join(input_folder, run, 'data')
            dataFold = os.path.join(input_folder, run)

            fq_filename = f"SQK-RBK114-96_barcode{barcode}.fastq.gz"
            fq_path = os.path.join(dataFold, fq_filename)
            
            if not os.path.isfile(fq_path):
                print(f"Warning: File {fq_path} does not exist. Skipping.", file=sys.stderr)
                continue
            
            #with open(fq_path, 'r') as infile:
            #    for line in infile:
            #        outfile.write(line)
            os.system(f'cat {fq_path} >> {output_file}')
        
                    
    print(f"Concatenation complete. Output written to {output_file}")

def main(input_folder, output_file, list_file, set='validation', passed_only=True):
    # Read the list file to get run and barcode information
    df_list = get_list(list_file, passed_only=passed_only)
    
    link_files(input_folder, output_file, df_list, set=set)
    #fq_cat(input_folder, output_file, df_list)




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Concatenate multiple FASTQ files into a single FASTQ file.")
    parser.add_argument('-i', '--input_folder', required=True, help="Input folder containing FASTQ files.")
    parser.add_argument('-o', '--output_file', required=True, help="Output FASTQ file.")
    parser.add_argument('-l', '--list_file', required=True, help="File contain the run and barcode of the files to concatenate.")
    parser.add_argument('-s', '--set', choices=['validation', 'derivation'], default='validation', help="Set type: validation or derivation.")
    parser.add_argument('-p', '--passed_only', action='store_true', help="Include only passed samples.")
    args = parser.parse_args()

    main(args.input_folder, args.output_file, args.list_file, set=args.set, passed_only=args.passed_only)