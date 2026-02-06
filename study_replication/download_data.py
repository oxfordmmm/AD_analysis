#!/usr/bin/env python3
import pandas as pd
import os


def read_ENA_accession(file_path):
    df = pd.read_csv(file_path)
    df['Run'] = df['alias'].str.split('_barcode').str[0]
    df['barcode'] = df['alias'].str.split('_barcode').str[1]
    return df

# ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR158/098/ERR15836098/ERR15836098_1.fastq.gz

def download_data(ena_data):
    # makes runs folders if they don't exist
    runs= ena_data['Run'].unique()
    for run in runs:
        run_folder = os.path.join('workflow_runs', run, 'data')
        os.makedirs(run_folder, exist_ok=True)

    for index, row in ena_data.iterrows():
        run = row['Run']
        barcode = row['barcode']
        accession = row['accession']
        url = f"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{accession[:6]}/{accession[-3:]}/{accession}/{accession}.fastq.gz"
        output_file = os.path.join('workflow_runs', run, 'data', f"barcode{barcode}.fastq.gz")
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        # download the file using wget
        os.system(f"wget -O {output_file} {url}")

if __name__ == "__main__":
    # Define the path to the CSV file
    csv_file_path = os.path.join('meta_data', 'ENA_accessions.csv')
    
    # Read the ENA accession data
    ena_data = read_ENA_accession(csv_file_path)
    
    # download the data into Run folders
    download_data(ena_data)
    
    