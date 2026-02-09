#!/usr/bin/env python3
import pandas as pd
import os


def read_ENA_accession(file_path):
    df = pd.read_csv(file_path)
    df['Run'] = df['run_name']
    #df['barcode'] = df['alias'].str.split('_barcode').str[1]
    return df

# https://ftp.sra.ebi.ac.uk/vol1/run/ERR158/ERR15818380/reads.fastq.gz

def download_data(ena_data):
    # makes runs folders if they don't exist
    runs= ena_data['Run'].unique()
    for run in runs:
        run_folder = os.path.join('workflow_runs', run, 'data')
        os.makedirs(run_folder, exist_ok=True)

    for index, row in ena_data.iterrows():
        run = row['Run']
        barcode = row['barcode']
        accession = row['run_accession']
        # get ERS accession from the ERR accession
        # construct the URL for the FASTQ file
        url = f"ftp://ftp.sra.ebi.ac.uk/vol1/run/{accession[:6]}/{accession}/reads.fastq.gz"
        output_file = os.path.join('workflow_runs', run, 'data', f"barcode{barcode:02d}.fastq.gz")
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        # download the file using wget
        os.system(f"wget -O {output_file} {url}")

if __name__ == "__main__":
    # Define the path to the CSV file
    csv_file_path = os.path.join('meta_data', 'ena_webin_reports.csv')
    
    # Read the ENA accession data
    ena_data = read_ENA_accession(csv_file_path)
    
    # download the data into Run folders
    download_data(ena_data)
    
    