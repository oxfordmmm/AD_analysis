#!/usr/bin/env python3
# From Nick's crumpit2
import pandas as pd
from argparse import ArgumentParser
import sys


def get_refUrl(df):
	'''return URL from dataframe which has already been filtered to taxid'''
	if len(df)==1:
		return df['ftp_path'].unique()[0]
	elif len(df)>1:
		if "reference genome" in df["refseq_category"].unique():
			df=df.query("refseq_category=='reference genome'")
			return df['ftp_path'].unique()[0]
		elif "representative genome" in df["refseq_category"].unique():
			df=df.query("refseq_category=='representative genome'")
			return df['ftp_path'].unique()[0]
		else:
			df=df.sample(1)
			return df['ftp_path'].unique()[0]


def getRsyncURL(url_path):
	'''change input ncbi FTP url to rsync compatible'''
	prefix=url_path.split('/')[-1]
	suffix='_genomic.fna.gz'
	newurl=url_path.replace('https','rsync')
	url_path=f"{newurl}/{prefix}{suffix}"
	return url_path.replace("\n",'')

def run_get_ref(taxid: str, bacteria_file: str, virus_file: str):
	'''Take a taxid and two assembly files, return the url to download'''
	df=pd.read_csv(bacteria_file, sep='\t', skiprows=1)
	df['species_taxid']=df['species_taxid'].map(str)
	df.dropna(subset='ftp_path',inplace=True)
	df=df.query('species_taxid==@taxid')
	if len(df)>=1:
		url=get_refUrl(df)
		if url == None: 
			print("No correct url, exciting")
			sys.exit()
		url=getRsyncURL(url)
		with open('ref_path.txt', 'wt') as outfile:
			outfile.write(url)
		return url
	
	df=pd.read_csv(virus_file, sep='\t', skiprows=1)
	df['species_taxid']=df['species_taxid'].map(str)
	df.dropna(subset='ftp_path',inplace=True)
	df=df.query('species_taxid==@taxid')
	if len(df)>=1:
		url=get_refUrl(df)
		if url == None: 
			print("No correct url, exciting")
			sys.exit()
		url=getRsyncURL(url)
		with open('ref_path.txt', 'wt') as outfile:
			outfile.write(url)
		return url


if __name__ == '__main__':
	parser = ArgumentParser(description='')
	parser.add_argument('-t', '--taxid', required=True, 
                             help='taxid to find reference for')
	parser.add_argument('-b', '--bacteria_file', required=True,
                             help='bacteria assembly meta data file')
	parser.add_argument('-v', '--viral_file', required=True, 
                             help='virus assembly meta data file')

	opts, unknown_args = parser.parse_known_args()

	run_get_ref(str(opts.taxid), opts.bacteria_file, opts.viral_file)