#!/usr/bin/env python3
# Scipt taken from Nick's crumpit2
import pandas as pd
import sys
from ete3 import NCBITaxa
ncbi = NCBITaxa()


def get_species(taxes):
	'''return species if above species taxon'''
	try:
		l=ncbi.get_lineage(taxes)
		r=ncbi.get_rank(l)
		reverse_dict={v: k for k, v in r.items()}
	except:
		return None
	if 'species' in reverse_dict:
		return str(reverse_dict['species'])
	else:
		return None


def run_batch_reads(read_tax):
	'''Take a list of reads and taxids and make files of reads in respective taxid.csv files.'''
	df=pd.read_csv(read_tax, sep=' ')
	df['taxID']=df['taxID'].map(str)
	df=df[df['taxID'].apply(lambda x: x.isnumeric())]
	df=df.query('taxID!="0"')
	df['species']=df['taxID'].map(get_species)
	df.dropna(inplace=True)
	species=df['species'].unique()
	for s in species:
		df2=df.query('species==@s')
		df2['readID'].to_csv(f'taxids/{s}.txt',index=False)
	return species


if __name__ == '__main__':
	run_batch_reads(sys.argv[1])