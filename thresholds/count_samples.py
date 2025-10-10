#!/usr/bin/env python3
import pandas as pd
import numpy as np
import sys
import os

df=pd.read_csv(sys.argv[1])

# Replace NULL with np.nan
df.replace('NULL', np.nan, inplace=True)
df.replace('NA', np.nan, inplace=True)
df.replace('nan', np.nan, inplace=True)


## explode colulmn 'iord_pcr_organisms' on :
#df['iord_pcr_organisms'] = df['iord_pcr_organisms'].str.split(':')
#df = df.explode('iord_pcr_organisms')
#print(df)

def split_tests_from_IORD(row):
    s=row['iord_pcr_organisms'].split(';') if pd.notna(row['iord_pcr_organisms']) else []
    for i in s:
        if ':' in i:
            i_split = i.split(':')
            if i_split[0] not in row:
                row[i_split[0]] = i_split[-1]
            elif row[i_split[0]] is np.nan:
                row[i_split[0]] = i_split[-1]
            else: # if there is already a value, append with ;
                if i_split[-1] not in row[i_split[0]].split(';'):
                    row[i_split[0]] = row[i_split[0]] + ';' + i_split[-1]
            #row[i_split[0]]=i_split[-1]
        else:
            row[i]=None
    return row

df = df.apply(split_tests_from_IORD, axis=1)
# remove empty columns
df = df.dropna(axis=1, how='all')

# add columns for each test
test_dict = {'FLRA': 'Alinity',	
             'FLRP': 'Ceipheid',
             'RPCR': 'BioFire'}

for col in ['FLRA', 'FLRP', 'RPCR']:
    df[test_dict[col]] = np.where(df[col].notna(), 1, 0)

# merge with metaDF.csv
metaDF = pd.read_csv('metaDF.csv')
metaDF=metaDF[metaDF['Biofire positive']==1]
metaDF['pathogen'] = metaDF['pathogen'].map(str)
metaDF=metaDF.groupby(['Run','barcode'])['pathogen'].apply(','.join).reset_index()
#metaDF.drop_duplicates(subset=['Run','barcode'], inplace=True)
df = pd.merge(metaDF, df, left_on=['Run', 'barcode'],right_on=['run_name', 'barcode'], how='left')

# count number of unique rows
print(f'Number of samples in dataset: {len(df)}')
# Count the number of unique combinations of accession_1 and accession_2
unique_samples = df[['accession_1', 'accession_2']].drop_duplicates()
print(f'Number of unique samples in dataset: {len(unique_samples)}')

# count the number of samples where pooled_with_accession_1 and pooled_with_accession_2 are not null
pooled_samples = df[~df['pooled_with_accession_1'].isna() & ~df['pooled_with_accession_2'].isna()]
print(f'Number of pooled samples in dataset: {len(pooled_samples)}')

# count the number of pathogens per sample
df['num_pathogens'] = df['pathogen'].str.split(',').apply(lambda x: len(x) if isinstance(x, list) else 0)
df['num_pathogens'] = np.where(df['pathogen']=='nan', 0, df['num_pathogens'])

print(f'Number of samples with no pathogens: {len(df[df["num_pathogens"] == 0])}')
print(f'Number of samples with one pathogen: {len(df[df["num_pathogens"] == 1])}')
print(f'Number of samples with more than one pathogen: {len(df[df["num_pathogens"] > 1])}')


# count the number of unique samples in pooled_with_accession_1 and pooled_with_accession_2
unique_pooled_samples = pooled_samples[['pooled_with_accession_1', 'pooled_with_accession_2']].drop_duplicates()
print(f'Number of unique pooled samples in dataset: {len(unique_pooled_samples)}')

# count number of Alinity, Cepheid, BioFire tests
num_Alinity_tests = df['Alinity'].sum()
num_Cepheid_tests = df['Ceipheid'].sum()
num_BioFire_tests = df['BioFire'].sum()
print(f'Number of Alinity tests in dataset: {num_Alinity_tests}')
print(f'Number of Cepheid tests in dataset: {num_Cepheid_tests}')
print(f'Number of BioFire tests in dataset: {num_BioFire_tests}')
num_had_tests = df[['Alinity', 'Ceipheid', 'BioFire']].any(axis=1).sum()
print(f'Number of samples with at least one test in dataset: {num_had_tests}')

# remove duplicate accessions
df2 = df.drop_duplicates(subset=['accession_1', 'accession_2'])
df2 = df2.copy()
df2 = df2[df2[['Alinity', 'Ceipheid', 'BioFire']].any(axis=1)]
num_unique_had_tests = len(df2)
print(f'Number of unique samples with at least one test in dataset: {num_unique_had_tests}')
num_Alinity_tests = df2['Alinity'].sum()
num_Cepheid_tests = df2['Ceipheid'].sum()
num_BioFire_tests = df2['BioFire'].sum()
print(f'Number of unique Alinity tests in dataset: {num_Alinity_tests}')
print(f'Number of unique Cepheid tests in dataset: {num_Cepheid_tests}')
print(f'Number of unique BioFire tests in dataset: {num_BioFire_tests}')

# print number of pathogens per sample
print(f'Number of unique samples with no pathogens: {len(df2[df2["num_pathogens"] == 0])}')
print(f'Number of unique samples with one pathogen: {len(df2[df2["num_pathogens"] == 1])}')
print(f'Number of unique samples with more than one pathogen: {len(df2[df2["num_pathogens"] > 1])}')

# caluclate number of tests per sample
df['num tests'] = df[['Alinity', 'Ceipheid', 'BioFire']].sum(axis=1)

df.to_csv('count_samples_expanded.csv', index=False)