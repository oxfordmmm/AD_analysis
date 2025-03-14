#!/usr/bin/env python3
import pandas as pd
import upsetplot as up
import matplotlib.pyplot as plt
import sys

df=df=pd.read_csv(sys.argv[1])
df=df[df['test_type'].isin(['BIOFIRE','ALINITY', 'CEPHEID'])]
df=df[(df['MS2_spike'] != 0) & (df['IC_virus_spike']!=0)]
df=df[df['PCs_passed']==1]
df=df[df['gold_standard']==1]
cols=['Run', 'barcode', 'OR pass',	'AND ratio pass',	'2 reads pass']
df=df[cols]
#df.drop_duplicates(inplace=True, keep='first')

# remove spiked == 0

# df long to wide with pathogen_reduced and prediction as columns and values
#df_wide=df.pivot(index=['Run', 'barcode'],columns='pathogen_reduced', values='prediction')

# group to generate counts
pass_types=[ 'OR pass',	'AND ratio pass',	'2 reads pass']
g=df.groupby(pass_types).size().reset_index().rename(columns={0:'count'})

# convert g to series
g=g.set_index(pass_types)['count']

# plot
up.plot(g, show_counts=True)
#plt.show()
plt.savefig('pass_upset_plot.pdf')

