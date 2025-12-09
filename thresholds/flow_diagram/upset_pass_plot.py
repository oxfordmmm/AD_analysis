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
df=df[df['Run']!='010524_Expt9_SISPA_Daisy']
cols=['Run', 'barcode', 'OR pass',	'AND ratio pass',	'2 reads pass']
df=df[cols]
#df.drop_duplicates(inplace=True, keep='first')

# remove spiked == 0

# df long to wide with pathogen_reduced and prediction as columns and values
#df_wide=df.pivot(index=['Run', 'barcode'],columns='pathogen_reduced', values='prediction')

# group to generate counts
pass_types=[ 'OR pass',	'AND ratio pass',	'2 reads pass']
# rename pass types for plotting
pt_rename={'OR pass':'Main criteria (Table 2)', 
           'AND ratio pass':'Ratio criteria (Table 2)', 
           '2 reads pass':'>=2 reads'}
df.rename(columns=pt_rename, inplace=True)
pass_types=[pt_rename[pt] for pt in pass_types]

g=df.groupby(pass_types).size().reset_index().rename(columns={0:'count'})




# convert g to series
g=g.set_index(pass_types)['count']

# plot
up.plot(g,sort_categories_by='input', show_counts=True)
#plt.show()
plt.savefig('pass_upset_plot.pdf')

