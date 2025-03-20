#!/usr/bin/env python3
import pandas as pd
import upsetplot as up
import matplotlib.pyplot as plt
import sys

df=df=pd.read_csv(sys.argv[1])
df=df[df['test_type'].isin(['BIOFIRE','ALINITY', 'CEPHEID'])]
cols=['Run', 'barcode','MS2_spike',	'IC_virus_spike', 'orthoreovirus passed',	'zika passed',	'MS2 passed',	'murine_respirovirus passed']
df=df[cols]
df.drop_duplicates(inplace=True, keep='first')

# remove spiked == 0
df=df[(df['MS2_spike'] != 0) & (df['IC_virus_spike']!=0)]

# df long to wide with pathogen_reduced and prediction as columns and values
#df_wide=df.pivot(index=['Run', 'barcode'],columns='pathogen_reduced', values='prediction')

# group to generate counts
pathogens=[ 'orthoreovirus passed',	'zika passed',	'MS2 passed',	'murine_respirovirus passed']
g=df.groupby(pathogens).size().reset_index().rename(columns={0:'count'})

# convert g to series
g=g.set_index(pathogens)['count']

# plot
up.plot(g, show_counts=True)
#plt.show()
plt.savefig('spikes_upset_plot.pdf')

