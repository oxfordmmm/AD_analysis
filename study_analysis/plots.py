#!/usr/bin/env python3
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import sys

df=pd.read_csv(sys.argv[1])

df['CT value'] = np.where(df['ct_1']!=0, df['ct_1'], df['ct_2'])
df2=df[df['CT value']!=0]
df2=df2.copy()

# bin CY values by ranges of 5
bins = [0,5,10,15,20,25,30]
df2['CT bins'] = pd.cut(df2['CT value'], bins, right=False)
#bin_order = ['[0, 5)', '[5, 10)', '[10, 15)', '[15, 20)', '[20, 25)', '[25, 30)']
#rev_bin_order = bin_order[::-1]

# plot reads over bins of CT values
df2['sample num reads'] = np.where(df2['sample num reads']==0, 0.9, df2['sample num reads'])
df2['pass'] = np.where(df2['pass']=='0', 'No reads', df2['pass'])
g1=sns.swarmplot(x='CT bins', y='sample num reads', data=df2, hue='pass')
g1.invert_xaxis()
plt.yscale('log')
plt.savefig('reads_vs_CT.pdf')
plt.clf()

# plot reads over bins of CT values with pathogen as hue
g1=sns.swarmplot(x='CT bins', y='sample num reads', data=df2, hue='pathogen')
g1.invert_xaxis()
plt.yscale('log')
plt.savefig('reads_vs_CT_pathogen.pdf')
plt.clf()

# plot reads over bins of CT values with pathogen as hue
g1=sns.scatterplot(x='CT value', y='sample num reads', data=df2, hue='pathogen')
g1.invert_xaxis()
plt.yscale('log')
plt.savefig('reads_vs_CT_pathogen_scatter.pdf')
plt.clf()

# plot meanDepth_trunc5 over bins of CT values with pathogen as hue
g1=sns.swarmplot(x='CT bins', y='meanDepth_trunc5', data=df2, hue='pathogen')
g1.invert_xaxis()
#plt.yscale('log')
plt.savefig('meanDepth_trunc5_vs_CT_pathogen.pdf')
plt.clf()

# plot meanDepth_trunc5 over bins of CT values with pathogen as hue
g1=sns.swarmplot(x='CT bins', y='Cov1_perc', data=df2, hue='pathogen')
g1.invert_xaxis()
#plt.yscale('log')
plt.savefig('Cov1_perc_vs_CT_pathogen.pdf')
plt.clf()