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

# plot reads over scatter of CT values with run name as hue
g1=sns.scatterplot(x='CT value', y='sample num reads', data=df2, hue='Run')
g1.invert_xaxis()
plt.yscale('log')
# move legend outside of plot
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.tight_layout()
plt.savefig('reads_vs_CT_run_scatter.pdf')
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

## Bar plots for each run/batch for sensitivity
# g1=df.groupby(['Run','Batch'])[['gold_standard']].sum().reset_index()
# df['pass'] = np.where(df['pass']=='0', 0, df['pass'])  # Convert 'pass' to numeric
# df['pass'] = np.where(df['pass']=='False', 0, df['pass'])
# df['pass'] = np.where(df['pass']=='True', 1, df['pass'])
# df['pass'] = df['pass'].astype(int)  # Ensure 'pass' is integer
# df['TP']= np.where((df['pass']==1) & (df['gold_standard']==1), 1, 0)  # Create a new column 'TP' for true positives
# g2=df.groupby(['Run','Batch'])[['TP']].sum().reset_index()
# # Merge the two dataframes on 'Run' and 'batch'
# g1 = g1.rename(columns={'gold_standard': 'gold_standard_count'})
# g2 = g2.rename(columns={'TP': 'pass_count'})
# g1 = g1.merge(g2, on=['Run', 'Batch'])
# g1['Total'] = g1['gold_standard_count'] - g1['pass_count']
# g1.drop(columns=['gold_standard_count'], inplace=True)
# print(g1)
# # plot stack bar plot of gold_standard and pass
# g1 = g1.set_index(['Run', 'Batch'])
# g1.plot(kind='bar', stacked=True, figsize=(10, 6))
# plt.title('Gold Standard vs Pass Counts by Run and Batch')
# plt.xlabel('Run and Batch')
# plt.ylabel('Count')
# plt.xticks(rotation=90)
# plt.tight_layout()
# plt.savefig('gold_standard_vs_pass_counts_by_run_and_batch.pdf')
# plt.clf()

## Bar plots for each run/batch for sensitivity, but only PCs_passed == 1
g0=df.groupby(['Run','Batch'])[['gold_standard']].sum().reset_index()
df_original = df.copy()  # Keep a copy of the original DataFrame
df = df[df['PCs_passed'] == 1]  # Filter for PCs_passed == 1
g1=df.groupby(['Run','Batch'])[['gold_standard']].sum().reset_index()
df['pass'] = np.where(df['pass']=='0', 0, df['pass'])  # Convert 'pass' to numeric
df['pass'] = np.where(df['pass']=='False', 0, df['pass'])
df['pass'] = np.where(df['pass']=='True', 1, df['pass'])
df['pass'] = df['pass'].astype(int)  # Ensure 'pass' is integer
df['TP']= np.where((df['pass']==1) & (df['gold_standard']==1), 1, 0)  # Create a new column 'TP' for true positives
g2=df.groupby(['Run','Batch'])[['TP']].sum().reset_index()
# Merge the two dataframes on 'Run' and 'batch'
g1 = g1.rename(columns={'gold_standard': 'gold_standard_count'})
g2 = g2.rename(columns={'TP': 'passed'})
g0 = g0.rename(columns={'gold_standard': 'all_gold_standard_count'})
g1 = g1.merge(g2, on=['Run', 'Batch'])
g1 = g1.merge(g0, on=['Run', 'Batch'])
g1['Failed criteria'] = g1['gold_standard_count'] - g1['passed']
g1['Failed PCs'] = g1['all_gold_standard_count'] - (g1['passed'] + g1['Failed criteria'])

print(g1)
g1.drop(columns=['gold_standard_count','all_gold_standard_count'], inplace=True)
print(g1)

# Sort the DataFrame by list of run names
run_order = {'AD_winter_study_201224':1,
             'AD_winter_study_220125':2,
             'AD_winter_study_070325':3,
             'AD_winter_study_170325':4,
             'AD_winter_study_100425':5,
             'AD_winter_study_160725':6,
             'AD_winter_study_220725':7,
             'AD_winter_study_240725':8,
             'AD_winter_study_300725':9,
             'AD_winter_study_010825':10,
             'AD_winter_study_130825_rpt050825':11
  }  # Define the order of runs

g1['Run_order'] = g1['Run'].map(run_order)  # Map the run names to their order
g1 = g1.sort_values(by=['Run_order', 'Batch'])  # Sort by Run order and Batch
g1 = g1.drop(columns=['Run_order'])  # Drop the temporary Run_order column

# plot stack bar plot of gold_standard and pass
g1 = g1.set_index(['Run', 'Batch'])
g1.plot(kind='bar', stacked=True, 
        figsize=(10, 6))
# make y scale full integers
#plt.ylim(0, g1['passed'].max())  # Adjust the y-axis limit for better visibility
plt.title('Gold Standard vs Pass Counts by Run and Batch')
plt.xlabel('Run and Batch')
plt.ylabel('Pathogen Count')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig('gold_standard_vs_pass_counts_by_run_and_batch_PCs_passed.pdf')
plt.clf()

print(g1['passed'].sum())
print(g1['Failed criteria'].sum())
print(g1['Failed PCs'].sum())

g1['sensitivity'] = g1['passed'] / (g1['passed'] + g1['Failed criteria'])
print(g1)

## Bar plots for each run/batch/pathogen for sensitivity, but only PCs_passed == 1
groups=['pathogen']
df = df_original.copy()  # Use the original DataFrame for grouping
g0=df.groupby(groups)[['gold_standard']].sum().reset_index()
df = df[df['PCs_passed'] == 1]  # Filter for PCs_passed == 1
g1=df.groupby(groups)[['gold_standard']].sum().reset_index()
df['pass'] = np.where(df['pass']=='0', 0, df['pass'])  # Convert 'pass' to numeric
df['pass'] = np.where(df['pass']=='False', 0, df['pass'])
df['pass'] = np.where(df['pass']=='True', 1, df['pass'])
df['pass'] = df['pass'].astype(int)  # Ensure 'pass' is integer
df['TP']= np.where((df['pass']==1) & (df['gold_standard']==1), 1, 0)  # Create a new column 'TP' for true positives
g2=df.groupby(groups)[['TP']].sum().reset_index()
# Merge the two dataframes on 'Run' and 'batch'
g1 = g1.rename(columns={'gold_standard': 'gold_standard_count'})
g2 = g2.rename(columns={'TP': 'passed'})
g0 = g0.rename(columns={'gold_standard': 'all_gold_standard_count'})
g1 = g1.merge(g2, on=groups)
g1 = g1.merge(g0, on=groups)
g1['Failed criteria'] = g1['gold_standard_count'] - g1['passed']
g1['Failed PCs'] = g1['all_gold_standard_count'] - (g1['passed'] + g1['Failed criteria'])

print(g1)
g1.drop(columns=['gold_standard_count','all_gold_standard_count'], inplace=True)
print(g1)
print(g1.index)

# remove 'unmapped' row if it exists
g1= g1[g1['pathogen'] != 'unmapped']

g1 = g1.set_index(groups)

# plot stack bar plot of gold_standard and pass
g1.plot(kind='bar', stacked=True, 
        figsize=(10, 6))
# make y scale full integers
#plt.ylim(0, g1['passed'].max())  # Adjust the y-axis limit for better visibility
plt.title('Gold Standard vs Pass Counts by pathogen')
plt.xlabel('pathogen')
plt.ylabel('Pathogen Count')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig('pathogoen_pass_fails.pdf')
plt.clf()

print(g1['passed'].sum())
print(g1['Failed criteria'].sum())
print(g1['Failed PCs'].sum())

g1['sensitivity'] = g1['passed'] / (g1['passed'] + g1['Failed criteria'])
print(g1)


# plot gold_standard pathogen counts by run and batch
g1 = df.groupby(['Run', 'Batch', 'pathogen'])[['gold_standard']].sum().reset_index()
g1 = g1.set_index(['Run', 'Batch', 'pathogen'])

g1 = g1.unstack(level='pathogen')

g1.columns = g1.columns.droplevel(0)  # Drop the top level of the MultiIndex
g1 = g1.fillna(0)  # Fill NaN values with 0
#g1.reset_index(inplace=True)  # Reset index to make it easier to plot
# sort dataframe by first index in MultiIndex by run_order
g1['Run_order'] = g1.index.get_level_values('Run').map(run_order)  # Map the run names to their order
g1 = g1.sort_values(by='Run_order')  # Sort by Run order
g1 = g1.drop(columns=['Run_order'])  # Drop the temporary Run_order column      
# drop 'unmapped' column if it exists
if 'unmapped' in g1.columns:
    g1 = g1.drop(columns=['unmapped'])
# drop empty columns
g1 = g1.loc[:, (g1 != 0).any(axis=0)]

print(g1)
g1.plot(kind='bar', stacked=True, 
        colormap='tab20',  # Use a colormap for better color differentiation
        figsize=(10, 6))
plt.title('Gold Standard Pathogen Counts by Run and Batch')
plt.xlabel('Run and Batch')
plt.ylabel('Pathogen Count')
plt.xticks(rotation=90)
# move legend outside of plot
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.tight_layout()
plt.savefig('gold_standard_pathogen_counts_by_run_and_batch.pdf')
plt.clf()