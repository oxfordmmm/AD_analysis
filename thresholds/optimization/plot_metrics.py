import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
sns.set_theme(style="darkgrid")


df=pd.read_csv(sys.argv[1])

cols=['Run',	'barcode',	'seq_name',	'pathogen',	'TP',	'FP',	'TN',	'FN',
       'Sample_reads_percent_of_run', 'Sample_reads_percent_of_refs',	'Sample_reads_percent_of_type_run', 'AuG_trunc10']
df=df[cols]

# melt dataframe by metrics
#df_melted = pd.melt(df, id_vars=['Run', 'barcode', 'seq_name', 'pathogen', 'TP', 'FP', 'TN', 'FN'], value_vars=['Sample_reads_percent_of_run', 'Sample_reads_percent_of_refs', 'Sample_reads_percent_of_type_run'], var_name='Metric', value_name='Metric Value')

# melt dataframe by 'TP', 'FP', 'TN', 'FN'
df_melted = pd.melt(df, id_vars=['Run', 'barcode', 'seq_name', 'pathogen','Sample_reads_percent_of_run', 'Sample_reads_percent_of_refs', 'Sample_reads_percent_of_type_run', 'AuG_trunc10'], value_vars=['TP', 'FP', 'TN', 'FN'], var_name='Test result', value_name='Test result Value')
df_melted=df_melted[df_melted['Test result Value']==1]

df_melted['AuG_trunc10/Sample_reads_percent_of_refs']=df_melted['Sample_reads_percent_of_refs']/df_melted['AuG_trunc10']

print(df_melted['AuG_trunc10/Sample_reads_percent_of_refs'].describe()) 
df_melted.to_csv('metrics_melted.csv', index=False)
# plot metrics

#g=sns.stripplot(x='Test result', y='Sample_reads_percent_of_refs', data=df_melted, hue='pathogen')
#plt.yscale('symlog')
#plt.show()

g2=sns.scatterplot(x='AuG_trunc10', y='Sample_reads_percent_of_refs', data=df_melted, hue='Test result')
g2.set_title('Derivation set, test_result does not incldue AND ratio')

# set axis limits
#g2.set(ylim=(0, 0.5))
#g2.set(xlim=(0, 1))
# draw a slope line with intercept 0 and slope 0.5
x = np.linspace(0,10,100)
y = 0.1*x
plt.plot(x, y, color='red')
plt.axhline(y=0.007, color='grey', linestyle='--')
plt.axvline(x=0.003, color='grey', linestyle='--')
plt.yscale('log')
plt.xscale('log')

plt.tight_layout()
plt.savefig('AuG_trunc10_vs_Sample_reads_percent_of_refs.pdf')
#plt.xscale('symlog')
plt.show()
plt.close()
plt.clf()

g3=sns.swarmplot(x='Test result', y='AuG_trunc10/Sample_reads_percent_of_refs', data=df_melted, hue='Run')
plt.yscale('symlog')
#plt.show()
plt.savefig('AuG_trunc10_vs_Sample_reads_percent_of_refs_swarm.pdf')
plt.close()
plt.clf()
