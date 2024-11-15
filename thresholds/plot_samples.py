#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df=pd.read_csv('metaDF_biofire_only.csv')
g=df.groupby(['pathogen'])[['Biofire positive']].sum()
g2=df.groupby(['pathogen'])[['Biofire positive']].count()
g['total']=g2['Biofire positive']
g['Biofire negative']=g['total']-g['Biofire positive']

g[['Biofire positive','Biofire negative']].plot(kind='bar',stacked=True)
plt.tight_layout()
plt.savefig('pathogens.pdf')