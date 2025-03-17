#!/usr/bin/env pythhon3
import pandas as pd
import numpy as np

biorifre_organisms=['Adenovirus',
'Coronavirus 229E',
'Coronavirus HKU',
'Coronavirus NL63',
'Coronavirus OC43',
'MERS',
'SARS-CoV-2',
'Influenza A',
'Influenza A/H1-2009',
'Influenza A/H1',
'Influenza A/H3',
'Influenza B',
'Rhinovirus/enterovirus',
'Metapneumovirus',
'Parainfluenza 1',
'Parainfluenza 2',
'Parainfluenza 3',
'Parainfluenza 4',
'RSV',
'Bordetella_parapertussis',
'Bordetella_pertussis',
'Chlamydia_pneumoniae',
'Mycoplasma_pneumoniae']

alinity_cephid_organisms=['Influenza A', 'Influenza B', 'RSV', 'SARS-CoV-2']

biofire_additional=[org for org in biorifre_organisms if org not in alinity_cephid_organisms]

df=pd.read_csv('biofire_results_merged.csv')
# repopulate df2 PC_passes with results from the same Run and barcode if PC_passes is 1
df['PCs_passed']=df.groupby(['Run', 'barcode'])['PCs_passed'].transform('max')

for col in ['orthoreovirus', 'zika',	'MS2',	'murine_respirovirus',	
            'orthoreovirus passed',	'zika passed',	'MS2 passed',	'murine_respirovirus passed',
            'total sample bases', 'total sample reads', 'PCs_passed', 'PC_PASSES']:
    df[col]=df.groupby(['Run', 'barcode'])[[col]].transform('max')

for col in ['total run bases',	'total run reads']:
    df[col]=df.groupby(['Run'])[[col]].transform('max')

df['sample_positive']=df.groupby(['Run', 'barcode'])['gold_standard'].transform('max')

# count pathogen columns by Run barcode 
df_PCR=pd.read_csv('pcr_organisms.csv')
df=df.merge(df_PCR, on=['Run','barcode'], how='left')
# if 'testcode_checked' == 'FLRA; C2VP' then remove extra biofire organisms
not_biofire_codes=['FLRA; C2VP', 'FLRP; C2VP']
# get potential additional yield 
df_ay=df[~((df['pathogen'].isin(biorifre_organisms)) & (df['testcode_checked'].isin(not_biofire_codes)))]
df_ay=df_ay[~df_ay['pathogen'].isin(alinity_cephid_organisms)]
df_ay.to_csv('additional_yield.csv', index=False)

# remove biofire tests from alinity cephid
df=df[~((df['pathogen'].isin(biofire_additional)) & (df['testcode_checked'].isin(not_biofire_codes)))]

df['pathogen_tests']=df.groupby(['Run', 'barcode'])[['pathogen']].transform('count')

df['test_type']=np.where(df['pathogen_tests']==24, 'Biofire', 'Alinity Cepheid')

# add in negative meta
negative_meta=pd.read_csv('negatives_meta.csv')
df=df.merge(negative_meta, on=['Run', 'barcode'], how='left')


# readjust spiked values for derivation set
df['PC_PASSES']=df['orthoreovirus passed']+df['zika passed']+df['murine_respirovirus passed']
df['PCs_passed']=np.where(df['PC_PASSES']>1, 1, 0)

df.to_csv('biofire_results_merged_adjusted.csv', index=False)

# remove negative controls
df=df[df['spiked']>=1]


# count number of tests by type
g=df.groupby(['test_type'])['pathogen_tests'].describe()
print(g)

# count number of samples by test type
bf=df[df['test_type']=='Biofire']
bf=bf.drop_duplicates(subset=['Run', 'barcode'])
ac=df[df['test_type']=='Alinity Cepheid']
ac=ac.drop_duplicates(subset=['Run', 'barcode'])
print('Number of Biofire samples:', bf.shape[0])
print('Number of Alinity Cepheid samples:', ac.shape[0])
total_samples=bf.shape[0]+ac.shape[0]
print(f'Total number of samples: {total_samples}')
unique_runs=df['Run'].nunique()
print('Total number of runs:', unique_runs )

# count number of samples that failed run/sample controls, higher than 426.6 MB for total run bases and 30,000 for total sample reads
dfF=df[df['total run bases']<400000000]
df=df[df['total run bases']>=400000000]
unique_runs=df['Run'].nunique()
failed_runs=dfF['Run'].nunique()
print(f'Number of runs with total run bases lower than 400 MB: {failed_runs}')
print(f'Number of runs with total run bases higher than 400 MB: {unique_runs}')
min_reads=25000
dfSF=df[df['total sample reads']<min_reads]
dfSunique=dfSF.drop_duplicates(subset=['Run', 'barcode'])
print(f'Number of samples with total sample reads lower than {min_reads:,}: {dfSunique.shape[0]}')
df=df[df['total sample reads']>=min_reads]
dfSunique=df.drop_duplicates(subset=['Run', 'barcode'])
print(f'Number of samples with total sample reads higher than {min_reads:,}: {dfSunique.shape[0]}')

print(f'Percentage of samples with total sample reads higher than {min_reads:,}: {dfSunique.shape[0]/total_samples*100:.2f}%')

# count number of samples that failed negative controls

# remove 010524_Expt9_SISPA_Daisy as barcode 74 contailed 7 FluA reads
df2=df[df['Run']!='010524_Expt9_SISPA_Daisy'] 

bf=df2[df2['test_type']=='Biofire']
bf=bf.drop_duplicates(subset=['Run', 'barcode'])
ac=df2[df2['test_type']=='Alinity Cepheid']
ac=ac.drop_duplicates(subset=['Run', 'barcode'])

total_samples_passing_BNC=bf.shape[0]+ac.shape[0]
print(f'Total number of samples passing batch negative controls: {total_samples_passing_BNC}')
print(f'Percentage of samples passing batch negative controls: {total_samples_passing_BNC/total_samples*100:.2f}%')

# count number of samples that failed positive controls
df3=df2[df2['PCs_passed']==0]
df3=df3.copy()
df3.to_csv('failed_positive_controls.csv', index=False)
df3.drop_duplicates(subset=['Run', 'barcode'], inplace=True)
print(f'Number of samples failing positive controls: {df3.shape[0]}')
print(f'Percentage of samples failing positive controls: {df3.shape[0]/total_samples*100:.2f}%')

# number of PCs that passed
print('PC_PASSES value count. 1=2 failed, 0=3 failed:')
print(df3['PC_PASSES'].value_counts())

bf=df3[df3['test_type']=='Biofire']
bf=bf.drop_duplicates(subset=['Run', 'barcode'])
ac=df3[df3['test_type']=='Alinity Cepheid']
ac=ac.drop_duplicates(subset=['Run', 'barcode'])
total_samples_failing_BNC=bf.shape[0]+ac.shape[0]

cMG_tests=(bf.shape[0]*len(biorifre_organisms))+(ac.shape[0]*len(alinity_cephid_organisms))
print(f'Total number of cMG tests for failed positive control samples: {cMG_tests}')

#print(f'Total number of samples failing batch negative controls: {total_samples_failing_BNC}')
#print(f'Percentage of samples failing batch negative controls: {total_samples_failing_BNC/total_samples*100:.2f}%')

# count number of samples that passed positive controls
df4=df2[df2['PCs_passed']==1]
df4=df4.copy()
df5=df4.drop_duplicates(subset=['Run', 'barcode'])
print(f'Number of samples passing positive controls: {df5.shape[0]}')
df5p=df5[df5['sample_positive']>=1]
print(f'Number of positive samples passing positive controls: {df5p.shape[0]}')
print(f'Percentage of samples passing positive controls of total samples: {df5.shape[0]/total_samples*100:.2f}%')
print(f'Percentage of samples passing positive controls of samples passing negative controls: {df5.shape[0]/total_samples_passing_BNC*100:.2f}%')


# count the number of samples pathogens detected in gold_standard
df5g=df4[df4['gold_standard']>=1]
df5g=df5g.copy()
df5=df5g[df5g['sample num reads']>=1]
#df5.drop_duplicates(subset=['Run', 'barcode'], inplace=True)
print(f'x pathogens identified by cMG: {df5.shape[0]}')
print(f'Percentage of samples with pathogens identified by cMG: {df5.shape[0]/df5g.shape[0]*100}%')

df5n=df4[df4['gold_standard']==0]
df5n0=df5n[df5n['sample num reads']==0]
print(f'Number of samples with no pathogens identified by cMG: {df5n0.shape[0]}')
print(f'Percentage of samples with no pathogens identified by cMG: {df5n0.shape[0]/df5n.shape[0]*100}%')


# count number of samples with > pathogen detected in gold_standard
df6=df4[df4['gold_standard']>=1]
df6=df6.copy()
df6['polymicrobial']=df6.duplicated(subset=['Run', 'barcode'], keep=False)
df7=df6[df6['polymicrobial']==True]
df7=df7.copy()
df7=df7[~df7['pathogen'].isin(['Influenza A','Influenza A/H1-2009'])]
print(f'Number of samples with > pathogen detected in gold_standard: {df7.shape[0]}')
df7.to_csv('more_than_1_pathogen_goldstandard.csv', index=False)

df7n=df6[df6['polymicrobial']==False]
print(f'Number of samples with 1 pathogen detected in gold_standard: {df7n.shape[0]}')


# count the number of samples where gold_standard is 1
df8=df4[df4['gold_standard']==1]
df8.to_csv('gold_standard_positives.csv', index=False)
df8flu=df8[df8['pathogen'].isin(['Influenza A','Influenza A/H1-2009'])]
df8flu=df8flu.copy()
df8flu['polyFlu_classification']=df8flu.duplicated(subset=['Run', 'barcode'], keep=False)
df8polyfly=df8flu[df8flu['polyFlu_classification']==True]
df8polyfly=df8polyfly.copy()
df8polyfly.drop_duplicates(subset=['Run', 'barcode'], inplace=True)
print(f'Total Z pathogens identified from routine laboratory testing: {df8.shape[0]}')
print(f'Number of samples with >1 Influenza detected in gold_standard: {df8polyfly.shape[0]}')

# count the number of samples where gold_standard is 0
df9=df4[df4['gold_standard']==0]
print(f'Total W pathogens not identified from routine laboratory testing: {df9.shape[0]}')

bf=df4[df4['test_type']=='Biofire']
bf=bf.drop_duplicates(subset=['Run', 'barcode'])
ac=df4[df4['test_type']=='Alinity Cepheid']
ac=ac.drop_duplicates(subset=['Run', 'barcode'])
total_passed_control_criteria=bf.shape[0]+ac.shape[0]
number_if_biofire=total_passed_control_criteria*23
actual_number_tested_for=(bf.shape[0]*23)+(ac.shape[0]*4)
X_not_tested_for=number_if_biofire-actual_number_tested_for
print(f'Total X pathogens not tested with routine laboratory testing: {X_not_tested_for}')

# count the number of samples where pass is 1
#print(df4['pass'].unique()) 
df10=df4[df4['pass']=='True']
print(f'Total X samples passing gold standard: {df10.shape[0]}')
print(f'Percentage of samples passing gold standard: {df10.shape[0]/df8.shape[0]*100:.2f}%')

# Count the number of samples where pass is 0 and gold_standard is 0
df11=df4[(df4['pass']!='True') & (df4['gold_standard']==0)]
print(f'Total x samples not identified by cMG: {df11.shape[0]}')
print(f'Percentage of samples not passing gold standard: {df11.shape[0]/df9.shape[0]*100}%')