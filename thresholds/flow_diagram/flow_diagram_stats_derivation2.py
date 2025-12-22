#!/usr/bin/env pythhon3
import pandas as pd
import numpy as np

def add_pass_criteria(df):
    df['Sample_reads_percent_of_refs_AuG_truc10']=df['Sample_reads_percent_of_refs']/df['AuG_trunc10']
    df['OR pass']=np.where((df['AuG_trunc10']>0.003) | (df['Cov1_perc']>0.25) | (df['Sample_reads_percent_of_refs']>0.007),True,False)
    df['AND ratio pass']=np.where(df['Sample_reads_percent_of_refs_AuG_truc10']>0.1,True,False)
    df['2 reads pass']=np.where(df['sample num reads']>=2,True,False)
    return df

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

# add in PCR meta data
df_PCR=pd.read_csv('pcr_organisms.csv')
df=df.merge(df_PCR, on=['Run','barcode'], how='left')
df['Batch']=1
# if 'testcode_checked' == 'FLRA; C2VP' then remove extra biofire organisms
not_biofire_codes=['FLRA; C2VP', 'FLRP; C2VP']
df['test']=np.where(df['testcode_checked'].isin(not_biofire_codes), 'ALINITY', None)
df['test']=np.where(df['testcode_checked'].isin(['RCPR; C2VP','RPCR; C2VP']), 'BIOFIRE', df['test'])
df['test']=np.where(df['testcode_inferred']=='RCPR', 'BIOFIRE', df['test'])

df_full=df.copy()
df=df[~((df['pathogen'].isin(biofire_additional)) & (df['testcode_checked'].isin(not_biofire_codes)))]
df['pathogen_tests']=df.groupby(['Run', 'barcode'])[['pathogen']].transform('count')

#df['test']=np.where(df['pathogen_tests']==24, 'BIOFIRE', 'ALINITY')
df['test_type']=df['test']

# readjust spiked values for derivation set
df['PC_PASSES']=df['orthoreovirus passed']+df['zika passed']+df['murine_respirovirus passed']
df['PCs_passed']=np.where((df['PC_PASSES']>=2 ), 1, 0)


# remove pathogens if test != BIOFIRE
# get potential additional yield 
#df_ay=df[~((df['pathogen'].isin(biorifre_organisms)) & (df['test']!='BIOFIRE'))]
#df_ay=df_ay[~df_ay['pathogen'].isin(alinity_cephid_organisms)]
#df_ay.to_csv('additional_yield.csv', index=False)


# remove biofire tests from alinity cephid
#df_full=df.copy()
#df=df[ ~((df['pathogen'].isin(biofire_additional)) & (df['test'].isin(['ALINITY', 'CEPHEID']))) ]

# add in negative meta
negative_meta=pd.read_csv('negatives_meta.csv')
df=df.merge(negative_meta, on=['Run', 'barcode'], how='left')
df['MS2_spike']=df['spiked']
df['IC_virus_spike']=df['spiked']

df=add_pass_criteria(df)
df.to_csv('biofire_results_merged_adjusted.csv', index=False)

# remove negative controls
#df=df[df['spiked']>=1]
#df=df[(df['MS2_spike']==1) & (df['IC_virus_spike']==1)]


# count number of tests by type
#g=df.groupby(['test_type'])['pathogen_tests'].describe()
#print(g)

# count number of samples by test type
bf=df[df['test_type']=='BIOFIRE']
bf=bf.drop_duplicates(subset=['Run', 'barcode'])
ac=df[df['test_type'].isin(['ALINITY', 'CEPHEID'])]
ac=ac.drop_duplicates(subset=['Run', 'barcode'])
cf=df[df['test_type'].isin(['CEPHEID'])]
cf=cf.drop_duplicates(subset=['Run', 'barcode'])
al=df[df['test_type'].isin(['ALINITY'])]
al=al.drop_duplicates(subset=['Run', 'barcode'])
total_samples=bf.shape[0]+ac.shape[0]
print(f'Total number of samples X: {total_samples}')
print('Number of Biofire samples Yb:', bf.shape[0])
print('Number of Alinity Cepheid samples Yc:', ac.shape[0])
print('Number of Cepheid samples:', cf.shape[0])
print('Number of Alinity samples:', al.shape[0])


unique_runs=df['Run'].nunique()
print('Total number of runs R:', unique_runs )
unique_batches=len(df.drop_duplicates(['Run','Batch']).index)
print('Total number of batches B:', unique_batches )

# count number of samples that failed run/sample controls, higher than 426.6 MB for total run bases and 30,000 for total sample reads
dfF=df[df['total run bases']<400000000]
df=df[df['total run bases']>=400000000]
unique_runs=df['Run'].nunique()
failed_runs=dfF['Run'].nunique()
print(f'Number of runs with total run bases lower than 400 MB xR: {failed_runs}')
print(f'Number of runs with total run bases higher than 400 MB: {unique_runs}')
min_reads=25000
dfSF=df[df['total sample reads']<min_reads]
dfSunique=dfSF.drop_duplicates(subset=['Run', 'barcode'])
dfSunique=dfSunique[dfSunique['test_type'].isin(['BIOFIRE', 'ALINITY', 'CEPHEID'])]
print(f'Number of samples with total sample reads lower than {min_reads:,} xS: {dfSunique.shape[0]}')
#df=df[df['total sample reads']>=min_reads]
dfSunique=df.drop_duplicates(subset=['Run', 'barcode'])
dfSunique=dfSunique[dfSunique['test_type'].isin(['BIOFIRE', 'ALINITY', 'CEPHEID'])]

print(f'Number of samples with total sample reads higher than {min_reads:,} X2: {dfSunique.shape[0]}')
print(f'Percentage of samples with total sample reads higher than {min_reads:,}: {dfSunique.shape[0]/total_samples*100:.2f}%')

# count number of samples that failed negative controls
print(df['IC_virus_spike'].unique())
df_negs=df[df['IC_virus_spike']==0]
df_negs2=df_negs.copy()
df_negs2=df_negs2[df_negs2['seq_name'].isin(['Negative control'])]
print(df_negs2[['Run','barcode']].drop_duplicates())
df_negs=df_negs.copy()
df_negs['neg_pass']=np.where((df_negs['pathogen']!='unmapped')&(df_negs['sample num reads']>=2), True, False)
df_negs_pass=df_negs[(df_negs['neg_pass']==True)|(df_negs['MS2 passed']==1)|(df_negs['zika passed']==1)|(df_negs['orthoreovirus passed']==1)|(df_negs['murine_respirovirus passed']==1)]
pd.options.display.max_columns = None
print(df_negs_pass[['Run','Batch','barcode','pathogen','neg_pass','sample num reads']])
#df_negs_pass=df_negs[(df_negs['pass']==True) | (df_negs['PCs_passed']==1)]
if len(df_negs_pass)>0:
    unique_batches=len(df_negs_pass.drop_duplicates(['Run','Batch'], keep='first').index)
    print(f'Batches that failed negative controls but passed run/sample controls:xB {unique_batches}')
    failedruns=list(df_negs_pass['Run'].unique())
    print('Run names failed negative controls but passed run/sample controls:', failedruns)
    df_failed_negs=df[df['Run'].isin(failedruns)]
    df_failed_negs=df_failed_negs.copy()
    df_failed_negs.drop_duplicates(subset=['Run', 'barcode'], keep='first', inplace=True)
    df_failed_negs=df_failed_negs[df_failed_negs['test_type'].isin(['BIOFIRE', 'ALINITY', 'CEPHEID'])] 
    print(f'Samples that failed negative controls but passed run/sample controls:Z {df_failed_negs.shape[0]}')
else:
    failedruns=[]
    print('Batches that failed negative controls but passed run/sample controls:xB 0')
    print('Samples that failed negative controls but passed run/sample controls:Z 0')


# samples that passed  negative controls passed
if len(failedruns)>0:
    df2=df[~df['Run'].isin(failedruns)]
    df2=df2.copy()
else:
    df2=df.copy()

bf=df2[df2['test_type']=='BIOFIRE']
bf=bf.drop_duplicates(subset=['Run', 'barcode'])
ac=df2[df2['test_type'].isin(['ALINITY', 'CEPHEID'])]
ac=ac.drop_duplicates(subset=['Run', 'barcode'])

total_samples_passing_BNC=bf.shape[0]+ac.shape[0]
print(f'Total number of samples passing batch negative controls: X3 {total_samples_passing_BNC} ({total_samples_passing_BNC/total_samples*100:.0f}%)')
print(f'Percentage of samples passing batch negative controls: {total_samples_passing_BNC/total_samples*100:.2f}%')

# count number of samples that failed batch ampification negative controls
df_amp_control=df[(df['MS2_spike']==0) & (df['IC_virus_spike']==1) ]#& (df['test']==0)]
# orthoreovirus passed	zika passed	MS2 passed	murine_respirovirus passed
df_amp_control_PCFAIL=df_amp_control[(df_amp_control['orthoreovirus passed']==0) & (df_amp_control['zika passed']==0) & (df_amp_control['murine_respirovirus passed']==0)]
#if len(df_amp_control_PCFAIL)>0:
#    # get Run and Batch tuples
#    failed_batches=list(df_amp_control_PCFAIL.drop_duplicates(['Run','Batch'])[['Run','Batch']].itertuples(index=False, name=None))
#    print(failed_batches)
#    unique_batches=len(df_amp_control_PCFAIL.drop_duplicates(['Run','Batch']).index)
#    print(f'Batches that failed amplification controls but passed run/sample controls:xB {unique_batches}')
#    # filter df by failed_batches 
#    df_ACF_samples=df_amp_control_PCFAIL[df_amp_control_PCFAIL[['Run','Batch']].apply(tuple, axis=1).isin(failed_batches)]

    #print('Samples that failed amplification controls but passed run/sample controls:Z {df_amp_control_PCFAIL.shape[0]}')
#else:
#    print('Batches that failed amplification controls but passed run/sample controls:xB 0')
#    print('Samples that failed amplification controls but passed run/sample controls:Z 0')


# samples that passed batch amplification controls
bf=df2[df2['test_type']=='BIOFIRE']
bf=bf.drop_duplicates(subset=['Run', 'barcode'])
ac=df2[df2['test_type'].isin(['ALINITY', 'CEPHEID'])]
ac=ac.drop_duplicates(subset=['Run', 'barcode'])

total_samples_passing_BNC=bf.shape[0]+ac.shape[0]
print(f'Total number of samples passing batch negative controls: X4 {total_samples_passing_BNC} ({total_samples_passing_BNC/total_samples*100:.0f}%)')
#print(f'Percentage of samples passing batch negative controls: {total_samples_passing_BNC/total_samples*100:.2f}%')



# count number of samples that failed positive controls
df3=df2[(df2['PCs_passed']==0) & (df2['test'].isin(['BIOFIRE', 'ALINITY', 'CEPHEID']))]
df3=df3.copy()
df3.to_csv('failed_positive_controls.csv', index=False)
df3.drop_duplicates(subset=['Run', 'barcode'], inplace=True)
print(f'Number of samples failing positive controls xF: {df3.shape[0]} ({df3.shape[0]/total_samples_passing_BNC*100:.0f}%)')
#print(f'Percentage of samples failing positive controls: {df3.shape[0]/total_samples*100:.2f}%')

# number of PCs that passed
failed_samples=df3.shape[0]
MS2fail=df3[df3['MS2 passed']==0]
print(f'Number of samples failing MS2: {MS2fail.shape[0]}({MS2fail.shape[0]/failed_samples*100:.0f}%)')
ICfail_1=df3[df3['PC_PASSES']==2]
ICfail_2=df3[df3['PC_PASSES']==1]
ICfail_3=df3[df3['PC_PASSES']==0]

print(f'Number of samples failing 1 IC virus: {ICfail_1.shape[0]}({ICfail_1.shape[0]/failed_samples*100:.0f}%)')
print(f'Number of samples failing 2 IC viruses: {ICfail_2.shape[0]}({ICfail_2.shape[0]/failed_samples*100:.0f}%)')
print(f'Number of samples failing 3 IC viruses: {ICfail_3.shape[0]}({ICfail_3.shape[0]/failed_samples*100:.0f}%)')

bf=df3[df3['test_type']=='BIOFIRE']
bf=bf.drop_duplicates(subset=['Run', 'barcode'])
ac=df3[df3['test_type'].isin( ['ALINITY', 'CEPHEID'])]
ac=ac.drop_duplicates(subset=['Run', 'barcode'])
total_samples_failing_BNC=bf.shape[0]+ac.shape[0]

cMG_tests=(bf.shape[0]*len(biorifre_organisms))+(ac.shape[0]*len(alinity_cephid_organisms))
df3=df2[(df2['PCs_passed']==0) & (df2['test'].isin(['BIOFIRE', 'ALINITY', 'CEPHEID']))]
df3=df3.copy()
df3g=df3[df3['gold_standard']>=1]
num_pathogens=df3g.shape[0]
df3fp=df3g[df3g['pass']=='True']
num_p_pass=df3fp.shape[0]
print(f'Number of samples failing positive controls with pathogens identified by cMG ix:{num_p_pass}/iZ:{num_pathogens} ({num_p_pass/num_pathogens*100:.0f}%)')
df3fp2=df3g[df3g['pass']==0].shape[0]
negs_not_identified=cMG_tests-df3fp2
print(f'Total number of cMG tests for failed positive control samples px:{negs_not_identified}/ pW:{cMG_tests} ({negs_not_identified/cMG_tests*100:.0f}%)') 

#print(f'Total number of samples failing batch negative controls: {total_samples_failing_BNC}')
#print(f'Percentage of samples failing batch negative controls: {total_samples_failing_BNC/total_samples*100:.2f}%')

# count number of samples that passed positive controls
df4=df2[(df2['PCs_passed']==1) & (df2['test'].isin(['BIOFIRE', 'ALINITY', 'CEPHEID']))]
df4=df4.copy()
df5=df4.drop_duplicates(subset=['Run', 'barcode'])
print(f'Number of samples passing positive controls X5: {df5.shape[0]}({df5.shape[0]/total_samples*100:.0f}% total, {df5.shape[0]/total_samples_passing_BNC*100:.0f}% of samples passing negative controls)')
#df5p=df5[df5['sample_positive']>=1]
#print(f'Number of positive samples passing positive controls : {df5p.shape[0]} ({df5.shape[0]/total_samples*100:.0f}% total, {df5p.shape[0]/df5.shape[0]*100:.0f}% of samples passing negative controls)')
#print(f'Percentage of samples passing positive controls of total samples: {df5.shape[0]/total_samples*100:.2f}%')
#print(f'Percentage of samples passing positive controls of samples passing negative controls: {df5.shape[0]/total_samples_passing_BNC*100:.2f}%')

# print number of samples per test type
bf=df4[df4['test_type']=='BIOFIRE']
bf=bf.drop_duplicates(subset=['Run', 'barcode'])
ac=df4[df4['test_type'].isin(['ALINITY', 'CEPHEID'])]
ac=ac.drop_duplicates(subset=['Run', 'barcode'])
total_samples=bf.shape[0]+ac.shape[0]
print(f'Number of Biofire samples passing controls: \t\t{bf.shape[0]}')
print(f'Number of Alinity Cepheid samples passsing controls: \t{ac.shape[0]}')

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

df7n=df6[df6['polymicrobial']==False]
print(f'Number of samples with 1 pathogen detected in gold_standard nX: {df7n.shape[0]}')
print(f'Number of samples with > pathogen detected in gold_standard nY: {df7.shape[0]/2} ({df7.shape[0]} pathogens)')
df7.to_csv('more_than_1_pathogen_goldstandard.csv', index=False)

# count the number of samples where gold_standard is 1
df8=df4[df4['gold_standard']==1]
df8.to_csv('gold_standard_positives.csv', index=False)
df8flu=df8[df8['pathogen'].isin(['Influenza A','Influenza A/H1-2009'])]
df8flu=df8flu.copy()
df8flu['polyFlu_classification']=df8flu.duplicated(subset=['Run', 'barcode'], keep=False)
df8polyfly=df8flu[df8flu['polyFlu_classification']==True]
df8polyfly=df8polyfly.copy()
df8polyfly.drop_duplicates(subset=['Run', 'barcode'], inplace=True)
print(f'Number of samples with >1 Influenza detected in gold_standard nW: {df8polyfly.shape[0]}')
print(f'Total Z pathogens identified from routine laboratory testing: {df8.shape[0]}')

# count number of samples where there are no pathogens detected in gold_standard
df8n=df4[df4['sample_positive']==0]
df8n=df8n.copy()
df8n.drop_duplicates(subset=['Run', 'barcode'], inplace=True, keep='first')
print(f'Samples negative on all targets assayed nZ: {df8n.shape[0]}')

# count the number of samples where gold_standard is 0
df9=df4[(df4['gold_standard']==0) & (df4['pathogen']!='unmapped')]
print(f'Total W pathogens not identified from routine laboratory testing: {df9.shape[0]}')

bf=df4[df4['test_type']=='BIOFIRE']
bf=bf.drop_duplicates(subset=['Run', 'barcode'])
ac=df4[df4['test_type'].isin(['ALINITY', 'CEPHEID'])]
ac=ac.drop_duplicates(subset=['Run', 'barcode'])
total_passed_control_criteria=bf.shape[0]+ac.shape[0]
number_if_biofire=total_passed_control_criteria*23
actual_number_tested_for=(bf.shape[0]*23)+(ac.shape[0]*4)
X_not_tested_for=number_if_biofire-actual_number_tested_for
print(f'Total XP pathogens not tested with routine laboratory testing XP: {X_not_tested_for}')

# count the number of samples where pass is 1
#print(df4['pass'].unique()) 
df10=df4[(df4['pass']=='True') & (df4['gold_standard']==1)]
df10.to_csv('xS_passing_samples.csv', index=False)
print(f'Total xS samples passing gold standard: {df10.shape[0]} ({df10.shape[0]/df8.shape[0]*100:.0f}%)')
#print(f'Percentage of samples passing gold standard: {df10.shape[0]/df8.shape[0]*100:.2f}%')

# Count the number of samples where pass is 0 and gold_standard is 0
df11=df4[(df4['pass']!='True') & (df4['gold_standard']==0) & (df4['pathogen']!='unmapped')]
fp=df4[(df4['pass']=='True') & (df4['gold_standard']==0) & (df4['pathogen']!='unmapped')]
fp.to_csv('false_positives.csv', index=False)
print(f'Total x samples not identified by cMG xf: {df11.shape[0]} ({df11.shape[0]/df9.shape[0]*100}%)')
#print(f'Percentage of samples not passing gold standard: {df11.shape[0]/df9.shape[0]*100}%')

# count number of additional yield
df_ay=df_full[((df_full['pathogen'].isin(biofire_additional)) & (df_full['test'].isin(['ALINITY', 'CEPHEID']))) ]
df_ay=df_ay[(df_ay['pass']=="True") & (df_ay['PCs_passed']==1) & (~df_ay['Run'].isin(failedruns))]
df_ay.to_csv('additional_yield_passed.csv', index=False)
print(f'Number of additional yield samples xN/XP:{df_ay.shape[0]}/{X_not_tested_for} ({df_ay.shape[0]/X_not_tested_for*100:.2f}%)')
