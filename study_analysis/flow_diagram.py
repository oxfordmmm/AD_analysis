#!/usr/bin/env pythhon3
import pandas as pd
import numpy as np
import argparse
import sys
import os

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

DNA_orgs=['Adenovirus',
'Bordetella_parapertussis',
'Bordetella_pertussis',
'Chlamydia_pneumoniae',
'Mycoplasma_pneumoniae']

HRE_orgs=['Rhinovirus/enterovirus']
          
alinity_cephid_organisms=['Influenza A', 'Influenza B', 'RSV', 'SARS-CoV-2']

biofire_additional=[org for org in biorifre_organisms if org not in alinity_cephid_organisms]

def readjust_pass(df, AND_ratio, AND_ratio_metric='Sample_reads_percent_of_refs_AuG_truc10'):
    """Change the pass colum in a each row based on the AND_ratio.
    
    Parameters
    ----------
    df : dataframe
        DataFrame containing the data to be processed.
    AND_ratio : float
        The ratio used to determine if a sample passes based on the number of reads and the total reads.
        
    Returns
    -------
    df : dataframe
        DataFrame with the updated 'pass' column.
    """
    df['AND ratio pass']=np.where(df[AND_ratio_metric]>AND_ratio, True, False)
    
    df['pass']=np.where((df['AuG_trunc10']>0.003) | (df['Cov1_perc']>0.25) | (df['Sample_reads_percent_of_refs']>0.007),'True','False')    
    df['pass']=np.where(df['AND ratio pass']==True, df['pass'],'False')
    df['pass']=np.where(df['sample num reads']>=2, df['pass'],'False')
    return df

def readjust_flu_pass(df):
    """Change pass to 'False' if pathogen is flu A genotype and TOP_FLU_A == FALSE"""
    flu_A_options=['Influenza A','Influenza A/H1-2009','Influenza A/H3','Influenza A/H1']

    df['pass']=np.where((df['pathogen'].isin(['Influenza A/H1-2009','Influenza A/H1','Influenza A/H3']) & (df['TOP_FLU_A']==False)), 'False', df['pass'])

    # change flu gold standard 
    df['gold_standard']=np.where( ( df['pathogen'].isin(flu_A_options) ) & ( df['FLU_A_POS']==1) & (df['TOP_FLU_A']==1), 1, df['gold_standard'] )
    df['gold_standard']=np.where( ( df['pathogen'].isin(flu_A_options) ) & ( df['FLU_A_POS']==1) & (df['TOP_FLU_A']==0), 0, df['gold_standard'] )
    df['gold_standard']=np.where( ( df['pathogen'].isin(['Influenza A']) ) & ( df['FLU_A_POS']==1) , 1, df['gold_standard'] )
    return df

def repopulate_columns(df):
    # repopulate df PC_passes with results from the same Run and barcode if PC_passes is 1
    df['PCs_passed']=df.groupby(['Run', 'barcode'])['PCs_passed'].transform('max')
    df['extraction_date']=pd.to_datetime(df['extraction_date'], format='%d.%m.%Y',errors='coerce')
    df['collection_date']=pd.to_datetime(df['collection_date'], format='%d.%m.%Y', errors='coerce') 
    for col in ['orthoreovirus', 'zika',	'MS2',	'murine_respirovirus',	
                'orthoreovirus passed',	'zika passed',	'MS2 passed',	'murine_respirovirus passed',
                'total sample bases', 'total sample reads', 'PCs_passed', 'PC_PASSES', 'FLU_A_POS', 
                'extraction_date','collection_date', 'amplification_control' ]:
        df[col]=df.groupby(['Run', 'barcode'])[[col]].transform('max')

    for col in ['total run bases',	'total run reads']:
        df[col]=df.groupby(['Run'])[[col]].transform('max')

    df['sample_positive']=df.groupby(['Run', 'barcode'])['gold_standard'].transform('max')
    return df

def readjust_spiked_values(df):
    # readjust spiked values for derivation set
    df['PC_PASSES']=df['orthoreovirus passed']+df['zika passed']+df['murine_respirovirus passed']
    df['PCs_passed']=np.where((df['PC_PASSES']>=2) & (df['MS2 passed']==1), 1, 0)
    return df

def get_additional_yield(df):
    # remove pathogens if test != BIOFIRE
    # get potential additional yield 
    df_ay=df[~((df['pathogen'].isin(biorifre_organisms)) & (~df['test'].isin(['BIOFIRE', 'BIOFIRE '])))]
    df_ay=df_ay[~df_ay['pathogen'].isin(alinity_cephid_organisms)]
    #df_ay.to_csv('additional_yield.csv', index=False)
    return df_ay

def remove_biofire_additional(df):
    # remove biofire tests from alinity cephid
    df_full=df.copy()
    df=df[ ~((df['pathogen'].isin(biofire_additional)) & (df['test'].isin(['ALINITY','ALINTY', 'CEPHEID']))) ]   
    df=df.copy()

    df['pathogen_tests']=df.groupby(['Run', 'barcode'])[['pathogen']].transform('count')

    df['test_type']=df['test']

    # add in negative meta
    #negative_meta=pd.read_csv('negatives_meta.csv')
    #df=df.merge(negative_meta, on=['Run', 'barcode'], how='left')

    #df.to_csv('biofire_results_merged_adjusted.csv', index=False)
    return df, df_full

def count_samples(df, metrics):
    # count number of samples by test type
    bf=df[df['test_type'].isin(['BIOFIRE','BIOFIRE '])]
    bf=bf.drop_duplicates(subset=['Run', 'barcode'])
    ac=df[df['test_type'].isin(['ALINITY','ALINTY', 'CEPHEID'])]
    ac=ac.drop_duplicates(subset=['Run', 'barcode'])
    cf=df[df['test_type'].isin(['CEPHEID'])]
    cf=cf.drop_duplicates(subset=['Run', 'barcode'])
    al=df[df['test_type'].isin(['ALINITY','ALINTY'])]
    al=al.drop_duplicates(subset=['Run', 'barcode'])
    total_samples=bf.shape[0]+ac.shape[0]
    print(f'Total number of samples X: {total_samples}')
    print('Number of Biofire samples Yb:', bf.shape[0])
    print('Number of Alinity Cepheid samples Yc:', ac.shape[0])
    print('Number of Cepheid samples Ycc:', cf.shape[0])
    print('Number of Alinity samples Yca:', al.shape[0])

    # count number of pathogens (gold_standard) identified by unique Run and barcode
    df['num_gold_standard'] = df.groupby(['Run', 'barcode'])['gold_standard'].transform('sum')
    df_unique_tested=df[df['test_type'].isin(['BIOFIRE','BIOFIRE ', 'ALINITY','ALINTY', 'CEPHEID'])]
    df_unique_tested=df_unique_tested.drop_duplicates(subset=['Run', 'barcode'])
    print(f'Number of samples with no pathogen identified by PCR: {len(df_unique_tested[df_unique_tested["num_gold_standard"]==0])}')
    print(f'Number of samples with one pathogen identified by PCR: {len(df_unique_tested[df_unique_tested["num_gold_standard"]==1])}')
    print(f'Number of samples with more than one pathogen identified by PCR: {len(df_unique_tested[df_unique_tested["num_gold_standard"]>1])}')

    unique_runs=df['Run'].nunique()
    print('Total number of runs R:', unique_runs )
    unique_batches=len(df.drop_duplicates(['Run','Batch']).index)
    print('Total number of batches B:', unique_batches )
    metrics2={'Total number of samples X': total_samples,
            'Number of Biofire samples Yb': bf.shape[0],
            'Number of Alinity Cepheid samples Yc': ac.shape[0],
            'Total number of runs R': unique_runs,
            'Total number of batches B': unique_batches}
    metrics.update(metrics2)
    return metrics

def count_failed_runs_samples(df, metrics):
    # count number of samples that failed run/sample controls, higher than 426.6 MB for total run bases and 30,000 for total sample reads
    total_samples=metrics['Total number of samples X']
    dfF=df[df['total run bases']<400000000]
    dfH=df[df['total run bases']>=400000000]
    unique_runs=dfH['Run'].nunique()
    failed_runs=dfF['Run'].nunique()
    print(f'Number of runs with total run bases lower than 400 MB xR: {failed_runs}')
    print(f'Number of runs with total run bases higher than 400 MB: {unique_runs}')
    min_reads=25000
    dfSF=df[df['total sample reads']<min_reads]
    dfSunique=dfSF.drop_duplicates(subset=['Run', 'barcode'])
    dfSunique=dfSunique[dfSunique['test_type'].isin(['BIOFIRE','BIOFIRE ', 'ALINITY','ALINTY', 'CEPHEID'])]
    print(f'Number of samples with total sample reads lower than {min_reads:,} xS: {dfSunique.shape[0]}')
    df2=df[(df['total sample reads']>=min_reads) | (df['amplification_control']==1) | (df['reverse_transcription_control']==1)] # to keep batch controls
    dfSunique=df2.drop_duplicates(subset=['Run', 'barcode'])
    dfSunique=dfSunique[dfSunique['test_type'].isin(['BIOFIRE','BIOFIRE ', 'ALINITY','ALINTY', 'CEPHEID'])]

    print(f'Number of samples with total sample reads higher than {min_reads:,} X2: {dfSunique.shape[0]}')
    print(f'Percentage of samples with total sample reads higher than {min_reads:,}: {dfSunique.shape[0]/total_samples*100:.2f}%')
    
    metrics['Number of runs with total run bases lower than 400 MB xR'] = failed_runs 
    metrics['Number of runs with total run bases higher than 400 MB'] = unique_runs
    metrics['Number of samples with total sample reads lower than 25,000 xS'] = dfSunique.shape[0]
    metrics['Number of samples with total sample reads higher than 25,000 X2'] = dfSunique.shape[0]
    metrics['Percentage of samples with total sample reads higher than 25,000'] = dfSunique.shape[0]/total_samples*100

    df['run_pass']=np.where((df['test_type'].isin(['BIOFIRE','BIOFIRE ', 'ALINITY','ALINTY', 'CEPHEID'])) & (df['total run bases']>=400000000), 'True', 'False')
    df['barcode_pass']=np.where((df['test_type'].isin(['BIOFIRE','BIOFIRE ', 'ALINITY','ALINTY', 'CEPHEID'])) & (df['total sample reads']>=min_reads), 'True', 'False')
    return df, metrics

def count_failed_negative_controls(df, metrics):
    # count number of samples that failed negative controls
    # spike negs
    #df_negs=df[(df['MS2_spike']==0) & (df['IC_virus_spike']==0)]
    # This is also called the batch negatice control in the flow diagram
    df_negs=df[df['amplification_control']==1]
    df_negs_pass=df_negs[(df_negs['pass']=='True') | (df_negs['PCs_passed']==1)]
    # clinical negs
    df_negs_clinical=df[(df['Negative control']==1) & (df['pass']=='True')]
    #df_negs_pass=pd.concat([df_negs_pass, df_negs_clinical], ignore_index=True)
    if len(df_negs_pass)>0:
        unique_batches=len(df_negs_pass.drop_duplicates(['Run','Batch'], keep='first').index)
        print(f'Batches that failed negative controls but passed run/sample controls:xB {unique_batches}')
        metrics['Batches that failed negative controls but passed run/sample controls:xB'] = unique_batches
        failedruns=list(df_negs_pass['Run'].unique())
        df_failed_negs=df[df['Run'].isin(failedruns)]
        df_failed_negs=df_failed_negs.copy()
        df_failed_negs.drop_duplicates(subset=['Run', 'barcode'], keep='first', inplace=True)
        df_failed_negs=df_failed_negs[df_failed_negs['test_type'].isin(['BIOFIRE', 'BIOFIRE ', 'ALINITY','ALINTY', 'CEPHEID']) & ((df_failed_negs['run_pass']=='True') | (df_failed_negs['barcode_pass']=='True'))] 
        print(f'Samples that failed negative controls but passed run/sample controls:Z {df_failed_negs.shape[0]}')
        metrics['Samples that failed negative controls but passed run/sample controls:Z'] = df_failed_negs.shape[0]
        
    else:
        failedruns=[]
        print('Batches that failed negative controls but passed run/sample controls:xB 0')
        print('Samples that failed negative controls but passed run/sample controls:Z 0')
        metrics['Batches that failed negative controls but passed run/sample controls:xB'] = 0
        metrics['Samples that failed negative controls but passed run/sample controls:Z'] = 0
    
    # samples that passed  negative controls passed
    if len(failedruns)>0:
        df2=df[~df['Run'].isin(failedruns)]
        df2=df2.copy()
    else:
        df2=df.copy()

    return df, df2, metrics

def count_passing_samples(df2, metrics):
    bf=df2[(df2['run_pass']=='True') & (df2['barcode_pass']=='True') & (df2['test_type'].isin(['BIOFIRE','BIOFIRE ']))]
    bf=bf.drop_duplicates(subset=['Run', 'barcode'])
    ac=df2[(df2['run_pass']=='True') & (df2['barcode_pass']=='True') & (df2['test_type'].isin(['ALINITY','ALINTY', 'CEPHEID']))]
    ac=ac.drop_duplicates(subset=['Run', 'barcode'])

    total_samples_passing_BNC=bf.shape[0]+ac.shape[0]
    total_samples=metrics['Total number of samples X']
    print(f'Total number of samples passing batch negative controls: X3 {total_samples_passing_BNC} ({total_samples_passing_BNC/total_samples*100:.0f}%)')
    print(f'Percentage of samples passing batch negative controls: {total_samples_passing_BNC/total_samples*100:.2f}%')

    metrics['Total number of samples passing batch negative controls: X3'] = total_samples_passing_BNC
    metrics['Percentage of samples passing batch negative controls'] = total_samples_passing_BNC/total_samples*100

    return metrics

def count_samples_failing_batch_controls(df, df2, metrics):
    # count number of samples that failed batch ampification negative controls
    df_RC_control=df[(df['reverse_transcription_control']==1) & (df['IC_virus_spike']==1)]
    # orthoreovirus passed	zika passed	MS2 passed	murine_respirovirus passed
    df_RC_control_PCFAIL=df_RC_control[(df_RC_control['orthoreovirus passed']==0) & (df_RC_control['zika passed']==0) & (df_RC_control['murine_respirovirus passed']==0) \
                                    | (df_RC_control['MS2 passed']==1) \
                                    | (df_RC_control['pass']=='True') ]
    if len(df_RC_control_PCFAIL)>0:
        # remove unmapped pathogens
        df2=df2[df2['pathogen']!='unmapped']
    #    # get Run and Batch tuples
        failed_batches=list(df_RC_control_PCFAIL.drop_duplicates(['Run','Batch'])[['Run','Batch']].itertuples(index=False, name=None))
        print(failed_batches)
        failed_unique_batches=len(df_RC_control_PCFAIL.drop_duplicates(['Run','Batch']).index)
        unique_batches=metrics['Total number of batches B']
        print(f'Batches that failed amplification controls but passed run/sample controls:xR/B {failed_unique_batches}/{unique_batches}')
        metrics['Batches that failed amplification controls but passed run/sample controls:xR/B'] = f'{failed_unique_batches} {failed_unique_batches}/{unique_batches}'
    #    # filter df by failed_batches 
        df_ACF_samples=df2[df2[['Run','Batch']].apply(tuple, axis=1).isin(failed_batches)]
        df_ACF_samples=df_ACF_samples.copy()
        df_ACF_samples=df_ACF_samples[(df_ACF_samples['test_type'].isin(['BIOFIRE', 'BIOFIRE ', 'ALINITY','ALINTY', 'CEPHEID'])) & ((df_ACF_samples['run_pass']=='True') & (df_ACF_samples['barcode_pass']=='True'))]
        df3r=df_ACF_samples.copy()
        df_ACF_samples.drop_duplicates(subset=['Run', 'barcode'], inplace=True)
        print(f'Samples in failed batch amplification controls:rZ {df_ACF_samples.shape[0]}')
        metrics['Samples in failed batch amplification controls:rZ'] = df_ACF_samples.shape[0]

        # count any samples that passed with a target pathogen
        df_pathogen_pass=df_RC_control_PCFAIL[(df_RC_control_PCFAIL['chrom']!='unmapped') & (df_RC_control_PCFAIL['pass']=='True')]
        if len(df_pathogen_pass)>0:
            print(f'Number of samples in batch amplification controls that passed with a target pathogen: {df_pathogen_pass.shape[0]} ({df_pathogen_pass.shape[0]/df_RC_control_PCFAIL.shape[0]*100:.0f}%)')
            metrics['Number of samples in batch amplification controls that passed with a target pathogen'] = df_pathogen_pass.shape[0]
        else:
            print(f'Number of samples in batch amplification controls that passed with a target pathogen: 0')
            metrics['Number of samples in batch amplification controls that passed with a target pathogen'] = 0

        df_RC_control_PCFAIL=df_RC_control_PCFAIL.drop_duplicates(subset=['Run', 'barcode'], keep='first')
        ms2_passed=df_RC_control_PCFAIL[df_RC_control_PCFAIL['MS2 passed']==1]
        print(f'Number of samples in failed batch amplification controls that passed MS2: {ms2_passed.shape[0]} ({ms2_passed.shape[0]/df_ACF_samples.shape[0]*100:.0f}%)')
        metrics['Number of samples in failed batch amplification controls that passed MS2'] = ms2_passed.shape[0]

        ICfail_1=df_RC_control_PCFAIL[df_RC_control_PCFAIL['PC_PASSES']==2]
        ICfail_2=df_RC_control_PCFAIL[df_RC_control_PCFAIL['PC_PASSES']==1]
        ICfail_3=df_RC_control_PCFAIL[df_RC_control_PCFAIL['PC_PASSES']==0]
        print(f'Number of samples in failed batch amplification controls that failed 1 IC virus: {ICfail_1.shape[0]} ({ICfail_1.shape[0]/df_RC_control_PCFAIL.shape[0]*100:.0f}%)')   
        metrics['Number of samples in failed batch amplification controls that failed 1 IC virus'] = ICfail_1.shape[0]
        print(f'Number of samples in failed batch amplification controls that failed 2 IC viruses: {ICfail_2.shape[0]} ({ICfail_2.shape[0]/df_RC_control_PCFAIL.shape[0]*100:.0f}%)')
        metrics['Number of samples in failed batch amplification controls that failed 2 IC viruses'] = ICfail_2.shape[0]
        print(f'Number of samples in failed batch amplification controls that failed 3 IC viruses: {ICfail_3.shape[0]} ({ICfail_3.shape[0]/df_RC_control_PCFAIL.shape[0]*100:.0f}%)')
        metrics['Number of samples in failed batch amplification controls that failed 3 IC viruses'] = ICfail_3.shape[0]

        # count pathogen tests that would have passed
        bf=df3r[df3r['test_type'].isin(['BIOFIRE','BIOFIRE '])]
        bf=bf.drop_duplicates(subset=['Run', 'barcode'])
        ac=df3r[df3r['test_type'].isin( ['ALINITY','ALINTY', 'CEPHEID'])]
        ac=ac.drop_duplicates(subset=['Run', 'barcode'])
        total_samples_failing_BNC=bf.shape[0]+ac.shape[0]

        cMG_tests=(bf.shape[0]*len(biorifre_organisms))+(ac.shape[0]*len(alinity_cephid_organisms))
        df3=df2[(df2['PCs_passed']==0) & (df2['test'].isin(['BIOFIRE','BIOFIRE ', 'ALINITY','ALINTY', 'CEPHEID']))]
        df3r=df3r.copy()
        df3g=df3r[df3r['gold_standard']>=1]
        num_pathogens=df3g.shape[0]
        df3fp=df3g[df3g['pass']=='True']
        num_p_pass=df3fp.shape[0]
        print(f'Number of samples failing positive batch with pathogens identified by cMG rx:{num_p_pass}/rZ:{num_pathogens} ({num_p_pass/num_pathogens*100:.0f}%)')
        metrics['Number of samples failing positive batch with pathogens identified by cMG'] = f'{num_p_pass}/ {num_pathogens} ({num_p_pass/num_pathogens*100:.0f}%)'

        df3n=df3r[df3r['gold_standard']==0]
        df3fp2=df3n[df3n['pass']!='True']
        #print(df3fp2['pass'].unique())
        negs_not_identified=df3fp2.shape[0]

        # count pathogen tests again that would have passed
        bf=df3n[df3n['test_type'].isin(['BIOFIRE','BIOFIRE '])]
        bf=bf.drop_duplicates(subset=['Run', 'barcode'])
        ac=df3n[df3n['test_type'].isin( ['ALINITY','ALINTY', 'CEPHEID'])]
        ac=ac.drop_duplicates(subset=['Run', 'barcode'])

        cMG_tests=df3n.shape[0]
        #cMG_tests=(bf.shape[0]*len(biorifre_organisms))+(ac.shape[0]*len(alinity_cephid_organisms))

        print(f'Total number of cMG tests for failed positive control samples rpx:{negs_not_identified}/ rW:{cMG_tests} ({negs_not_identified/cMG_tests*100:.0f}%)')
        metrics['Total number of cMG tests for failed positive control samples rpx/rW'] = f'{negs_not_identified}/ {cMG_tests} ({negs_not_identified/cMG_tests*100:.0f}%)'

        # remove failed batches from the rest of the analysis
        df2=df2[~df2[['Run','Batch']].apply(tuple, axis=1).isin(failed_batches)]
    else:
        print('Batches that failed amplification controls but passed run/sample controls:xR 0')
        metrics['Batches that failed amplification controls but passed run/sample controls:xR'] = 0
        print('Samples that failed amplification controls but passed run/sample controls:Z 0')
        metrics['Samples that failed amplification controls but passed run/sample controls:Z'] = 0

    return df2, metrics

def count_samples_passing_batch_controls(df2, metrics):
    # samples that passed batch amplification controls
    bf=df2[(df2['test_type'].isin(['BIOFIRE', 'BIOFIRE '])) & ( (df2['run_pass']=='True') & (df2['barcode_pass']=='True') )]
    bf=bf.drop_duplicates(subset=['Run', 'barcode'])
    ac=df2[(df2['test_type'].isin(['ALINITY','ALINTY', 'CEPHEID'])) & ( (df2['run_pass']=='True') & (df2['barcode_pass']=='True') )]
    ac=ac.drop_duplicates(subset=['Run', 'barcode'])

    total_samples_passing_BNC=bf.shape[0]+ac.shape[0]
    total_samples=metrics['Total number of samples X']
    print(f'Total number of samples passing batch negative controls: X4 {total_samples_passing_BNC} ({total_samples_passing_BNC/total_samples*100:.0f}%)')
    metrics['Total number of samples passing batch negative controls: X4'] = total_samples_passing_BNC
    #print(f'Percentage of samples passing batch negative controls: {total_samples_passing_BNC/total_samples*100:.2f}%')
    return metrics

def count_samples_failing_positive_controls(df2, metrics):
    # count number of samples that failed positive controls
    df3=df2[(df2['PCs_passed']==0) & (df2['test'].isin(['BIOFIRE', 'BIOFIRE ', 'ALINITY','ALINTY', 'CEPHEID']))]
    df3=df3.copy()
    df3=df3[df3['pathogen']!='unmapped']
    #df3.to_csv('failed_positive_controls.csv', index=False)
    df3.drop_duplicates(subset=['Run', 'barcode'], inplace=True)
    total_samples=metrics['Total number of samples X']
    print(f'Number of samples failing positive controls xF: {df3.shape[0]} ({df3.shape[0]/total_samples*100:.0f}%)')
    metrics['Number of samples failing positive controls xF'] = df3.shape[0]
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
    metrics['Number of samples failing MS2'] = MS2fail.shape[0]
    metrics['Number of samples failing 1 IC virus'] = ICfail_1.shape[0]
    metrics['Number of samples failing 2 IC viruses'] = ICfail_2.shape[0]
    metrics['Number of samples failing 3 IC viruses'] = ICfail_3.shape[0]


    #bf=df3[df3['test_type']=='BIOFIRE']
    #bf=bf.drop_duplicates(subset=['Run', 'barcode'])
    #ac=df3[df3['test_type'].isin( ['ALINITY','ALINTY', 'CEPHEID'])]
    #ac=ac.drop_duplicates(subset=['Run', 'barcode'])
    #total_samples_failing_BNC=bf.shape[0]+ac.shape[0]
    #cMG_tests=(bf.shape[0]*len(biorifre_organisms))+(ac.shape[0]*len(alinity_cephid_organisms))
    
    df3=df2[(df2['PCs_passed']==0) & (df2['test'].isin(['BIOFIRE', 'BIOFIRE ', 'ALINITY','ALINTY', 'CEPHEID']))]
    df3=df3.copy()
    df3g=df3[df3['gold_standard']>=1]
    num_pathogens=df3g.shape[0]
    df3fp=df3g[df3g['pass']=='True']
    num_p_pass=df3fp.shape[0]
    print(f'Number of samples failing positive controls with pathogens identified by cMG ix:{num_p_pass}/iZ:{num_pathogens} ({num_p_pass/num_pathogens*100:.0f}%)')
    metrics['Number of samples failing positive controls with pathogens identified by cMG'] = f'{num_p_pass}/ {num_pathogens} ({num_p_pass/num_pathogens*100:.0f}%)'
    
    df3n=df3[df3['gold_standard']==0]
    #df3n.to_csv('df3n.csv',index=False)
    df3fp2=df3n[df3n['pass']!='True']
    #print(df3fp2['pass'].unique())
    #df3fp2.to_csv('df3fp2.csv', index=False)
    negs_not_identified=df3fp2.shape[0]

    # calculate number of cMG tests for failed positive control samples
    #bf=df3n[df3n['test_type']=='BIOFIRE']
    #bf=bf.drop_duplicates(subset=['Run', 'barcode'])
    #ac=df3n[df3n['test_type'].isin( ['ALINITY','ALINTY', 'CEPHEID'])]
    #ac=ac.drop_duplicates(subset=['Run', 'barcode'])
    #print('alinity,biofire samples: ', ac.shape[0], bf.shape[0])
    
    cMG_tests=df3n.shape[0] 
    #cMG_tests=(bf.shape[0]*len(biorifre_organisms))+(ac.shape[0]*len(alinity_cephid_organisms))

    print(f'Total number of cMG tests for failed positive control samples px:{negs_not_identified}/ pW:{cMG_tests} ({negs_not_identified/cMG_tests*100:.0f}%)') 
    metrics['Total number of cMG tests for failed positive control samples'] = f'{negs_not_identified}/ {cMG_tests} ({negs_not_identified/cMG_tests*100:.0f}%)'

    #print(f'Total number of samples failing batch negative controls: {total_samples_failing_BNC}')
    #print(f'Percentage of samples failing batch negative controls: {total_samples_failing_BNC/total_samples*100:.2f}%')
    return metrics

def count_samples_passing_positive_controls(df2, metrics):
    # count number of samples that passed positive controls
    df4=df2[(df2['PCs_passed']==1) & (df2['test'].isin(['BIOFIRE','BIOFIRE ', 'ALINITY','ALINTY', 'CEPHEID'])) & (df2['run_pass']=='True') & (df2['barcode_pass']=='True')]
    df4=df4.copy()
    df5=df4.drop_duplicates(subset=['Run', 'barcode'])
    total_samples_passing_BNC= metrics['Total number of samples passing batch negative controls: X4']
    total_samples=metrics['Total number of samples X']
    print(f'Number of samples passing positive controls X5: {df5.shape[0]}({df5.shape[0]/total_samples*100:.0f}% total, {df5.shape[0]/total_samples_passing_BNC*100:.0f}% of samples passing negative controls)')
    metrics['Number of samples passing positive controls X5'] = df5.shape[0]
    #df5p=df5[df5['sample_positive']>=1]
    #print(f'Number of positive samples passing positive controls : {df5p.shape[0]} ({df5.shape[0]/total_samples*100:.0f}% total, {df5p.shape[0]/df5.shape[0]*100:.0f}% of samples passing negative controls)')
    #print(f'Percentage of samples passing positive controls of total samples: {df5.shape[0]/total_samples*100:.2f}%')
    #print(f'Percentage of samples passing positive controls of samples passing negative controls: {df5.shape[0]/total_samples_passing_BNC*100:.2f}%')

    # print number of samples per test type
    bf=df4[df4['test_type'].isin(['BIOFIRE','BIOFIRE '])]
    bf=bf.drop_duplicates(subset=['Run', 'barcode'])
    ac=df4[df4['test_type'].isin(['ALINITY','ALINTY', 'CEPHEID'])]
    ac=ac.drop_duplicates(subset=['Run', 'barcode'])
    alinSamples=df4[df4['test_type'].isin(['ALINITY','ALINTY'])]
    cephSamples=df4[df4['test_type'].isin(['CEPHEID'])]
    alinSamples=alinSamples.drop_duplicates(subset=['Run', 'barcode'])
    cephSamples=cephSamples.drop_duplicates(subset=['Run', 'barcode'])
    total_samples=bf.shape[0]+ac.shape[0]
    print(f'Number of Biofire samples passing controls: \t\t{bf.shape[0]}')
    metrics['Number of Biofire samples passing controls'] = bf.shape[0]
    print(f'Number of Alinity Cepheid samples passsing controls: \t{ac.shape[0]}')
    metrics['Number of Alinity Cepheid samples passsing controls'] = ac.shape[0]
    print(f'Number of Alinity samples passing controls: \t\t{alinSamples.shape[0]}')
    metrics['Number of Alinity samples passing controls'] = alinSamples.shape[0]
    print(f'Number of Cepheid samples passing controls: \t\t{cephSamples.shape[0]}')
    metrics['Number of Cepheid samples passing controls'] = cephSamples.shape[0]

    # count the number of samples pathogens detected in gold_standard
    df5g=df4[df4['gold_standard']>=1]
    df5g=df5g.copy()
    df5=df5g[df5g['sample num reads']>=1]
    #df5.drop_duplicates(subset=['Run', 'barcode'], inplace=True)
    print(f'x pathogens identified by cMG: {df5.shape[0]}')
    metrics['x pathogens identified by cMG'] = df5.shape[0]
    print(f'Percentage of samples with pathogens identified by cMG: {df5.shape[0]/df5g.shape[0]*100}%')
    metrics['Percentage of samples with pathogens identified by cMG'] = df5.shape[0]/df5g.shape[0]*100

    df5n=df4[df4['gold_standard']==0]
    df5n0=df5n[df5n['sample num reads']==0]
    print(f'Number of samples with no pathogens identified by cMG: {df5n0.shape[0]}')
    metrics['Number of samples with no pathogens identified by cMG'] = df5n0.shape[0]
    print(f'Percentage of samples with no pathogens identified by cMG: {df5n0.shape[0]/df5n.shape[0]*100}%')
    metrics['Percentage of samples with no pathogens identified by cMG'] = df5n0.shape[0]/df5n.shape[0]*100
    return df4, metrics

def count_multipathogen_samples(df4, metrics):
    # count number of samples with > pathogen detected in gold_standard
    df6=df4[df4['gold_standard']>=1]
    df6=df6.copy()
    df6['polymicrobial']=df6.duplicated(subset=['Run', 'barcode'], keep=False)
    df7=df6[df6['polymicrobial']==True]
    df7=df7.copy()
    #df7=df7[~df7['pathogen'].isin(['Influenza A','Influenza A/H1-2009','Influenza A/H3', 'Influenza A/H1'])]
    df7_uniq=df7.drop_duplicates(subset=['Run', 'barcode'])

    df7n=df6[df6['polymicrobial']==False]
    print(f'Number of samples with 1 pathogen detected in gold_standard nX: {df7n.shape[0]}')
    metrics['Number of samples with 1 pathogen detected in gold_standard nX'] = df7n.shape[0]
    print(f'Number of samples with > pathogen detected in gold_standard nY: {df7_uniq.shape[0]} ({df7.shape[0]} pathogens)')
    metrics['Number of samples with > pathogen detected in gold_standard nY'] =df7_uniq.shape[0]
    #df7.to_csv('more_than_1_pathogen_goldstandard.csv', index=False)

    return metrics

def count_multi_gold_standard_samples(df4, metrics):
    # count the number of samples where gold_standard is 1
    df8=df4[df4['gold_standard']==1]
    #df8['polymicrobial']=df8.duplicated(subset=['Run', 'barcode'], keep=False)
    #df8.to_csv('tmp/gold_standard_positives.csv', index=False)
    df8flu=df8[df8['pathogen'].isin(['Influenza A','Influenza A/H1-2009', 'Influenza A/H3', 'Influenza A/H1'])]
    df8flu=df8flu.copy()
    df8flu['polyFlu_classification']=df8flu.duplicated(subset=['Run', 'barcode'], keep=False)
    df8polyfly=df8flu[df8flu['polyFlu_classification']==True]
    df8polyfly=df8polyfly.copy()
    df8polyfly.drop_duplicates(subset=['Run', 'barcode'], inplace=True)
    print(f'Number of samples with >1 Influenza detected in gold_standard nW: {df8polyfly.shape[0]}')
    metrics['Number of samples with >1 Influenza detected in gold_standard nW'] = df8polyfly.shape[0]
    print(f'Total Z pathogens identified from routine laboratory testing: {df8.shape[0]}')
    metrics['Total Z pathogens identified from routine laboratory testing'] = df8.shape[0]

    return metrics

def count_full_passing_samples(df4,df_full, metrics, AND_ratio):
    # count number of samples where there are no pathogens detected in gold_standard
    df8=df4[df4['gold_standard']==1]
    df8n=df4[df4['sample_positive']==0]
    df8n=df8n.copy()
    df8n=df8n[df8n['test_type'].isin(['BIOFIRE','BIOFIRE ', 'ALINITY','ALINTY', 'CEPHEID'])]
    df8n.drop_duplicates(subset=['Run', 'barcode'], inplace=True, keep='first')
    print(f'Samples negative on all targets assayed nZ: {df8n.shape[0]}')
    metrics['Samples negative on all targets assayed nZ'] = df8n.shape[0]

    # count the number of samples where gold_standard is 0
    df9=df4[(df4['gold_standard']==0) & (df4['pathogen']!='unmapped') & (df4['test_type'].isin(['BIOFIRE','BIOFIRE ', 'ALINITY','ALINTY', 'CEPHEID'])) & ((df4['run_pass']=='True') & (df4['barcode_pass']=='True'))]
    print(f'Total W pathogens not identified from routine laboratory testing: {df9.shape[0]}')
    metrics['Total W pathogens not identified from routine laboratory testing'] = df9.shape[0]

    bf=df4[df4['test_type'].isin(['BIOFIRE','BIOFIRE '])]
    bf=bf.drop_duplicates(subset=['Run', 'barcode'])
    ac=df4[df4['test_type'].isin(['ALINITY','ALINTY', 'CEPHEID'])]
    ac=ac.drop_duplicates(subset=['Run', 'barcode'])
    total_passed_control_criteria=bf.shape[0]+ac.shape[0]
    number_if_biofire=total_passed_control_criteria*23
    actual_number_tested_for=(bf.shape[0]*23)+(ac.shape[0]*4)
    X_not_tested_for=number_if_biofire-actual_number_tested_for
    print(f'Total XP pathogens not tested with routine laboratory testing XP: {X_not_tested_for}')
    metrics['Total XP pathogens not tested with routine laboratory testing XP'] = X_not_tested_for

    # count the number of samples where pass is 1
    #print(df4['pass'].unique()) 
    df10=df4[(df4['pass']=='True') & (df4['gold_standard']>=1)]
    print(f'Total xS samples passing gold standard: {df10.shape[0]} ({df10.shape[0]/df8.shape[0]*100:.0f}%)')
    metrics['Total xS samples passing gold standard'] = f'{df10.shape[0]} ({df10.shape[0]/df8.shape[0]*100:.0f}%)'
    #print(f'Percentage of samples passing gold standard: {df10.shape[0]/df8.shape[0]*100:.2f}%')

    # count the number of samples where pass is 1 not including the AND fails
    df10a=df4[(df4['OR pass']==True) & (df4['2 reads pass']==True) & (df4['gold_standard']>=1)]
    #print(f'Total grey-xS samples passing gold standard: {df10a.shape[0]} ({df10a.shape[0]/df8.shape[0]*100:.0f}%)')

    # Count the number of samples where pass is 0 and gold_standard is 0
    df11=df4[(df4['pass']!='True') & (df4['gold_standard']==0) & (df4['pathogen']!='unmapped') & (df4['test_type'].isin(['BIOFIRE','BIOFIRE ', 'ALINITY','ALINTY', 'CEPHEID'])) & ((df4['run_pass']=='True') & (df4['barcode_pass']=='True'))]
    print(f'Total x samples not identified by cMG xf: {df11.shape[0]} ({df11.shape[0]/df9.shape[0]*100}%)')
    metrics['Total x samples not identified by cMG xf'] = f'{df11.shape[0]} ({df11.shape[0]/df9.shape[0]*100:.0f}%)'
    df11.to_csv('xf.csv', index=False)

    # Count the number of samples where pass is 0 and gold_standard is 1 not including the AND fails
    df11a=df4[(df4['OR pass']!=True) & (df4['2 reads pass']!=True) & (df4['gold_standard']==0) & (df4['pathogen']!='unmapped')]

    #print(f'Total grey-x samples not identified by cMG xf: {df11a.shape[0]} ({df11a.shape[0]/df9.shape[0]*100}%)')

    # count number of additional yield
    df_ay=df_full[((df_full['pathogen'].isin(biofire_additional)) & (df_full['test'].isin(['ALINITY','ALINTY', 'CEPHEID']))) ]
    df_ay=readjust_pass(df_ay, AND_ratio=AND_ratio, AND_ratio_metric='Sample_reads_percent_of_refs_AuG_truc10')
    
    # readjust TOP_FLU_A for additional yield
    flu_A_options=['Influenza A/H1-2009','Influenza A/H3','Influenza A/H1']
    
    # FLU_A_POS is 1 if biofire positive for Influenza A for the run and barcode
    df_ay['condition'] = np.where(df_ay['pathogen'].isin(flu_A_options), 1, 0)
    df_ay['FLU_A_POS'] = df_ay.groupby(['Run', 'barcode'])['condition'].transform('any')
    filtered_df = df_ay[df_ay['pathogen'].isin(flu_A_options)]

    # Step 2: Group by 'Run' and 'barcode' and find the index of the row with the highest 'num_reads'
    filtered_df['TOP_FLU_A'] = filtered_df.groupby(['Run', 'barcode'])['sample num reads'].transform('idxmax')

    # Step 3: Set the 'top_flu' column to 1 for the row with the highest 'num_reads', and 0 for others
    df_ay['TOP_FLU_A'] = df_ay.index.isin(filtered_df['TOP_FLU_A'])



    # readjust pass criteria for additional yield
    df_ay['pass']=np.where((df_ay['pathogen'].isin(['Influenza A/H1-2009','Influenza A/H1','Influenza A/H3']) & (df_ay['TOP_FLU_A']==False)), 'False', df_ay['pass'])
    # remove runs and batches in list of tuple [('AD_winter_study_130825_rpt050825', 1.0), ('AD_winter_study_220125', 2.0)]
    df_ay=df_ay[~df_ay[['Run','Batch']].apply(tuple, axis=1).isin([('AD_winter_study_130825_rpt050825', 1.0), ('AD_winter_study_220125', 2.0)])]

    df_ay=df_ay[(df_ay['pass']=="True") & (df_ay['PCs_passed']==1) ]

    df_ay.to_csv('tmp/df_ay_pre_filter.csv', index=False)
    
    df_ay.to_csv(f'additional_yield/additional_yield_pass_{AND_ratio}.csv', index=False)
    
    X_not_tested_for=metrics['Total XP pathogens not tested with routine laboratory testing XP']
    print(f'Number of additional yield samples xN/XP:{df_ay.shape[0]}/{X_not_tested_for} ({df_ay.shape[0]/X_not_tested_for*100:.2f}%)')
    metrics['Number of additional yield samples xN/XP'] = f'{df_ay.shape[0]}/{X_not_tested_for} ({df_ay.shape[0]/X_not_tested_for*100:.2f}%)'

    return metrics

def appy_test_types(df, test_type_meta_file):
    """Read in test type meta data and apply to main dataframe
    Check test_to_use column only contains BIOFIRE, ALINITY, CEPHEID
    After merging, check column test_to_use that has missing values are negative controls, with amplification_control or reverse_transcription_control columns == 1
    Except for 'Negative Sequencing control 1' which is also a negative control but wasn't spiked correctly.
    """
    test_type_meta=pd.read_csv(test_type_meta_file)
    test_type_meta.dropna(subset=['sample_name'], inplace=True)
    test_type_meta=test_type_meta[['sample_name', 'test_to_use']]
    test_type_meta['test_to_use']=test_type_meta['test_to_use'].str.strip()
    df=df.merge(test_type_meta, left_on=['seq_name'], right_on=['sample_name'], how='left')
    
    # check test_to_use column only contains BIOFIRE, ALINITY, CEPHEID or is NaN
    valid_test_types=['BIOFIRE', 'ALINITY', 'CEPHEID']
    invalid_test_types=df[~df['test_to_use'].isin(valid_test_types) & (df['test_to_use'].notna())]
    if len(invalid_test_types)>0:
        print('Invalid test types found in test_type_meta_file:')
        print(invalid_test_types[['seq_name', 'test_to_use']].drop_duplicates())
        sys.exit('Please correct test_type_meta_file to only contain valid test types: BIOFIRE, ALINITY, CEPHEID')

    # check that rows with missing test_to_use are negative controls
    missing_test_types=df[df['test_to_use'].isna()]
    non_negative_controls=missing_test_types[(missing_test_types['amplification_control']!=1) & (missing_test_types['reverse_transcription_control']!=1) & (missing_test_types['seq_name']!='Negative Sequencing control 1')]
    if len(non_negative_controls)>0:
        print('Rows with missing test types that are not negative controls found:')
        print(non_negative_controls[['seq_name', 'amplification_control', 'reverse_transcription_control']].drop_duplicates())
        sys.exit('Please correct test_type_meta_file to include test types for all non-negative control samples')

    # replace missing test_to_use with 'Negative Control'
    df['test_to_use']=df['test_to_use'].fillna('Negative Control')

    # save the test type information for scrutiny
    df[['seq_name','Run','barcode','test', 'test_to_use', 'amplification_control', 'reverse_transcription_control']].drop_duplicates().to_csv('test_type_assignment.csv', index=False)

    # replace test with test_to_use
    df['test_type']=df['test_to_use']
    df['test']=df['test_to_use']

    df.drop(columns=['sample_name', 'test_to_use'], inplace=True)
    return df

def run_analysis(AND_ratio=0.1, AND_ratio_metric='Sample_reads_percent_of_refs_AuG_truc10', organisms='biofire_set'):
    df=pd.read_csv('biofire_results_merged.csv')

    #test_meta_file='extracted_samples_storage_and_processing_info_pcr_organisms_test_code_check.csv'
    test_meta_file='meta_data/validation_PCR_results_anonymised.csv'
    df=appy_test_types(df, test_meta_file)

    # remove DNA pathogens 
    if 'no_DNA' in organisms:
        df=df[df['pathogen'].isin(DNA_orgs)==False]
    # remove HRE pathogens
    if 'no_HRE' in organisms:
        df=df[df['pathogen'].isin(HRE_orgs)==False]
    
    metrics={'AND_ratio': AND_ratio}

    df=readjust_pass(df, AND_ratio, AND_ratio_metric=AND_ratio_metric)

    df=repopulate_columns(df)

    df=readjust_spiked_values(df)

    df=readjust_flu_pass(df)

    df_ay=get_additional_yield(df)

    df, df_full=remove_biofire_additional(df)

    metrics=count_samples(df, metrics)

    df, metrics=count_failed_runs_samples(df, metrics)

    df,df2, metrics=count_failed_negative_controls(df, metrics)

    df.to_csv(f'AND_ratios/biofire_results_merged_adjusted_{AND_ratio}.csv', index=False)
    
    metrics=count_passing_samples(df2, metrics)
    
    df2, metrics=count_samples_failing_batch_controls(df, df2, metrics)
    
    metrics=count_samples_passing_batch_controls(df2, metrics)
    
    metrics=count_samples_failing_positive_controls(df2, metrics)
    
    df4, metrics=count_samples_passing_positive_controls(df2, metrics)

    metrics=count_multipathogen_samples(df4, metrics)
    
    metrics=count_multi_gold_standard_samples(df4, metrics)

    df4.to_csv(f'passing_samples/passing_samples_{AND_ratio}.csv', index=False)

    metrics=count_full_passing_samples(df4,df_full, metrics, AND_ratio)
       
    return metrics

if __name__ == '__main__':
    argparse=argparse.ArgumentParser(description='Run Biofire analysis with different AND ratios')
    argparse.add_argument('--AND_ratio_metric', type=str, default='Sample_reads_percent_of_refs_AuG_truc10',
                          help='metric used for calculating AND ratio')
    argparse.add_argument('--organisms', type=str, nargs='+', default='biofire_set',
                          help='set of organisms to include in the analysis: biofire_set, no_DNA, no_HRE')
    args=argparse.parse_args()

    metrics=[]
    ratios=[0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 
            0.07, 0.071, 0.072, 0.073, 0.074, 0.075,
            0.076, 0.077, 0.078, 0.079,
            0.08, 0.081, 0.082, 0.083, 0.084, 0.085, 0.086, 0.087, 0.088, 0.089,
            0.09, 0.091, 0.092, 0.093, 0.094, 0.095, 0.096, 0.097, 0.098, 0.099,  0.1]
    ratios=[0.1, 0.03]
    for ratio in ratios:
        print(f'Running analysis with AND_ratio={ratio}')
        metrics.append(run_analysis(AND_ratio=ratio, AND_ratio_metric=args.AND_ratio_metric, organisms=args.organisms))
    #metrics.append(run_analysis(AND_ratio=0.1))
    #metrics.append(run_analysis(AND_ratio=0.0))

    df=pd.DataFrame(metrics)
    df.to_csv('biofire_analysis_metrics.csv', index=False)
