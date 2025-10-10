from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, auc, roc_auc_score
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import itertools
import os
import argparse
import sys


# metrics to use in model
metrics=['AuG_trunc10', 'sample num reads', 'Cov1_perc',
          'Sample_reads_percent_of_run', 'Sample_reads_percent_of_refs', 
          'Sample_reads_percent_of_type_run', 'Sample_reads_percent_of_type_sample',
          'mean_read_length','mean_aligned_length']

# removed 'sample num reads' from metrics
metrics=['AuG_trunc10', 'Cov1_perc',
          'Sample_reads_percent_of_run', 'Sample_reads_percent_of_refs', 
          'Sample_reads_percent_of_type_run', 'Sample_reads_percent_of_type_sample',
          'mean_read_length','mean_aligned_length']

global_metrics=metrics.copy()

# all other metrics
all_metrics=['meanDepth', 'meanDepth_trunc5', 'meanDepth_trunc10','AuG','AuG_trunc5','AuG_trunc10',
          'bases', 'bases_perc', 
          'Cov1', 'Cov3', 'Cov5', 'Cov10', 
          'Cov1_perc','Cov3_perc','Cov5_perc','Cov10_perc',
          'sample num reads','total run reads mapped', 'total run reads inc unmapped', 
          'Sample_reads_percent_of_run', 'Sample_reads_percent_of_refs', 'Sample_reads_percent_of_type_run', 'Sample_reads_percent_of_type_sample',
          'median_read_length', 'median_aligned_length', 'mean_read_length', 'mean_aligned_length',
          'Sample_num_reads_200', 'Sample_num_reads_300', 'Sample_num_reads_400']

# all_other metrics without metrics
all_other_metrics=[m for m in all_metrics if m not in metrics]
    
# step in range for each metric when mapping threshold to value
range_values={'AuG_trunc10': 0.001, 
             'sample num reads': 0.1, 
             'Cov1_perc': 0.05, 
             'Sample_reads_percent_of_run': 0.00001, 
             'Sample_reads_percent_of_refs': 0.001, 
             'Sample_reads_percent_of_type_run': 0.001, 
             'Sample_reads_percent_of_type_sample': 0.1, 
             'mean_read_length': 1, 
             'mean_aligned_length': 1,
             'AuG_trunc10_Sample_reads_percent_of_refs': 0.01,
             'meanDepth': 0.1,
             'meanDepth_trunc5': 0.1,
             'meanDepth_trunc10': 0.1,
             'AuG': 0.1,
             'AuG_trunc5': 0.1,
             'bases': 10,
             'bases_perc': 0.01,
             'Cov1': 1,
             'Cov3': 1,
             'Cov5': 1,
             'Cov10': 1,
             'Cov3_perc': 0.05,
             'Cov5_perc': 0.05,
             'Cov10_perc': 0.05,
             'total run reads mapped': 10000,
             'total run reads inc unmapped': 10000,
             'median_read_length': 1}

# add all other metrics with a default step of 0.001
for m in all_other_metrics:
    if m not in range_values:
        range_values[m]=0.1

# round values to nearest value for each metric
round_values={'AuG_trunc10': 0.001, 
             'sample num reads': 0.1, 
             'Cov1_perc': 0.05, 
             'Sample_reads_percent_of_run': 0.00001, 
             'Sample_reads_percent_of_refs': 0.001, 
             'Sample_reads_percent_of_type_run': 0.001, 
             'Sample_reads_percent_of_type_sample': 0.1, 
             'mean_read_length': 5, 
             'mean_aligned_length': 5,
             'AuG_trunc10_Sample_reads_percent_of_refs': 0.01,
             'meanDepth': 0.1,
             'meanDepth_trunc5': 0.1,
             'meanDepth_trunc10': 0.1,
             'AuG': 0.1,
             'AuG_trunc5': 0.1,
             'bases': 10,
             'bases_perc': 0.01,
             'Cov1': 1,
             'Cov3': 1,
             'Cov5': 1,
             'Cov10': 1,
             'Cov3_perc': 0.05,
             'Cov5_perc': 0.05,
             'Cov10_perc': 0.05,
             'total run reads mapped': 10000,
             'total run reads inc unmapped': 100000,
             'median_read_length': 1}

# add all other metrics with a default round of 0.001
for m in all_other_metrics:
    if m not in round_values:
        round_values[m]=0.001

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

def get_combinations(metrics):
    ''' Return a tuple of all combinations of metrics'''
    combinations=[]
    for L in range(0, len(metrics)+1): # for length of combinations
    #for L in range(0, 4): # for length of combinations 4
        for subset in itertools.combinations(metrics, L):
            if len(subset)>0:
                combinations.append(subset)
    return combinations

def get_permutations(metrics):
    ''' Return a tuple of all permutations of metrics'''
    combinations=[]
    for L in range(0, len(metrics)+1):
        for subset in itertools.permutations(metrics, L):
            if len(subset)>0:
                combinations.append(subset)
    return combinations

def get_value(maxVal, model, threshold, metric):
    '''Get the minimum value for a given metric that corresponds to a given threshold'''
    if maxVal<1:
        maxVal=1
    if metric in range_values:
        val_range=np.arange(0, maxVal, range_values[metric])
    else:
        val_range=np.arange(0, maxVal, maxVal/1000)
    vals=np.array(val_range).reshape(-1,1)
    probs=model.predict_proba(vals)
    preds=model.predict(vals) 

    df=pd.DataFrame({'value':list(val_range), 'probs':probs[:,1], 'preds':preds})
    df.sort_values(by='value', ascending=True, inplace=True)
    #print(df)
    if df['probs'].max()>=threshold:
        df2=df[df['probs']>=threshold]
        df2=df2.copy()
        minVal=df2['value'].values[0]
    else:
        minVal=df['value'].max()
    
    if metric not in round_values:
        return minVal
    if type(round_values[metric]) == int:
        minVal=round_values[metric] * round(minVal/round_values[metric])
    return minVal

def run_logit(df, metric):
    '''Runs logistic regression on a given metric and returns the ROC curve data'''

    # remove rows with gold standard = 1 and metric = 0
    #df=df_original[~((df_original['gold_standard']==1) & (df_original[metric]==0))]

    # get X and y for metric and gold standard
    X=df[[metric]].values
    y=df['gold_standard']

    # fit logistic regression
    logreg = LogisticRegression(max_iter=1000)
    logreg.fit(X, y)
    y_pred = logreg.predict(X)
    y_probs = logreg.predict_proba(X)[:,1]
    df['logit_probs']=y_probs
    df['logit_pred']=y_pred

    # ROC curve for logistic regression
    #logit_roc_auc = roc_auc_score(y, logreg.predict(X))
    logit_roc_auc = roc_auc_score(y, logreg.predict_proba(X)[:,1])
    fpr2, tpr2, thresholds2 = roc_curve(y, logreg.predict_proba(X)[:,1])
    df4=pd.DataFrame({'fpr':fpr2, 'tpr':tpr2})
    df4[f'thresholds'] = thresholds2
    df4['sensitivity']=df4['tpr']
    df4['specificity']=1-df4['fpr']
    df4['1-specificity']=df4['fpr']
    df4['Youden J statistic']=df4['sensitivity']+df4['specificity']-1
    df4['metric']=metric
    df4.to_csv(f'logreg_roc_data/roc_data_{metric}.csv', index=False)

    # calculate threshold for max Youden J statistic
    max_J = df4[df4['Youden J statistic']==df4['Youden J statistic'].max()]
    threshold2=max_J['thresholds'].values[0]
    df5=df4[df4['thresholds']==threshold2]
    df5=df5.copy()
    
    # get min value for threshold
    minVal=get_value(df[metric].max(), logreg, threshold2, metric)  
    
    ## Standard ROC curve to check absalute values correspond
    # this part just runs ROC directly on the metric so the thresholds are the input values
    # I'm just using this to check that the values from get_value are correct
    fpr, tpr, thresholds = roc_curve(y, X)
    roc_auc = auc(fpr, tpr)
    #roc_auc_score = roc_auc_score(y, X)

    # make data frame for fpr, tpr, thresholds
    df2=pd.DataFrame({'fpr':fpr, 'tpr':tpr})
    df2[f'thresholds'] = thresholds
    df2['sensitivity']=df2['tpr']
    df2['specificity']=1-df2['fpr']
    df2['1-specificity']=df2['fpr']
    df2['Youden J statistic']=df2['sensitivity']+df2['specificity']-1
    df2['metric']=metric
    df2['ROC AUC']=roc_auc
    df2['logit ROC AUC']=logit_roc_auc
    df2.to_csv(f'roc_data/roc_data_{metric}.csv', index=False)

    # calculate threshold for max Youden J statistic
    max_J = df2[df2['Youden J statistic']==df2['Youden J statistic'].max()]
    threshold=max_J['thresholds'].values[0]
    df3=df2[df2['thresholds']==threshold]
    df3=df3.copy()

    # Add min val to output table
    df3['minVal threshold']=minVal
    df3['logit_roc_auc']=logit_roc_auc

    # plot ROC and curve
    plots(df, metric, fpr, tpr, roc_auc)
    
    return df3

def plots(df, metric, fpr, tpr, roc_auc):
    # plot ROC
    plt.figure()
    plt.plot(fpr, tpr, label='Logistic Regression (area = %0.2f)' % roc_auc)
    plt.plot([0, 1], [0, 1],'r--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(f'Receiver operating characteristic for {metric}')
    plt.legend(loc="lower right")
    plt.savefig(f'rocs/Log_ROC_{metric}.pdf')
    plt.clf()
    plt.close()

    # plot curve
    g=sns.regplot(x=metric, y='gold_standard', data=df, logistic=True, ci=None ,scatter_kws={'color': 'black'}, line_kws={'color': 'red'})
    plt.xscale('log')
    plt.savefig(f'curves/Log_curve_{metric}.pdf')
    plt.clf()
    plt.clf()

def test_combination(combination, df, thresholds, combination_ID):
    '''Test a combination of metrics and return sensitivity, specificity, PPV, NPV'''
    num_metrics=len(combination)
    df['prediction']=0
    for metric in combination:
        # predict based on threshold, if value is greater than threshold, predict 1 OR keep previous prediction
        df['prediction']=np.where((df[metric]>=thresholds[metric]) & (df['sample num reads']>=2), 1, df['prediction'])

    df['TP']=np.where((df['gold_standard']==1) & (df['prediction']==1), 1, 0)
    df['FP']=np.where((df['gold_standard']==0) & (df['prediction']==1), 1, 0)
    df['TN']=np.where((df['gold_standard']==0) & (df['prediction']==0), 1, 0)
    df['FN']=np.where((df['gold_standard']==1) & (df['prediction']==0), 1, 0)

    TP=df['TP'].sum()
    FP=df['FP'].sum()
    TN=df['TN'].sum()
    FN=df['FN'].sum()

    sensitivity=TP/(TP+FN)
    specificity=TN/(TN+FP)
    PPV=TP/(TP+FP)
    NPV=TN/(TN+FN)

    youden=sensitivity+specificity-1

    cols=['Run', 'barcode',	'seq_name',	'pathogen',	
        'Biofire positive', 'gold_standard', 'PCs_passed',
        'prediction', 'TP',	'FP',	'TN',	'FN',
        'AuG_trunc10', 'sample num reads', 'Cov1_perc',
        'Sample_reads_percent_of_run', 'Sample_reads_percent_of_refs', 
        'Sample_reads_percent_of_type_run', 
        'Sample_reads_percent_of_type_sample', 
        'mean_read_length', 
        'mean_aligned_length',
        'orthoreovirus', 'zika', 'MS2', 'murine_respirovirus',
        'orthoreovirus passed', 'zika passed', 'MS2 passed', 'murine_respirovirus passed']
    df2=df[cols]

    df2.to_csv(f'predictions/predictions_{combination_ID}_metrics.csv', index=False)

    d={'num_metrics': num_metrics,
        'combination':combination, 
        'sensitivity':sensitivity, 
        'specificity':specificity, 'PPV':PPV, 'NPV':NPV, 'Youden':youden}
    return d
        
def plot_combaintion_results(df):
    # facet scatter plot by mapQ    
    #g=sns.FacetGrid(df, col='mapQ', col_wrap=4,  hue='num_metrics')
    #g.map(sns.scatterplot, 'sensitivity', 'specificity')
    g=sns.scatterplot(x='specificity', y='sensitivity', data=df, hue='num_metrics')
    # add legend
    plt.legend(loc='lower right')
    plt.savefig('plots/sensitivity_specificity.pdf')
    plt.clf()
    plt.close()

    # facet scatter plot by mapQ
    #g2=sns.FacetGrid(df, col='mapQ', col_wrap=4,  hue='num_metrics')
    #g2.map(sns.scatterplot, 'PPV', 'NPV')
    g2=sns.scatterplot(x='NPV', y='PPV', data=df, hue='num_metrics')
    plt.legend(loc='lower right')

    #g2=sns.scatterplot(x='NPV', y='PPV', data=df, hue='num_metrics')
    plt.savefig('plots/NPV_PPV.pdf')
    plt.clf()
    plt.close()

def makeFolders():
    '''Make folders for output files'''
    folders=['rocs', 'curves', 'predictions', 'logreg_roc_data', 'roc_data', 'plots']
    for folder in folders:
        if not os.path.exists(folder):
            os.makedirs(folder)
            print(f'Created folder {folder}')
        else:
            print(f'Folder {folder} already exists')

def run(input_file, set_type, remove_no_reads, use_metrics=False):
    # make some folders
    makeFolders()
    # Step 1: Load data
    df=pd.read_csv(input_file)

    
    if set_type=='derivation':
        # repopulate metrics with value from sample run,barcode groupby
        for col in ['orthoreovirus passed', 'zika passed', 'MS2 passed', 'murine_respirovirus passed',
                   'orthoreovirus', 'zika', 'MS2', 'murine_respirovirus']:

            df[col]=df.groupby(['Run', 'barcode'])[[col]].transform('max')

        # Step 1.1: Remove Flu A with no biofire from gold standard
        flu_A_types=['Influenza A/H1', 'Influenza A/H1-2009', 'Influenza A/H3']
        df=df[~((df['pathogen'].isin(flu_A_types)) & (df['Biofire positive']==0))]
        
        # Step 1.2: remove rows where biofire wasn't run.    
        df_PCR=pd.read_csv('pcr_organisms.csv')
        df=df.merge(df_PCR, on=['Run','barcode'], how='left')

        not_biofire_codes=['FLRA; C2VP', 'FLRP; C2VP']
        df['test']=np.where(df['testcode_checked'].isin(not_biofire_codes), 'ALINITY', None)
        df['test']=np.where(df['testcode_checked'].isin(['RCPR; C2VP','RPCR; C2VP']), 'BIOFIRE', df['test'])
        df['test']=np.where(df['testcode_inferred']=='RCPR', 'BIOFIRE', df['test'])
        df['test_type']=df['test']

        # if 'testcode_checked' == 'FLRA; C2VP' then remove extra biofire organisms
        not_biofire_codes=['FLRA; C2VP', 'FLRP; C2VP']
        df=df[~((df['pathogen'].isin(biofire_additional)) & (df['testcode_checked'].isin(not_biofire_codes)))]

        # Step 1.3: load positive control results and remove rows where positive controls = 0
        #df_positive=pd.read_csv('control_results.csv')
        #spike_reads=df_positive[['Run', 'barcode', 'seq_name', 'pathogen_reduced', 'sample num reads']]
        #df_sr_T=spike_reads.pivot(index=['Run', 'barcode', 'seq_name'], columns='pathogen_reduced', values='sample num reads')
        #print(df_sr_T)

        #df=df.merge(df_sr_T, on=['Run', 'barcode'], how='left')
        #print(df)
        spikes=['orthoreovirus', 'zika', 'MS2', 'murine_respirovirus']
        for spike in spikes:
            df[f'{spike} passed']=np.where(df[spike]>=2, 1, 0)
        #df['PCs_passed']=np.where(df['orthoreovirus passed']+df['MS2 passed']+df['zika passed']+df['murine_respirovirus passed']>=2, 1, 0)
        df['PCs_passed']=np.where(df['orthoreovirus passed']+df['zika passed']+df['murine_respirovirus passed']>=2, 1, 0)

    # choose metrics
    if use_metrics == False:
        metrics=global_metrics
        print('Using normal metrics')
    elif use_metrics == True:
        print('Using specific metrics')
        metrics=all_other_metrics
        
    if use_metrics == False:
        # Add in all the permuations of metrics to create ratio
        print('Getting lists of metric permutations')
        perms=get_permutations(metrics)
        perms2=[p for p in perms if len(p)==2]
        for perm in perms2:
            metric1=perm[0]
            metric2=perm[1]
            print(f'Creating ratio of {metric1} and {metric2}')
            #df[metric1].fillna(0, inplace=True)
            #df[metric2].fillna(0, inplace=True)
            df[f'{metric1}_{metric2}_ratio']=df[metric1]/df[metric2]
            #df[f'{metric1}_{metric2}_ratio'].fillna(0, inplace=True)
            df.fillna({f'{metric1}_{metric2}_ratio': 0}, inplace=True)
            # replace inf with 0
            #df[f'{metric1}_{metric2}_ratio'].replace([np.inf, -np.inf], 0, inplace=True)
            df[f'{metric1}_{metric2}_ratio'] = df[f'{metric1}_{metric2}_ratio'].replace([np.inf, -np.inf], 0)
            #df.replace({f'{metric1}_{metric2}_ratio': [np.inf, -np.inf], 0}, inplace=True)
        #    metrics.append(f'{metric1}_{metric2}_ratio')

    # save data
    #df.to_csv('train_data.csv', index=False)

    # remove failed batches if validation set
    if set_type=='validation':
        print(df['pass'].unique())
        # count number of samples that failed batch ampification negative controls
        df_RC_control=df[(df['reverse_transcription_control']==1) & (df['IC_virus_spike']==1)]
        # orthoreovirus passed	zika passed	MS2 passed	murine_respirovirus passed
        df_RC_control_PCFAIL=df_RC_control[(df_RC_control['orthoreovirus passed']==0) & (df_RC_control['zika passed']==0) & (df_RC_control['murine_respirovirus passed']==0) \
                                    | (df_RC_control['MS2 passed']==1) \
                                    | (df_RC_control['pass']==True) ]
        failed_batches=list(df_RC_control_PCFAIL.drop_duplicates(['Run','Batch'])[['Run','Batch']].itertuples(index=False, name=None))
        print(f'Failed batches: {failed_batches}')
        df=df[~df[['Run','Batch']].apply(tuple, axis=1).isin(failed_batches)]
        nsamples=df.groupby(['Run', 'barcode'])[['pathogen']].nunique()
        print(f'Number of unique samples after removing failed batches: {len(nsamples)}')

    # print number of unique samples
    nsamples=df.groupby(['Run', 'barcode'])[['pathogen']].nunique()
    print(f'Number of unique samples inc negs: {len(nsamples)}')
    # remove negs with no test_type
    print(df['test_type'].unique())
    df=df[df['test_type'].isin(['BIOFIRE', 'BIOFIRE ', 'ALINITY', 'ALINTY', 'CEPHEID'])]
    # print number of unique samples
    nsamples=df.groupby(['Run', 'barcode'])[['pathogen']].nunique()
    print(f'Number of unique samples not inc negs: {len(nsamples)}')

    # remove 010524_Expt9_SISPA_Daisy that failed negative controls but passed run/sample controls
    df=df[~(df['Run']=='010524_Expt9_SISPA_Daisy')]
    # count number of samples that failed PCs
    n_failed_PCs=df[df['PCs_passed']==0][['Run', 'barcode']].drop_duplicates()
    print(f'Number of samples that failed PCs: {len(n_failed_PCs)}')
    #print(n_failed_PCs)
    # remove rows where PCs failed
    df=df[df['PCs_passed']>0]
    df.to_csv('train_data_PCs_passed.csv', index=False)

    # if remove_no_reads is True then remove rows where sample num reads = 0
    if remove_no_reads:
        df=df[df['sample num reads']>0]
        print(f'Removed samples with no reads, new number of samples: {len(df.groupby(["Run", "barcode"])["pathogen"].nunique())}')

    # print number of unique samples
    nsamples=df.groupby(['Run', 'barcode'])[['pathogen']].nunique()
    print(f'Number of unique samples: {len(nsamples)}')
    # print number of unique samples with gold standard = 1
    nsamples_gs=len(df[df['gold_standard']==1])
    print(f'Number of unique samples with gold standard = 1: {nsamples_gs}')

    # Step 2: Run logistic regression for each metric and get thresholds for max Youden J statistic
    dfs=[]

    for metric in metrics:
        print('running logit for', metric)  
        df2=run_logit(df, metric)      
        dfs.append(df2)

    df4=pd.concat(dfs)
    df4.to_csv('roc_data_all.csv', index=False)

    # stop if use_metrics is True
    if use_metrics == True:
        print('use_metrics is True, stopping after calculating metrics')
        sys.exit()

    # Step 3: create a dictionary of thresholds from df4 with metric as key and threshold as value
    print('creating dictionary of thresholds')
    df5=df4[['metric', 'minVal threshold']]
    df5.set_index('metric', inplace=True)
    thresholds=df5.to_dict()['minVal threshold']
    
    # Step 4: Run through all combiantions of metrics
    print('getting all combinations')
    combinations=get_combinations(metrics)
    combination_results=[]
    print('running through all combinations')
    for i,combination in enumerate(combinations):
        #for mapQ in mapQs:
        #    df3=df[(df['mapQ']==mapQ) & (df['pathogen'].isin(bacteria_chroms))]
        #    df4=df[(df['mapQ']==0) & (~df['pathogen'].isin(bacteria_chroms))]
        #    df2=pd.concat([df3, df4])
        #print('testing combination', i, combination)
        d=test_combination(combination, df, thresholds, i)
        d['combination number']=i
        #d['mapQ']=mapQ
        combination_results.append(d)
        #print(combination)

        # test combination ratio


    df6=pd.DataFrame(combination_results)
    df6.sort_values(by='Youden', ascending=False, inplace=True)
    df6.to_csv('combination_results.csv', index=False)

    # Step 5: Plot results
    plot_combaintion_results(df6)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run optimization models")
    parser.add_argument("--input_file", type=str, default="biofire_results_merged.csv", help="Path to input CSV file")
    parser.add_argument("--set", type=str, default="derivation", help="derivation or validation set")
    parser.add_argument("--remove_no_reads", type=bool, default=False, help="Remove samples with no reads")
    parser.add_argument("--use_metrics", type=bool, default=False, help="Use specific metrics")

    args = parser.parse_args()

    run(args.input_file, args.set, args.remove_no_reads, use_metrics=args.use_metrics)