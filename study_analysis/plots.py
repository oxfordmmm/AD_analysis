#!/usr/bin/env python3
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import sys
from mizani.formatters import label_number, label_comma
from scipy.stats import spearmanr
from scipy.interpolate import make_interp_spline
from statsmodels.nonparametric.smoothers_lowess import lowess
from plotnine import (
    ggplot,
    aes,
    after_stat,
    stage,
    geom_bar,
    geom_point,
    geom_jitter,
    geom_boxplot,
    geom_hline,
    geom_smooth,
    geom_text,
    facet_wrap,
    facet_grid,
    theme,
    theme_minimal,
    element_text,
    scale_fill_manual,
    scale_x_discrete,
    scale_y_discrete,
    scale_y_log10,
    scale_x_reverse,
    scale_color_discrete,
    geom_bin_2d,
    stat_bin_2d,
    labs,
    scale_shape_manual,
    scale_color_manual,
    scale_fill_manual,
    geom_line,
    geom_hline,
    geom_vline,
    scale_x_log10,
    scale_y_log10,
    annotate,
    element_blank)

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

test_type_mapping = {
        'Biofire': 'Biofire\nFilmArray Respiratory\n2.1 Panel\n(23 pathogen targets)',
        'Biofire ': 'Biofire\nFilmArray Respiratory\n2.1 Panel\n(23 pathogen targets)',
        'BIOFIRE': 'Biofire\nFilmArray Respiratory\n2.1 Panel\n(23 pathogen targets)',
        'BIOFIRE ': 'Biofire\nFilmArray Respiratory\n2.1 Panel\n(23 pathogen targets)',
        'Alinity': 'Alinity\nm Resp-4-plex assay\n(4 pathogen targets)',
        'ALINITY': 'Alinity\nm Resp-4-plex assay\n(4 pathogen targets)',
        'ALINTY': 'Alinity\nm Resp-4-plex assay\n(4 pathogen targets)',
        'Cepheid': 'Cepheid\nXpert Xpress\nCoV-2/Flu/RSV assay\n(4 pathogen targets)',
        'CEPHEID': 'Cepheid\nXpert Xpress\nCoV-2/Flu/RSV assay\n(4 pathogen targets)',
    }

def plot_cts_vs_reads():
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

def plot_Ct_reads_plotnine(df):
    df=df.copy()
    df['CT value'] = np.where(df['ct_1']!=0, df['ct_1'], df['ct_2'])
    df2=df[df['CT value']!=0]
    df2=df2.copy()
    df2=df2[df2['pathogen']!='Influenza B']

    p=(ggplot(df2, aes(x='CT value', y='sample num reads', color='pathogen')) +
       geom_point() +
       #geom_smooth(method='lm', se=True) + 
       geom_smooth(se=True) + 
       scale_y_log10(labels=label_comma()) +  # <- plain numbers
       scale_x_reverse() + # reverse x axis
       labs(x='Ct Value', y='Reads') +
       theme(legend_title=element_blank())
    )
    p.save('Figure04.pdf', width=6, height=4)

    # same again but with reads per 10k sample reads
    df2['reads_per_10k'] = (df2['sample num reads'] / df2['total sample reads']) * 10000
    p2=(ggplot(df2, aes(x='CT value', y='reads_per_10k', color='pathogen')) +
         geom_point() +
         geom_smooth(se=True) +
         scale_y_log10() + # log y axis
         scale_x_reverse() + # reverse x axis
         labs(x='Ct Value', y='Reads per 10k total sample reads'))
    p2.save('Figure04_per10kreads.pdf', width=6, height=4)

    # just RSV
    df3=df2[df2['pathogen']=='RSV']
    p3=(ggplot(df3, aes(x='CT value', y='reads_per_10k', color='chrom')) +
        geom_point() +
        geom_smooth(se=True) +
        scale_y_log10() + # log y axis
        scale_x_reverse() + # reverse x axis
        labs(x='Ct Value', y='Reads per 10k total sample reads'))
    p3.save('Figure04_RSV.pdf', width=6, height=4)

def plot_sensitivity_by_run(df):
    pass
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

def plot_sensitivity_by_run_and_pathogen(df):
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

def plot_sensitivity_by_pathogen(df_original):
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
    plt.savefig('pathogen_pass_fails.pdf')
    plt.clf()

    print(g1['passed'].sum())
    print(g1['Failed criteria'].sum())
    print(g1['Failed PCs'].sum())

    g1['sensitivity'] = g1['passed'] / (g1['passed'] + g1['Failed criteria'])
    print(g1)


def get_CI(df, set=None):
    # calculate 95% confidence intervals for sensitivity
    from statsmodels.stats.proportion import proportion_confint
    df = df.copy()
    df['sensitivity'] = (df['Pathogen detected'] / (df['Pathogen detected'] + df['Failed criteria']))*100
    # using normal approximation method
    methods=['beta','normal','agresti_coull','wilson','jeffreys','binom_test']
    for method in methods:
        ci_lower, ci_upper = proportion_confint(df['Pathogen detected'], df['Pathogen detected'] + df['Failed criteria'], 
                                                alpha=0.05, method=method)
        df[f'CI lower {method}'] = ci_lower * 100
        df[f'CI upper {method}'] = ci_upper * 100
    #ci_lower, ci_upper = proportion_confint(df['Pathogen detected'], df['Pathogen detected'] + df['Failed criteria'], 
    #                                        alpha=0.05, method='beta')
    #df['CI lower'] = ci_lower
    #df['CI upper'] = ci_upper
    print(df)
    if set:
        df.to_csv(f'pathogen_sensitivity_CI_{set}.csv')
    else:
        df.to_csv('pathogen_sensitivity_CI.csv')
    #return df

def get_pathogen_sensistivity_data(df_original, set=None):
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
    g2 = g2.rename(columns={'TP': 'Pathogen detected'})
    g0 = g0.rename(columns={'gold_standard': 'all_gold_standard_count'})
    g1 = g1.merge(g2, on=groups)
    g1 = g1.merge(g0, on=groups)
    g1['Failed criteria'] = g1['gold_standard_count'] - g1['Pathogen detected']
    g1['Failed PCs'] = g1['all_gold_standard_count'] - (g1['Pathogen detected'] + g1['Failed criteria'])

    print(g1)
    g1.drop(columns=['gold_standard_count','all_gold_standard_count'], inplace=True)
    print(g1)
    print(g1.index)

    # remove 'unmapped' row if it exists
    g1= g1[g1['pathogen'] != 'unmapped']

    g1 = g1.set_index(groups)
    print(g1)

    # save to csv
    get_CI(g1, set=set)

    # convert to long format for plotnine
    g1 = g1.reset_index().melt(id_vars=['pathogen'], var_name='Result', value_name='count')
    return g1

def plot_sensitivity_by_pathogen_plotnine(df_original, dfd):
    df_original = df_original.copy()  # Use the original DataFrame for grouping
    df_original= df_original[~((df_original['Run']=='AD_winter_study_220125') & (df_original['Batch']==2))]
    g1=get_pathogen_sensistivity_data(df_original, set='Validation')
    g1['Set']='Validation'
    print(g1)

    # repeat with derivation set
    dfd.to_csv('dfd.csv', index=False)
    g2=get_pathogen_sensistivity_data(dfd, set='Derivation')
    g2['Set']='Derivation'
    g1=pd.concat([g1, g2])

    # remove _ from pathogen names
    g1['pathogen'] = g1['pathogen'].str.replace('_', ' ')

    g1.to_csv('pathogen_sensitivity_data.csv', index=False)
    # Red, blue, green colourbline safe
    result_colors = ["#393b79", '#ff7f0e', '#2ca02c']
    #result_colors = ['#1f77b4', '#ff7f0e', '#2ca02c'] 
    
    # stacked bar plot of pathogen, with pass, gold_standard - pass, all_gold_standard - gold_standard as different colours
    p=(ggplot(g1, aes(x='pathogen', y='count', fill='Result')) +
       geom_bar(stat='identity', position='stack') +
       scale_fill_manual(values=result_colors, name='Pathogen') +
       facet_wrap('~Set') +
       labs(x='Pathogen',
            y='Count',
            fill='Result') +
       theme(axis_text_x=element_text(rotation=90, hjust=1))
       )
    p.save('Figure_S6.pdf', width=10, height=6)

def plot_gold_standard_pathogen_counts(df, dfd):
    # plot gold_standard pathogen counts by run and batch
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
    # remove AD_winter_study_220125 batch 2
    df= df[~((df['Run']=='AD_winter_study_220125') & (df['Batch']==2))]
    df['Run_order'] = df['Run'].map(run_order)  # Map the run names to their order
    df.sort_values(by=['Run_order', 'Batch'], inplace=True)  # Sort by Run order and Batch
    df['Run']= 'Run ' + df['Run_order'].astype(str).str.zfill(2) + ' batch ' + df['Batch'].astype(int).astype(str)  # Combine Run and Batch into a single column, add leading zero to run number
    run_ordered = df['Run'].unique().tolist()
    run_order = {name: i for i, name in enumerate(run_ordered)}  # Create a mapping of run names to their order

    # correct df test_type
    df['test_type'] = df['test_type'].map(test_type_mapping)
    dfd['test_type'] = dfd['test_type'].map(test_type_mapping)
    dfd.to_csv('dfd.csv', index=False)
    g1 = df.groupby(['Run', 'pathogen', 'test_type'])[['gold_standard']].sum().reset_index()
    g1 = g1.set_index(['Run', 'pathogen', 'test_type'])
    g1 = g1.unstack(level='pathogen')
    g1.columns = g1.columns.droplevel(0)  # Drop the top level of the MultiIndex
    g1 = g1.fillna(0)  # Fill NaN values with 0
    g1.to_csv('validation_gold_standard_counts.csv', index=True)

    g2 = df.groupby(['Run', 'pathogen','test_type'])[['barcode']].count().reset_index()
    g2 = g2.set_index(['Run', 'pathogen', 'test_type'])
    g2 = g2.unstack(level='pathogen')
    g2.columns = g2.columns.droplevel(0)  # Drop the top level of the MultiIndex
    g2 = g2.fillna(0)  #
    g2.to_csv('validation_test_counts.csv')
    
    # sort dataframe by first index in MultiIndex by run_order
    g1['Run_order'] = g1.index.get_level_values('Run').map(run_order)  # Map the run names to their order
    g1 = g1.sort_values(by='Run_order')  # Sort by Run order
    g1 = g1.drop(columns=['Run_order'])  # Drop the temporary Run_order column      

    # drop 'unmapped' column if it exists
    if 'unmapped' in g1.columns:
        g1 = g1.drop(columns=['unmapped'])
    # drop empty columns
    g1 = g1.loc[:, (g1 != 0).any(axis=0)]
    ax=g1.plot(kind='bar', stacked=True, 
            colormap='tab20',  # Use a colormap for better color differentiation
            figsize=(10, 6))
    # make y scale full integers
    ax.yaxis.get_major_locator().set_params(integer=True)

    plt.title('Gold Standard Pathogen Counts by Run and Batch')
    plt.xlabel('Run and Batch')
    plt.ylabel('Pathogen Count')
    plt.xticks(rotation=90)
    # move legend outside of plot
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    plt.tight_layout()
    plt.savefig('gold_standard_pathogen_counts_by_run_and_batch.pdf')
    plt.clf()

    vgi=g1.copy()
    vg2=g2.copy()

    g1 = dfd.groupby(['Run', 'pathogen','test_type'])[['gold_standard']].sum().reset_index()
    g1 = g1.set_index(['Run', 'pathogen', 'test_type'])
    g1 = g1.unstack(level='pathogen')
    g1.columns = g1.columns.droplevel(0)  # Drop the top level of the MultiIndex
    g1 = g1.fillna(0)  # Fill NaN values with 0
    g1.to_csv('derivation_gold_standard_counts.csv', index=True)

    # total test counts
    g2 = dfd.groupby(['Run', 'pathogen','test_type'])[['barcode']].count().reset_index()
    g2 = g2.set_index(['Run', 'pathogen', 'test_type'])
    g2 = g2.unstack(level='pathogen')
    g2.columns = g2.columns.droplevel(0)  # Drop the top level of the MultiIndex
    g2 = g2.fillna(0)  #

    # sort dataframe by first index in MultiIndex by run_order
    #g1['Run_order'] = g1.index.get_level_values('Run').map(run_order)  # Map the run names to their order
    #g1 = g1.sort_values(by='Run_order')  # Sort by Run order
    #g1 = g1.drop(columns=['Run_order'])  # Drop the temporary Run_order column      

    # drop 'unmapped' column if it exists
    if 'unmapped' in g1.columns:
        g1 = g1.drop(columns=['unmapped'])

    # drop empty columns
    g1 = g1.loc[:, (g1 != 0).any(axis=0)]
    print(g1)
    g1.to_csv('derivation_gold_standard_counts.csv', index=True)

    # get columns from vg1
    cols=vgi.columns
    for col in cols:
        if col not in g1.columns:
            g1[col]=0
    #g1=g1[cols]  # reorder columns to match vgi
    print(g1)
    g1.to_csv('derivation_gold_standard_counts.csv', index=True)

    ax=g1.plot(kind='bar', stacked=True, 
            colormap='tab20',  # Use a colormap for better color differentiation
            figsize=(10, 6))

    # make y scale full integers
    ax.yaxis.get_major_locator().set_params(integer=True)

    plt.title('Gold Standard Pathogen Counts by Run and Batch')
    plt.xlabel('Run and Batch')
    plt.ylabel('Pathogen Count')
    plt.xticks(rotation=90)
    # move legend outside of plot
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    plt.tight_layout()
    plt.savefig('gold_standard_pathogen_counts_by_run_and_batch_derivation.pdf')
    plt.clf()

    # combine derivation and validation set plots
    vgi['Dataset']='Validation'
    vg2['Dataset']='Validation'
    
    g1['Dataset']='Derivation'
    g2['Dataset']='Derivation'
    #g1.reset_index(inplace=True,drop=True)
    #vgi.reset_index(inplace=True,drop=True)
    combined=pd.concat([g1, vgi])
    combined_test_counts=pd.concat([g2, vg2])
    combined_test_counts.to_csv('combined_test_counts.csv', index=True)

    ax2=combined.plot(kind='bar', stacked=True, column='Dataset',
            colormap='tab20',  # Use a colormap for better color differentiation
            figsize=(10, 6))        

    # make y scale full integers
    ax2.yaxis.get_major_locator().set_params(integer=True)

    plt.title('Gold Standard Pathogen Counts by Run and Batch')
    plt.xlabel('Run and Batch')
    plt.ylabel('Pathogen Count')
    plt.xticks(rotation=90)
    # move legend outside of plot
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    plt.tight_layout()
    plt.savefig('gold_standard_pathogen_counts_by_run_and_batch_combined.pdf')
    plt.clf()
    return combined, combined_test_counts


def plot_pathogen_counts(combined, combined_test_counts):
    # try with plotnine
    # convert combined to long format
    combined.reset_index(inplace=True)
    combined_long = combined.melt(id_vars=['Run', 'Dataset', 'test_type'], var_name='Pathogen', value_name='Count')
    combined_long = combined_long[combined_long['Run'] != 0]  # Filter out zero Runs
    combined_test_counts.reset_index(inplace=True)
    combined_test_counts_long = combined_test_counts.melt(id_vars=['Run', 'Dataset', 'test_type'],
                                                          var_name='Pathogen', value_name='Total_Count')
    combined_test_counts_long = combined_test_counts_long[combined_test_counts_long['Run'] != 0]  #Filter out zero Runs
    run_ordered = combined_long['Run'].unique().tolist()
    # remove _ from pathogen names for better display
    combined_long['Pathogen'] = combined_long['Pathogen'].str.replace('_', ' ')
    combined_test_counts_long['Pathogen'] = combined_test_counts_long['Pathogen'].str.replace('_', ' ')
    # map test_type to different names
    
    test_types=combined_long['test_type'].unique()
    for test_type in test_types:
        if test_type not in test_type_mapping:
            print(f'Warning: test type "{test_type}" not found in mapping, using original value')
            test_type_mapping[test_type] = test_type  # use original value if not found in mapping
    combined_long=combined_long[combined_long['test_type']!='0']
    combined_long['test_type'] = combined_long['test_type'].map(test_type_mapping)
    combined_long.to_csv('gold_standard_pathogen_counts_by_run_and_batch_combined_long.csv', index=False)
    combined_test_counts_long=combined_test_counts_long[combined_test_counts_long['test_type']!='0']
    combined_test_counts_long['test_type'] = combined_test_counts_long['test_type'].map(test_type_mapping)
    combined_test_counts_long.to_csv('total_test_counts_by_run_and_batch_combined_long.csv', index=False)

    pathogen_colors = [
        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
        "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#a55194", "#393b79",
        "#637939", "#8c6d31", "#843c39", "#7b4173", "#ad494a", "#3182bd"
    ]

    test_colours = [
        "#e377c2", "#17becf", "#2ca02c"
    ]    
 
    p = (ggplot(combined_long, aes(x='Run', y='Count', fill='Pathogen'))
        + geom_bar(stat='identity', position='stack')
        + theme(axis_text_x=element_text(rotation=90, hjust=1))
        + scale_fill_manual(values=pathogen_colors, name='Pathogen')
        + facet_wrap('~Dataset', scales='free_x')
        )

    p.save('gold_standard_pathogen_counts_by_run_and_batch_combined_plotnine.pdf', width=10, height=6)
    p.save('Figure02_OG.pdf', width=10, height=6)

    # try plotting pathogen on the x axis and counts on the y axis, facet by Dataset
    p2 = (ggplot(combined_long, aes(x='Pathogen', y='Count', fill='test_type'))
        + geom_bar(stat='identity', position='stack')
        + theme(axis_text_x=element_text(rotation=90, hjust=1),
               legend_key_spacing_y = 20, # adjust lengend spacing to avoid overlapping text
               legend_text=element_text(ha='left', ma='left'),
               legend_title=element_text(ha='left', ma='left', size=10, weight='bold', linespacing=3)
        )    
        + scale_fill_manual(values=pathogen_colors, name='Test platform')
        + facet_wrap('~Dataset', scales='free_x')
        ) 
    p2.save('gold_standard_pathogen_counts_by_pathogen_and_batch_combined_plotnine.pdf', width=10, height=6)
    p2.save('Figure02.pdf', width=10, height=6)

    ## Add negative counts
    combined_test_counts_long=combined_test_counts_long[combined_test_counts_long['Pathogen']!='unmapped']
    print(len(combined_long[combined_long['Pathogen']=='Rhinovirus/enterovirus']))
    combined_long=combined_long.merge(combined_test_counts_long, on=['Run', 'Dataset', 'test_type', 'Pathogen'],how='outer')
    
    combined_long=combined_long.groupby(['Dataset', 'test_type', 'Pathogen'])[['Count', 'Total_Count']].sum().reset_index()
    combined_long['Negative tests']= combined_long['Total_Count'] - combined_long['Count']
    combined_long.to_csv('combined_long_with_negatives.csv')
    

    # plot negative counts by pathogen
    p3 = (ggplot(combined_long, aes(x='Pathogen', y='Negative tests', fill='test_type'))
        + geom_bar(stat='identity', position='stack')
        + theme(axis_text_x=element_text(rotation=90, hjust=1),
               legend_key_spacing_y = 20, # adjust lengend spacing to avoid overlapping text
               legend_text=element_text(ha='left', ma='left'),
               legend_title=element_text(ha='left', ma='left', size=10, weight='bold', linespacing=3)
        )    
        + scale_fill_manual(values=pathogen_colors, name='Test platform')
        + facet_wrap('~Dataset', scales='free_x')
        ) 
    #p3.save('gold_standard_pathogen_counts_by_pathogen_and_batch_combined_plotnine.pdf', width=10, height=6)
    p3.save('Figure02_negative_counts.pdf', width=10, height=6)

    # melt combined_long to long format for plotting positives and negatives together
    combined_long['Positive tests'] = combined_long['Count']
    combined_melted = combined_long.melt(id_vars=['Dataset', 'test_type', 'Pathogen'],
                                        value_vars=['Positive tests', 'Negative tests'], var_name='Result', value_name='Number of tests')
    # order Result so that Positive tests comes before Negative tests
    combined_melted['Result'] = pd.Categorical(combined_melted['Result'], categories=['Positive tests', 'Negative tests'], ordered=True)
    
    p4 = (ggplot(combined_melted, aes(x='Pathogen', y='Number of tests', fill='test_type'))
        + geom_bar(stat='identity', position='stack')
        + theme(axis_text_x=element_text(rotation=90, hjust=1),
               legend_key_spacing_y = 20, # adjust lengend spacing to avoid overlapping text
               legend_text=element_text(ha='left', ma='left'),
               legend_title=element_text(ha='left', ma='left', size=10, weight='bold', linespacing=3)
        )
        + scale_fill_manual(values=pathogen_colors, name='Test platform')
        + facet_grid(rows='Result', cols='Dataset', scales='free_y')
        ) 
    p4.save('Figure02_positives_negatives.pdf', width=10, height=6)


def plot_reads_bases(df_original, dfd):
    ## Plot number of reads and bases per sample
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
    df4=df_original.copy()
    df4.dropna(subset=['total run bases','total run reads'], inplace=True)
    df4=df4[df4['total run reads']!=0]
    df4=df4[df4['total run bases']!=0]
    df4.drop_duplicates(subset=['Run', 'barcode'], inplace=True)
    df5=df4.copy()
    df5['Run_order'] = df5['Run'].map(run_order)  # Map the run names to their order
    df5.sort_values(by=['Run_order'], inplace=True)  # Sort by Run order 
    df5['Run']= 'Run ' + df5['Run_order'].astype(str).str.zfill(2)
    dfd['Dataset']='Derivation'
    dfd.dropna(subset=['total run bases','total run reads'], inplace=True)
    dfd=dfd[dfd['total sample reads']!=0]
    dfd=dfd[dfd['total sample bases']!=0]
    dfd.drop_duplicates(subset=['Run', 'barcode'], inplace=True)
    df5['Dataset']='Validation'
    df5=pd.concat([df5, dfd])
    df5['Million Reads'] = df5['total sample reads'] / 1e6
    df5['Mega Bases'] = df5['total sample bases'] / 1e6
    df5=df5[['Run', 'barcode','Dataset','Million Reads', 'Mega Bases']]
    df5=df5.drop_duplicates(subset=['Run', 'barcode', 'Dataset'], keep='last')
    df5['seqID'] = list(range(1, len(df5) + 1))
    print(df5)

    # melt df5 to long format
    df5_long = df5.melt(id_vars=['Run', 'barcode', 'seqID', 'Dataset'], var_name='Metric', value_name='Value')
    print(df5_long)

    # plot bar plot of total run reads inc unmapped by run
    p3 = (ggplot(df5_long, aes(x='Run', y='Value', fill='Metric'))
        + geom_jitter(stat='identity', position='jitter')
        + geom_boxplot(aes(group='Run'), outlier_shape=None, alpha=0.5)
        + theme(axis_text_x=element_text(rotation=90, hjust=1))
        + facet_grid(rows='Metric', cols='Dataset', scales='free')
        # draw horizontal line at y=0.025 for million reads and y=0.4 for mega bases
        #+ geom_hline(yintercept=0.025, linetype='dashed', color='red')
        #+ geom_hline(yintercept=0.4, linetype='dashed', color='blue')
        # + scale_y_discrete(name='Total run reads inc unmapped')
        )
    p3.save('total_sample_reads_plotnine.pdf', width=10, height=6)

def plot_reads_bases_per_run(df_original, dfd):
    ## Plot number of reads per run
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
    df4=df_original.copy()
    df4.dropna(subset=['total run bases','total run reads'], inplace=True)
    df4=df4[df4['total run reads']!=0]
    df4=df4[df4['total run bases']!=0]
    df4.drop_duplicates(subset=['Run', 'barcode'], inplace=True)
    df5=df4.copy()
    df5['Run_order'] = df5['Run'].map(run_order)  # Map the run names to their order
    df5.sort_values(by=['Run_order'], inplace=True)  # Sort by Run order 
    df5['Run']= 'Run ' + df5['Run_order'].astype(str).str.zfill(2)
    dfd['Dataset']='Derivation'
    dfd.dropna(subset=['total run bases','total run reads'], inplace=True)
    dfd=dfd[dfd['total run reads']!=0]
    dfd=dfd[dfd['total run bases']!=0]
    dfd.drop_duplicates(subset=['Run', 'barcode'], inplace=True)
    df5['Dataset']='Validation'
    df5=pd.concat([df5, dfd])
    df5['Million Reads'] = df5['total run reads'] / 1e6
    df5['Mega Bases'] = df5['total run bases'] / 1e6
    df5=df5[['Run', 'Dataset','Million Reads', 'Mega Bases']]
    df5=df5.drop_duplicates(subset=['Run', 'Dataset'], keep='last')
    print(df5)

    # melt df5 to long format
    df5_long = df5.melt(id_vars=['Run', 'Dataset'], var_name='Metric', value_name='Value')
    print(df5_long)

    # plot bar plot of total run reads inc unmapped by run
    p3 = (ggplot(df5_long, aes(x='Run', y='Value', fill='Metric'))
        + geom_bar(stat='identity', position='dodge')
        + theme(axis_text_x=element_text(rotation=90, hjust=1))
        + facet_grid(rows='Metric', cols='Dataset', scales='free')
        # + scale_y_discrete(name='Total run reads inc unmapped')
        )
    p3.save('total_run_reads_inc_unmapped_by_run_plotnine.pdf', width=10, height=6)


def plot_figure_S2(df_original, dfd):
    ## Plot number of reads and bases per sample
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
    df4=df_original.copy()
    df4.dropna(subset=['total run bases','total run reads'], inplace=True)
    df4=df4[df4['total run reads']!=0]
    df4=df4[df4['total run bases']!=0]
    df4.drop_duplicates(subset=['Run', 'barcode'], inplace=True)
    df5=df4.copy()
    df5['Run_order'] = df5['Run'].map(run_order)  # Map the run names to their order
    df5.sort_values(by=['Run_order'], inplace=True)  # Sort by Run order 
    df5['Run']= 'Run ' + df5['Run_order'].astype(str).str.zfill(2)
    dfd['Dataset']='Derivation'
    dfd.dropna(subset=['total run bases','total run reads'], inplace=True)
    dfd=dfd[dfd['total sample reads']!=0]
    dfd=dfd[dfd['total sample bases']!=0]
    dfd.drop_duplicates(subset=['Run', 'barcode'], inplace=True)
    df5['Dataset']='Validation'
    df5=pd.concat([df5, dfd])
    df5['Thousand Reads'] = df5['total sample reads'] / 1e3
    df5['Gigabases'] = df5['total run bases'] / 1e9
    df5=df5[['Run', 'seq_name','barcode','Dataset','Thousand Reads', 'Gigabases']]
    df5=df5.drop_duplicates(subset=['Run', 'barcode', 'Dataset'], keep='last')
    df5['Run_Sample'] = 'Sample'
    df5['Threshold']=0.025
    df6=df5.copy()
    df6['Run_Sample'] = 'Run'
    df6['Threshold']=0.4
    df6=df6.drop_duplicates(subset=['Run', 'Dataset' ], keep='last')
    df5=pd.concat([df5, df6])
    #df5['seqID'] = list(range(1, len(df5) + 1))
    print(df5[df5['Thousand Reads']<25])

    # melt df5 to long format
    df5_long = df5.melt(id_vars=['Run','seq_name','barcode', 'Dataset', 'Run_Sample', 'Threshold'], var_name='Metric', value_name='Value')
    #print(df5_long[df5_long['Thousand Reads']<0.025])

    # plot bar plot of total run reads inc unmapped by run
    p3 = (ggplot(df5_long, aes(x='Run', y='Value', fill='Metric'))
        + facet_grid(rows='Metric', cols='Dataset', scales='free')
        #+ facet_wrap('~Run_Sample', scales='free_x')
        + geom_jitter(data=df5_long[(df5_long['Run_Sample'] == 'Sample') & (df5_long['Metric'] == 'Thousand Reads')], stat='identity', position='jitter')
        + geom_boxplot(aes(group='Run'), outlier_shape=None, alpha=0.5, data=df5_long[(df5_long['Run_Sample'] == 'Sample') & (df5_long['Metric'] == 'Thousand Reads')])
        + geom_bar(data=df5_long[(df5_long['Run_Sample'] == 'Run') & (df5_long['Metric'] == 'Gigabases')], stat='identity', position='dodge')
        + theme(axis_text_x=element_text(rotation=90, hjust=1))
        # draw horizontal line at y=0.025 for thousand reads and y=0.4 for gigabases
        + geom_hline(aes(yintercept=0.4), linetype='dashed', color='red', data=df5_long[(df5_long['Run_Sample'] == 'Run') & (df5_long['Metric'] == 'Gigabases')])
        + geom_hline(aes(yintercept=25), linetype='dashed', color='red', data=df5_long[(df5_long['Run_Sample'] == 'Run') & (df5_long['Metric'] == 'Thousand Reads')])
        #+ geom_hline(yintercept=0.4, linetype='dashed', color='blue')
        # + scale_y_discrete(name='Total run reads inc unmapped')
        )
    p3.save('Figure_S2.pdf', width=10, height=6)
    p3.save('Figure_S2.svg', width=10, height=6)

def plot_ratios(df, dfd):
    # plot AuG_trunc10 on x axis and Sample_reads_percent_of_refs on y axis with 
    df['Set']='Validation'
    df['batch positive amplification control']=df['reverse_transcription_control'].astype(str)
    df['batch PCR negative control']=df['amplification_control']
    dfd['Set']='Derivation'
    # change derivation set pass to boolean
    dfd['pass']=dfd['pass'].astype(bool)
    dfd['batch positive amplification control']='0.0'
    dfd['batch PCR negative control']=np.where((dfd['MS2_spike']==0) & (dfd['IC_virus_spike']==0), 1, 0)
    df=pd.concat([df, dfd])
    df['sample positive'] = df.groupby(['Run','barcode'])['gold_standard'].transform('max')
    # filter out samples with 0 reads
    df=df[df['sample num reads']>0]
    df['full_pass']=np.where((df['pass']==True) & (df['PCs_passed']== 1), True, False)
    # filter out unmapped pathogens
    df=df[~df['pathogen'].isin(['unmapped','unmapped_2'])]

    print('Unique values in pass column:', df['pass'].unique())
    print('Unique values in PCs_passed column:', df['PCs_passed'].unique())

    # create sample type of sample, batch positive amplification control, batch PCR negative control
    df['sample_type'] = np.where(df['sample positive'] == 1, 'Positive sample', 'Negative sample')
    df['sample_type'] = np.where(df['batch positive amplification control'] == '1.0', 'Batch positive amplification control', df['sample_type'])
    df['sample_type'] = np.where(df['batch PCR negative control'] == 1, 'Batch PCR negative control', df['sample_type'])


    # create classification results column
    conditions = [
        (df['gold_standard'] == 1) & (df['full_pass'] == True),  # True Positive
        (df['gold_standard'] == 1) & (df['full_pass'] == False),  # False Negative
        (df['gold_standard'] == 0) & (df['full_pass'] == False),  # True Negative
        (df['gold_standard'] == 0) & (df['full_pass'] == True)]  # False Positive
    choices = ['TP', 'FN', 'TN', 'FP']

    df["Test result"] = np.select(conditions, choices, default='Unknown')

    test_result_colors = {"TP": "#009E73", "FP": "#CC79A7", "TN": "#D55E00", "FN": "#0072B2"}

    # line_data of line y = 0.1 * x
    line_data = pd.DataFrame({
        'x': [0.001, 10],
        'y': [0.0001, 1]
    })

    # Create the fill_group column (use 'None' string instead of np.nan)
    df['fill_group'] = np.where(df['sample num reads'] >= 2, 
                                 df['Test result'], 
                                 'Failed 2 Reads')

    # Create the color column
    df['color_group'] = np.where(df['PCs_passed'], 
                                  df['Test result'], 
                                  'Failed Positive Controls')
    
    # Convert reverse_transcription_control to string and ensure only 'True'/'False' values
    df['batch positive amplification control'] = df['batch positive amplification control'].astype(str)

    df.to_csv('figure_03_input_data.csv', index=False)
    
    # add sample QC pass/fail based on PCs_passed, and 2 reads pass
    df['Sample_QC_pass'] = np.where((df['PCs_passed'] == True) & (df['sample num reads'] >= 2), 'Sample QC Pass', 'Sample QC Fail')

    # order Sample_QC_pass so that 'Sample QC Fail' comes after 'Sample QC Pass'
    df['Sample_QC_pass'] = pd.Categorical(df['Sample_QC_pass'], categories=['Sample QC Pass', 'Sample QC Fail'], ordered=True)

    # Update color mappings to include 'None'
    fill_colors = {**test_result_colors, 'Failed 2 Reads': 'white'}
    color_colors = {**test_result_colors, 'Failed Positive Controls': 'black'}

    # sample type shapes
    sample_type_shapes = {'Batch positive amplification control': '^', 
                                'Batch PCR negative control': 'v', 
                                'Positive sample': 'o', 
                                'Negative sample': 's'}

    # Plot using plotnine
    p = (ggplot(df, aes(x='AuG_trunc10', y='Sample_reads_percent_of_refs')) +
     geom_point(aes(fill='fill_group', 
                    color='color_group',
                    shape='sample_type'),
                size=2,
                stroke=0.5) +
     scale_shape_manual(values=sample_type_shapes, 
                        name='Sample type') +
     scale_color_manual(values=color_colors,
                        name='Outline color') +
     scale_fill_manual(values=fill_colors,
                       name='Fill color (only if 2 reads pass)') +
     geom_line(data=line_data, mapping=aes(x='x', y='y'), 
               color='red', inherit_aes=False) +
     geom_hline(yintercept=0.007, linetype='dashed', color='grey') +
     geom_vline(xintercept=0.003, linetype='dashed', color='grey') +
     scale_x_log10(limits=(0.001, 10)) +
     scale_y_log10(limits=(0.0001, 10)) +
     facet_grid('Sample_QC_pass ~ Set') +
     labs(x='Area under the genome truncated at depth 10 ',
          y='Number of reads mapping to the target reference as a percentage\nof the number of reads mapping to any reference in the run'))
    
    # Facet labels a/b/c/d
    label_df = pd.DataFrame({
        'Set': ['Derivation', 'Validation', 'Derivation', 'Validation'],
        'Sample_QC_pass': ['Sample QC Pass', 'Sample QC Pass', 'Sample QC Fail', 'Sample QC Fail'],
        'x': [0.002, 0.002, 0.002, 0.002],
        'y': [5, 5, 5, 5],
        'label': ['a', 'b', 'c', 'd'],
    })

    # order label_df to match facet order
    label_df['Sample_QC_pass'] = pd.Categorical(label_df['Sample_QC_pass'], categories=['Sample QC Pass', 'Sample QC Fail'], ordered=True)
    label_df = label_df.sort_values(by=['Sample_QC_pass', 'Set'])

    p = p + geom_text(
        data=label_df,
        mapping=aes(x='x', y='y', label='label'),
        inherit_aes=False,
        size=12,
        color='black',
        fontweight='bold'
    )

    p.save('Figure03.pdf', width=10, height=6)

def plot_alt_ratios(df, dfd):
    # plot AuG_trunc10 on x axis and Sample_reads_percent_of_refs on y axis with 
    df['Set']='Validation'
    df['batch positive amplification control']=df['reverse_transcription_control'].astype(str)
    df['batch PCR negative control']=df['amplification_control']
    # readjust pass based on Sample_reads_percent_of_type_run > 0.03
    df['pass'] = np.where((df['Sample_reads_percent_of_type_run'] > 0.03) 
                          & (df['OR pass']==True) & (df['2 reads pass']==True),
                           True, False)
    # refilter the Flu that aren't the top flu found
    df=df[~(df['FLU_A_POS']==True)&(df['TOP_FLU_A']==False)]
    #df['pass'] = np.where((df['FLU_A_POS']==True)&(df['TOP_FLU_A']==False), False, df['pass'])
    dfd['Set']='Derivation'
    # change derivation set pass to boolean
    dfd['pass']=dfd['pass'].astype(bool)
    dfd['batch positive amplification control']='0.0'
    dfd['batch PCR negative control']=np.where((dfd['MS2_spike']==0) & (dfd['IC_virus_spike']==0), 1, 0)
    df=pd.concat([df, dfd])
    df['sample positive'] = df.groupby(['Run','barcode'])['gold_standard'].transform('max')
    # filter out samples with 0 reads
    df=df[df['sample num reads']>0]
    df['full_pass']=np.where((df['pass']==True) & (df['PCs_passed']== 1), True, False)
    # filter out unmapped pathogens
    df=df[~df['pathogen'].isin(['unmapped','unmapped_2'])]

    # create sample type of sample, batch positive amplification control, batch PCR negative control
    df['sample_type'] = np.where(df['sample positive'] == 1, 'Positive sample', 'Negative sample')
    df['sample_type'] = np.where(df['batch positive amplification control'] == '1.0', 'Batch positive amplification control', df['sample_type'])
    df['sample_type'] = np.where(df['batch PCR negative control'] == 1, 'Batch PCR negative control', df['sample_type'])


    # create classification results column
    conditions = [
        (df['gold_standard'] == 1) & (df['full_pass'] == True),  # True Positive
        (df['gold_standard'] == 1) & (df['full_pass'] == False),  # False Negative
        (df['gold_standard'] == 0) & (df['full_pass'] == False),  # True Negative
        (df['gold_standard'] == 0) & (df['full_pass'] == True)]  # False Positive
    choices = ['TP', 'FN', 'TN', 'FP']

    df["Test result"] = np.select(conditions, choices, default='Unknown')

    test_result_colors = {"TP": "#009E73", "FP": "#CC79A7", "TN": "#D55E00", "FN": "#0072B2"}


    # Create the fill_group column (use 'None' string instead of np.nan)
    df['fill_group'] = np.where(df['sample num reads'] >= 2, 
                                 df['Test result'], 
                                 'Failed 2 Reads')

    # Create the color column
    df['color_group'] = np.where(df['PCs_passed'], 
                                  df['Test result'], 
                                  'Failed Positive Controls')
    
    # Convert reverse_transcription_control to string and ensure only 'True'/'False' values
    df['batch positive amplification control'] = df['batch positive amplification control'].astype(str)
    
    # add sample QC pass/fail based on PCs_passed, and 2 reads pass
    df['Sample_QC_pass'] = np.where((df['PCs_passed'] == True) & (df['sample num reads'] >= 2), 'Sample QC Pass', 'Sample QC Fail')

    # order Sample_QC_pass so that 'Sample QC Fail' comes after 'Sample QC Pass'
    df['Sample_QC_pass'] = pd.Categorical(df['Sample_QC_pass'], categories=['Sample QC Pass', 'Sample QC Fail'], ordered=True)

    df.to_csv('figure_S6_input_data.csv', index=False)

    # Update color mappings to include 'None'
    fill_colors = {**test_result_colors, 'Failed 2 Reads': 'white'}
    color_colors = {**test_result_colors, 'Failed Positive Controls': 'black'}

    # sample type shapes
    sample_type_shapes = {'Batch positive amplification control': '^', 
                                'Batch PCR negative control': 'v', 
                                'Positive sample': 'o', 
                                'Negative sample': 's'}

    # Plot using plotnine
    p = (ggplot(df, aes(x='AuG_trunc10', y='Sample_reads_percent_of_type_run')) +
     geom_point(aes(fill='fill_group', 
                    color='color_group',
                    shape='sample_type'),
                size=2,
                stroke=0.5) +
     scale_shape_manual(values=sample_type_shapes, 
                        name='Sample type') +
     scale_color_manual(values=color_colors,
                        name='Outline color') +
     scale_fill_manual(values=fill_colors,
                       name='Fill color (only if 2 reads pass)') +
     #geom_hline(yintercept=0.007, linetype='dashed', color='grey') +
     geom_vline(xintercept=0.003, linetype='dashed', color='grey') +
     geom_hline(yintercept=0.03,  color='red') +
     scale_x_log10(limits=(0.001, 10)) +
     scale_y_log10(limits=(0.0001, 100)) +
     facet_grid('Sample_QC_pass ~ Set') +
     labs(x='Area under the genome truncated at depth 10 ',
          y='Reads mapping to target reference as percentage of reads\nmapping to references of same type in run'))
    
    # Facet labels a/b/c/d
    label_df = pd.DataFrame({
        'Set': ['Derivation', 'Validation', 'Derivation', 'Validation'],
        'Sample_QC_pass': ['Sample QC Pass', 'Sample QC Pass', 'Sample QC Fail', 'Sample QC Fail'],
        'x': [0.0002, 0.0002, 0.0002, 0.0002],
        'y': [100, 100, 100, 100],
        'label': ['a', 'b', 'c', 'd'],
    })

    # order label_df to match facet order
    label_df['Sample_QC_pass'] = pd.Categorical(label_df['Sample_QC_pass'], categories=['Sample QC Pass', 'Sample QC Fail'], ordered=True)
    label_df = label_df.sort_values(by=['Sample_QC_pass', 'Set'])

    p = p + geom_text(
        data=label_df,
        mapping=aes(x='x', y='y', label='label'),
        inherit_aes=False,
        size=12,
        color='black',
        fontweight='bold'
    )

    p.save('FigureS6.pdf', width=10, height=6)

def lowess_with_confidence_bounds(
    x, y, eval_x, N=200, conf_interval=0.95, lowess_kw=None
):
    """
    Perform Lowess regression and determine a confidence interval by bootstrap resampling
    """
    # Lowess smoothing
    smoothed = lowess(exog=x, endog=y, xvals=eval_x, frac=0.1)  

    # Perform bootstrap resamplings of the data
    # and  evaluate the smoothing at a fixed set of points
    smoothed_values = np.empty((N, len(eval_x)))
    for i in range(N):
        sample = np.random.choice(len(x), len(x), replace=True)
        sampled_x = x[sample]
        sampled_y = y[sample]

        smoothed_values[i] = lowess(
            exog=sampled_x, endog=sampled_y, xvals=eval_x, frac=0.1
        )

    # Get the confidence interval
    sorted_values = np.sort(smoothed_values, axis=0)
    bound = int(N * (1 - conf_interval) / 2)
    bottom = sorted_values[bound - 1]
    top = sorted_values[-bound]

    return smoothed, bottom, top


def correlation_plot(dfd):
    # Correlation plot with spearman correlations

    # Load data
    dat2 = dfd[dfd['sample num reads']>0].copy()
    
    # Select metrics of interest
    cols = ['AuG_trunc10', 'sample num reads', 'Cov1_perc',
            'Sample_reads_percent_of_run', 'Sample_reads_percent_of_refs', 
            'Sample_reads_percent_of_type_run', 'Sample_reads_percent_of_type_sample',
            'mean_read_length', 'mean_aligned_length']

    alt_labels = {
        'AuG_trunc10': 'Area under genome\ntruncated at depth 10',
        'sample num reads': 'Number of reads\nmapping to target reference',
        'Cov1_perc': 'Percentage of target\nreference covered at depth 1',
        'Sample_reads_percent_of_run': 'Reads mapping to target\nreference as percentage of\ntotal reads in run',
        'Sample_reads_percent_of_refs': 'Reads mapping to target\nreference as percentage of all reads\nmapping to references in run',
        'Sample_reads_percent_of_type_run': 'Reads mapping to target\nreference as percentage of reads\nmapping to references of\nsame type in run',
        'Sample_reads_percent_of_type_sample': 'Reads mapping to target\nreference as percentageof reads\nmapping to references of\nsame type in sample',
        'mean_read_length': 'Mean read length (bp)',
        'mean_aligned_length': 'Mean aligned length (bp)'
    }
    
    dat = dat2[cols].copy()
    
    # Replace 0 with NaN
    dat = dat.replace(0, np.nan)
    
    n = len(cols)
    labels = cols
    
    # Step 1: Compute all correlations
    cor_data = []
    for i in range(n-1):
        for j in range(i+1, n):
            x = dat.iloc[:, j].dropna()
            y = dat.iloc[:, i].dropna()
            
            # Get common indices
            common_idx = x.index.intersection(y.index)
            x_clean = x.loc[common_idx]
            y_clean = y.loc[common_idx]
            
            if len(x_clean) > 0 and len(y_clean) > 0:
                correlation, p_value = spearmanr(x_clean, y_clean)
                
                print(f"{cols[j]}|{cols[i]}|{correlation}")
                
                cor_data.append({
                    'x_var': cols[j],
                    'y_var': cols[i],
                    'estimate': correlation,
                    'p_value': p_value,
                    'i': i,
                    'j': j
                })
    
    cor_df = pd.DataFrame(cor_data)
    
    # Categorize correlations
    cor_df['corr_group'] = pd.cut(cor_df['estimate'], 
                                   bins=[-np.inf, 0, np.inf],
                                   labels=['negative', 'positive'])
    
    # Categorize p-values
    cor_df['p_value_group'] = pd.cut(cor_df['p_value'],
                                      bins=[-np.inf, 0.001, 0.01, 0.05, 0.1, np.inf],
                                      labels=['***', '**', '*', '.', ''])
    
    # Create correlation labels
    cor_df['cor'] = cor_df.apply(lambda row: f"r: {row['estimate']:.2f}{row['p_value_group']}", axis=1)
    
    min_cor = cor_df['estimate'].abs().min()
    max_cor = cor_df['estimate'].abs().max()
    
    # Create figure with GridSpec
    fig = plt.figure(figsize=(14, 14))
    gs = GridSpec(n, n, figure=fig, hspace=0.05, wspace=0.05)
    
    # RdBu color palette (reversed)
    colors_list = ['#053061', '#2166ac', '#4393c3', '#92c5de', '#d1e5f0',
                   '#f7f7f7', '#fddbc7', '#f4a582', '#d6604d', '#b2182b', '#67001f']
    colors_list = colors_list[::-1]  # Reverse
    rdbu_cmap = LinearSegmentedColormap.from_list('RdBu', colors_list)
    
    # Plot each cell
    for i in range(n):
        for j in range(n):
            ax = fig.add_subplot(gs[i, j])
            #print (f'Plotting cell ({i}, {j})')
            if i == j:
                #print('Diagonal cell - grey tile')
                # Diagonal: grey tiles
                ax.add_patch(mpatches.Rectangle((0, 0), 1, 1, 
                                    facecolor='darkgrey', edgecolor='none'))
                ax.set_xlim(0, 1)
                ax.set_ylim(0, 1)
                #ax.axis('off')
                ax.tick_params(axis='x', labelbottom=False, labelleft=False, bottom=False, left=False)
                ax.tick_params(axis='y', labelbottom=False, labelleft=False, bottom=False, left=False)
                # add row and column labels
                if i == n-1:
                    ax.set_xlabel(alt_labels[labels[i]], fontsize=12, rotation=90, ha='right', va='top')
                if j == 0:
                    ax.set_ylabel(alt_labels[labels[j]], fontsize=12, rotation=360, ha='right', va='center')
                
            elif i < j:
                #print('Upper triangle cell - correlation circle')
                # Upper triangle: correlation circles
                cor_row = cor_df[(cor_df['i'] == i) & (cor_df['j'] == j)]
                
                if not cor_row.empty:
                    est = cor_row.iloc[0]['estimate']
                    cor_label = cor_row.iloc[0]['cor']
                    
                    # Normalize size
                    size = (abs(est) - min_cor) / (max_cor - min_cor) * 1400 + 100
                    size = abs(est) / 1 * 0.4  # Scale to [0, 0.4] for axes coordinates
                    
                    # Plot circle
                    circle = plt.Circle((0.5, 0.5), size, 
                                       color=rdbu_cmap((est + 1) / 2),
                                       transform=ax.transAxes)
                    ax.add_patch(circle)
                    
                    # Add text
                    ax.text(0.5, 0.1, cor_label, ha='center', va='center',
                           fontsize=10, backgroundcolor='white', 
                           transform=ax.transAxes)
                
                ax.set_xlim(0, 1)
                ax.set_ylim(0, 1)
                ax.axis('off')
                
            elif i > j:
                #print('Lower triangle cell - scatter plot')
                # Lower triangle: scatter plots
                x_data = dat.iloc[:, j].replace(0, np.nan).dropna()
                y_data = dat.iloc[:, i].replace(0, np.nan).dropna()
                
                common_idx = x_data.index.intersection(y_data.index)
                df_point = pd.DataFrame({
                    'x': x_data.loc[common_idx],
                    'y': y_data.loc[common_idx]
                }).dropna()
                
                if len(df_point) > 0:
                    ax.scatter(df_point['x'], df_point['y'], c='k', s=1, alpha=0.5)
                    
                    # Add loess smooth
                    if len(df_point) > 3:
                        sorted_df = df_point.sort_values('x')
                        smoothed = lowess(sorted_df['y'], sorted_df['x'], frac=0.2)
                        #eval_x = np.linspace(0, 4 * np.pi, 31)
                        #smoothed, bottom, top = lowess_with_confidence_bounds(sorted_df['y'], sorted_df['x'], eval_x, lowess_kw={"frac": 0.1})
                        x_smooth = smoothed[:, 0]
                        y_smooth = smoothed[:, 1]
                        # add smoothed line with confidence interval
                        ax.plot(x_smooth, y_smooth, 'r-', linewidth=1, alpha=0.7)
                        #ax.plot(eval_x, smoothed, c="k")
                        #ax.fill_between(eval_x, bottom, top, alpha=0.5, color="b")
                        # add black line
                        #ax.plot(x_smooth, y_smooth, 'k-', linewidth=1)
                        
                
                ax.tick_params(axis='x', labelbottom=False, labelleft=False, bottom=False, left=False)
                ax.tick_params(axis='y', labelbottom=False, labelleft=False, bottom=False, left=False)
                #ax.set_xticks([], [])
                #ax.set_yticks([], [])
                ax.grid(False)
                #ax.set_axis_off()

                # add row and column labels
                if i == n - 1:
                    ax.set_xlabel(alt_labels[labels[j]], fontsize=12, rotation=90, ha='center', va='top')
                if j == 0:
                    ax.set_ylabel(alt_labels[labels[i]], fontsize=12, rotation=360, ha='right', va='center')
            else:
                print('Unexpected cell position')

    #plt.tight_layout()
    #plt.savefig('correlation_plot.pdf', dpi=300, bbox_inches='tight')
    plt.savefig('Figure_S5.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    


if __name__ == "__main__":
    df=pd.read_csv(sys.argv[1])
    #total_samples=df.shape[0]
    #print(f'Total number of samples: {total_samples}')

    #biorifre_organisms=['Adenovirus', 'Coronavirus HKU1', 'Coronavirus NL63', 'Coronavirus OC43', 'Coronavirus 229E', 'Human Metapneumovirus', 'Influenza A', 'Influenza B', 'Parainfluenza 1', 'Parainfluenza 2', 'Parainfluenza 3', 'Parainfluenza 4', 'Respiratory Syncytial Virus', 'Rhinovirus/Enterovirus']
    #alinity_cephid_organisms=['Adenovirus', 'Coronavirus HKU1', 'Coronavirus NL63', 'Coronavirus OC43', 'Coronavirus 229E', 'Human Metapneumovirus', 'Influenza A', 'Influenza B', 'Parainfluenza Virus 1-4', 'Respiratory Syncytial Virus', 'Rhinovirus/Enterovirus']

    #plot_cts_vs_reads()

    plot_Ct_reads_plotnine(df)

    #plot_sensitivity_by_run_and_pathogen(df)

    #plot_sensitivity_by_pathogen(df)

    # read in derivation set
    if sys.argv[2]:
        dfd=pd.read_csv(sys.argv[2])
        if sys.argv[3]: # add in test type info
            dfd_tt=pd.read_csv(sys.argv[3])
            dfd=dfd.merge(dfd_tt, on=['Run', 'barcode'], how='left')
            dfd['test_type']=dfd['test_type_y']

        run_order={'Expt9A_SISPA_Katie':1,
            '010524_Expt9_SISPA_Daisy':2,   
            '010524_Expt9_SISPA_Kate':3,
            'expt10_03072024':4,
            'expt10A_17072024':5,
            'expt11_150824':6}

        dfd['Run_order'] = dfd['Run'].map(run_order)  # Map the run names to their order
        dfd.sort_values(by=['Run_order'], inplace=True)  # Sort by Run order and Batch
        dfd['Run']= 'Run ' + dfd['Run_order'].astype(str).str.zfill(2)
        dfd.to_csv('derivation_set_with_run_names_test_types.csv', index=False)
    else:
        sys.exit('Please provide derivation set csv as second argument to continue plotting')

    #plot_reads_bases(df, dfd)

    #plot_figure_S2(df, dfd)

    #plot_reads_bases_per_run(df, dfd)

    #combined, combined_test_counts = plot_gold_standard_pathogen_counts(df, dfd)

    #plot_pathogen_counts(combined, combined_test_counts)

    #plot_sensitivity_by_pathogen_plotnine(df, dfd)

    #plot_ratios(df, dfd)

    plot_alt_ratios(df, dfd)

    #correlation_plot(dfd)
