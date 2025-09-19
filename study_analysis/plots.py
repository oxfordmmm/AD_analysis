#!/usr/bin/env python3
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import sys
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
    geom_text,
    facet_wrap,
    facet_grid,
    theme,
    element_text,
    scale_fill_manual,
    scale_x_discrete,
    scale_y_discrete,
    scale_color_discrete,
    geom_bin_2d,
    stat_bin_2d,
    labs
    )

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


def get_pathogen_sensistivity_data(df_original):
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
    print(g1)

    # convert to long format for plotnine
    g1 = g1.reset_index().melt(id_vars=['pathogen'], var_name='Result', value_name='count')
    return g1

def plot_sensitivity_by_pathogen_plotnine(df_original, dfd):
    df_original = df_original.copy()  # Use the original DataFrame for grouping
    df_original= df_original[~((df_original['Run']=='AD_winter_study_220125') & (df_original['Batch']==2))]
    g1=get_pathogen_sensistivity_data(df_original)
    g1['Set']='Validation'
    print(g1)

    # repeat with derivation set
    g2=get_pathogen_sensistivity_data(dfd)
    g2['Set']='Derivation'
    g1=pd.concat([g1, g2])

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
    g1 = df.groupby(['Run', 'pathogen'])[['gold_standard']].sum().reset_index()


    g1 = g1.set_index(['Run', 'pathogen'])

    g1 = g1.unstack(level='pathogen')

    g1.columns = g1.columns.droplevel(0)  # Drop the top level of the MultiIndex
    g1 = g1.fillna(0)  # Fill NaN values with 0
    #g1.reset_index(inplace=True)  # Reset index to make it easier to plot
    # sort dataframe by first index in MultiIndex by run_order
    g1['Run_order'] = g1.index.get_level_values('Run').map(run_order)  # Map the run names to their order
    #g1['Run']= 'Run ' + g1['Run_order'].astype(str)  # Combine Run and Batch into a single column
    g1 = g1.sort_values(by='Run_order')  # Sort by Run order
    g1 = g1.drop(columns=['Run_order'])  # Drop the temporary Run_order column      

    # drop 'unmapped' column if it exists
    if 'unmapped' in g1.columns:
        g1 = g1.drop(columns=['unmapped'])
    # drop empty columns
    g1 = g1.loc[:, (g1 != 0).any(axis=0)]

    print(g1)
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

    # read in derivation set
    # if sys.argv[2]:
    #     dfd=pd.read_csv(sys.argv[2])
    # else:
    #     sys.exit('Please provide derivation set csv as second argument to continue plotting')

    # run_order={'Expt9A_SISPA_Katie':1,
    #         '010524_Expt9_SISPA_Daisy':2,   
    #         '010524_Expt9_SISPA_Kate':3,
    #         'expt10_03072024':4,
    #         'expt10A_17072024':5,
    #         'expt11_150824':6}

    # dfd['Run_order'] = dfd['Run'].map(run_order)  # Map the run names to their order
    # dfd.sort_values(by=['Run_order'], inplace=True)  # Sort by Run order and Batch
    # dfd['Run']= 'Run ' + dfd['Run_order'].astype(str).str.zfill(2)  # Combine Run and Batch into a single column, add leading zero to run number 
    # run_ordered = dfd['Run'].unique().tolist()
    # run_order = {name: i for i, name in enumerate(run_ordered)}  # Create a mapping of run names to their order

    g1 = dfd.groupby(['Run', 'pathogen'])[['gold_standard']].sum().reset_index()
    g1 = g1.set_index(['Run', 'pathogen'])
    g1 = g1.unstack(level='pathogen')
    g1.columns = g1.columns.droplevel(0)  # Drop the top level of the MultiIndex
    g1 = g1.fillna(0)  # Fill NaN values with 0

    # sort dataframe by first index in MultiIndex by run_order
    g1['Run_order'] = g1.index.get_level_values('Run').map(run_order)  # Map the run names to their order
    #g1['Run']= 'Run ' + g1['Run_order'].astype(str)  # Combine Run and Batch into a single column
    g1 = g1.sort_values(by='Run_order')  # Sort by Run order
    g1 = g1.drop(columns=['Run_order'])  # Drop the temporary Run_order column      

    # drop 'unmapped' column if it exists
    if 'unmapped' in g1.columns:
        g1 = g1.drop(columns=['unmapped'])

    # drop empty columns
    g1 = g1.loc[:, (g1 != 0).any(axis=0)]
    print(g1)

    # get columns from vg1
    cols=vgi.columns
    for col in cols:
        if col not in g1.columns:
            g1[col]=0
    g1=g1[cols]  # reorder columns to match vgi
    print(g1)

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
    g1['Dataset']='Derivation'
    #g1.reset_index(inplace=True,drop=True)
    #vgi.reset_index(inplace=True,drop=True)
    combined=pd.concat([g1, vgi])


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
    return combined


def plot_pathogen_counts(combined):
    # try with plotnine
    # convert combined to long format
    combined.reset_index(inplace=True)
    print(combined)
    combined_long = combined.melt(id_vars=['Run', 'Dataset'], var_name='Pathogen', value_name='Count')
    combined_long = combined_long[combined_long['Run'] != 0]  # Filter out zero Runs
    print(combined_long)
    run_ordered = combined_long['Run'].unique().tolist()
    combined_long.to_csv('gold_standard_pathogen_counts_by_run_and_batch_combined_long.csv', index=False)

    pathogen_colors = [
        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
        "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#a55194", "#393b79",
        "#637939", "#8c6d31", "#843c39", "#7b4173", "#ad494a", "#3182bd"
    ]

    p = (ggplot(combined_long, aes(x='Run', y='Count', fill='Pathogen'))
        + geom_bar(stat='identity', position='stack')
        + theme(axis_text_x=element_text(rotation=90, hjust=1))
        + scale_fill_manual(values=pathogen_colors, name='Pathogen')
        + facet_wrap('~Dataset', scales='free_x')
        )

    p.save('gold_standard_pathogen_counts_by_run_and_batch_combined_plotnine.pdf', width=10, height=6)
    p.save('Figure_1.pdf', width=10, height=6)

    # try plotting pathogen on the x axis and counts on the y axis, facet by Dataset
    p2 = (ggplot(combined_long, aes(x='Pathogen', y='Count', fill='Pathogen'))
        + geom_bar(stat='identity', position='stack')
        + theme(axis_text_x=element_text(rotation=90, hjust=1))
        + scale_fill_manual(values=pathogen_colors, name='Pathogen')
        + facet_wrap('~Dataset', scales='free_x')
        ) 
    p2.save('gold_standard_pathogen_counts_by_pathogen_and_batch_combined_plotnine.pdf', width=10, height=6)
    p2.save('Figure_1_alt.pdf', width=10, height=6)


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

if __name__ == "__main__":
    df=pd.read_csv(sys.argv[1])
    #total_samples=df.shape[0]
    #print(f'Total number of samples: {total_samples}')

    #biorifre_organisms=['Adenovirus', 'Coronavirus HKU1', 'Coronavirus NL63', 'Coronavirus OC43', 'Coronavirus 229E', 'Human Metapneumovirus', 'Influenza A', 'Influenza B', 'Parainfluenza 1', 'Parainfluenza 2', 'Parainfluenza 3', 'Parainfluenza 4', 'Respiratory Syncytial Virus', 'Rhinovirus/Enterovirus']
    #alinity_cephid_organisms=['Adenovirus', 'Coronavirus HKU1', 'Coronavirus NL63', 'Coronavirus OC43', 'Coronavirus 229E', 'Human Metapneumovirus', 'Influenza A', 'Influenza B', 'Parainfluenza Virus 1-4', 'Respiratory Syncytial Virus', 'Rhinovirus/Enterovirus']

    #plot_cts_vs_reads()

    #plot_sensitivity_by_run_and_pathogen(df)

    #plot_sensitivity_by_pathogen(df)

    # read in derivation set
    if sys.argv[2]:
        dfd=pd.read_csv(sys.argv[2])
        run_order={'Expt9A_SISPA_Katie':1,
            '010524_Expt9_SISPA_Daisy':2,   
            '010524_Expt9_SISPA_Kate':3,
            'expt10_03072024':4,
            'expt10A_17072024':5,
            'expt11_150824':6}

        dfd['Run_order'] = dfd['Run'].map(run_order)  # Map the run names to their order
        dfd.sort_values(by=['Run_order'], inplace=True)  # Sort by Run order and Batch
        dfd['Run']= 'Run ' + dfd['Run_order'].astype(str).str.zfill(2)
    else:
        sys.exit('Please provide derivation set csv as second argument to continue plotting')

    #plot_reads_bases(df, dfd)

    #plot_figure_S2(df, dfd)

    #plot_reads_bases_per_run(df, dfd)

    combined=plot_gold_standard_pathogen_counts(df, dfd)

    plot_pathogen_counts(combined)


    plot_sensitivity_by_pathogen_plotnine(df, dfd)
