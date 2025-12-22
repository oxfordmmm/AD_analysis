#!/usr/bin/env python3
import pandas as pd
import numpy as np
#import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import upsetplot as up
import sys
import argparse
from mizani.formatters import label_number, label_comma
from scipy.stats import spearmanr
from scipy.interpolate import make_interp_spline
from statsmodels.nonparametric.smoothers_lowess import lowess
from statsmodels.stats.proportion import proportion_confint
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
    geom_errorbar,
    geom_histogram,
    facet_wrap,
    facet_grid,
    theme,
    theme_minimal,
    element_text,
    scale_fill_manual,
    scale_color_brewer,
    scale_x_discrete,
    scale_y_discrete,
    scale_y_log10,
    scale_x_reverse,
    scale_x_continuous,
    scale_y_continuous,
    position_dodge,
    position_jitter,
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
    geom_abline,
    annotate,
    element_blank)

ward_generalisation={'C-DC Renal':	'Other_inpatient_adult',
'C-WD CICU':	'Critical_care',
'C-WD EPCTU':	'Haematology_oncology',
'C-WD Haem':	'Haematology_oncology',
'C-WD OncHaemAmb':	'Haematology_oncology',
'C-WD OncHTriage':	'Haematology_oncology',
'C-WD Oncology':	'Haematology_oncology',
'C-WD Renal':	'Other_inpatient_adult',
'C-WD Transplant':	'Other_inpatient_adult',
'C-WD Urology':	'Other_inpatient_adult',
'H-WD Childrens':	'Other_inpatient_paeds',
'H-WD EAU':	'Other_inpatient_adult',
'H-WD Juniper':	'Other_inpatient_adult',
'H-WD Laburnum':	'Other_inpatient_adult',
'H-WD Oak HCU':	'Other_inpatient_adult',
'H-WD Rowan AU':	'Ambulatory_and_ED',
'H-WD Trauma F':	'Other_inpatient_adult',
'J-WD 5A':	'Other_inpatient_adult',
'J-WD 5C-5D':	'Other_inpatient_adult',
'J-WD Bell-Dray':	'Other_inpatient_paeds',
'J-WD Child CDU':	'Ambulatory_and_ED',
'J-WD CMU-A':	'Other_inpatient_adult',
'J-WD CMU-B':	'Other_inpatient_adult',
'J-WD CMU-C':	'Other_inpatient_adult',
'J-WD CritCareL2':	'Critical_care',
'J-WD EAU':	'Other_inpatient_adult',
'J-WD ED ORU':	'Ambulatory_and_ED',
'J-WD KamDaycare':	'Haematology_oncology',
'J-WD Kamrans':	'Haematology_oncology',
'J-WD L4 AAU':	'Ambulatory_and_ED',
'J-WD Maty AU':	'Other_inpatient_adult',
'J-WD Maty L6':	'Other_inpatient_adult',
'J-WD Melanies':	'Other_inpatient_paeds',
'J-WD Neuro Red':	'Other_inpatient_adult',
'J-WD NeuroPurpl':	'Other_inpatient_adult',
'J-WD Newborn IC':	'Critical_care',
'J-WD PCC':	'Critical_care',
'J-WD Robins':	'Other_inpatient_paeds',
'J-WD SEU E':	'Other_inpatient_adult',
'J-WD Toms':	'Other_inpatient_paeds',
'J-WD Trauma 2A':	'Other_inpatient_adult',
'J-WD Trauma 3A':	'Other_inpatient_adult',
'J-WD Ward 5E-5F':	'Other_inpatient_adult',
'O-WD AHatHome':	'Ambulatory_and_ED',
'O-WD CHIM':	'Ambulatory_and_ED',
'V-WD HAHCentral':	'Ambulatory_and_ED',
'V-WD HAHNorth':	'Ambulatory_and_ED'}

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
        'CEPHEID': 'Cepheid\nXpert Xpress\nCoV-2/Flu/RSV assay\n(4 pathogen targets)'
    }

test_type_normalisation = {
        'Biofire': 'BIOFIRE',
        'Biofire ': 'BIOFIRE',
        'BIOFIRE': 'BIOFIRE',
        'BIOFIRE ': 'BIOFIRE',
        'Alinity': 'ALINITY',
        'ALINITY': 'ALINITY',
        'ALINTY': 'ALINITY',
        'Cepheid': 'CEPHEID',
        'CEPHEID': 'CEPHEID'
    }

dna_pathogens=['Adenovirus', 'Bordetella_pertussis', 'Bordetella_parapertussis', 'Mycoplasma_pneumoniae', 'Chlamydia_pneumoniae',
               'Bordetella parapertussis',
               'Bordetella pertussis',
               'Chlamydia pneumoniae',
               'Mycoplasma pneumoniae']

flu_A_pathogens=['Influenza A/H1', 'Influenza A/H3', 'Influenza A/H1N1pdm09','Influenza A/H1-2009']

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
    return df2

def readjust_pass(df):
    df=df.copy()
    # adjust pass column to be 'True' if sample num reads > 50 else 'False'
    df['pass'] = np.where(df['Sample_reads_percent_of_type_run']>=0.037, True, False)
    df['pass'] = np.where(df['sample num reads']<2, False, df['pass'])
    df['pass'] = np.where(df['pathogen']=='unmapped', False, df['pass'])

    # readjust pass for flu that are not TOP_FLU_A
    df['pass'] = np.where((df['pathogen'].isin(flu_A_pathogens)) & (df['TOP_FLU_A']==False), False, df['pass'])
    return df

    
def remove_failed_runs(df):
    # count number of samples that failed negative controls
    # spike negs
    # This is also called the batch negatice control in the flow diagram
    df['pass'] = np.where(df['pathogen']=='unmapped', False, df['pass'])
    df_negs=df[(df['amplification_control']==1) & (df['pathogen']!='unmapped')]
    df_negs_pass=df_negs[(df_negs['pass']==True) | (df_negs['PCs_passed']==1)]
    if len(df_negs_pass)>0:
        unique_batches=len(df_negs_pass.drop_duplicates(['Run','Batch'], keep='first').index)
        print(f'Batches that failed negative controls but passed run/sample controls:xB {unique_batches}')
        failedruns=list(df_negs_pass['Run'].unique())
        print(failedruns)
        df2=df[~df['Run'].isin(failedruns)]
        df2=df2.copy()
    else:
        df2 = df.copy()

    # count number of samples that failed batch ampification negative controls
    df_RC_control=df2[(df2['reverse_transcription_control']==1) & (df2['IC_virus_spike']==1)]
    #print(df_RC_control[['Run','Batch','seq_name','barcode']].drop_duplicates())
    test=df_RC_control[(df_RC_control['Run']=='AD_winter_study_130825_rpt050825')& (df_RC_control['barcode']==89)]
    print(test['pass'].unique())
    # orthoreovirus passed	zika passed	MS2 passed	murine_respirovirus passed
    df_RC_control_PCFAIL=df_RC_control[(df_RC_control['orthoreovirus passed']==0) & (df_RC_control['zika passed']==0) & (df_RC_control['murine_respirovirus passed']==0) \
                                    | (df_RC_control['MS2 passed']==1) \
                                    | (df_RC_control['pass']==True) ]
    print(df_RC_control_PCFAIL[['Run','Batch','seq_name','pathogen','orthoreovirus passed','zika passed','MS2 passed','murine_respirovirus passed','pass']])
    if len(df_RC_control_PCFAIL)>0:
        df2=df2[df2['pathogen']!='unmapped']
        failed_batches=list(df_RC_control_PCFAIL.drop_duplicates(['Run','Batch'])[['Run','Batch']].itertuples(index=False, name=None))
        print(f'Batches that failed reverse transcription controls but passed run/sample controls: {len(failed_batches)}')
        print(failed_batches)
        df2=df2[~df2[['Run','Batch']].apply(tuple, axis=1).isin(failed_batches)]

    df2=df2[df2['total sample reads']>25_000]
    df2=df2[df2['total run bases']>400_000]
    return df2

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
    return g1, df_original

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
    return g1, df_original

def get_CI(df, set=None):
    # calculate 95% confidence intervals for sensitivity
    df = df.copy()
    df['sensitivity'] = (df['Pathogen detected\nby sequencing'] / (df['Pathogen detected\nby sequencing'] + df['Not detected\nby sequencing']))*100
    # using normal approximation method
    methods=['beta','normal','agresti_coull','wilson','jeffreys','binom_test']
    for method in methods:
        ci_lower, ci_upper = proportion_confint(df['Pathogen detected\nby sequencing'], df['Pathogen detected\nby sequencing'] + df['Not detected\nby sequencing'], 
                                                alpha=0.05, method=method)
        df[f'CI lower {method}'] = ci_lower * 100
        df[f'CI upper {method}'] = ci_upper * 100
    #ci_lower, ci_upper = proportion_confint(df['Pathogen detected'], df['Pathogen detected'] + df['Failed criteria'], 
    #                                        alpha=0.05, method='beta')
    #df['CI lower'] = ci_lower
    #df['CI upper'] = ci_upper
    #print(df)
    if set:
        df.to_csv(f'pathogen_sensitivity_CI_{set}.csv')
    else:
        df.to_csv('pathogen_sensitivity_CI.csv')
    return df

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
    g2 = g2.rename(columns={'TP': 'Pathogen detected\nby sequencing'})
    g0 = g0.rename(columns={'gold_standard': 'all_gold_standard_count'})
    g1 = g1.merge(g2, on=groups)
    g1 = g1.merge(g0, on=groups)
    g1['Not detected\nby sequencing'] = g1['gold_standard_count'] - g1['Pathogen detected\nby sequencing']
    g1['Failed QC checks'] = g1['all_gold_standard_count'] - (g1['Pathogen detected\nby sequencing'] + g1['Not detected\nby sequencing'])

    #print(g1)
    g1.drop(columns=['gold_standard_count','all_gold_standard_count'], inplace=True)
    #print(g1)
    #print(g1.index)

    # remove 'unmapped' row if it exists
    g1= g1[g1['pathogen'] != 'unmapped']

    g1 = g1.set_index(groups)
    #print(g1)

    # Get sensitivity and CI and save to csv
    sensitivity=get_CI(g1, set=set)
    sensitivity['Set']=set

    # convert to long format for plotnine
    g1 = g1.reset_index().melt(id_vars=['pathogen'], var_name='Result', value_name='count')
    return g1, sensitivity

def sens_pretty(row):
    detections=row['Pathogen detected\nby sequencing']
    s=f'{row["sensitivity"]:.0f}%[{row["ci_lower"]:.0f}%-{row["ci_upper"]:.0f}%] ({detections:.0f}/{row["gold_standard_count"]:.0f})'
    return s

def specs_pretty(row):
    s=f'{row["specificity"]:.0f}%[{row["spec_ci_lower"]:.0f}%-{row["spec_ci_upper"]:.0f}%] ({row["True Negatives"]:.0f}/{row["Negative Tests"]:.0f})'
    return s

def sens_pretty_per_sample(row):
    s=f'{row["Sensitivity (PS)"]:.0f}%[{row["ci_lower_PS"]:.0f}%-{row["ci_upper_PS"]:.0f}%] ({row["Sample True Positives"]:.0f}/{row["Positive samples"]:.0f})'
    return s

def specs_pretty_per_sample(row):
    TN=row['Number Samples']-row['Sample False Positives']
    s=f'{row["Specificity (PS)"]:.0f}%[{row["spec_ci_lower_PS"]:.0f}%-{row["spec_ci_upper_PS"]:.0f}%] ({TN:.0f}/{row["Number Samples"]:.0f})'
    return s

def get_set_sensistivity_data(df_original, set=None):
    groups=['set']
    df = df_original.copy()  # Use the original DataFrame for grouping
    df = df[df['test_type'].isin(test_type_mapping.keys())]
    # remove 'unmapped' row if it exists
    df = df[df['pathogen'] != 'unmapped']
    g0=df.groupby(groups)[['gold_standard']].sum().reset_index()
    df = df[df['PCs_passed'] == 1]  # Filter for PCs_passed == 1
    g1=df.groupby(groups)[['gold_standard']].sum().reset_index()
    df['pass'] = np.where(df['pass']=='0', 0, df['pass'])  # Convert 'pass' to numeric
    df['pass'] = np.where(df['pass']=='False', 0, df['pass'])
    df['pass'] = np.where(df['pass']=='True', 1, df['pass'])
    df['pass'] = df['pass'].astype(int)  # Ensure 'pass' is integer
    df['TP']= np.where((df['pass']==1) & (df['gold_standard']==1), 1, 0)  # Create a new column 'TP' for true positives
    df['FP']= np.where((df['pass']==1) & (df['gold_standard']==0), 1, 0)  # Create a new column 'FP' for false positives
    df['TN']= np.where((df['pass']==0) & (df['gold_standard']==0), 1, 0)  # Create a new column 'TN' for true negatives
    df['FN']= np.where((df['pass']==0) & (df['gold_standard']==1), 1, 0)  # Create a new column 'FN' for false negatives
    df['negative test']= np.where(df['gold_standard']==0, 1, 0)  # Create a new column 'negative test' for false negatives
    df['sample TP']=df.groupby(['set','Run','barcode'])['TP'].transform(max)
    df['sample FP']=df.groupby(['set','Run','barcode'])['FP'].transform(max)
    # save input data
    df.to_csv('table_3_input_data.csv', index=False)

    df2=df.drop_duplicates(subset=['set','Run','barcode'])

    # groupby various metrics
    g2=df.groupby(groups)[['TP']].sum().reset_index()
    g3=df.groupby(groups)[['FP']].sum().reset_index()
    g4=df.groupby(groups)[['TN']].sum().reset_index()
    g5=df.groupby(groups)[['negative test']].sum().reset_index()
    g6=df.groupby(groups)[['FN']].sum().reset_index()
    g7=df2.groupby(groups)[['sample TP','sample FP','sample_positive']].sum().reset_index()
    g8=df2.groupby(groups)[['sample_name']].nunique().reset_index()
    # Merge the two dataframes on 'Run' and 'batch'
    g1 = g1.rename(columns={'gold_standard': 'gold_standard_count'})
    g2 = g2.rename(columns={'TP': 'Pathogen detected\nby sequencing'})
    g0 = g0.rename(columns={'gold_standard': 'all_gold_standard_count'})
    g3 = g3.rename(columns={'FP': 'False-positive detections'})
    g4 = g4.rename(columns={'TN': 'True Negatives'})
    g5 = g5.rename(columns={'negative test': 'Negative Tests'})
    g6 = g6.rename(columns={'FN': 'False Negatives'})
    g7 = g7.rename(columns={'sample TP': 'Sample True Positives', 'sample FP': 'Sample False Positives', 'sample_positive': 'Positive samples'})
    g8 = g8.rename(columns={'sample_name': 'Number Samples'})
    g1 = g1.merge(g2, on=groups)
    g1 = g1.merge(g0, on=groups)
    g1 = g1.merge(g3, on=groups)
    g1 = g1.merge(g4, on=groups)
    g1 = g1.merge(g5, on=groups)
    g1 = g1.merge(g6, on=groups)
    g1 = g1.merge(g7, on=groups)
    g1 = g1.merge(g8, on=groups)
    g1['Not detected\nby sequencing'] = g1['gold_standard_count'] - g1['Pathogen detected\nby sequencing']
    g1['Failed QC checks'] = g1['all_gold_standard_count'] - (g1['Pathogen detected\nby sequencing'] + g1['Not detected\nby sequencing'])

    #g1.drop(columns=['gold_standard_count','all_gold_standard_count'], inplace=True)
    g1 = g1.set_index(groups)

    g1['sensitivity'] = (g1['Pathogen detected\nby sequencing'] / (g1['gold_standard_count'] ))*100
    g1['specificity'] = (g1['True Negatives'] / (g1['Negative Tests']))*100
    ci_lower, ci_upper = proportion_confint(g1['Pathogen detected\nby sequencing'], g1['Pathogen detected\nby sequencing'] + g1['Not detected\nby sequencing'], 
                                                alpha=0.05, method='beta')

    g1['ci_lower'] = ci_lower * 100
    g1['ci_upper'] = ci_upper * 100

    specs_ci_lower, specs_ci_upper = proportion_confint(g1['True Negatives'], g1['Negative Tests'], 
                                                alpha=0.05, method='beta')
    g1['spec_ci_lower'] = specs_ci_lower * 100
    g1['spec_ci_upper'] = specs_ci_upper * 100

    g1['Sensitivity (per target)'] = g1.apply(sens_pretty, axis=1)
    g1['Specificity (per target)'] = g1.apply(specs_pretty, axis=1)

    g1['Sensitivity (PS)'] = (g1['Sample True Positives'] / g1['Positive samples']) * 100
    g1['Specificity (PS)'] = ((g1['Number Samples'] - g1['Sample False Positives']) / g1['Number Samples']) * 100

    ci_lower_PS, ci_upper_PS = proportion_confint(g1['Sample True Positives'], g1['Positive samples'], 
                                                alpha=0.05, method='beta')
    g1['ci_lower_PS'] = ci_lower_PS * 100
    g1['ci_upper_PS'] = ci_upper_PS * 100

    specs_ci_lower_PS, specs_ci_upper_PS = proportion_confint(g1['Number Samples'] - g1['Sample False Positives'], g1['Number Samples'], 
                                                alpha=0.05, method='beta')
    g1['spec_ci_lower_PS'] = specs_ci_lower_PS * 100
    g1['spec_ci_upper_PS'] = specs_ci_upper_PS * 100

    g1['Sensitivity (per sample)'] = g1.apply(sens_pretty_per_sample, axis=1)
    g1['Specificity (per sample)'] = g1.apply(specs_pretty_per_sample, axis=1)

    return g1

def plot_sensitivity_by_pathogen_plotnine(df, dfd):
    df_original = df.copy()  # Use the original DataFrame for grouping
    #df_original= df_original[~((df_original['Run']=='AD_winter_study_220125') & (df_original['Batch']==2))]
    df_original = remove_failed_runs(df_original)
    g1,s1=get_pathogen_sensistivity_data(df_original, set='Validation - main criteria')
    g1['Set']='Validation - main criteria'
    #print(g1)
    # repeat with readjusted pass
    df_adjusted = readjust_pass(df)
    df_adjusted = remove_failed_runs(df_adjusted)
    g3,s3=get_pathogen_sensistivity_data(df_adjusted, set='Validation - alternative criteria')
    g3['Set']='Validation - alternative criteria'
    # repeat with derivation set
    dfd.to_csv('dfd.csv', index=False)
    g2,s2=get_pathogen_sensistivity_data(dfd, set='Derivation - main criteria')
    g2['Set']='Derivation - main criteria'
    g1=pd.concat([g1, g2, g3])
    s1=pd.concat([s1, s3])
    s1.reset_index(inplace=True)
    # remove _ from pathogen names
    g1['pathogen'] = g1['pathogen'].str.replace('_', ' ')
    s1['pathogen'] = s1['pathogen'].str.replace('_', ' ')
    g1.to_csv('pathogen_sensitivity_data.csv', index=False)
    # Red, blue, green colourbline safe
    result_colors = ["#393b79", '#ff7f0e', '#2ca02c']
    #result_colors = ['#1f77b4', '#ff7f0e', '#2ca02c'] 
    # stacked bar plot of pathogen, with pass, gold_standard - pass, all_gold_standard - gold_standard as different colours

    # order facet by Derivation - main criteria, Validation - main criteria, Validation - alternative criteria
    g1['Set'] = pd.Categorical(g1['Set'], 
                               categories=['Derivation - main criteria', 'Validation - main criteria', 'Validation - alternative criteria'], 
                               ordered=True)
    p=(ggplot(g1, aes(x='pathogen', y='count', fill='Result')) +
       geom_bar(stat='identity', position='stack') +
       scale_fill_manual(values=result_colors, name='Pathogen') +
       facet_wrap('~Set',ncol=1) +
       # put legend in the top right corner inside plot
       theme(legend_position=(0.05, 0.9),
             legend_justification='left',
             legend_title=element_blank(),
             legend_key_size=8) +
       labs(x='',
            y='Count',
            fill='Result',
            tag='A') +
       theme(axis_text_x=element_text(rotation=90, hjust=1, size=6),
             legend_text=element_text(size=6)
       )
    )

    # Facet labels a/b/c
    label_df = pd.DataFrame({
        'Set': ['Derivation - main criteria', 'Validation - main criteria', 
                'Validation - alternative criteria'],
        #'Sample_QC_pass': ['Sample QC Pass', 'Sample QC Pass', 'Sample QC Fail', 'Sample QC Fail'],
        'x': ['Adenovirus', 'Adenovirus', 'Adenovirus' ],
        'y': [57, 57, 57],
        'label': ['A', 'B', 'C'],
    })

    # order label_df to match facet order
    label_df['Set'] = pd.Categorical(label_df['Set'], 
                                     categories=['Derivation - main criteria', 'Validation - main criteria', 'Validation - alternative criteria'], 
                                     ordered=True)
    label_df = label_df.sort_values(by=['Set'])

    #p = p + geom_text(
    #    data=label_df,
    #    mapping=aes(x='x', y='y', label='label'),
    #    inherit_aes=False,
    #    size=8,
    #    color='black',
    #    fontweight='bold'
    #)
    
    #p.save('Figure_S7ab.pdf')
    # plot sensitivity with error bars
    s1['Total tests'] = s1['Pathogen detected\nby sequencing'] + s1['Not detected\nby sequencing']
    s1['>20 tests'] = np.where(s1['Total tests']>=20, '20 or more tests', 'Less than 20 tests')
    s1['Sensitivity (%)'] = np.where(s1['>20 tests']=='20 or more tests', s1['sensitivity'], np.nan)
    CI=['CI lower beta', 'CI upper beta']
    for ci in CI:
        s1[ci] = np.where(s1['>20 tests']=='20 or more tests', s1[ci], np.nan)
    # test colours for Total tests black and light gray
    #test_colours= ["#393b79", "#7f7f7f"]
    s1.dropna(subset=['Sensitivity (%)'], inplace=True)
    # Order Set by Derivation - main criteria, Validation - main criteria, Validation - alternative criteria
    s1['Set'] = pd.Categorical(s1['Set'], categories=['Validation - main criteria', 'Validation - alternative criteria'], ordered=True)
    p2=(ggplot(s1, aes(x='pathogen', y='Sensitivity (%)', color='Set')) +
        geom_point(size=3, position=position_dodge(width=0.5)) +
        geom_errorbar(aes(ymin='CI lower beta', ymax='CI upper beta'), width=0.2, position=position_dodge(width=0.5)) +
        scale_color_manual(values=["#1f77b4", '#ff7f0e', '#2ca02c']) +
        labs(x='',
             y='Sensitivity (%)',
             tag='B') +
        theme(axis_text_x=element_text(rotation=90, hjust=1, size=6),
              legend_text=element_text(size=5),
              legend_title=element_blank(),
              legend_position='top'
        )
    )
    #p2.save('Figure_S7cd.pdf')
    

    px = p | (p2)
    #px + theme(figure_size=(8, 3))
    px.save('supplemental/Figure_S6_counts_sensitivity.pdf')
    return g1, s1

def table_3_sensitivity_specificity(df):
    """Calculte sensitivity and specificity for each data validation set with different criteria"""
    # only use clinical samples
    df_original = df.copy()  # Use the original DataFrame for grouping
    df_original['pass'] = np.where((df_original['pathogen']!='unmapped') & (df_original['2 reads pass']==True) & (df_original['OR pass']==True) & (df_original['AND ratio pass']==True), 1, 0)
    df_original['pass'] = np.where((df_original['pathogen'].isin(flu_A_pathogens)) & (df_original['TOP_FLU_A']==False), 0, df_original['pass'])
    df_original = remove_failed_runs(df_original)
    df_original['set'] = 'Validation - main criteria'
    
    # repeat with > 0 sample_num_reads
    df_positive_reads = df_original[df_original['sample num reads'] > 0]
    df_positive_reads['set'] = 'Validation - main criteria > 0 reads'

    # repeat with non null CT values
    df_original_ct = df_original[df_original['qpcr_ct_value'].notnull()]
    df_original_ct['set'] = 'Validation - main criteria non null CT'

    # repeat with Ct < 35
    df_original_ct35 = df_original[df_original['qpcr_ct_value'] < 35]
    df_original_ct35['set'] = 'Validation - main criteria CT < 35'
    
    # repeat with readjusted pass
    df_adjusted = readjust_pass(df)
    df_adjusted = remove_failed_runs(df_adjusted)
    df_adjusted['set'] = 'Validation - alternative criteria'

    # repeat with non null CT values
    df_adjusted_ct = df_adjusted[df_adjusted['qpcr_ct_value'].notnull()]
    df_adjusted_ct['set'] = 'Validation - alternative criteria non null CT'

    # repeat with Ct < 35
    df_adjusted_ct35 = df_adjusted[df_adjusted['qpcr_ct_value'] < 35]
    df_adjusted_ct35['set'] = 'Validation - alternative criteria CT < 35'
    
    # repeat with > 0 sample_num_reads
    df_adjusted_positive_reads = df_adjusted[df_adjusted['sample num reads'] > 0]
    df_adjusted_positive_reads['set'] = 'Validation - alternative criteria > 0 reads'

    # repeat with only RNA pathogens
    df_adjusted_rna = df_adjusted[~df_adjusted['pathogen'].isin(dna_pathogens)]
    df_adjusted_rna['set'] = 'Validation - alternative criteria RNA only'

    # repeat with non null CT values
    df_adjusted_rna_ct = df_adjusted_rna[df_adjusted_rna['qpcr_ct_value'].notnull()]
    df_adjusted_rna_ct['set'] = 'Validation - alternative criteria RNA only non null CT'

    # repeat with Ct < 35
    df_adjusted_rna_ct35 = df_adjusted_rna[df_adjusted_rna['qpcr_ct_value'] < 35]
    df_adjusted_rna_ct35['set'] = 'Validation - alternative criteria RNA only CT < 35'

    # repeat with only RNA pathogens > 0 reads
    df_adjusted_rna_positive_reads = df_adjusted_rna[df_adjusted_rna['sample num reads'] > 0]
    df_adjusted_rna_positive_reads['set'] = 'Validation - alternative criteria RNA only > 0 reads'

    # repeat with only RNA pathogens without HRE
    df_adjusted_rna_no_hre = df_adjusted_rna[~df_adjusted_rna['pathogen'].isin(['Rhinovirus/enterovirus'])]
    df_adjusted_rna_no_hre['set'] = 'Validation - alternative criteria RNA only without HRE'

    # repeat with non null CT values
    df_adjusted_rna_no_hre_ct = df_adjusted_rna_no_hre[df_adjusted_rna_no_hre['qpcr_ct_value'].notnull()]
    df_adjusted_rna_no_hre_ct['set'] = 'Validation - alternative criteria RNA only without HRE non null CT'

    # repeat with Ct < 35
    df_adjusted_rna_no_hre_ct35 = df_adjusted_rna_no_hre[df_adjusted_rna_no_hre['qpcr_ct_value'] < 35]
    df_adjusted_rna_no_hre_ct35['set'] = 'Validation - alternative criteria RNA only without HRE CT < 35'

    # repeat with only RNA pathogenes without HRE > 0 reads
    df_adjusted_rna_no_hre_positive_reads = df_adjusted_rna_no_hre[df_adjusted_rna_no_hre['sample num reads'] > 0]
    df_adjusted_rna_no_hre_positive_reads['set'] = 'Validation - alternative criteria RNA only without HRE > 0 reads'

    # combine datasets
    df2=pd.concat([df_original, df_positive_reads, df_original_ct, df_original_ct35, 
                   df_adjusted, df_adjusted_ct, df_adjusted_positive_reads, df_adjusted_ct35,
                   df_adjusted_rna, df_adjusted_rna_ct, df_adjusted_rna_positive_reads, 
                   df_adjusted_rna_no_hre, df_adjusted_rna_no_hre_ct, df_adjusted_rna_no_hre_positive_reads,
                   df_adjusted_rna_ct35, df_adjusted_rna_no_hre_ct35])
    
    g1=get_set_sensistivity_data(df2)
    g1.sort_values(by=['sensitivity'], ascending=True, inplace=True)

    g1.to_csv('table_3_sensitivity_specificity_summary.csv')

    # restructure for table 3
    g1=g1.reset_index()
    main_sets=['Validation - main criteria','Validation - alternative criteria', 'Validation - alternative criteria RNA only', 'Validation - alternative criteria RNA only without HRE']
    df3=g1[g1['set'].isin(main_sets)]
    
    # Get selection of colimns for >0 reads
    no_reads_sets=['Validation - main criteria > 0 reads', 'Validation - alternative criteria > 0 reads', 'Validation - alternative criteria RNA only > 0 reads', 'Validation - alternative criteria RNA only without HRE > 0 reads']
    df3_no_reads=g1[g1['set'].isin(no_reads_sets)]
    df3_no_reads=df3_no_reads[['set', 'Sensitivity (per target)', 'gold_standard_count']]
    df3_no_reads.rename(columns={'Sensitivity (per target)': 'Sensitivity in true-positive targets with >0 reads', 'gold_standard_count': 'gold_standard_count_>0'}, inplace=True)

    set_mapper={'Validation - main criteria > 0 reads': 'Validation - main criteria',
                'Validation - alternative criteria > 0 reads': 'Validation - alternative criteria',
                'Validation - alternative criteria RNA only > 0 reads': 'Validation - alternative criteria RNA only',
                'Validation - alternative criteria RNA only without HRE > 0 reads': 'Validation - alternative criteria RNA only without HRE'}
    
    df3_no_reads['set'] = df3_no_reads['set'].map(set_mapper)

    df3=df3.merge(df3_no_reads, on='set')
    df3['Pathogen targets with 0 mapped reads']=df3['gold_standard_count'] - df3['gold_standard_count_>0']

    # get selection of columns for CT values
    ct_sets=['Validation - main criteria non null CT', 'Validation - alternative criteria non null CT', 'Validation - alternative criteria RNA only non null CT', 'Validation - alternative criteria RNA only without HRE non null CT']
    df3_ct=g1[g1['set'].isin(ct_sets)]
    df3_ct=df3_ct[['set', 'Sensitivity (per target)', 'gold_standard_count']]
    df3_ct.rename(columns={'Sensitivity (per target)': 'Sensitivity in true-positive targets with non null CT', 'gold_standard_count': 'gold_standard_count_non_null_CT'}, inplace=True)
    set_mapper={'Validation - main criteria non null CT': 'Validation - main criteria',
                'Validation - alternative criteria non null CT': 'Validation - alternative criteria',
                'Validation - alternative criteria RNA only non null CT': 'Validation - alternative criteria RNA only',
                'Validation - alternative criteria RNA only without HRE non null CT': 'Validation - alternative criteria RNA only without HRE'}
    
    df3_ct['set'] = df3_ct['set'].map(set_mapper)
    
    df3=df3.merge(df3_ct, on='set')
    df3['Pathogen targets not detected on qPCR*'] = df3['gold_standard_count'] - df3['gold_standard_count_non_null_CT']

    # get selection of columns for CT < 35 values
    ct35_sets=['Validation - main criteria CT < 35', 'Validation - alternative criteria CT < 35', 'Validation - alternative criteria RNA only CT < 35', 'Validation - alternative criteria RNA only without HRE CT < 35']
    df3_ct35=g1[g1['set'].isin(ct35_sets)]
    df3_ct35=df3_ct35[['set', 'Sensitivity (per target)']]
    df3_ct35.rename(columns={'Sensitivity (per target)': 'Sensitivity in true-positive targets with CT < 35'}, inplace=True)
    set_mapper={'Validation - main criteria CT < 35': 'Validation - main criteria',
                'Validation - alternative criteria CT < 35': 'Validation - alternative criteria',
                'Validation - alternative criteria RNA only CT < 35': 'Validation - alternative criteria RNA only',
                'Validation - alternative criteria RNA only without HRE CT < 35': 'Validation - alternative criteria RNA only without HRE'}
    df3_ct35['set'] = df3_ct35['set'].map(set_mapper)
    df3=df3.merge(df3_ct35, on='set')

    cols=['set', 'Sensitivity (per target)','Pathogen targets not detected on qPCR*','Sensitivity in true-positive targets with non null CT','Sensitivity in true-positive targets with CT < 35', 'Pathogen targets with 0 mapped reads', 'Sensitivity in true-positive targets with >0 reads', 
    'Sensitivity (per sample)', 'Specificity (per target)',  'False-positive detections', 'Specificity (per sample)']
    df3=df3[cols]
    df3.to_csv('table_3.csv')


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
    p.save('figures/Figure02_OG.pdf', width=10, height=6)

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
    p2.save('figures/Figure02.pdf', width=10, height=6)

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
    p3.save('figures/Figure02_negative_counts.pdf', width=10, height=6)

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
    # add A, B, C, D labels to facets
    label_df = pd.DataFrame({
        'Dataset': ['Derivation', 'Derivation', 'Validation', 'Validation'],
        'Result': ['Positive tests', 'Negative tests', 'Positive tests', 'Negative tests'],
        'x': ['Adenovirus', 'Adenovirus', 'Adenovirus', 'Adenovirus'],
        'y': [57, 300, 57, 300],
        'label': ['A', 'C', 'B', 'D'],
    })
    label_df['Result'] = pd.Categorical(label_df['Result'], categories=['Positive tests', 'Negative tests'], ordered=True)
    p4 = p4 + geom_text(
        data=label_df,
        mapping=aes(x='x', y='y', label='label'),
        inherit_aes=False,
        size=12,
        color='black',
        fontweight='bold'
    )
    p4.save('figures/Figure02_positives_negatives.pdf', width=10, height=6)

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

def plot_figure_S1(df_original, dfd):
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
    p3.save('supplemental/Figure_S1.pdf', width=10, height=6)
    p3.save('supplemental/Figure_S1.svg', width=10, height=6)

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
    df['Sample_QC_pass'] = np.where((df['PCs_passed'] == True) & (df['sample num reads'] >= 2) & (df['run_pass']==1), 'Sample QC Pass', 'Sample QC Fail')

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
     #geom_hline(yintercept=0.007, linetype='dashed', color='grey') +
     geom_hline(yintercept=0.006, linetype='dashed', color='grey') +
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
        'x': [0.001, 0.001, 0.001, 0.001],
        'y': [10, 10, 10, 10],
        'label': ['A', 'B', 'C', 'D'],
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

    p.save('figures/Figure03.pdf', width=10, height=6)

def plot_alt_ratios(df, dfd):
    # plot AuG_trunc10 on x axis and Sample_reads_percent_of_refs on y axis with 
    df['Set']='Validation'
    df['batch positive amplification control']=df['reverse_transcription_control'].astype(str)
    df['batch PCR negative control']=df['amplification_control']
    # readjust pass based on Sample_reads_percent_of_type_run > 0.03
    #df['pass'] = np.where((df['Sample_reads_percent_of_type_run'] > 0.03) 
    #                      & (df['OR pass']==True) & (df['2 reads pass']==True),
    #                       True, False)
    df['pass'] = np.where(df['Sample_reads_percent_of_type_run'] > 0.037, True, False)
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

    # add sample QC pass/fail based on PCs_passed, and 2 reads pass
    df['Sample_QC_pass'] = np.where((df['PCs_passed'] == True) & (df['sample num reads'] >= 2) & (df['run_pass']==1), 'Sample QC Pass', 'Sample QC Fail')

    # order Sample_QC_pass so that 'Sample QC Fail' comes after 'Sample QC Pass'
    df['Sample_QC_pass'] = pd.Categorical(df['Sample_QC_pass'], categories=['Sample QC Pass', 'Sample QC Fail'], ordered=True)

    df.to_csv('figure_S6_input_data.csv', index=False)

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
     geom_hline(yintercept=0.037,  color='red') +
     #geom_hline(yintercept=0.003,  color='pink') +
     scale_x_log10(limits=(0.001, 10)) +
     scale_y_log10(limits=(0.0001, 100)) +
     facet_grid('Sample_QC_pass ~ Set') +
     labs(x='Area under the genome truncated at depth 10 ',
          y='Reads mapping to target reference as percentage of reads\nmapping to references of same type in run'))
    
    # Facet labels a/b/c/d
    label_df = pd.DataFrame({
        'Set': ['Derivation', 'Validation', 'Derivation', 'Validation'],
        'Sample_QC_pass': ['Sample QC Pass', 'Sample QC Pass', 'Sample QC Fail', 'Sample QC Fail'],
        'x': [0.001, 0.001, 0.001, 0.001],
        'y': [100, 100, 100, 100],
        'label': ['A', 'B', 'C', 'D'],
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

    p.save('supplemental/Figure_S5.pdf', width=10, height=6)

    # plot swarm plot of the same data with Sample_reads_percent_of_type_run on y axis
    p2 = (ggplot(df, aes(x='Sample_QC_pass', y='Sample_reads_percent_of_type_run')) +
     geom_jitter(aes(fill='fill_group', 
                    color='color_group',
                    shape='sample_type'),
                size=2,
                stroke=0.5,
                position=position_jitter(width=0.2, height=0)) +
     scale_shape_manual(values=sample_type_shapes, 
                        name='Sample type') +
     scale_color_manual(values=color_colors,
                        name='Outline color') +
     scale_fill_manual(values=fill_colors,
                       name='Fill color (only if 2 reads pass)') +
     geom_hline(yintercept=0.037,  color='red') +
     #geom_hline(yintercept=0.003,  color='pink') +
     scale_y_log10(limits=(0.0001, 100)) +
     facet_grid('. ~ Set') +
     labs(x='Sample QC pass/fail',
          y='Reads mapping to target reference as percentage of reads\nmapping to references of same type in run'))
    
    p2.save('supplemental/Figure_S5_swarm.pdf', width=8, height=6)

def upset_pass_plot(df, set='Validation'):
    df=df[(df['MS2_spike'] != 0) & (df['IC_virus_spike']!=0)]
    df=df[df['PCs_passed']==1]
    #if set=='Validation':
    df=df[df['run_pass']==1]
    df=df[df['gold_standard']==1]

    #df['QC_pass']=np.where((df['PCs_passed']==1) & (df['run_pass']==1), 1, 0)
    #cols=['Run', 'barcode', 'OR pass',	'AND ratio pass',	'2 reads pass', 'QC_pass' ]
    cols=['Run', 'barcode', 'OR pass',	'AND ratio pass',	'2 reads pass' ]
    df=df[cols]

    # group to generate counts
    pass_types=[ 'OR pass',	'AND ratio pass',	'2 reads pass' ]
    # rename pass types for plotting
    pt_rename={'OR pass':'Main criteria (Table 2)', 
            'AND ratio pass':'Ratio criteria (Table 2)', 
            '2 reads pass':'>=2 reads',
            'PCs_passed':'Passed internal positive controls',
            'QC_pass':'QC pass'}    
    df.rename(columns=pt_rename, inplace=True)
    pass_types=[pt_rename[pt] for pt in pass_types]

    g=df.groupby(pass_types).size().reset_index().rename(columns={0:'count'})

    # convert g to series
    g=g.set_index(pass_types)['count']

    # plot
    up.plot(g, sort_categories_by='input', show_counts=True)
    if set=='Validation':
        subplot='d'
    if set=='Derivation':
        subplot='c'
    plt.savefig(f'supplemental/Figure_S3{subplot}_pass_upset_plot.pdf')
    plt.savefig(f'supplemental/Figure_S3{subplot}_pass_upset_plot.svg')

def upset_spiked_plot(df, set='Validation'):
    df=df[df['test_type'].isin(test_type_normalisation.keys())]
    # remove run failures
    df=df[df['run_pass']==1]
    df=df[(df['run_bases_pass']==True) & (df['sample_reads_pass']==True)]

    cols=['Run', 'barcode','MS2_spike',	'IC_virus_spike', 'orthoreovirus passed',	'zika passed',	'MS2 passed',	'murine_respirovirus passed']
    df=df[cols]
    df.drop_duplicates(inplace=True, keep='first')

    # remove spiked == 0
    df=df[(df['MS2_spike'] != 0) & (df['IC_virus_spike']!=0)]

    # group to generate counts
    pathogens=[ 'orthoreovirus passed',	'zika passed',	'MS2 passed',	'murine_respirovirus passed']
    # rename pathogens for plotting
    pt_rename={'orthoreovirus passed':'Orthoreovirus',
            'zika passed':'Zika virus',
            'MS2 passed':'MS2',
            'murine_respirovirus passed':'Murine respirovirus'}
    df.rename(columns=pt_rename, inplace=True)
    pathogens=[pt_rename[pt] for pt in pathogens]
    g=df.groupby(pathogens).size().reset_index().rename(columns={0:'count'})

    # convert g to series
    g=g.set_index(pathogens)['count']

    # plot
    up.plot(g, sort_categories_by='input', show_counts=True)
    if set=='Validation':
        subplot='b'
    if set=='Derivation':
        subplot='a'
    plt.savefig(f'supplemental/Figure_S3{subplot}_spikes_upset_plot.pdf')
    plt.savefig(f'supplemental/Figure_S3{subplot}_spikes_upset_plot.svg')


def upset_spiked_stacked_plot(df, set='Validation'):
    df=df[df['test_type'].isin(test_type_normalisation.keys())]

    df['run_pass']=df['run_pass'].map({True:'QC run pass', False:'QC run fail'})

    cols=['Run', 'barcode','MS2_spike',	'IC_virus_spike', 'orthoreovirus passed',	'zika passed',	'MS2 passed',	'murine_respirovirus passed', 'run_pass']
    df=df[cols]
    df.drop_duplicates(inplace=True, keep='first')

    # remove spiked == 0
    df=df[(df['MS2_spike'] != 0) & (df['IC_virus_spike']!=0)]

    # group to generate counts
    pathogens=[ 'orthoreovirus passed',	'zika passed',	'MS2 passed',	'murine_respirovirus passed']
    # rename pathogens for plotting
    pt_rename={'orthoreovirus passed':'Orthoreovirus',
            'zika passed':'Zika virus',
            'MS2 passed':'MS2',
            'murine_respirovirus passed':'Murine respirovirus'}
    df.rename(columns=pt_rename, inplace=True)
    pathogens=[pt_rename[pt] for pt in pathogens]
    
    df=df.set_index(df['MS2']==1).set_index(df['Orthoreovirus']==1, append=True).set_index(df['Zika virus']==1, append=True).set_index(df['Murine respirovirus']==1, append=True)
    
    # plot
    upset=up.UpSet(df, sort_categories_by='input', intersection_plot_elements=0, show_counts=True, )
    upset.add_stacked_bars('run_pass', colors=['black', 'lightgray'], title='Intersection size')
    upset.plot()
    if set=='Validation':
        subplot='b'
    if set=='Derivation':
        subplot='a'
    plt.savefig(f'supplemental/Figure_S3{subplot}_spikes_upset_stacked_plot.pdf')
    plt.savefig(f'supplemental/Figure_S3{subplot}_spikes_upset_stacked_plot.svg')

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
    plt.savefig('Figure_S4.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    
def plot_qPCR(df, qPCR):
    qPCR_name_mapping = {'PARAINFLUENZA VIRUS 3': 'Parainfluenza 3',
                         'PARAINFLUENZA VIRUS 4': 'Parainfluenza 4',
                          'CORONAVIRUS HKU': 'Coronavirus HKU',
                          'CORONAVIRUS OC43': 'Coronavirus OC43',
                          'CORONAVIRUS NL63': 'Coronavirus NL63',
                          'METAPNEUMOVIRUS': 'Metapneumovirus',
                          'INFLUENZA A': 'Influenza A',
                          'INFLUENZA B': 'Influenza B',
                          'RHINOVIRUS': 'Rhinovirus/enterovirus',
                          'ENTEROVIRUS': 'Rhinovirus/enterovirus',
                          'PARAINFLUENZA VIRUS 1': 'Parainfluenza 1',
                          'PARAINFLUENZA VIRUS 2': 'Parainfluenza 2',
                          'SARSCOV2': 'SARS-CoV-2',
                          'RSVA': 'RSV',
                          'RSVB': 'RSV'}
    qPCR['pathogen'] = qPCR['target'].map(qPCR_name_mapping)
    qPCR=qPCR[qPCR['qpcr_ct'].notna()]
    flu_qPCR=qPCR[qPCR['pathogen']=='Influenza A']
    fluAH1=flu_qPCR.copy()
    fluAH1['pathogen'] = 'Influenza A/H1'
    fluAH3=flu_qPCR.copy()
    fluAH3['pathogen'] = 'Influenza A/H3'
    fluA12009=flu_qPCR.copy()
    fluA12009['pathogen'] = 'Influenza A/H1-2009'
    qPCR=pd.concat([qPCR, fluAH1, fluAH3, fluA12009])
    qPCR['gold_standard']=1
    # replace _ from pathogen names in df
    df['pathogen'] = df['pathogen'].str.replace('_', ' ')
    print('Unique pathogens in qPCR data:', qPCR['target'].unique())
    # merge qPCR data with df on Run and barcode
    df_merged = df.merge(qPCR, on=['Run', 'barcode', 'pathogen', 'gold_standard'], how='left')
    df_merged['pQCR_status'] = np.where(df_merged['qpcr_ct'].notna(), 'qPCR performed', 'qPCR not performed')
    df_merged['qpcr undetermined'] = np.where(df_merged['qpcr_ct']=='UNDETERMINED', 1, 0)
    df_merged['qpcr_ct_value'] = np.where((df_merged['qpcr_ct']!='UNDETERMINED')&(df_merged['pQCR_status']=='qPCR performed'), df_merged['qpcr_ct'], np.nan)
    df_merged['qpcr_ct_value']=df_merged['qpcr_ct_value'].astype(float)
    df_merged.drop_duplicates(subset=['Run', 'barcode', 'pathogen'], keep='first', inplace=True)

    df_merged['ct_1']=df_merged['ct_1'].replace(0, np.nan)
    df_merged['ct_2']=df_merged['ct_2'].replace(0, np.nan)
    # create routine_CT_value column that uses ct_1 if present, else ct_2
    df_merged['routine_CT_value'] = df_merged['ct_1']
    df_merged['routine_CT_value'] = np.where(df_merged['ct_1'].isna(), df_merged['ct_2'], df_merged['routine_CT_value'])
    
    df_merged_orginal=df_merged.copy()
    print(f'Merged dataframe shape: {df_merged.shape}')

    df_merged=df_merged[df_merged['gold_standard']==1]
    df_merged['Test platform']=df_merged['test'].map(test_type_mapping)
    df_merged['chrom']=np.where(df_merged['sample num reads']>0, df_merged['chrom'], 'No reads detected')

    

    df2=df_merged[(df_merged['routine_CT_value'].notna()) & (df_merged['routine_CT_value']>0) & (df_merged['chrom']!='0') & (df_merged['qpcr_ct'].notna()) & (df_merged['qpcr_ct']!='UNDETERMINED')]
    df2.to_csv('FigureS8_input.csv')

    # plot Ct value vs Sample_reads_percent_of_refs
    df2['Flu genotype']=df2['chrom']
    df2['Reads']=df2['sample num reads']
    df2['qpcr_ct'] = df2['qpcr_ct'].astype(float)
    
    p = (ggplot(df2, aes(y='qpcr_ct', x='routine_CT_value')) +
         geom_point(aes(color='pathogen', shape='Test platform'),  alpha=0.7) +
         # add straight line y=x
        geom_abline(slope=1, intercept=0, linetype='dashed', color='red') +
        labs(y='qPCR Ct Value', x='Clinical Ct Value', tag='A') +
        # limit axes
        scale_x_continuous(limits=(15, 40)) +
        scale_y_continuous(limits=(15, 40)) +
        theme(legend_position='right',
              legend_text=element_text(size=10),
              legend_title=element_blank())
        # add facet by pathogen
        #+ facet_wrap('~pathogen')
    )

    p.save('supplemental/Figure_S7a.pdf', width=8, height=6)

    # box plot of Ct values by method
    # melt cts to long format
    df_long = pd.melt(df_merged, id_vars=['Run', 'barcode', 'pathogen','chrom','sample num reads'], 
                      value_vars=['qpcr_ct', 'routine_CT_value'], var_name='Method', value_name='Ct Value')

    df_long['Method'] = df_long['Method'].map({'qpcr_ct': 'Assay qPCR', 'routine_CT_value': 'Clinical qPCR'}) 
    # order pathogens by list
    pathogen_order = ['Influenza A', 'Influenza B','RSV','SARS-CoV-2',
        'Adenovirus', 
        'Coronavirus HKU', 'Coronavirus NL63', 'Coronavirus OC43', 'Coronavirus 229E', 
        'Metapneumovirus',  
        'Parainfluenza 1', 'Parainfluenza 2', 'Parainfluenza 3', 'Parainfluenza 4',  
        'Rhinovirus/enterovirus']
    df_long['pathogen'] = pd.Categorical(df_long['pathogen'], categories=pathogen_order, ordered=True)
    df_long=df_long[df_long['Ct Value'].notna()]
    df_long=df_long[df_long['Ct Value']!= 'UNDETERMINED']
    #df_long['Ct Value'] = df_long['Ct Value'].astype(float)
    df_long['Ct Value'] = df_long['Ct Value'].astype(float)
    df_long=df_long[df_long['Ct Value']>0]
    df_long.to_csv('FigureS7b_input.csv')
     # boxplot of Ct values by pathogen and method
     # save as FigureS8b.pdf
    p2 = (ggplot(df_long, aes(x='pathogen', y='Ct Value', fill='Method')) +
         geom_boxplot(alpha=0.7) +
         #geom_jitter(width=0.2, alpha=0.5) +
         labs(x='Pathogen target', y='Ct Value', tag='B') +
         theme(axis_text_x=element_text(rotation=90, hjust=1),
               legend_text=element_text(size=10),
               legend_title=element_blank())
    )
    p2.save('supplemental/Figure_S7b.pdf', width=8, height=6)

    px=p/p2
    px.save('supplemental/Figure_S7_combined.pdf', width=10, height=15)


    # Ct vs reads plot
    df3=df_merged[df_merged['sample num reads']>0]
    df3=df3[df3['qpcr_ct'].notna()]
    df3=df3[df3['qpcr_ct']!='UNDETERMINED']
    df3['qpcr_ct'] = df3['qpcr_ct'].astype(float)
    df3=df3[df3['qpcr_ct']>0]
    p=(ggplot(df3, aes(x='qpcr_ct', y='sample num reads')) +
       geom_point() +
       #geom_smooth(method='lm', se=True) + 
       #geom_smooth(se=True) + 
       scale_y_log10(labels=label_comma()) +  # <- plain numbers
       scale_x_reverse() + # reverse x axis
       labs(x='Ct Value', y='Reads') +
       facet_wrap('~pathogen') +
       theme(legend_title=element_blank(), panel_spacing = 0.025)
    )
    p.save('figures/Figure05_assayCt.pdf', width=8, height=4)

    # clinical Ct vs reads plot
    df4=df_merged[df_merged['sample num reads']>0]
    df4=df4[df4['routine_CT_value'].notna()]
    #df4=df4[df4['ct_1']!='UNDETERMINED']
    df4['routine_CT_value'] = df4['routine_CT_value'].astype(float)
    df4=df4[df4['routine_CT_value']>0]
    p=(ggplot(df4, aes(x='routine_CT_value', y='sample num reads')) +
       geom_point() +
       #geom_smooth(method='lm', se=True) + 
       #geom_smooth(se=True) + 
       scale_y_log10(labels=label_comma()) +  # <- plain numbers
       scale_x_reverse() + # reverse x axis
       labs(x='Ct Value', y='Reads') +
       facet_wrap('~pathogen') +
       theme(legend_title=element_blank(), panel_spacing = 0.025)
    )
    p.save('figures/Figure05_clinicalCt.pdf', width=8, height=4)

    # Quant vs reads plot
    #p=(ggplot(df3, aes(x='qpcr_quantity', y='sample num reads')) +
    #   geom_point() +
    #   #geom_smooth(method='lm', se=True) + 
    #   #geom_smooth(se=True) + 
    #   scale_y_log10(labels=label_comma()) +  # <- plain numbers
    #   scale_x_log10(labels=label_comma()) +
    #   #scale_x_reverse() + # reverse x axis
    #   labs(x='Assay Quantity', y='Reads') +
    #   facet_wrap('~pathogen') +
    #   theme(legend_title=element_blank())
    #)
    #p.save('figures/Figure05_assayQuantity.pdf', width=6, height=4)

    # plot both types of Cts

    p4=(ggplot(df_long, aes(x='Ct Value', y='sample num reads', color='pathogen')) +
       geom_point() +
       #geom_smooth(method='lm', se=True) + 
       #geom_smooth(se=True) + 
       scale_y_log10(labels=label_comma()) +  # <- plain numbers
       #scale_x_log10(labels=label_comma()) +
       scale_x_reverse() + # reverse x axis
       labs(x='Ct Value', y='Reads') +
       facet_wrap('~Method') +
       theme(legend_title=element_blank())
       # change colour scale to divergent
         + scale_color_brewer(type='qual', palette='Set1')
    )
    # Add A/B tags
    label_df = pd.DataFrame({
        'Method': ['Assay qPCR', 'Clinical qPCR'],
        'x': [40, 40],
        'y': [max(df_long['sample num reads']), max(df_long['sample num reads'])],
        'label': ['A', 'B'],
    })
    p4 = p4 + geom_text(
        data=label_df,
        mapping=aes(x='x', y='y', label='label'),
        inherit_aes=False,
        size=12,
        color='black',
        fontweight='bold'
    )
    
    p4.save('figures/Figure05_assayCt_combined.pdf', width=12, height=4)
    return df_merged_orginal


def plot_mapQ(df, dfd):
    dfd=dfd[['Run', 'barcode','pathogen','gold_standard','FLU_A_POS', 'TOP_FLU_A']]
    # readjust gold_standard for flu A
    print(dfd['pathogen'].unique())
    dfd['gold_standard'] = np.where((dfd['FLU_A_POS']==True)&(dfd['TOP_FLU_A']==True), 1, dfd['gold_standard'])
    df['pathogen'] = df['pathogen_reduced']
    # convert barcodes to int
    df['barcode'] = df['barcode'].astype(int)
    dfd['barcode'] = dfd['barcode'].astype(int)
    df=df.merge(dfd, on=['Run', 'barcode','pathogen'], how='left')
    df=df[df['gold_standard'].notna()]
    df['pathogen'] = df['pathogen'].str.replace('_', ' ')
    print(df['gold_standard'].unique())
    dfu=df[['pathogen']][df['gold_standard'].isna()]
    print(dfu['pathogen'].unique())
    df['Species correct'] = np.where(df['gold_standard']==1, True, False)
    # Replace 0,1 with False, True in 'True species hit' column
    df['True species hit'] = df['True species hit'].map({0: 'False', 1: 'True', np.nan: 'Unknown'})
    # remove duplicate readIDs
    df=df[df['readID_dup']==False]
    # bacteria species list
    bacteria_species=['Bordetella pertussis', 'Bordetella parapertussis', 'Chlamydia pneumoniae', 'Mycoplasma pneumoniae']
    # plot the edit distance by mapQ score for each pathogen, colour by 'True species hit'
    p=(ggplot(df, aes(x='edit_percentage', y='mapQ')) +
         geom_point(aes(color='Species correct'), alpha=0.5) +
            facet_wrap('~pathogen',ncol=3) +
            labs(x='Read difference to reference genome (%)', 
                 y='Mapping quality score (mapQ)') +
            theme(legend_position='bottom', 
                  legend_title=element_text(text='Species positive by gold standard'))
            # draw red line at mapQ=40 only for Bordetta pertussis
            + geom_hline(aes(yintercept=40), linetype="dashed", color = "red",
                          data=df[df['pathogen'].isin(bacteria_species)])
    )
    p.save('supplemental/Figure_S8_mapQ_vs_edit_distance.pdf', width=8, height=8)

def split_tests_from_IORD(row):
    s=row['iord_pcr_organisms'].split(';') if pd.notna(row['iord_pcr_organisms']) else []
    for i in s:
        if ':' in i:
            i_split = i.split(':')
            if i_split[0] not in row:
                row[i_split[0]] = i_split[-1]
            elif row[i_split[0]] is np.nan:
                row[i_split[0]] = i_split[-1]
            else: # if there is already a value, append with ;
                if i_split[-1] not in row[i_split[0]].split(';'):
                    row[i_split[0]] = row[i_split[0]] + ';' + i_split[-1]
            #row[i_split[0]]=i_split[-1]
        else:
            row[i]=None
    return row

def plot_ages(df):
    # plot age distribution of samples in df
    df2=df.drop_duplicates(subset=['Run', 'barcode'])
    p=(ggplot(df2, aes(x='AgeInYearsAtCollection', fill='LinkedSex')) +
       geom_histogram(binwidth=5, alpha=0.7) +
       labs(x='Age (years)', y='Number of samples') 
    )
    p.save('figures/age_distribution.pdf', width=6, height=4)

def derivation_set_characteristics(pcr_types, metaDF, dfd):
    df=pd.read_csv(pcr_types)
        # Replace NULL with np.nan
    df.replace('NULL', np.nan, inplace=True)
    df.replace('NA', np.nan, inplace=True)
    df.replace('nan', np.nan, inplace=True)

    ## explode colulmn 'iord_pcr_organisms' on :
    #df['iord_pcr_organisms'] = df['iord_pcr_organisms'].str.split(':')
    #df = df.explode('iord_pcr_organisms')
    #print(df)

    df = df.apply(split_tests_from_IORD, axis=1)
    # remove empty columns
    df = df.dropna(axis=1, how='all')

    # add columns for each test
    test_dict = {'FLRA': 'Alinity',	
                'FLRP': 'Ceipheid',
                'RPCR': 'BioFire'}

    for col in ['FLRA', 'FLRP', 'RPCR']:
        df[test_dict[col]] = np.where(df[col].notna(), 1, 0)

    # merge with metaDF.csv
    metaDF = pd.read_csv(metaDF)
    metaDF=metaDF[metaDF['Biofire positive']==1]
    metaDF['pathogen'] = metaDF['pathogen'].map(str)
    metaDF=metaDF.groupby(['Run','barcode'])['pathogen'].apply(','.join).reset_index()
    #metaDF.drop_duplicates(subset=['Run','barcode'], inplace=True)
    df = pd.merge(metaDF, df, left_on=['Run', 'barcode'],right_on=['run_name', 'barcode'], how='left')

    # remove rows with all False in test columns
    df = df[df[['Alinity', 'Ceipheid', 'BioFire']].any(axis=1)]

    # count number of unique rows
    print(f'Number of samples in dataset: {len(df)}')
    # Count the number of unique combinations of accession_1 and accession_2
    unique_samples = df[['accession_1', 'accession_2']].drop_duplicates()
    print(f'Number of unique samples in dataset: {len(unique_samples)}')

    # count the number of samples where pooled_with_accession_1 and pooled_with_accession_2 are not null
    pooled_samples = df[~df['pooled_with_accession_1'].isna() & ~df['pooled_with_accession_2'].isna()]
    print(f'Number of pooled samples in dataset: {len(pooled_samples)}')

    # count the number of pathogens per sample
    df['num_pathogens'] = df['pathogen'].str.split(',').apply(lambda x: len(x) if isinstance(x, list) else 0)
    df['num_pathogens'] = np.where(df['pathogen']=='nan', 0, df['num_pathogens'])

    print(f'Number of samples with no pathogens: {len(df[df["num_pathogens"] == 0])}')
    print(f'Number of samples with one pathogen: {len(df[df["num_pathogens"] == 1])}')
    print(f'Number of samples with more than one pathogen: {len(df[df["num_pathogens"] > 1])}')


    # count the number of unique samples in pooled_with_accession_1 and pooled_with_accession_2
    unique_pooled_samples = pooled_samples[['pooled_with_accession_1', 'pooled_with_accession_2']].drop_duplicates()
    print(f'Number of unique pooled samples in dataset: {len(unique_pooled_samples)}')

    # count number of Alinity, Cepheid, BioFire tests
    num_Alinity_tests = df['Alinity'].sum()
    num_Cepheid_tests = df['Ceipheid'].sum()
    num_BioFire_tests = df['BioFire'].sum()
    print(f'Number of Alinity tests in dataset: {num_Alinity_tests}')
    print(f'Number of Cepheid tests in dataset: {num_Cepheid_tests}')
    print(f'Number of BioFire tests in dataset: {num_BioFire_tests}')
    num_had_tests = df[['Alinity', 'Ceipheid', 'BioFire']].any(axis=1).sum()
    print(f'Number of samples with at least one test in dataset: {num_had_tests}')
    # calculate bases and reads median and IQR
    #dfd=pd.read_csv(derivation_set)
    #dfd['pathogen'] = dfd['pathogen_reduced']
    #dfd=dfd[['Run', 'barcode','total run bases','total run reads','run_pass','PCs_passed', 'pathogen','sample num reads']]
    dfd['total run bases']=dfd.groupby(['Run', 'barcode'])['total run bases'].transform('max')
    dfd['total run reads']=dfd.groupby(['Run', 'barcode'])['total run reads'].transform('max')
    dfdd=dfd.drop_duplicates(subset=['Run', 'barcode'])
    df=df.merge(dfdd, on=['Run', 'barcode'], how='left')
    bases=df.drop_duplicates(subset=['Run', 'barcode'])['total run bases'].median()
    gb=bases/1_000_000_000
    q75, q25 = np.percentile(df.drop_duplicates(subset=['Run', 'barcode'])['total run bases'],  [75 ,25])
    bases_IQR = f'{q25/1_000_000_000:.2f} - {q75/1_000_000_000:.2f}'
    reads=df.drop_duplicates(subset=['Run', 'barcode'])['total run reads'].median()
    q75, q25 = np.percentile(df.drop_duplicates(subset=['Run', 'barcode'])['total run reads'],  [75 ,25])
    reads_IQR = f'{q25/1_000_000:.3f} - {q75/1_000_000:.3f}'

    # Any target pathogen reads identified = sum of all the samples with any reads for target pathogens (median IQR)
    dfd['test normalized'] = dfd['test'].map(test_type_normalisation)
    print(dfd['test normalized'].unique())
    dfd['reads_to_any_reference']=dfd[dfd['pathogen']!='unmapped'].groupby(['Run', 'barcode'])['sample num reads'].transform('sum')
    median_atpri = dfd[(dfd['test normalized'].isin(test_type_normalisation.values())) ].drop_duplicates(subset=['Run', 'barcode'])['reads_to_any_reference'].median()
    q75, q25 = np.percentile(dfd[(dfd['test normalized'].isin(test_type_normalisation.values()))].drop_duplicates(subset=['Run', 'barcode'])['reads_to_any_reference'],  [75 ,25])
    atpri_N = dfd[(dfd['test normalized'].isin(test_type_normalisation.values()))].drop_duplicates(subset=['Run', 'barcode']).shape[0]
    atpri_IQR = f'{q25} - {q75}'
    Any_target_pathogen_reads_identified = f'{median_atpri:.0f} (IQR: {q25:.0f} - {q75:.0f}) n={atpri_N:.0f}'

    # same again but only for samples with reads_to_any_reference > 0
    median_atpri_gs = dfd[(dfd['test normalized'].isin(test_type_normalisation.values())) & (dfd['reads_to_any_reference']>0)]['reads_to_any_reference'].median()
    q75_gs, q25_gs = np.percentile(dfd[(dfd['test normalized'].isin(test_type_normalisation.values())) & (dfd['reads_to_any_reference']>0)]['reads_to_any_reference'],  [75 ,25])
    atpri_N_gs = dfd[(dfd['test normalized'].isin(test_type_normalisation.values())) & (dfd['reads_to_any_reference']>0)].drop_duplicates(subset=['Run', 'barcode']).shape[0]
    atpri_IQR_gs = f'{q25_gs} - {q75_gs}'
    Any_target_pathogen_reads_identified_gs = f'{median_atpri_gs:.0f} (IQR: {q25_gs:.0f} - {q75_gs:.0f}) n={atpri_N_gs:.0f}'

    # Gold standard pathogen reads identified 
    median_gspri =  dfd[(dfd['test normalized'].isin(test_type_normalisation.values())) & (dfd['gold_standard']==1) ]['sample num reads'].median()
    q75, q25 = np.percentile(dfd[(dfd['test normalized'].isin(test_type_normalisation.values())) & (dfd['gold_standard']==1) ]['sample num reads'],  [75 ,25])
    gspri_N = dfd[(dfd['test normalized'].isin(test_type_normalisation.values())) & (dfd['gold_standard']==1)].shape[0]
    gspri_IQR = f'{q25} - {q75}'
    Gold_standard_pathogen_reads_identified = f'{median_gspri:.0f} (IQR: {q25:.0f} - {q75:.0f}) n={gspri_N:.0f}'

    # Gold standard pathogen reads identified where reads > 0
    median_gspri_gs = dfd[(dfd['test normalized'].isin(test_type_normalisation.values())) & (dfd['gold_standard']==1) & (dfd['sample num reads']>0) ]['sample num reads'].median()
    q75_gs, q25_gs = np.percentile(dfd[(dfd['test normalized'].isin(test_type_normalisation.values())) & (dfd['gold_standard']==1) & (dfd['sample num reads']>0) ]['sample num reads'],  [75 ,25])
    gspri_N_gs = dfd[(dfd['test normalized'].isin(test_type_normalisation.values())) & (dfd['gold_standard']==1) & (dfd['sample num reads']>0)].shape[0]
    gspri_IQR_gs = f'{q25_gs} - {q75_gs}'
    Gold_standard_pathogen_reads_identified_gs = f'{median_gspri_gs:.0f} (IQR: {q25_gs:.0f} - {q75_gs:.0f}) n={gspri_N_gs:.0f}'

    # Count number of FPs with single read
    dfd['PT_single_read_FP'] = np.where((dfd['test normalized'].isin(test_type_normalisation.values())) & (dfd['pathogen']!='unmapped') & (dfd['gold_standard']==0) & (dfd['sample num reads']==1), 1, 0)
    dfd['PS_single_read_FPs'] = dfd.groupby(['Run', 'barcode'])['PT_single_read_FP'].transform('sum')
    dfd['PS_single_read_FP'] = np.where(dfd['PS_single_read_FPs']>0, 1, 0)
    # Count number of FPs with two or more reads
    dfd['PT_twoMore_reads_FP'] = np.where((dfd['test normalized'].isin(test_type_normalisation.values())) & (dfd['pathogen']!='unmapped') & (dfd['gold_standard']==0) & (dfd['sample num reads']>1), 1, 0)
    dfd['PS_twoMore_reads_FPs'] = dfd.groupby(['Run', 'barcode'])['PT_twoMore_reads_FP'].transform('sum')
    dfd['PS_twoMore_reads_FP'] = np.where(dfd['PS_twoMore_reads_FPs']>0, 1, 0)
    # calculate per sample FP rates
    num_fps_single_read = dfd[dfd['PS_single_read_FP']>0].drop_duplicates(subset=['Run', 'barcode'])['PS_single_read_FPs'].shape[0]
    num_fps_twoMore_reads = dfd[dfd['PS_twoMore_reads_FP']>0].drop_duplicates(subset=['Run', 'barcode'])['PS_twoMore_reads_FPs'].shape[0]
    n_denominator = dfd[(dfd['test normalized'].isin(test_type_normalisation.values()))].drop_duplicates(subset=['Run', 'barcode']).shape[0]
    fp_reads=f'1={num_fps_single_read}, 2+={num_fps_twoMore_reads}, n={n_denominator}'

    # Count number of TPs with single read
    dfd['PT_single_read_TP'] = np.where((dfd['test normalized'].isin(test_type_normalisation.values())) & (dfd['pathogen']!='unmapped') & (dfd['gold_standard']==1) & (dfd['sample num reads']==1), 1, 0)
    dfd['PS_single_read_TPs'] = dfd.groupby(['Run', 'barcode'])['PT_single_read_TP'].transform('sum')
    dfd['PS_single_read_TP'] = np.where(dfd['PS_single_read_TPs']>0, 1, 0)
    # Count number of TPs with two or more reads
    dfd['PT_twoMore_reads_TP'] = np.where((dfd['test normalized'].isin(test_type_normalisation.values())) & (dfd['pathogen']!='unmapped') & (dfd['gold_standard']==1) & (dfd['sample num reads']>1), 1, 0)
    dfd['PS_twoMore_reads_TPs'] = dfd.groupby(['Run', 'barcode'])['PT_twoMore_reads_TP'].transform('sum')
    dfd['PS_twoMore_reads_TP'] = np.where(dfd['PS_twoMore_reads_TPs']>0, 1, 0)
    # calculate per sample TP rates
    num_tps_single_read = dfd[dfd['PS_single_read_TP']>0].drop_duplicates(subset=['Run', 'barcode'])['PS_single_read_TPs'].shape[0]
    num_tps_twoMore_reads = dfd[dfd['PS_twoMore_reads_TP']>0].drop_duplicates(subset=['Run', 'barcode'])['PS_twoMore_reads_TPs'].shape[0]
    n_denominator = dfd[(dfd['test normalized'].isin(test_type_normalisation.values()))].drop_duplicates(subset=['Run', 'barcode']).shape[0]
    tp_reads=f'1={num_tps_single_read}, 2+={num_tps_twoMore_reads}, n={n_denominator}'

    d1={
        'Set': 'derivation set sequences',
        'Total sequences': len(df),
        'total samples': len(unique_samples),
        'pooled samples': len(unique_pooled_samples),
        'samples with at least one test': num_had_tests,
        'Alinity m Respiratory Panel samples': num_Alinity_tests,
        'BioFire Respiratory Panel 2.1 samples': num_BioFire_tests,
        'Cepheid Xpert Xpress SARS-CoV-2/Flu/RSV samples': num_Cepheid_tests,
        '0 pathogens': len(df[df['num_pathogens'] == 0]),
        '1 pathogens': len(df[df['num_pathogens'] == 1]),
        '2 pathogens': len(df[df['num_pathogens'] == 2]),
        '3 pathogens': len(df[df['num_pathogens'] == 3]),
        'Total pathogens detected': df['num_pathogens'].sum(),
        'Any target pathogen reads identified': Any_target_pathogen_reads_identified,
        'Any target pathogen reads identified (if reads>0)': Any_target_pathogen_reads_identified_gs,
        '(target) Gold standard pathogen reads identified': Gold_standard_pathogen_reads_identified,
        '(target) Gold standard pathogen reads identified (if reads>0)': Gold_standard_pathogen_reads_identified_gs,
        'Gigabases': f'{gb:.2f} (IQR: {bases_IQR})',
        'Million Reads': f'{reads/1_000_000:.3f} (IQR: {reads_IQR})',
        '(Per sample) False positives by reads': fp_reads,
        '(Per sample) True positives by reads': tp_reads,
    }

    # remove duplicate accessions
    df2 = df.drop_duplicates(subset=['accession_1', 'accession_2'])
    df2 = df2.copy()
    df2 = df2[df2[['Alinity', 'Ceipheid', 'BioFire']].any(axis=1)]
    num_unique_had_tests = len(df2)
    print(f'Number of unique samples with at least one test in dataset: {num_unique_had_tests}')
    num_Alinity_tests = df2['Alinity'].sum()
    num_Cepheid_tests = df2['Ceipheid'].sum()
    num_BioFire_tests = df2['BioFire'].sum()
    print(f'Number of unique Alinity tests in dataset: {num_Alinity_tests}')
    print(f'Number of unique Cepheid tests in dataset: {num_Cepheid_tests}')
    print(f'Number of unique BioFire tests in dataset: {num_BioFire_tests}')

    # print number of pathogens per sample
    print(f'Number of unique samples with no pathogens: {len(df2[df2["num_pathogens"] == 0])}')
    print(f'Number of unique samples with one pathogen: {len(df2[df2["num_pathogens"] == 1])}')
    print(f'Number of unique samples with more than one pathogen: {len(df2[df2["num_pathogens"] > 1])}')

    # caluclate number of tests per sample
    df['num tests'] = df[['Alinity', 'Ceipheid', 'BioFire']].sum(axis=1)

    # Add test_type column
    df['test_type'] = np.where(df['Alinity'] == 1, 'ALINITY',None)
    df['test_type'] = np.where(df['Ceipheid'] == 1, 'CEPHEID', df['test_type'])
    df['test_type'] = np.where(df['BioFire'] == 1, 'BIOFIRE', df['test_type'])

    d2={
        'Set': 'derivation set samples (unique accessions)',
        'Total sequences': len(df2),
        'total samples': len(df2),
        'samples with at least one test': num_unique_had_tests,
        'Alinity m Respiratory Panel samples': num_Alinity_tests,
        'BioFire Respiratory Panel 2.1 samples': num_BioFire_tests,
        'Cepheid Xpert Xpress SARS-CoV-2/Flu/RSV samples': num_Cepheid_tests,
        '0 pathogens': len(df2[df2['num_pathogens'] == 0]),
        '1 pathogens': len(df2[df2['num_pathogens'] == 1]),
        '2 pathogens': len(df2[df2['num_pathogens'] == 2]),
        '3 pathogens': len(df2[df2['num_pathogens'] == 3]),
        'Total pathogens detected': df2['num_pathogens'].sum(),
        'Any target pathogen reads identified': None,
        'Any target pathogen reads identified (if reads>0)': None,
        '(target) Gold standard pathogen reads identified': None,
        '(target) Gold standard pathogen reads identified (if reads>0)': None,
        'Gigabases': f'{gb:.2f} (IQR: {bases_IQR})',
        'Million Reads': f'{reads/1_000_000:.3f} (IQR: {reads_IQR})',
        '(Per sample) False positives by reads': None,
    }

    # Remove samples that failed QC
    df3=df[(df['run_pass']==True) & (df['PCs_passed']==True)]
    df3=df3.copy()
    unique_samples_qc = df3[['accession_1', 'accession_2']].drop_duplicates()
    print(f'Number of unique samples in dataset after QC: {len(unique_samples_qc)}')
     # count number of Alinity, Cepheid, BioFire tests
    df3 = df3[df3[['Alinity', 'Ceipheid', 'BioFire']].any(axis=1)]
    num_unique_had_tests = len(df3)
    print(f'Number of unique samples with at least one test in dataset: {num_unique_had_tests}')
    num_Alinity_tests = df3['Alinity'].sum()
    num_Cepheid_tests = df3['Ceipheid'].sum()
    num_BioFire_tests = df3['BioFire'].sum()
    print(f'Number of unique Alinity tests in dataset: {num_Alinity_tests}')
    print(f'Number of unique Cepheid tests in dataset: {num_Cepheid_tests}')
    print(f'Number of unique BioFire tests in dataset: {num_BioFire_tests}')

    # print number of pathogens per sample
    print(f'Number of unique samples with no pathogens: {len(df3[df3["num_pathogens"] == 0])}')
    print(f'Number of unique samples with one pathogen: {len(df3[df3["num_pathogens"] == 1])}')
    print(f'Number of unique samples with more than one pathogen: {len(df3[df3["num_pathogens"] > 1])}')

    # caluclate number of tests per sample
    df3['num tests'] = df3[['Alinity', 'Ceipheid', 'BioFire']].sum(axis=1)

    # Add test_type column
    df3['test_type'] = np.where(df3['Alinity'] == 1, 'ALINITY',None)
    df3['test_type'] = np.where(df3['Ceipheid'] == 1, 'CEPHEID', df3['test_type'])
    df3['test_type'] = np.where(df3['BioFire'] == 1, 'BIOFIRE', df3['test_type'])

    # Any target pathogen reads identified = sum of all the samples with any reads for target pathogens (median IQR)
    df2=dfd[(dfd['run_pass']==True) & (dfd['PCs_passed']==True)]
    df2['reads_to_any_reference']=df2[df2['pathogen']!='unmapped'].groupby(['Run', 'barcode'])['sample num reads'].transform('sum')
    median_atpri = df2[(df2['test normalized'].isin(test_type_normalisation.values())) ].drop_duplicates(subset=['Run', 'barcode'])['reads_to_any_reference'].median()
    q75, q25 = np.percentile(df2[(df2['test normalized'].isin(test_type_normalisation.values()))].drop_duplicates(subset=['Run', 'barcode'])['reads_to_any_reference'],  [75 ,25])
    atpri_N = df2[(df2['test normalized'].isin(test_type_normalisation.values()))].drop_duplicates(subset=['Run', 'barcode']).shape[0]
    atpri_IQR = f'{q25} - {q75}'
    Any_target_pathogen_reads_identified = f'{median_atpri:.0f} (IQR: {q25:.0f} - {q75:.0f}) n={atpri_N:.0f}'

    # same again but only for samples with reads_to_any_reference > 0
    median_atpri_gs = df2[(df2['test normalized'].isin(test_type_normalisation.values())) & (df2['reads_to_any_reference']>0)]['reads_to_any_reference'].median()
    q75_gs, q25_gs = np.percentile(df2[(df2['test normalized'].isin(test_type_normalisation.values())) & (df2['reads_to_any_reference']>0)]['reads_to_any_reference'],  [75 ,25])
    atpri_N_gs = df2[(df2['test normalized'].isin(test_type_normalisation.values())) & (df2['reads_to_any_reference']>0)].drop_duplicates(subset=['Run', 'barcode']).shape[0]
    atpri_IQR_gs = f'{q25_gs} - {q75_gs}'
    Any_target_pathogen_reads_identified_gs = f'{median_atpri_gs:.0f} (IQR: {q25_gs:.0f} - {q75_gs:.0f}) n={atpri_N_gs:.0f}'

    # Gold standard pathogen reads identified 
    median_gspri = df2[(df2['test normalized'].isin(test_type_normalisation.values())) & (df2['gold_standard']==1) ]['sample num reads'].median()
    q75, q25 = np.percentile(df2[(df2['test normalized'].isin(test_type_normalisation.values())) & (df2['gold_standard']==1) ]['sample num reads'],  [75 ,25])
    gspri_N = df2[(df2['test normalized'].isin(test_type_normalisation.values())) & (df2['gold_standard']==1)].shape[0]
    gspri_IQR = f'{q25} - {q75}'
    Gold_standard_pathogen_reads_identified = f'{median_gspri:.0f} (IQR: {q25:.0f} - {q75:.0f}) n={gspri_N:.0f}'

    # Gold standard pathogen reads identified where reads > 0
    median_gspri_gs = df2[(df2['test normalized'].isin(test_type_normalisation.values())) & (df2['gold_standard']==1) & (df2['sample num reads']>0) ]['sample num reads'].median()
    q75_gs, q25_gs = np.percentile(df2[(df2['test normalized'].isin(test_type_normalisation.values())) & (df2['gold_standard']==1) & (df2['sample num reads']>0) ]['sample num reads'],  [75 ,25])
    gspri_N_gs = df2[(df2['test normalized'].isin(test_type_normalisation.values())) & (df2['gold_standard']==1) & (df2['sample num reads']>0)].shape[0]
    gspri_IQR_gs = f'{q25_gs} - {q75_gs}'
    Gold_standard_pathogen_reads_identified_gs = f'{median_gspri_gs:.0f} (IQR: {q25_gs:.0f} - {q75_gs:.0f}) n={gspri_N_gs:.0f}'

    # Count number of FPs with single read
    df2['PT_single_read_FP'] = np.where((df2['test normalized'].isin(test_type_normalisation.values())) & (df2['pathogen']!='unmapped') & (df2['gold_standard']==0) & (df2['sample num reads']==1), 1, 0)
    df2['PS_single_read_FPs'] = df2.groupby(['Run', 'barcode'])['PT_single_read_FP'].transform('sum')
    df2['PS_single_read_FP'] = np.where(df2['PS_single_read_FPs']>0, 1, 0)
    # Count number of FPs with two or more reads
    df2['PT_twoMore_reads_FP'] = np.where((df2['test normalized'].isin(test_type_normalisation.values())) & (df2['pathogen']!='unmapped') & (df2['gold_standard']==0) & (df2['sample num reads']>1), 1, 0)
    df2['PS_twoMore_reads_FPs'] = df2.groupby(['Run', 'barcode'])['PT_twoMore_reads_FP'].transform('sum')
    df2['PS_twoMore_reads_FP'] = np.where(df2['PS_twoMore_reads_FPs']>0, 1, 0)
    # calculate per sample FP rates
    num_fps_single_read = df2[df2['PS_single_read_FP']>0].drop_duplicates(subset=['Run', 'barcode'])['PS_single_read_FPs'].shape[0]
    num_fps_twoMore_reads = df2[df2['PS_twoMore_reads_FP']>0].drop_duplicates(subset=['Run', 'barcode'])['PS_twoMore_reads_FPs'].shape[0]
    n_denominator = df2[(df2['test normalized'].isin(test_type_normalisation.values()))].drop_duplicates(subset=['Run', 'barcode']).shape[0]
    fp_reads=f'1={num_fps_single_read}, 2+={num_fps_twoMore_reads}, n={n_denominator}'

    # Count number of TPs with single read
    df2['PT_single_read_TP'] = np.where((df2['test normalized'].isin(test_type_normalisation.values())) & (df2['pathogen']!='unmapped') & (df2['gold_standard']==1) & (df2['sample num reads']==1), 1, 0)
    df2['PS_single_read_TPs'] = df2.groupby(['Run', 'barcode'])['PT_single_read_TP'].transform('sum')
    df2['PS_single_read_TP'] = np.where(df2['PS_single_read_TPs']>0, 1, 0)
    # Count number of TPs with two or more reads
    df2['PT_twoMore_reads_TP'] = np.where((df2['test normalized'].isin(test_type_normalisation.values())) & (df2['pathogen']!='unmapped') & (df2['gold_standard']==1) & (df2['sample num reads']>1), 1, 0)
    df2['PS_twoMore_reads_TPs'] = df2.groupby(['Run', 'barcode'])['PT_twoMore_reads_TP'].transform('sum')
    df2['PS_twoMore_reads_TP'] = np.where(df2['PS_twoMore_reads_TPs']>0, 1, 0)
    # calculate per sample TP rates
    num_tps_single_read = df2[df2['PS_single_read_TP']>0].drop_duplicates(subset=['Run', 'barcode'])['PS_single_read_TPs'].shape[0]
    num_tps_twoMore_reads = df2[df2['PS_twoMore_reads_TP']>0].drop_duplicates(subset=['Run', 'barcode'])['PS_twoMore_reads_TPs'].shape[0]
    n_denominator = df2[(df2['test normalized'].isin(test_type_normalisation.values()))].drop_duplicates(subset=['Run', 'barcode']).shape[0]
    tp_reads=f'1={num_tps_single_read}, 2+={num_tps_twoMore_reads}, n={n_denominator}'

    df3={
        'Set': 'Derivation set sequences (QC passed)',
        'Total sequences': len(df3),
        'total samples': len(unique_samples_qc),
        'samples with at least one test': num_unique_had_tests,
        'Alinity m Respiratory Panel samples': num_Alinity_tests,
        'BioFire Respiratory Panel 2.1 samples': num_BioFire_tests,
        'Cepheid Xpert Xpress SARS-CoV-2/Flu/RSV samples': num_Cepheid_tests,
        '0 pathogens': len(df3[df3['num_pathogens'] == 0]),
        '1 pathogens': len(df3[df3['num_pathogens'] == 1]),
        '2 pathogens': len(df3[df3['num_pathogens'] == 2]),
        '3 pathogens': len(df3[df3['num_pathogens'] == 3]),
        'Total pathogens detected': df3['num_pathogens'].sum(),
        'Any target pathogen reads identified': Any_target_pathogen_reads_identified,
        'Any target pathogen reads identified (if reads>0)': Any_target_pathogen_reads_identified_gs,
        '(target) Gold standard pathogen reads identified': Gold_standard_pathogen_reads_identified,
        '(target) Gold standard pathogen reads identified (if reads>0)': Gold_standard_pathogen_reads_identified_gs,
        'Gigabases': f'{gb:.2f} (IQR: {bases_IQR})',
        'Million Reads': f'{reads/1_000_000:.3f} (IQR: {reads_IQR})',
        '(Per sample) False positives by reads': fp_reads,
        '(Per sample) True positives by reads': tp_reads,
    }

    df[df["num_pathogens"] == 0].to_csv('derivation_set_zero_pathogen_samples.csv', index=False)
    df_table = pd.DataFrame([d2, d1, df3])
    df_table=df_table.transpose().reset_index()
    df_table.columns=['Characteristic', 'derivation set samples (unique accessions)', 'Derivation set sequences', 'Derivation set sequences (QC passed)']

    return df_table, df[['Run', 'barcode', 'test_type']]

def table_1_characteristics(df, IORD, dfd_table):
    l=[]
    df['sample_name'] = df['seq_name'].astype(str)#.astype(str)
    df['sample_name'] = df['sample_name'].str.replace('.0', '', regex=False)
    #print(df['sample_name'].unique())
    #print(IORD['sample_name'].unique())
    df=df.merge(IORD, on=['sample_name'], how='left') 
    df['test normalized'] = df['test'].map(test_type_normalisation)
    df_not_tested=df[~df['test normalized'].isin(test_type_normalisation.values())]
    df_not_tested.to_csv('not_tested_samples_table1.csv', index=False)
    df['num_pathogens_detected'] = df.groupby(['Run', 'barcode'])['gold_standard'].transform('sum')
    # Count number of samples notna for collection_date and extraction_date
    df['collection_date_notna'] = df['collection_date'].notna()
    df['extraction_date_notna'] = df['extraction_date'].notna()
    collection_date_count = df[df['test normalized'].isin(test_type_normalisation.values())].drop_duplicates(subset=['Run', 'barcode'])['collection_date_notna'].sum()
    extraction_date_count = df[df['test normalized'].isin(test_type_normalisation.values())].drop_duplicates(subset=['Run', 'barcode'])['extraction_date_notna'].sum()
    df['days_between'] = (pd.to_datetime(df['extraction_date']) - pd.to_datetime(df['collection_date'])).dt.days
    df['days_between_notna'] = df['days_between'].notna()
    days_between_count = df[df['test normalized'].isin(test_type_normalisation.values())].drop_duplicates(subset=['Run', 'barcode'])['days_between_notna'].sum()
    days_between_median = df[df['test normalized'].isin(test_type_normalisation.values())].drop_duplicates(subset=['Run', 'barcode'])['days_between'].median()
    days_between_75, days_between_25 = np.percentile(df[(df['test normalized'].isin(test_type_normalisation.values()) & (df['days_between_notna']==True) )].drop_duplicates(subset=['Run', 'barcode'])['days_between'], [75, 25])
    days_between_median_IQR_N = f'{days_between_median} (IQR: {days_between_25:.0f} - {days_between_75:.0f}) n={days_between_count}'

    # Any target pathogen reads identified = sum of all the samples with any reads for target pathogens (median IQR)
    df['reads_to_any_reference']=df[df['pathogen']!='unmapped'].groupby(['Run', 'barcode'])['sample num reads'].transform('sum')
    median_atpri = df[(df['test normalized'].isin(test_type_normalisation.values())) ].drop_duplicates(subset=['Run', 'barcode'])['reads_to_any_reference'].median()
    q75, q25 = np.percentile(df[(df['test normalized'].isin(test_type_normalisation.values()))].drop_duplicates(subset=['Run', 'barcode'])['reads_to_any_reference'],  [75 ,25])
    atpri_N = df[(df['test normalized'].isin(test_type_normalisation.values()))].drop_duplicates(subset=['Run', 'barcode']).shape[0]
    atpri_IQR = f'{q25} - {q75}'
    Any_target_pathogen_reads_identified = f'{median_atpri:.0f} (IQR: {q25:.0f} - {q75:.0f}) n={atpri_N:.0f}'

    # same again but only for samples with reads_to_any_reference > 0
    median_atpri_gs = df[(df['test normalized'].isin(test_type_normalisation.values())) & (df['reads_to_any_reference']>0)]['reads_to_any_reference'].median()
    q75_gs, q25_gs = np.percentile(df[(df['test normalized'].isin(test_type_normalisation.values())) & (df['reads_to_any_reference']>0)]['reads_to_any_reference'],  [75 ,25])
    atpri_N_gs = df[(df['test normalized'].isin(test_type_normalisation.values())) & (df['reads_to_any_reference']>0)].drop_duplicates(subset=['Run', 'barcode']).shape[0]
    atpri_IQR_gs = f'{q25_gs} - {q75_gs}'
    Any_target_pathogen_reads_identified_gs = f'{median_atpri_gs:.0f} (IQR: {q25_gs:.0f} - {q75_gs:.0f}) n={atpri_N_gs:.0f}'

    # Gold standard pathogen reads identified 
    median_gspri = df[(df['test normalized'].isin(test_type_normalisation.values())) & (df['gold_standard']==1) ]['sample num reads'].median()
    q75, q25 = np.percentile(df[(df['test normalized'].isin(test_type_normalisation.values())) & (df['gold_standard']==1) ]['sample num reads'],  [75 ,25])
    gspri_N = df[(df['test normalized'].isin(test_type_normalisation.values())) & (df['gold_standard']==1)].shape[0]
    gspri_IQR = f'{q25} - {q75}'
    Gold_standard_pathogen_reads_identified = f'{median_gspri:.0f} (IQR: {q25:.0f} - {q75:.0f}) n={gspri_N:.0f}'

    # Gold standard pathogen reads identified where reads > 0
    median_gspri_gs = df[(df['test normalized'].isin(test_type_normalisation.values())) & (df['gold_standard']==1) & (df['sample num reads']>0) ]['sample num reads'].median()
    q75_gs, q25_gs = np.percentile(df[(df['test normalized'].isin(test_type_normalisation.values())) & (df['gold_standard']==1) & (df['sample num reads']>0) ]['sample num reads'],  [75 ,25])
    gspri_N_gs = df[(df['test normalized'].isin(test_type_normalisation.values())) & (df['gold_standard']==1) & (df['sample num reads']>0)].shape[0]
    gspri_IQR_gs = f'{q25_gs} - {q75_gs}'
    Gold_standard_pathogen_reads_identified_gs = f'{median_gspri_gs:.0f} (IQR: {q25_gs:.0f} - {q75_gs:.0f}) n={gspri_N_gs:.0f}'

    # CT values available from clinical PCR
    clinical_CT_values_N = df[(df['routine_CT_value'].notna()) & (df['routine_CT_value']!=0) & (df['test normalized']=='ALINITY')].drop_duplicates(subset=['Run', 'barcode']).shape[0]
    clinical_CT_values_median = df[(df['routine_CT_value'].notna()) & (df['routine_CT_value']!=0) & (df['test normalized']=='ALINITY')].drop_duplicates(subset=['Run', 'barcode'])['routine_CT_value'].median()
    q75, q25 = np.percentile(df[(df['routine_CT_value'].notna()) & (df['routine_CT_value']!=0) & (df['test normalized']=='ALINITY')].drop_duplicates(subset=['Run', 'barcode'])['routine_CT_value'],  [75 ,25])
    clinical_CT_values_IQR = f'{q25:.2f} - {q75:.2f}'
    clinical_CT_values = f'{clinical_CT_values_median:.2f} (IQR: {clinical_CT_values_IQR}) n={clinical_CT_values_N}'

    clinical_CT_values_N_cepheid = df[(df['routine_CT_value'].notna()) & (df['routine_CT_value']!=0) & (df['test normalized']=='CEPHEID')].drop_duplicates(subset=['Run', 'barcode']).shape[0]
    clinical_CT_values_median_cepheid = df[(df['routine_CT_value'].notna()) & (df['routine_CT_value']!=0) & (df['test normalized']=='CEPHEID')].drop_duplicates(subset=['Run', 'barcode'])['routine_CT_value'].median()
    q75_c, q25_c = np.percentile(df[(df['routine_CT_value'].notna()) & (df['routine_CT_value']!=0) & (df['test normalized']=='CEPHEID')].drop_duplicates(subset=['Run', 'barcode'])['routine_CT_value'],  [75 ,25])
    clinical_CT_values_IQR_cepheid = f'{q25_c:.2f} - {q75_c:.2f}'
    clinical_CT_values_cepheid = f'{clinical_CT_values_median_cepheid:.2f} (IQR: {clinical_CT_values_IQR_cepheid}) n={clinical_CT_values_N_cepheid}'

    # CT values available from study qPCR
    study_CT_values_N = df[(df['qpcr_ct'].notna()) & (df['qpcr_ct']!=0)].drop_duplicates(subset=['Run', 'barcode']).shape[0]
    undetermined_count = df[(df['qpcr_ct']=='UNDETERMINED')].shape[0]
    df['qpcr_ct2'] = df['qpcr_ct'].replace('UNDETERMINED', np.nan)
    df['qpcr_ct2'] = df['qpcr_ct2'].astype(float)
    study_CT_values_median = df[(df['qpcr_ct2'].notna()) & (df['qpcr_ct2']!=0)].drop_duplicates(subset=['Run', 'barcode'])['qpcr_ct2'].median()
    q75, q25 = np.percentile(df[(df['qpcr_ct2'].notna()) & (df['qpcr_ct2']!=0)].drop_duplicates(subset=['Run', 'barcode'])['qpcr_ct2'],  [75 ,25])
    study_CT_values_IQR = f'{q25:.2f} - {q75:.2f}'
    study_CT_values = f'{study_CT_values_median:.2f} (IQR: {study_CT_values_IQR}) n={study_CT_values_N} U={undetermined_count}'

    # Ct values available from study qPCR per target pathogen
    study_CT_values_N = df[(df['qpcr_ct'].notna()) & (df['qpcr_ct']!=0)].shape[0]
    missing_CT= df[(df['qpcr_ct'].isna()) & (df['gold_standard']==1)].shape[0]
    undetermined_count = df[(df['qpcr_ct']=='UNDETERMINED')].shape[0]
    study_CT_values_median = df[(df['qpcr_ct2'].notna()) & (df['qpcr_ct2']!=0)]['qpcr_ct2'].median()
    q75, q25 = np.percentile(df[(df['qpcr_ct2'].notna()) & (df['qpcr_ct2']!=0)]['qpcr_ct2'],  [75 ,25])
    study_CT_values_IQR = f'{q25:.2f} - {q75:.2f}'
    study_CT_values_per_pathogen = f'{study_CT_values_median:.2f} (IQR: {study_CT_values_IQR}) n={study_CT_values_N} U={undetermined_count} M={missing_CT}'

    # Age median and IQR
    df['AgeInYearsAtCollection']=df['AgeInYearsAtCollection'].replace('NULL',np.nan)
    age_median = df[(df['AgeInYearsAtCollection'].notna()) & df['test normalized'].isin(test_type_normalisation.values())].drop_duplicates(subset=['Run', 'barcode'])['AgeInYearsAtCollection'].median()
    q75, q25 = np.percentile(df[(df['AgeInYearsAtCollection'].notna()) & df['test normalized'].isin(test_type_normalisation.values())].drop_duplicates(subset=['Run', 'barcode'])['AgeInYearsAtCollection'],  [75 ,25])
    age_IQR = f'{q25:.1f} - {q75:.1f}'
    age_N = df[(df['AgeInYearsAtCollection'].notna()) & df['test normalized'].isin(test_type_normalisation.values())].drop_duplicates(subset=['Run', 'barcode']).shape[0]
    age_summary = f'{age_median:.1f} (IQR: {age_IQR}) n={age_N}'

    # Patient sex split
    F=df[(df['LinkedSex']=='F') & (df['test normalized'].isin(test_type_normalisation.values()))].drop_duplicates(subset=['Run', 'barcode']).shape[0]
    M=df[(df['LinkedSex']=='M') & (df['test normalized'].isin(test_type_normalisation.values()))].drop_duplicates(subset=['Run', 'barcode']).shape[0]
    U=df[(df['LinkedSex']=='U') & (df['test normalized'].isin(test_type_normalisation.values()))].drop_duplicates(subset=['Run', 'barcode']).shape[0]
    patient_sex = f'F: {F}, M: {M}, U: {U}, n={F+M+U}'

    # location tablulation
    location_counts = df[df['test normalized'].isin(test_type_normalisation.values())].drop_duplicates(subset=['Run', 'barcode'])['LastKnownLocation_micro'].shape[0]
    locationDesc = df[df['test normalized'].isin(test_type_normalisation.values())].drop_duplicates(subset=['Run', 'barcode']).groupby(['LastKnownLocation_micro'])['LastKnownLocation_micro'].count()
    locationDesc = pd.DataFrame({'LastKnownLocation_micro':locationDesc.index, 'Count':locationDesc.values})
    locationDesc.sort_values(by='Count',ascending=False, inplace=True)
    #locationDesc.to_csv('table_1_tabs/location_description.csv')

    # ward tabulation
    df['WardName_linked_inpat']=df['WardName_linked_inpat'].replace('NULL',np.nan)
    ward_counts = df[df['test normalized'].isin(test_type_normalisation.values()) & (df['WardName_linked_inpat'].notna())].drop_duplicates(subset=['Run', 'barcode'])['WardName_linked_inpat'].shape[0]
    wardDesc = df[df['test normalized'].isin(test_type_normalisation.values())].drop_duplicates(subset=['Run', 'barcode']).groupby(['WardName_linked_inpat'])['WardName_linked_inpat'].count()
    wardDesc = pd.DataFrame({'WardName_linked_inpat':wardDesc.index, 'Count':wardDesc.values})
    wardDesc.sort_values(by='Count',ascending=False, inplace=True)
    wardDesc['Ward_generalized']=wardDesc['WardName_linked_inpat'].map(ward_generalisation)
    wardDesc.to_csv('table_1_tabs/ward_description.csv')

    # Ward generalized
    df['Ward']=df['WardName_linked_inpat'].map(ward_generalisation)


    # facility tabulation
    facility_counts = df[df['test normalized'].isin(test_type_normalisation.values())].drop_duplicates(subset=['Run', 'barcode'])['Facility_linked_inpat'].shape[0]
    facilityDesc = df[df['test normalized'].isin(test_type_normalisation.values())].drop_duplicates(subset=['Run', 'barcode']).groupby(['Facility_linked_inpat'])['Facility_linked_inpat'].count()
    facilityDesc = pd.DataFrame({'Facility_linked_inpat':facilityDesc.index, 'Count':facilityDesc.values})
    facilityDesc.sort_values(by='Count',ascending=False, inplace=True)
    #facilityDesc.to_csv('table_1_tabs/facility_description.csv')


    # calculate bases and reads median and IQR
    bases=df.drop_duplicates(subset=['Run', 'barcode'])['total run bases'].median()
    gb=bases/1_000_000_000
    q75, q25 = np.percentile(df.drop_duplicates(subset=['Run', 'barcode'])['total run bases'],  [75 ,25])
    bases_IQR = f'{q25/1_000_000_000:.2f} - {q75/1_000_000_000:.2f}'
    reads=df.drop_duplicates(subset=['Run', 'barcode'])['total run reads'].median()
    q75, q25 = np.percentile(df.drop_duplicates(subset=['Run', 'barcode'])['total run reads'],  [75 ,25])
    reads_IQR = f'{q25/1_000_000:.3f} - {q75/1_000_000:.3f}'

    # Count number of FPs with single read
    df['PT_single_read_FP'] = np.where((df['test normalized'].isin(test_type_normalisation.values())) & (df['pathogen']!='unmapped') & (df['gold_standard']==0) & (df['sample num reads']==1), 1, 0)
    df['PS_single_read_FPs'] = df.groupby(['Run', 'barcode'])['PT_single_read_FP'].transform('sum')
    df['PS_single_read_FP'] = np.where(df['PS_single_read_FPs']>0, 1, 0)
    # Count number of FPs with two or more reads
    df['PT_twoMore_reads_FP'] = np.where((df['test normalized'].isin(test_type_normalisation.values())) & (df['pathogen']!='unmapped') & (df['gold_standard']==0) & (df['sample num reads']>1), 1, 0)
    df['PS_twoMore_reads_FPs'] = df.groupby(['Run', 'barcode'])['PT_twoMore_reads_FP'].transform('sum')
    df['PS_twoMore_reads_FP'] = np.where(df['PS_twoMore_reads_FPs']>0, 1, 0)
    # calculate per sample FP rates
    num_fps_single_read = df[df['PS_single_read_FP']>0].drop_duplicates(subset=['Run', 'barcode'])['PS_single_read_FPs'].shape[0]
    num_fps_twoMore_reads = df[df['PS_twoMore_reads_FP']>0].drop_duplicates(subset=['Run', 'barcode'])['PS_twoMore_reads_FPs'].shape[0]
    n_denominator = df[(df['test normalized'].isin(test_type_normalisation.values()))].drop_duplicates(subset=['Run', 'barcode']).shape[0]
    fp_reads=f'1={num_fps_single_read}, 2+={num_fps_twoMore_reads}, n={n_denominator}'

    # Count number of TPs with single read
    df['PT_single_read_TP'] = np.where((df['test normalized'].isin(test_type_normalisation.values())) & (df['pathogen']!='unmapped') & (df['gold_standard']==1) & (df['sample num reads']==1), 1, 0)
    df['PS_single_read_TPs'] = df.groupby(['Run', 'barcode'])['PT_single_read_TP'].transform('sum')
    df['PS_single_read_TP'] = np.where(df['PS_single_read_TPs']>0, 1, 0)
    # Count number of TPs with two or more reads
    df['PT_twoMore_reads_TP'] = np.where((df['test normalized'].isin(test_type_normalisation.values())) & (df['pathogen']!='unmapped') & (df['gold_standard']==1) & (df['sample num reads']>1), 1, 0)
    df['PS_twoMore_reads_TPs'] = df.groupby(['Run', 'barcode'])['PT_twoMore_reads_TP'].transform('sum')
    df['PS_twoMore_reads_TP'] = np.where(df['PS_twoMore_reads_TPs']>0, 1, 0)
    # calculate per sample TP rates
    num_tps_single_read = df[df['PS_single_read_TP']>0].drop_duplicates(subset=['Run', 'barcode'])['PS_single_read_TPs'].shape[0]
    num_tps_twoMore_reads = df[df['PS_twoMore_reads_TP']>0].drop_duplicates(subset=['Run', 'barcode'])['PS_twoMore_reads_TPs'].shape[0]
    n_denominator = df[(df['test normalized'].isin(test_type_normalisation.values()))].drop_duplicates(subset=['Run', 'barcode']).shape[0]
    tp_reads=f'1={num_tps_single_read}, 2+={num_tps_twoMore_reads}, n={n_denominator}'

    # alinity samples
    n=df[df['test normalized']=='ALINITY'].drop_duplicates(subset=['Run', 'barcode']).shape[0]
    tp_alinity=df[(df['test normalized']=='ALINITY') & (df['gold_standard']==1)].shape[0]
    neg_alinity=df[(df['test normalized']=='ALINITY') & (df['sample_positive']==0)].drop_duplicates(subset=['Run', 'barcode']).shape[0]
    alinity_summary=f'n={n}, negs={neg_alinity}, p={tp_alinity}'

    # biofire samples
    n=df[df['test normalized']=='BIOFIRE'].drop_duplicates(subset=['Run', 'barcode']).shape[0]
    tp_biofire=df[(df['test normalized']=='BIOFIRE') & (df['gold_standard']==1)].shape[0]
    neg_biofire=df[(df['test normalized']=='BIOFIRE') & (df['sample_positive']==0)].drop_duplicates(subset=['Run', 'barcode']).shape[0]
    biofire_summary=f'n={n}, negs={neg_biofire}, p={tp_biofire}'

    # cepheid samples
    n=df[df['test normalized']=='CEPHEID'].drop_duplicates(subset=['Run', 'barcode']).shape[0]
    tp_cepheid=df[(df['test normalized']=='CEPHEID') & (df['gold_standard']==1)].shape[0]
    neg_cepheid=df[(df['test normalized']=='CEPHEID') & (df['sample_positive']==0)].drop_duplicates(subset=['Run', 'barcode']).shape[0]
    cepheid_summary=f'n={n}, negs={neg_cepheid}, p={tp_cepheid}'

    # samples with no reads or Ct == UNDETERMINED but gold standard positive
    no_reads=df[(df['gold_standard']==1) & (df['sample num reads']==0) & (df['test normalized'].isin(test_type_normalisation.values()))].shape[0]
    no_ct=df[(df['qpcr_ct'].isin(['UNDETERMINED', np.nan])) & (df['gold_standard']==1) & (df['test normalized'].isin(test_type_normalisation.values()))].shape[0]
    no_reads_or_undetermined_Ct_positive = df[(df['qpcr_ct'].isin(['UNDETERMINED', np.nan])) & (df['gold_standard']==1) & (df['sample num reads']==0) & (df['test normalized'].isin(test_type_normalisation.values()))].shape[0]
    total_pathogens_detected = df[(df['test normalized'].isin(test_type_normalisation.values()))].drop_duplicates(subset=['Run', 'barcode'])['num_pathogens_detected'].sum()
    no_reads_or_undetermined_Ct_positive_summary = f'no reads={no_reads}, no Ct={no_ct}, n={no_reads_or_undetermined_Ct_positive}({no_reads_or_undetermined_Ct_positive / total_pathogens_detected * 100}%)'

    d1={
        'Set': 'validation set samples',
        'Total sequences': df.drop_duplicates(subset=['Run', 'barcode']).shape[0],
        'total samples': df[df['test normalized'].isin(test_type_normalisation.values())].drop_duplicates(subset=['Run', 'barcode']).shape[0],
        'Alinity m Respiratory Panel samples': alinity_summary,
        'BioFire Respiratory Panel 2.1 samples': biofire_summary,
        'Cepheid Xpert Xpress SARS-CoV-2/Flu/RSV samples': cepheid_summary,
        '0 pathogens': df[(df['num_pathogens_detected']==0) & (df['test normalized'].isin(test_type_normalisation.values()))].drop_duplicates(subset=['Run', 'barcode']).shape[0],
        '1 pathogens': df[(df['num_pathogens_detected']==1) & (df['test normalized'].isin(test_type_normalisation.values()))].drop_duplicates(subset=['Run', 'barcode']).shape[0],
        '2 pathogens': df[(df['num_pathogens_detected']==2) & (df['test normalized'].isin(test_type_normalisation.values()))].drop_duplicates(subset=['Run', 'barcode']).shape[0],
        '3 pathogens': df[(df['num_pathogens_detected']==3) & (df['test normalized'].isin(test_type_normalisation.values()))].drop_duplicates(subset=['Run', 'barcode']).shape[0],
        'Total pathogens detected': df[(df['test normalized'].isin(test_type_normalisation.values()))].drop_duplicates(subset=['Run', 'barcode'])['num_pathogens_detected'].sum(),
        'Samples with collection date': collection_date_count,
        'Samples with extraction date': extraction_date_count,
        'Samples with days between collection and extraction date': days_between_median_IQR_N,
        'Ct (Alinity) available': clinical_CT_values,
        'Ct (Xpert) available': clinical_CT_values_cepheid,
        'Ct (study qPCR) available': study_CT_values,
        'Ct (study qPCR) available per target pathogen': study_CT_values_per_pathogen,
        #'Any target pathogen reads identified': df[(df['gold_standard']==1) & (df['sample num reads']>0)].drop_duplicates(subset=['Run', 'barcode']).shape[0],
        'Any target pathogen reads identified': Any_target_pathogen_reads_identified,
        #'Number of reads (if>0)': df[(df['sample num reads']>0) & (df['test normalized'].isin(test_type_normalisation.values())) & (df['pathogen']!='unmapped')].shape[0],
        'Any target pathogen reads identified (if reads>0)': Any_target_pathogen_reads_identified_gs,
        '(target) Gold standard pathogen reads identified': Gold_standard_pathogen_reads_identified,
        '(target) Gold standard pathogen reads identified (if reads>0)': Gold_standard_pathogen_reads_identified_gs,
        'Samples with patient age': age_summary,
        'Age range <1': len(df[(df['AgeInYearsAtCollection']<1) & (df['test normalized'].isin(test_type_normalisation.values()))].drop_duplicates(subset=['Run', 'barcode'])),
        'Age range 1-12': len(df[(df['AgeInYearsAtCollection']>=1) & (df['AgeInYearsAtCollection']<=12) & (df['test normalized'].isin(test_type_normalisation.values()))].drop_duplicates(subset=['Run', 'barcode'])),
        'Age range 13-49': len(df[(df['AgeInYearsAtCollection']>=13) & (df['AgeInYearsAtCollection']<=49) & (df['test normalized'].isin(test_type_normalisation.values()))].drop_duplicates(subset=['Run', 'barcode'])), 
        'Age range 50-70': len(df[(df['AgeInYearsAtCollection']>=50) & (df['AgeInYearsAtCollection']<=70) & (df['test normalized'].isin(test_type_normalisation.values()))].drop_duplicates(subset=['Run', 'barcode'])),
        'Age range 71-84': len(df[(df['AgeInYearsAtCollection']>=71) & (df['AgeInYearsAtCollection']<=84) & (df['test normalized'].isin(test_type_normalisation.values()))].drop_duplicates(subset=['Run', 'barcode'])),
        'Age range 85+': len(df[(df['AgeInYearsAtCollection']>=85) & (df['test normalized'].isin(test_type_normalisation.values()))].drop_duplicates(subset=['Run', 'barcode'])),
        'Samples with patient sex': patient_sex,
        'Samples with patient location': location_counts,
        'Samples with patient ward': ward_counts,
        'Haematology_oncology': df[df['Ward']=='Haematology_oncology'].drop_duplicates(subset=['Run', 'barcode']).shape[0],
        'Other_inpatient_adult': df[df['Ward']=='Other_inpatient_adult'].drop_duplicates(subset=['Run', 'barcode']).shape[0],
        'Other_inpatient_paeds': df[df['Ward']=='Other_inpatient_paeds'].drop_duplicates(subset=['Run', 'barcode']).shape[0],
        'Ambulatory_and_ED': df[df['Ward']=='Ambulatory_and_ED'].drop_duplicates(subset=['Run', 'barcode']).shape[0],
        'Critical_care': df[df['Ward']=='Critical_care'].drop_duplicates(subset=['Run', 'barcode']).shape[0],        
        'Samples with patient facility': facility_counts,
        'Samples with sample type NTS': df[df['sample_type']=='NTS'].drop_duplicates(subset=['Run', 'barcode']).shape[0],
        'Samples with sample type NS': df[df['sample_type']=='NS'].drop_duplicates(subset=['Run', 'barcode']).shape[0],
        'Gigabases': f'{gb:.2f} (IQR: {bases_IQR})',
        'Million Reads': f'{reads/1_000_000:.3f} (IQR: {reads_IQR})',
        '(Per sample) False positives by reads': fp_reads,
        '(Per sample) True positives by reads': tp_reads,
        '(per target pathogen) Samples with no reads or Ct UNDETERMINED but gold standard positive': no_reads_or_undetermined_Ct_positive_summary
    }

    # remove QC failed samples and repeat
    df2=df[(df['PCs_passed']==True) & (df['run_pass']==True)]
    collection_date_count = df2[df2['test normalized'].isin(test_type_normalisation.values())].drop_duplicates(subset=['Run', 'barcode'])['collection_date_notna'].sum()
    extraction_date_count = df2[df2['test normalized'].isin(test_type_normalisation.values())].drop_duplicates(subset=['Run', 'barcode'])['extraction_date_notna'].sum()
    days_between_count = df2[df2['test normalized'].isin(test_type_normalisation.values())].drop_duplicates(subset=['Run', 'barcode'])['days_between_notna'].sum()
    days_between_median = df2[df2['test normalized'].isin(test_type_normalisation.values())].drop_duplicates(subset=['Run', 'barcode'])['days_between'].median()
    days_between_25, days_between_75 = np.percentile(df2[(df2['test normalized'].isin(test_type_normalisation.values()) & (df2['days_between_notna']==True))].drop_duplicates(subset=['Run', 'barcode'])['days_between'], [25, 75])
    days_between_median_IQR_N = f'{days_between_median} (IQR: {days_between_25:.0f} - {days_between_75:.0f}) n={days_between_count}'

    # Any target pathogen reads identified = sum of all the samples with any reads for target pathogens (median IQR)
    median_atpri = df2[(df2['test normalized'].isin(test_type_normalisation.values()))]['reads_to_any_reference'].median()
    q75, q25 = np.percentile(df2[(df2['test normalized'].isin(test_type_normalisation.values())) & (df2['reads_to_any_reference']>0)]['reads_to_any_reference'],  [75 ,25])
    atpri_N = df2[(df2['test normalized'].isin(test_type_normalisation.values()))].drop_duplicates(subset=['Run', 'barcode']).shape[0]
    atpri_IQR = f'{q25:.0f} - {q75:.0f}'
    Any_target_pathogen_reads_identified = f'{median_atpri:.0f} (IQR: {atpri_IQR}) n={atpri_N:.0f}'

    # Same again but only for samples with reads_to_any_reference > 0
    median_atpri_gs = df2[(df2['test normalized'].isin(test_type_normalisation.values())) & (df2['reads_to_any_reference']>0)]['reads_to_any_reference'].median()
    q75_gs, q25_gs = np.percentile(df2[(df2['test normalized'].isin(test_type_normalisation.values())) & (df2['reads_to_any_reference']>0)]['reads_to_any_reference'],  [75 ,25])
    atpri_N_gs = df2[(df2['test normalized'].isin(test_type_normalisation.values())) & (df2['reads_to_any_reference']>0)].drop_duplicates(subset=['Run', 'barcode']).shape[0]
    atpri_IQR_gs = f'{q25_gs:.0f} - {q75_gs:.0f}'
    Any_target_pathogen_reads_identified_gs = f'{median_atpri_gs:.0f} (IQR: {atpri_IQR_gs}) n={atpri_N_gs:.0f}'

    # Gold standard pathogen reads identified 
    median_gspri = df2[(df2['test normalized'].isin(test_type_normalisation.values())) & (df2['gold_standard']==1) ]['sample num reads'].median()
    q75, q25 = np.percentile(df2[(df2['test normalized'].isin(test_type_normalisation.values())) & (df2['gold_standard']==1) ]['sample num reads'],  [75 ,25])
    gspri_N = df2[(df2['test normalized'].isin(test_type_normalisation.values())) & (df2['gold_standard']==1)].shape[0]
    gspri_IQR = f'{q25} - {q75}'
    Gold_standard_pathogen_reads_identified = f'{median_gspri:.0f} (IQR: {q25:.0f} - {q75:.0f}) n={gspri_N:.0f}'

    # Gold standard pathogen reads identified where reads > 0
    median_gspri_gs = df2[(df2['test normalized'].isin(test_type_normalisation.values())) & (df2['gold_standard']==1) & (df2['sample num reads']>0) ]['sample num reads'].median()
    q75_gs, q25_gs = np.percentile(df2[(df2['test normalized'].isin(test_type_normalisation.values())) & (df2['gold_standard']==1) & (df2['sample num reads']>0) ]['sample num reads'],  [75 ,25])
    gspri_N_gs = df2[(df2['test normalized'].isin(test_type_normalisation.values())) & (df2['gold_standard']==1) & (df2['sample num reads']>0)].shape[0]
    gspri_IQR_gs = f'{q25_gs} - {q75_gs}'
    Gold_standard_pathogen_reads_identified_gs = f'{median_gspri_gs:.0f} (IQR: {q25_gs:.0f} - {q75_gs:.0f}) n={gspri_N_gs:.0f}'

    # Ct values available from clinical PCR
    clinical_CT_values_N = df2[(df2['routine_CT_value'].notna()) & (df2['routine_CT_value']!=0) & (df2['test normalized']=='ALINITY')].drop_duplicates(subset=['Run', 'barcode']).shape[0]
    clinical_CT_values_median = df2[(df2['routine_CT_value'].notna()) & (df2['routine_CT_value']!=0) & (df2['test normalized']=='ALINITY')].drop_duplicates(subset=['Run', 'barcode'])['routine_CT_value'].median()
    q75, q25 = np.percentile(df2[(df2['routine_CT_value'].notna()) & (df2['routine_CT_value']!=0) & (df2['test normalized']=='ALINITY')].drop_duplicates(subset=['Run', 'barcode'])['routine_CT_value'],  [75 ,25])
    clinical_CT_values_IQR = f'{q25:.2f} - {q75:.2f}'
    clinical_CT_values = f'{clinical_CT_values_median:.2f} (IQR: {clinical_CT_values_IQR}) n={clinical_CT_values_N}'

    clinical_CT_values_N_cepheid = df2[(df2['routine_CT_value'].notna()) & (df2['routine_CT_value']!=0) & (df2['test normalized']=='CEPHEID')].drop_duplicates(subset=['Run', 'barcode']).shape[0]
    clinical_CT_values_median_cepheid = df2[(df2['routine_CT_value'].notna()) & (df2['routine_CT_value']!=0) & (df2['test normalized']=='CEPHEID')].drop_duplicates(subset=['Run', 'barcode'])['routine_CT_value'].median()
    q75_c, q25_c = np.percentile(df2[(df2['routine_CT_value'].notna()) & (df2['routine_CT_value']!=0) & (df2['test normalized']=='CEPHEID')].drop_duplicates(subset=['Run', 'barcode'])['routine_CT_value'],  [75 ,25])
    clinical_CT_values_IQR_cepheid = f'{q25_c:.2f} - {q75_c:.2f}'
    clinical_CT_values_cepheid = f'{clinical_CT_values_median_cepheid:.2f} (IQR: {clinical_CT_values_IQR_cepheid}) n={clinical_CT_values_N_cepheid}'

    # Ct values available from study qPCR
    study_CT_values_N = df2[(df2['qpcr_ct'].notna()) & (df2['qpcr_ct']!=0)].drop_duplicates(subset=['Run', 'barcode']).shape[0]
    undetermined_count = df2[(df2['qpcr_ct']=='UNDETERMINED')].shape[0]
    study_CT_values_median = df2[(df2['qpcr_ct2'].notna()) & (df2['qpcr_ct2']!=0)].drop_duplicates(subset=['Run', 'barcode'])['qpcr_ct2'].median()
    q75, q25 = np.percentile(df2[(df2['qpcr_ct2'].notna()) & (df2['qpcr_ct2']!=0)].drop_duplicates(subset=['Run', 'barcode'])['qpcr_ct2'],  [75 ,25])
    study_CT_values_IQR = f'{q25:.2f} - {q75:.2f}'
    study_CT_values = f'{study_CT_values_median:.2f} (IQR: {study_CT_values_IQR}) n={study_CT_values_N} U={undetermined_count}'

    # Ct values available from study qPCR per target pathogen
    study_CT_values_N = df2[(df2['qpcr_ct'].notna()) & (df2['qpcr_ct']!=0)].shape[0]
    missing_CT= df2[(df2['qpcr_ct'].isna()) & (df2['gold_standard']==1)].shape[0]
    undetermined_count = df2[(df2['qpcr_ct']=='UNDETERMINED')].shape[0]
    study_CT_values_median = df2[(df2['qpcr_ct2'].notna()) & (df2['qpcr_ct2']!=0)]['qpcr_ct2'].median()
    q75, q25 = np.percentile(df2[(df2['qpcr_ct2'].notna()) & (df2['qpcr_ct2']!=0)]['qpcr_ct2'],  [75 ,25])
    study_CT_values_IQR = f'{q25:.2f} - {q75:.2f}'
    study_CT_values_per_pathogen = f'{study_CT_values_median:.2f} (IQR: {study_CT_values_IQR}) n={study_CT_values_N} U={undetermined_count} M={missing_CT}'

    # Age median and IQR
    age_median = df2[(df2['AgeInYearsAtCollection'].notna()) & (df2['test normalized'].isin(test_type_normalisation.values()))].drop_duplicates(subset=['Run', 'barcode'])['AgeInYearsAtCollection'].median()
    q75, q25 = np.percentile(df2[(df2['AgeInYearsAtCollection'].notna()) & (df2['test normalized'].isin(test_type_normalisation.values()))].drop_duplicates(subset=['Run', 'barcode'])['AgeInYearsAtCollection'],  [75 , 25])
    age_IQR = f'{q25:.1f} - {q75:.1f}'
    age_N = df2[(df2['AgeInYearsAtCollection'].notna()) & (df2['test normalized'].isin(test_type_normalisation.values()))].drop_duplicates(subset=['Run', 'barcode']).shape[0]
    age_summary = f'{age_median:.1f} (IQR: {age_IQR}) n={age_N}'

    # Patient sex split
    F=df2[(df2['LinkedSex']=='F') & (df2['test normalized'].isin(test_type_normalisation.values()))].drop_duplicates(subset=['Run', 'barcode']).shape[0]
    M=df2[(df2['LinkedSex']=='M') & (df2['test normalized'].isin(test_type_normalisation.values()))].drop_duplicates(subset=['Run', 'barcode']).shape[0]
    U=df2[(df2['LinkedSex']=='U') & (df2['test normalized'].isin(test_type_normalisation.values()))].drop_duplicates(subset=['Run', 'barcode']).shape[0]
    patient_sex = f'F: {F}, M: {M}, U: {U}, n={F+M+U}'

    # location tablulation
    location_counts2 = df2[df2['test normalized'].isin(test_type_normalisation.values())].drop_duplicates(subset=['Run', 'barcode'])['LastKnownLocation_micro'].value_counts()
    locationDesc2 = df2[df2['test normalized'].isin(test_type_normalisation.values())].drop_duplicates(subset=['Run', 'barcode']).groupby(['LastKnownLocation_micro'])['LastKnownLocation_micro'].count()
    locationDesc2 = pd.DataFrame({'LastKnownLocation_micro':locationDesc2.index, 'Count pass QC':locationDesc2.values})
    locationDesc=locationDesc.merge(locationDesc2,on='LastKnownLocation_micro',how='left')
    locationDesc.sort_values(by='Count',ascending=False, inplace=True)
    locationDesc.to_csv('table_1_tabs/location_description.csv')

    # ward tabulation
    ward_counts = df2[df2['test normalized'].isin(test_type_normalisation.values())].drop_duplicates(subset=['Run', 'barcode'])['WardName_linked_inpat'].value_counts()
    wardDesc2 = df2[df2['test normalized'].isin(test_type_normalisation.values())].drop_duplicates(subset=['Run', 'barcode']).groupby(['WardName_linked_inpat'])['WardName_linked_inpat'].count()
    wardDesc2 = pd.DataFrame({'WardName_linked_inpat':wardDesc2.index, 'Count pass QC':wardDesc2.values})
    wardDesc=wardDesc.merge(wardDesc2,on='WardName_linked_inpat',how='left')
    wardDesc.sort_values(by='Count',ascending=False, inplace=True)
    wardDesc.to_csv('table_1_tabs/ward_description.csv')

    # facility tabulation
    facility_counts = df2[df2['test normalized'].isin(test_type_normalisation.values())].drop_duplicates(subset=['Run', 'barcode'])['Facility_linked_inpat'].value_counts()
    facilityDesc2 = df2[df2['test normalized'].isin(test_type_normalisation.values())].drop_duplicates(subset=['Run', 'barcode']).groupby(['Facility_linked_inpat'])['Facility_linked_inpat'].count()
    facilityDesc2 = pd.DataFrame({'Facility_linked_inpat':facilityDesc2.index, 'Count pass QC':facilityDesc2.values})
    facilityDesc=facilityDesc.merge(facilityDesc2,on='Facility_linked_inpat',how='left')
    facilityDesc.sort_values(by='Count',ascending=False, inplace=True)
    facilityDesc.to_csv('table_1_tabs/facility_description.csv')


    # calculate bases and reads median and IQR
    bases=df2.drop_duplicates(subset=['Run', 'barcode'])['total run bases'].median()
    gb=bases/1_000_000_000
    q75, q25 = np.percentile(df2.drop_duplicates(subset=['Run', 'barcode'])['total run bases'],  [75 ,25])
    bases_IQR = f'{q25/1_000_000_000:.2f} - {q75/1_000_000_000:.2f}'
    reads=df2.drop_duplicates(subset=['Run', 'barcode'])['total run reads'].median()
    q75, q25 = np.percentile(df2.drop_duplicates(subset=['Run', 'barcode'])['total run reads'],  [75 ,25])
    reads_IQR = f'{q25/1_000_000:.3f} - {q75/1_000_000:.3f}'

    # Count number of FPs with single read
    df2['PT_single_read_FP'] = np.where((df2['test normalized'].isin(test_type_normalisation.values())) & (df2['pathogen']!='unmapped') & (df2['gold_standard']==0) & (df2['sample num reads']==1), 1, 0)
    df2['PS_single_read_FPs'] = df2.groupby(['Run', 'barcode'])['PT_single_read_FP'].transform('sum')
    df2['PS_single_read_FP'] = np.where(df2['PS_single_read_FPs']>0, 1, 0)
    # Count number of FPs with two or more reads
    df2['PT_twoMore_reads_FP'] = np.where((df2['test normalized'].isin(test_type_normalisation.values())) & (df2['pathogen']!='unmapped') & (df2['gold_standard']==0) & (df2['sample num reads']>1), 1, 0)
    df2['PS_twoMore_reads_FPs'] = df2.groupby(['Run', 'barcode'])['PT_twoMore_reads_FP'].transform('sum')
    df2['PS_twoMore_reads_FP'] = np.where(df2['PS_twoMore_reads_FPs']>0, 1, 0)
    # calculate per sample FP rates
    num_fps_single_read = df2[df2['PS_single_read_FP']>0].drop_duplicates(subset=['Run', 'barcode'])['PS_single_read_FPs'].shape[0]
    num_fps_twoMore_reads = df2[df2['PS_twoMore_reads_FP']>0].drop_duplicates(subset=['Run', 'barcode'])['PS_twoMore_reads_FPs'].shape[0]
    n_denominator = df2[(df2['test normalized'].isin(test_type_normalisation.values()))].drop_duplicates(subset=['Run', 'barcode']).shape[0]
    fp_reads=f'1={num_fps_single_read}, 2+={num_fps_twoMore_reads}, n={n_denominator}'

    # Count number of TPs with single read
    df2['PT_single_read_TP'] = np.where((df2['test normalized'].isin(test_type_normalisation.values())) & (df2['pathogen']!='unmapped') & (df2['gold_standard']==1) & (df2['sample num reads']==1), 1, 0)
    df2['PS_single_read_TPs'] = df2.groupby(['Run', 'barcode'])['PT_single_read_TP'].transform('sum')
    df2['PS_single_read_TP'] = np.where(df2['PS_single_read_TPs']>0, 1, 0)
    # Count number of TPs with two or more reads
    df2['PT_twoMore_reads_TP'] = np.where((df2['test normalized'].isin(test_type_normalisation.values())) & (df2['pathogen']!='unmapped') & (df2['gold_standard']==1) & (df2['sample num reads']>1), 1, 0)
    df2['PS_twoMore_reads_TPs'] = df2.groupby(['Run', 'barcode'])['PT_twoMore_reads_TP'].transform('sum')
    df2['PS_twoMore_reads_TP'] = np.where(df2['PS_twoMore_reads_TPs']>0, 1, 0)
    # calculate per sample TP rates
    num_tps_single_read = df2[df2['PS_single_read_TP']>0].drop_duplicates(subset=['Run', 'barcode'])['PS_single_read_TPs'].shape[0]
    num_tps_twoMore_reads = df2[df2['PS_twoMore_reads_TP']>0].drop_duplicates(subset=['Run', 'barcode'])['PS_twoMore_reads_TPs'].shape[0]
    n_denominator = df2[(df2['test normalized'].isin(test_type_normalisation.values()))].drop_duplicates(subset=['Run', 'barcode']).shape[0]
    tp_reads=f'1={num_tps_single_read}, 2+={num_tps_twoMore_reads}, n={n_denominator}'
    
    # alinity samples
    n=df2[df2['test normalized']=='ALINITY'].drop_duplicates(subset=['Run', 'barcode']).shape[0]
    tp_alinity=df2[(df2['test normalized']=='ALINITY') & (df2['gold_standard']==1)].shape[0]
    neg_alinity=df2[(df2['test normalized']=='ALINITY') & (df2['sample_positive']==0)].drop_duplicates(subset=['Run', 'barcode']).shape[0]
    alinity_summary=f'n={n}, negs={neg_alinity}, p={tp_alinity}'

    # biofire samples
    n=df2[df2['test normalized']=='BIOFIRE'].drop_duplicates(subset=['Run', 'barcode']).shape[0]
    tp_biofire=df2[(df2['test normalized']=='BIOFIRE') & (df2['gold_standard']==1)].shape[0]
    neg_biofire=df2[(df2['test normalized']=='BIOFIRE') & (df2['sample_positive']==0)].drop_duplicates(subset=['Run', 'barcode']).shape[0]
    biofire_summary=f'n={n}, negs={neg_biofire}, p={tp_biofire}'

    # cepheid samples
    n=df2[df2['test normalized']=='CEPHEID'].drop_duplicates(subset=['Run', 'barcode']).shape[0]
    tp_cepheid=df2[(df2['test normalized']=='CEPHEID') & (df2['gold_standard']==1)].shape[0]
    neg_cepheid=df2[(df2['test normalized']=='CEPHEID') & (df2['sample_positive']==0)].drop_duplicates(subset=['Run', 'barcode']).shape[0]
    cepheid_summary=f'n={n}, negs={neg_cepheid}, p={tp_cepheid}' 

    # samples with no reads or Ct == UNDETERMINED but gold standard positive
    no_reads=df2[(df2['gold_standard']==1) & (df2['sample num reads']==0) & (df2['test normalized'].isin(test_type_normalisation.values()))].shape[0]
    no_ct=df2[(df2['qpcr_ct'].isin(['UNDETERMINED', np.nan])) & (df2['gold_standard']==1) & (df2['test normalized'].isin(test_type_normalisation.values()))].shape[0]
    no_reads_or_undetermined_Ct_positive = df2[(df2['qpcr_ct'].isin(['UNDETERMINED', np.nan])) & (df2['gold_standard']==1) & (df2['sample num reads']==0) & (df2['test normalized'].isin(test_type_normalisation.values()))].shape[0]
    total_pathogens_detected = df2[(df2['test normalized'].isin(test_type_normalisation.values()))].drop_duplicates(subset=['Run', 'barcode'])['num_pathogens_detected'].sum()
    no_reads_or_undetermined_Ct_positive_summary = f'no reads={no_reads}, no Ct={no_ct}, n={no_reads_or_undetermined_Ct_positive}({no_reads_or_undetermined_Ct_positive / total_pathogens_detected * 100}%)'
    

    d2={
        'Set': 'validation set samples after QC',
        'Total sequences': df2.drop_duplicates(subset=['Run', 'barcode']).shape[0],
        'total samples': df2[df2['test normalized'].isin(test_type_normalisation.values())].drop_duplicates(subset=['Run', 'barcode']).shape[0],
        'Alinity m Respiratory Panel samples': alinity_summary,
        'BioFire Respiratory Panel 2.1 samples': biofire_summary,
        'Cepheid Xpert Xpress SARS-CoV-2/Flu/RSV samples': cepheid_summary,
        '0 pathogens': df2[(df2['num_pathogens_detected']==0) & (df2['test normalized'].isin(test_type_normalisation.values()))].drop_duplicates(subset=['Run', 'barcode']).shape[0],
        '1 pathogens': df2[(df2['num_pathogens_detected']==1) & (df2['test normalized'].isin(test_type_normalisation.values()))].drop_duplicates(subset=['Run', 'barcode']).shape[0],
        '2 pathogens': df2[(df2['num_pathogens_detected']==2) & (df2['test normalized'].isin(test_type_normalisation.values()))].drop_duplicates(subset=['Run', 'barcode']).shape[0],
        '3 pathogens': df2[(df2['num_pathogens_detected']==3) & (df2['test normalized'].isin(test_type_normalisation.values()))].drop_duplicates(subset=['Run', 'barcode']).shape[0],
        'Total pathogens detected': df2[(df2['test normalized'].isin(test_type_normalisation.values()))].drop_duplicates(subset=['Run', 'barcode'])['num_pathogens_detected'].sum(),
        'Samples with collection date': collection_date_count,
        'Samples with extraction date': extraction_date_count,
        'Samples with days between collection and extraction date': days_between_median_IQR_N,
        'Ct (Alinity) available': clinical_CT_values,
        'Ct (Xpert) available': clinical_CT_values_cepheid,
        'Ct (study qPCR) available': study_CT_values,
        'Ct (study qPCR) available per target pathogen': study_CT_values_per_pathogen,
        #'Any target pathogen reads identified': df2[(df2['gold_standard']==1) & (df2['sample num reads']>0)].drop_duplicates(subset=['Run', 'barcode']).shape[0],
        'Any target pathogen reads identified': Any_target_pathogen_reads_identified,
        #'Number of reads (if>0)': df2[(df2['sample num reads']>0) & (df2['test normalized'].isin(test_type_normalisation.values()))& (df2['pathogen']!='unmapped')].shape[0],
        'Any target pathogen reads identified (if reads>0)': Any_target_pathogen_reads_identified_gs,
        '(target) Gold standard pathogen reads identified': Gold_standard_pathogen_reads_identified,
        '(target) Gold standard pathogen reads identified (if reads>0)': Gold_standard_pathogen_reads_identified_gs,
        'Samples with patient age': age_summary,
        'Age range <1': len(df2[(df2['AgeInYearsAtCollection']<1) & (df2['test normalized'].isin(test_type_normalisation.values()))].drop_duplicates(subset=['Run', 'barcode'])),
        'Age range 1-12': len(df2[(df2['AgeInYearsAtCollection']>=1) & (df2['AgeInYearsAtCollection']<=12) & (df2['test normalized'].isin(test_type_normalisation.values()))].drop_duplicates(subset=['Run', 'barcode'])),
        'Age range 13-49': len(df2[(df2['AgeInYearsAtCollection']>=13) & (df2['AgeInYearsAtCollection']<=49) & (df2['test normalized'].isin(test_type_normalisation.values  ()))].drop_duplicates(subset=['Run', 'barcode'])), 
        'Age range 50-70': len(df2[(df2['AgeInYearsAtCollection']>=50) & (df2['AgeInYearsAtCollection']<=70) & (df2['test normalized'].isin(test_type_normalisation.values()))].drop_duplicates(subset=['Run', 'barcode'])),
        'Age range 71-84': len(df2[(df2['AgeInYearsAtCollection']>=71) & (df2['AgeInYearsAtCollection']<=84) & (df2['test normalized'].isin(test_type_normalisation.values()))].drop_duplicates(subset=['Run', 'barcode'])),
        'Age range 85+': len(df2[(df2['AgeInYearsAtCollection']>=85) & (df2['test normalized'].isin(test_type_normalisation.values()))].drop_duplicates(subset=['Run', 'barcode'])),
        'Samples with patient sex': patient_sex,
        'Samples with patient location': df2[df2['LastKnownLocation_micro'].notna()].drop_duplicates(subset=['Run', 'barcode']).shape[0],
        'Samples with patient ward': df2[df2['WardName_linked_inpat'].notna()].drop_duplicates(subset=['Run', 'barcode']).shape[0],
        'Haematology_oncology': df2[df2['Ward']=='Haematology_oncology'].drop_duplicates(subset=['Run', 'barcode']).shape[0],
        'Other_inpatient_adult': df2[df2['Ward']=='Other_inpatient_adult'].drop_duplicates(subset=['Run', 'barcode']).shape[0],
        'Other_inpatient_paeds': df2[df2['Ward']=='Other_inpatient_paeds'].drop_duplicates(subset=['Run', 'barcode']).shape[0],
        'Ambulatory_and_ED': df2[df2['Ward']=='Ambulatory_and_ED'].drop_duplicates(subset=['Run', 'barcode']).shape[0],
        'Critical_care': df2[df2['Ward']=='Critical_care'].drop_duplicates(subset=['Run', 'barcode']).shape[0],
        'Samples with patient facility': df2[df2['Facility_linked_inpat'].notna()].drop_duplicates(subset=['Run', 'barcode']).shape[0],
        'Samples with sample type NTS': df2[df2['sample_type']=='NTS'].drop_duplicates(subset=['Run', 'barcode']).shape[0],
        'Samples with sample type NS': df2[df2['sample_type']=='NS'].drop_duplicates(subset=['Run', 'barcode']).shape[0],
        'Gigabases': f'{gb:.2f} (IQR: {bases_IQR})',
        'Million Reads': f'{reads/1_000_000:.3f} (IQR: {reads_IQR})',
        '(Per sample) False positives by reads': fp_reads,
        '(Per sample) True positives by reads': tp_reads,
        '(per target pathogen) Samples with no reads or Ct UNDETERMINED but gold standard positive': no_reads_or_undetermined_Ct_positive_summary
    }
    df_table=pd.DataFrame([d1, d2])
    # pivot table so that 'Set' is the index and all the other columns are columns
    df_table=df_table.transpose().reset_index()
    df_table.columns=['Characteristic', 'All validation set samples', 'Validation set samples after QC']
    df_table=df_table.merge(dfd_table, on='Characteristic', how='left')
    cols=['Characteristic', 'derivation set samples (unique accessions)', 'Derivation set sequences','Derivation set sequences (QC passed)', 
          'All validation set samples', 'Validation set samples after QC']
    df_table[cols].to_csv('table_1_sample_characteristics.csv', index=False)
    df.to_csv('table_1_input_data.csv', index=False)
    plot_ages(df)
                                                                    
def run_pass(df):
    # add column if run passes QC
    # spike negs
    # Add column if Total run bases < 0.4 million bases
    df['run_bases_pass']=np.where(df['total run bases']<400_000, False, True)
    # Add column if Total sample reads < 25,000 reads
    df['sample_reads_pass']=np.where(df['total sample reads']<25_000, False, True)
    # This is also called the batch negatice control in the flow diagram
    df_negs=df[df['amplification_control']==1]
    df_negs_pass=df_negs[(df_negs['pass']==True) | (df_negs['PCs_passed']==1)]
    failed_amp_batches=list(df_negs_pass.drop_duplicates(['Run','Batch'])[['Run','Batch']].itertuples(index=False, name=None))
    # Use apply to create tuple and check if in failed_amp_batches
    df['batch_neg_control_fail']=df[['Run','Batch']].apply(tuple, axis=1).isin(failed_amp_batches)

    # count number of samples that failed batch ampification negative controls
    df_RC_control=df[(df['reverse_transcription_control']==1) & (df['IC_virus_spike']==1)]
    # orthoreovirus passed	zika passed	MS2 passed	murine_respirovirus passed
    df_RC_control_PCFAIL=df_RC_control[(df_RC_control['orthoreovirus passed']==0) & (df_RC_control['zika passed']==0) & (df_RC_control['murine_respirovirus passed']==0) \
                                    | (df_RC_control['MS2 passed']==1) \
                                    | (df_RC_control['pass']==True) ]
    failed_batches=list(df_RC_control_PCFAIL.drop_duplicates(['Run','Batch'])[['Run','Batch']].itertuples(index=False, name=None))
    # Use apply to create tuple and check if in failed_batches
    df['batch_RT_control_fail']=df[['Run','Batch']].apply(tuple, axis=1).isin(failed_batches)

    # add run_pass column
    df['run_pass']=np.where((df['batch_neg_control_fail']==False) & (df['batch_RT_control_fail']==False) & (df['run_bases_pass']==True) & (df['sample_reads_pass']==True), True, False)
    # print runs and batches that failed
    failed_runs_batches=df[df['run_pass']==False][['Run','Batch', 'barcode', 'seq_name']].drop_duplicates()
    print('Runs and Batches that failed QC:')
    print(failed_runs_batches)
    return df

def run_pass_derivation(df):
    # add column if run passes QC
    # spike negs
    # Add column if Total run bases < 0.4 million bases
    df['run_bases_pass']=np.where(df['total run bases']<400_000, False, True)
    # Add column if Total sample reads < 25,000 reads
    df['sample_reads_pass']=np.where(df['total sample reads']<25_000, False, True)
    # This is also called the batch negatice control in the flow diagram
    # count number of samples that failed batch ampification negative controls
    df_negs=df[df['spiked']==0]
    df_negs=df_negs[df_negs['pathogen']!='unmapped']
    df_negs_pass=df_negs[df_negs['sample num reads']>1]
    if len(df_negs_pass)>0:
        failedruns=list(df_negs_pass['Run'].unique())
        print(f'Runs that failed negative controls: {failedruns}')
        # remove rows where negative controls failed
        df['batch_RT_control_fail']=df['Run'].isin(failedruns)
        #df=df[~df['Run'].isin(failedruns)]

    # add run_pass column
    df['run_pass']=np.where((df['batch_RT_control_fail']==False) & (df['run_bases_pass']==True) & (df['sample_reads_pass']==True), True, False)
    # print runs and batches that failed
    failed_runs_batches=df[df['run_pass']==False][['Run','Batch', 'barcode', 'seq_name']].drop_duplicates()
    print('Runs and Batches that failed QC:')
    print(failed_runs_batches)
    return df

def table_S2_roc_all_data(roc_all_data_files):
    # read in all roc all data files, match to expected prefix and concatenate
    prefixes={'derivation_roc_data_all': {'set': 'Derivation', 'zero_reads': True},
              'validation_roc_data_all': {'set': 'Validation', 'zero_reads': True},
              'derivation_non_zero_roc_data_all': {'set': 'Derivation', 'zero_reads': False},
              'validation_non_zero_roc_data_all': {'set': 'Validation', 'zero_reads': False}}


    df_list=[]
    for file in roc_all_data_files:
        dft=pd.read_csv(file)
        p=prefixes[file.split('/')[-1].split('.')[0]]
        dft['set']=p['set']
        dft['zero_reads']=p['zero_reads']
        if p['zero_reads']:
            suf=''
        else:
            suf=' excluding true-positives with 0 reads'
        dft.rename(columns={'sensitivity': f'{p["set"]} Sensitivity{suf}', 
                            'specificity': f'{p["set"]} Specificity{suf}',
                            'ROC AUC': f'{p["set"]} AUROC{suf}',
                            'minVal threshold': f'{p["set"]} Youdens optimal threshold{suf}'}, inplace=True)
        df_list.append(dft)
    # merge all dataframes on metric
    df_roc_all=df_list[0]
    for dft in df_list[1:]:
        df_roc_all=df_roc_all.merge(dft, on='metric', how='left', suffixes=('', '_drop'))

    cols=['metric', 'Derivation Sensitivity', 'Derivation Sensitivity excluding true-positives with 0 reads',
          'Derivation Specificity', 'Derivation AUROC', 'Derivation Youdens optimal threshold',
          'Validation Sensitivity', 'Validation Sensitivity excluding true-positives with 0 reads',
          'Validation Specificity', 'Validation AUROC', 'Validation Youdens optimal threshold']
    df=df_roc_all[cols].copy()
    # round to 3 decimal places
    df['Derivation Sensitivity']=df['Derivation Sensitivity'].round(3)
    df['Derivation Sensitivity excluding true-positives with 0 reads']=df['Derivation Sensitivity excluding true-positives with 0 reads'].round(3)
    df['Derivation Specificity']=df['Derivation Specificity'].round(3)
    df['Derivation AUROC']=df['Derivation AUROC'].round(3)
    df['Validation Sensitivity']=df['Validation Sensitivity'].round(3)
    df['Validation Sensitivity excluding true-positives with 0 reads']=df['Validation Sensitivity excluding true-positives with 0 reads'].round(3)
    df['Validation Specificity']=df['Validation Specificity'].round(3)
    df['Validation AUROC']=df['Validation AUROC'].round(3)
    
    #df_roc_all=pd.concat(df_list.values(), ignore_index=True)
    # save to csv
    df.to_csv('table_S2_roc_all_data.csv', index=False)

def get_files():
    args=argparse.ArgumentParser(description='Plotting script for mNGS respiratory validation data')
    args.add_argument('-v', '--validation_csv', help='Validation (study) CSV from flow diagram with AND_ratio applied', required=True)
    args.add_argument('-d', '--derivation_set_csv', help='Derivation set CSV from flow diagram and adjusted', required=False)
    args.add_argument('-p', '--pcr_types', help='Derivation set PCR types for each sample and sequence CSV file', required=False)
    args.add_argument('-q', '--qPCR_data', help='qPCR data from study assay CSV file', required=False)
    args.add_argument('-m', '--mapQ_data', help='mapQ from derivation mapped bam files data CSV file', required=False)
    args.add_argument('-I', '--IORD_meta', help='IORD meta data for patient level information CSV file', required=True)
    args.add_argument('-s', '--sample_derivation_meta', help='Sample derivation meta data CSV file', required=False)
    args.add_argument('-o', '--output_prefix', help='Output prefix for plots', required=False)
    args.add_argument('-r', '--roc_all_data', help='ROC all data CSV files', nargs='+', required=False)
    return args.parse_args()

if __name__ == "__main__":
    args=get_files()
    df=pd.read_csv(args.validation_csv)
    df=run_pass(df)
    upset_spiked_plot(df, set='Validation')
    upset_spiked_stacked_plot(df, set='Validation')
    upset_pass_plot(df, set='Validation')


    # generate table S2 roc all data if provided
    #if args.roc_all_data:
    #    table_S2_roc_all_data(args.roc_all_data)
    #total_samples=df.shape[0]
    #print(f'Total number of samples: {total_samples}')

    #biorifre_organisms=['Adenovirus', 'Coronavirus HKU1', 'Coronavirus NL63', 'Coronavirus OC43', 'Coronavirus 229E', 'Human Metapneumovirus', 'Influenza A', 'Influenza B', 'Parainfluenza 1', 'Parainfluenza 2', 'Parainfluenza 3', 'Parainfluenza 4', 'Respiratory Syncytial Virus', 'Rhinovirus/Enterovirus']
    #alinity_cephid_organisms=['Adenovirus', 'Coronavirus HKU1', 'Coronavirus NL63', 'Coronavirus OC43', 'Coronavirus 229E', 'Human Metapneumovirus', 'Influenza A', 'Influenza B', 'Parainfluenza Virus 1-4', 'Respiratory Syncytial Virus', 'Rhinovirus/Enterovirus']

    #plot_cts_vs_reads()

    #plot_Ct_reads_plotnine(df)

    #plot_sensitivity_by_run_and_pathogen(df)

    #plot_sensitivity_by_pathogen(df)

    # read in derivation set
    if args.derivation_set_csv: # derivation set csv
        dfd=pd.read_csv(args.derivation_set_csv)
        dfd=run_pass_derivation(dfd)
        if args.pcr_types: # add in test type info
            dfd_table, dfd_tt=derivation_set_characteristics(args.pcr_types, args.sample_derivation_meta, dfd)
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
        #dfd=run_pass_derivation(dfd)
        dfd.to_csv('derivation_set_with_run_names_test_types.csv', index=False)
    else:
        sys.exit('Please provide derivation set csv as second argument to continue plotting')

    #upset_spiked_plot(dfd, set='Derivation')
    #upset_spiked_stacked_plot(dfd, set='Derivation')
    #upset_pass_plot(dfd, set='Derivation')

    #plot_reads_bases(df, dfd)

    #plot_figure_S1(df, dfd)

    #plot_reads_bases_per_run(df, dfd)

    #combined, combined_test_counts = plot_gold_standard_pathogen_counts(df, dfd)

    #plot_pathogen_counts(combined, combined_test_counts)

    #plot_sensitivity_by_pathogen_plotnine(df, dfd)

    #plot_ratios(df, dfd)

    plot_alt_ratios(df, dfd)

    #correlation_plot(dfd)

    if args.qPCR_data: # qPCR data 
        qPCR=pd.read_csv(args.qPCR_data)
        qPCR['Run']=qPCR['seq_run_name']
        
        df_merged=plot_qPCR(df, qPCR)
        # IORD meta data
        dfI=pd.read_csv(args.IORD_meta)
        table_1_characteristics(df_merged, dfI, dfd_table)
        table_3_sensitivity_specificity(df_merged)

    #if sys.argv[5]: # mapQ data
    #    mapQ_df=pd.read_csv(sys.argv[5])
    #    mapQ_df['Run_order'] = mapQ_df['Run'].map(run_order)  # Map the run names to their order
    #    mapQ_df.dropna(subset=['Run_order'], inplace=True)  # Drop rows where Run_order is NaN
    #    mapQ_df.sort_values(by=['Run_order'], inplace=True)  # Sort by Run order and Batch
    #    mapQ_df['Run_order']=mapQ_df['Run_order'].astype(int)
    #    mapQ_df['Run']= 'Run ' + mapQ_df['Run_order'].astype(str).str.zfill(2)
    #    plot_mapQ(mapQ_df, dfd)
