# Winter virus study 

Here is how to replicate the results from our winter virus study using the workflow in this repository and the raw fastq data uploaded to the ENA. 

## Download and prepare data

The fastq files need to be downloaded from the ENA and sorted into the runs they were sequenced in. This is because the run totals are used to caluculate thresholds independantly of other runs. 

```bash
git clone https://github.com/oxfordmmm/AD_analysis.git
cd AD_analysis/study_replication
python3 download_data.py 
```

## Run the workflow for each run

The nextflow workflow must be run for each sequencing flow cell run. This will take a long time and has some steps that could be optimised. 

```bash
bash run_all.bash
```

## Run the analysis

This scripts collate all of the output results from the workflow and compile it into output files for analysis and plotting.

```bash
data='workflow_data/'
repo_location='~/soft/'
python3  ${repo_location}AD_analysis/study_analysis/analysis.py -i \
    ${data}AD_winter_study_220125 \
    ${data}AD_winter_study_201224 \
    ${data}AD_winter_study_070325 \
    ${data}AD_winter_study_170325 \
    ${data}AD_winter_study_100425 \
    ${data}AD_winter_study_160725 \
    ${data}AD_winter_study_220725 \
    ${data}AD_winter_study_240725 \
    ${data}AD_winter_study_300725 \
    ${data}AD_winter_study_010825 \
    ${data}AD_winter_study_130825_rpt050825 \
    -m ${repo_location}AD_analysis/study_analysis/meta.csv \
    -p ${repo_location}AD_analysis/thresholds/pathogens.csv \
    -pr ${repo_location}AD_analysis/thresholds/pathogen_reduced.csv \
    -bf ${repo_location}AD_analysis/thresholds/biofire_pathogens.csv \
    -o results
```

This scripts generates the reuslts for flow diagram and an input file for the plots.

```bash
repo_location='~/soft/'
python3 ${repo_location}AD_analysis/study_analysis/flow_diagram.py | tee flow_diagram_results.txt
```

This script generates the plots and tables.

```bash
repo_location='~/soft/'
python3 ${repo_location}AD_analysis/study_analysis/plots.py \
        -v AND_ratios/biofire_results_merged_adjusted_0.1.csv \
        -d ${repo_location}AD_analysis/thresholds/flow_diagram/biofire_results_merged_adjusted.csv \
        -p meta_data/expt_9to11_Nick_pcr_organisms_plus_machine_extra_daisy.csv \
        -q ${repo_location}AD_analysis/study_analysis/qpcr_info.csv \
        -m meta_data/mapQ.csv \
        -I meta_data/IORD2.csv \
        -t meta_data/extracted_samples_storage_and_processing_info_pcr_organisms_test_code_check.csv \
        -s meta_data/metaDF.csv \
        -r model_data/derivation_roc_data_all.csv \
    model_data/derivation_non_zero_roc_data_all.csv \
    ${repo_location}/AD_analysis/thresholds/optimization/validation/validation_roc_data_all.csv \
    ${repo_location}/AD_analysis/thresholds/optimization/validation/no_zero/validation_non_zero_roc_data_all.csv \
    -k seqkit_stats/derivation_all_fqs.txt seqkit_stats/derivation_passed_fqs.txt \
    seqkit_stats/validation_all_fqs.txt seqkit_stats/validation_passed_fqs.txt | tee plots_output.txt
```