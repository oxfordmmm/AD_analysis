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

The nextflow workflow must be run for each sequencing flow cell run. This will take a long time and has some steps that could be optimised. The workflow was originally written to deal with Phi generated cicular reads, that are cut into sub-reads. The process of finding the SISPA primer and trimming it was mainted for this workflow, but it is expensive. Work could be done to test whether it is needed, with reads directly mapped without any pre trimming or filtering. However, the derivation thresholds were set using this workflow, so the study validation needs to use the same method. This script (`run_all.bash`) runs the `workflow_commands.bash` script for each sequencing run. It uses the `standard` profile in the nextflow.config files. You might need to edit this profile (or create your own) to something that suits your computing environment. 

```bash
repo_location='/home/nick/soft/'
bash run_all.bash ${repo_location}
```

## Run the analysis

This scripts collate all of the output results from the workflow and compile it into output files for analysis and plotting.

```bash
data='workflow_runs/'
repo_location='/home/nick/soft/'
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
mkdir -p AND_ratios additional_yield figures passing_samples supplemental table_1_tabs tmp
repo_location='/home/nick/soft/'
python3 ${repo_location}AD_analysis/study_analysis/flow_diagram.py | tee flow_diagram_results.txt
```

This script generates the plots and tables.

```bash
repo_location='/home/nick/soft/'
python3 ${repo_location}AD_analysis/study_analysis/plots.py \
        -v AND_ratios/biofire_results_merged_adjusted_0.1.csv \
        -d ${repo_location}AD_analysis/thresholds/flow_diagram/biofire_results_merged_adjusted.csv \
        -p meta_data/derivation_PCR_results_anonymised.csv \
        -q ${repo_location}AD_analysis/study_analysis/qpcr_info.csv \
        -m meta_data/mapQ.csv \
        -I meta_data/IORD2.csv \
        -t meta_data/validation_PCR_results_anonymised.csv \
        -s meta_data/derivation_sample_meta.csv \
        -E meta_data/ENA_accessions.csv \
        -r model_data/derivation_roc_data_all.csv \
    model_data/derivation_non_zero_roc_data_all.csv \
    ${repo_location}/AD_analysis/thresholds/optimization/validation/validation_roc_data_all.csv \
    ${repo_location}/AD_analysis/thresholds/optimization/validation/no_zero/validation_non_zero_roc_data_all.csv \
    -k seqkit_stats/derivation_all_fqs.txt seqkit_stats/derivation_passed_fqs.txt \
    seqkit_stats/validation_all_fqs.txt seqkit_stats/validation_passed_fqs.txt | tee plots_output.txt
```

The tables and plots are generated in the working directory, `figures`, and `supplemental` folders. There is a known discrepency between the published results and the results that should be generated here. For extra vigilence, a further round of human read removal was conducted before data was uploaded to the ENA using Deacon. This slightly changed the overall read numbers in the run and the denominator for some of the metrics. 2 samples subsequently passed thresholds that had previously failed. Increasing the sensitivity from 51% (133/260) to 52% (135/260) (Sample numbers 198 and 211, both FluA). There were also 3 more false positive detections. The results using the alternative criteria do not appear to be affected.