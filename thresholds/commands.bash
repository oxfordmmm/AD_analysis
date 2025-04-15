
data='/mnt/data/analysis/nick/agnostic_diagnostic/'

python3  ~/soft/AD_analysis/thresholds/analysis.py -i \
    ${data}expt10A_17072024 \
    ${data}expt10_03072024 \
    ${data}expt11_150824 \
    ${data}Expt9A_SISPA_Katie \
    ${data}010524_Expt9_SISPA_Daisy \
    ${data}010524_Expt9_SISPA_Kate \
    -m ~/soft/AD_analysis/thresholds/meta.csv \
    -p ~/soft/AD_analysis/thresholds/pathogens.csv \
    -pr ~/soft/AD_analysis/thresholds/pathogen_reduced.csv \
    -bf ~/soft/AD_analysis/thresholds/biofire_pathogens.csv \
    -o results

