python3 ../models.py --input_file ../../biofire_results_merged.csv --use_metrics True --negatives_meta_file ../../flow_diagram/negatives_meta.csv

mv roc_data_all.csv roc_data_all_other_metrics.csv
