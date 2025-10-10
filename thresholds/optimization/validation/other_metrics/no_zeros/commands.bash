python3 ../../../models.py \
	--input_file ../../biofire_results_merged_adjusted_0.1.csv \
       	--use_metrics True \
	--remove_no_reads True \
	--set validation

mv roc_data_all.csv roc_data_all_other_metrics_no_zeros.csv
