python3 ../../models.py \
	--input_file ../biofire_results_merged_adjusted_0.1.csv \
	--set validation \
	--remove_no_reads True

mv roc_data_all.csv validation_non_zero_roc_data_all.csv
