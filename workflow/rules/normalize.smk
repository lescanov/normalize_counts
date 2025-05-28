# Consult comments in script for params description
rule normalize_counts:
	input:
		raw_counts=config['raw_counts'],
		clinical=config['clinical']
	params:
		gene_column='ensembl_id',
		id_column='patient_id',
		is_batch_corrected=False,
		batch_column=None,
		is_threshold_filtered=False,
		count_threshold=None,
		is_subtype_filtered=False,
		subtype_column=None,
		subtype_value=None,
		is_mutation_filtered=False,
		mutation_column=None,
		mutation_value=None,
		is_counts_log_transformed=True
	output:
		normalized_counts='results/normalized_test_counts.csv'
	conda:
		'../envs/normalization.yaml'
	script:
		'../scripts/normalize_counts.R'

