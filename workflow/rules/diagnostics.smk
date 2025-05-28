rule plot_pca:
	input:
		normalized_counts=rules.normalize_counts.output,
		clinical=config['clinical']
	params:
		is_plotting_grouping=True,
		grouping_column='subtype',
		gene_column='gene_id',
		id_column='patient_id'
	output:
		plot='results/pca.pdf'
	conda:
		'../envs/plotting.yaml'
	script:
		'../scripts/generate_diagnostic_plot.R'
