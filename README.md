## Purpose:
To normalize bulk gene based miRNA and RNA sequencing counts for downstream statistical and visual analysis, with cancer datasets specifically in mind.

## Features:
- edgeR TMM library normalization of raw gene based bulk counts for cross sample comparisons
- cpm count quantification
- snakemake for reproducible results
- Can generate either UMAP or PCA plots, with samples grouped according to a column in the metadata, to assess batch effects
- Can select a subset of samples in raw counts for normalization, based on:
1. A subtype of patients that is denoted in the a column of the metadata
2. Binary status of a certain variable as specified by a column in the metadata, intended for a mutations of a certain gene
3. Both 1 and 2, allowing for normalization based on patient subtypes that harbour a particular mutation

## Methods:
Will use edgeR's TMM normalization technique to normalize gene based bulk reads for cross sample comparison. If there is a known batch effect, the user can correct for this by specifying a column in the metadata (termed clinical in this workflow) corresponding to the batch. Corrections will be done using limma's removeBatchEffect function. The workflow will be written using snakemake and thus requires snakemake to run. Before normalization, the raw counts are filtered by default using edgeR's filterByExpr function. The user has the option to instead filter genes by setting a threshold for what will be considered lowly expressed. Once TMM normalization is done, expression is calculated using edgeR's cpm function. This workflow will only use unique samples in the raw counts, selecting the first sample if the same sample has replicates.

Further reading on edgeR can be done [here](https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf)

## Output:
A csv file of normalized counts. The resulting csv is in a column major format, with a column representing genes.

## Usage:
The two main components for this workflow is the metadata (clinical) and the raw bulk RNA-sequencing counts that assumes a column major format (specifically where each column represents a sample, with the exception for a column denoting the genes). The metadata must contain a column that has the same identifiers that are present in the raw counts.
- Modify the the main Snakefile to specify the outputs of the pipeline.
- Modify individual rules in the workflow folder to adjust parameters on how normalization is done
- Modify the config to specify raw counts and clinical

