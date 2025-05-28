# Purpose:
# Normalize RNA seq counts using edgeR's TMM technique and then
# quantify counts using counts per million.

# WARNING:
# This script implicity filters raw counts for samples found in
# the clinical. It is assumed that the clinical contains identifiers
# for samples that you want to be normalized and that you have
# done your due deligence to ensure this.

library(edgeR)
library(limma)
library(dplyr)
library(readr)
library(stringr)
library(tibble)

# Input
# -----
raw_counts <- read_csv(snakemake@input[["raw_counts"]])
clinical <- read_csv(snakemake@input[["clinical"]])

# Pipeline parameters
# -----
gene_column <- snakemake@params[["gene_column"]]

# Id column is the list of identifiers that map clinical to raw counts
id_column <- snakemake@params[["id_column"]]

# Batch column is only necessary if is_batch_corrected is TRUE
# Otherwise it can be left as NA
is_batch_corrected <- snakemake@params[["is_batch_corrected"]]
batch_column <- snakemake@params[["batch_column"]]

# Threshold filtering prevents genes below a certain threshold from
# being filtered out prior to normalization. If is_threshold_filtered
# is FALSE then will use edgeR's filterByExpr. If FALSE, can keep
# count_threshold as NA
is_threshold_filtered <- snakemake@params[["is_threshold_filtered"]]
count_threshold <- snakemake@params[["count_threshold"]]

# Option to filter raw counts for samples in a given subtype as specified
# in clinical, with subtype_value corresponding to the subtype you
# wish to filter for.
# Both subtype_column and subtype_value may be left as NA if
# is_subtype_filtered is FALSE
is_subtype_filtered <- snakemake@params[["is_subtype_filtered"]]
subtype_column <- snakemake@params[["subtype_column"]]
subtype_value <- snakemake@params[["subtype_value"]]

# Option to filter raw counts for samples that have a particular status
# concerning a mutation. mutation_value corresponds to the status you
# wish to filter for.
# Both mutation_column and mutation_value can be left as NA if
# is_mutation_filtered is FALSE
is_mutation_filtered <- snakemake@params[["is_mutation_filtered"]]
mutation_column <- snakemake@params[["mutation_column"]]
mutation_value <- snakemake@params[["mutation_value"]]

# Whether or not the user wants log2 transformed cpm counts
is_counts_log_transformed <- snakemake@params[["is_counts_log_transformed"]]

# Output
# -----
normalized_counts_path <- snakemake@output[["normalized_counts"]]

stopifnot(
  is.logical(is_batch_corrected),
  is.logical(is_threshold_filtered),
  is.logical(is_subtype_filtered),
  is.logical(is_mutation_filtered),
  is.logical(is_counts_log_transformed),
  id_column %in% colnames(clinical),
  gene_column %in% colnames(raw_counts)
)

# Filtering for specific subtypes
if (is_subtype_filtered == TRUE) {
  stopifnot(subtype_value %in% clinical[[subtype_column]])
  clinical <- clinical %>%
    filter(!!sym(subtype_column) == subtype_value)
}

# Fitlering for specific mutation status
if (is_mutation_filtered == TRUE) {
  stopifnot(mutation_value %in% clinical[[mutation_column]])
  clinical <- clinical %>%
    filter(!!sym(mutation_column) == mutation_value)
}

# Ensuring that raw counts and clinical share same sample pool
# Also that there are no duplicated samples in either clinical or raw counts
clinical <- clinical %>%
  filter(duplicated(!!sym(id_column)) == FALSE)
shared_samples <- intersect(clinical[[id_column]], colnames(raw_counts))

if (length(test) < 1) {
  stop("There are no overlapping samples between clinical and raw counts")
}

clinical <- clinical %>%
  filter(!!sym(id_column) %in% shared_samples)

raw_counts <- raw_counts %>%
  select(!!sym(gene_column), all_of(shared_samples)) %>%
  column_to_rownames(var = gene_column)

# Now preparing and normalizing counts
counts <- as.matrix(raw_counts)
counts <- DGEList(counts)

if (is_threshold_filtered == FALSE) {
  keep <- filterByExpr(counts)
} else {
  stopifnot(is.numeric(count_threshold))
  keep <- rowSums(counts$counts) > count_threshold
}

counts <- counts[keep, , keep.lib.sizes = FALSE]
counts <- calcNormFactors(counts, method = "TMM")
counts <- cpm(counts, log = is_counts_log_transformed)

# Correcting counts if needed
if (is_batch_corrected == TRUE) {
  batch <- clinical %>%
    arrange(match(!!sym(id_column), colnames(counts))) %>%
    mutate(batch_factor = as.factor(!!sym(batch_column))) %>%
    select(batch_factor) %>%
    deframe()

  counts <- counts %>%
    removeBatchEffect(batch)
}

counts <- counts %>%
  as.data.frame() %>%
  rownames_to_column(var = gene_column)

write_csv(counts, normalized_counts_path)
