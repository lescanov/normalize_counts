# Purpose:
# Perform dimensionality reduction on normalized counts
# and assess if any outliers are present or if further correction needed

# Normalized counts are plotting using PCA that are scaled

library(dplyr)
library(readr)
library(ggplot2)
library(tibble)

# Input
# -----
normalized_counts <- read_csv(snakemake@input[["normalized_counts"]])
clinical <- read_csv(snakemake@input[["clinical"]])

# Pipeline parameters
# -----
# Can visualize how patients clusters in accordance with a grouping variable
# If is_plotting_grouping is FALSE, grouping_column can be NA
is_plotting_grouping <- snakemake@params[["is_plotting_grouping"]]
grouping_column <- snakemake@params[["grouping_column"]]

gene_column <- snakemake@params[["gene_column"]]
id_column <- snakemake@params[["id_column"]]

# Output
# -----
plot_outpath <- snakemake@output[["plot"]]

stopifnot(
  is.logical(is_plotting_grouping),
  gene_column %in% colnames(normalized_counts)
)

# Formatting normalized counts
# Counts would be in column major and need to convert to row-major
counts <- normalized_counts %>%
  column_to_rownames(var = gene_column) %>%
  t()

pca_counts <- prcomp(counts, scale = TRUE)[["x"]] %>%
  as.data.frame() %>%
  select(PC1, PC2) %>%
  rownames_to_column(var = id_column)

if (is_plotting_grouping == TRUE) {

  # Formatting clinical to match patient order as normalized counts
  clinical <- clinical %>%
    filter(duplicated(!!sym(id_column)) == FALSE) %>%
    filter(!!sym(id_column) %in% colnames(normalized_counts)) %>%
    arrange(match(!!sym(id_column), colnames(normalized_counts))) %>%
    select(all_of(id_column), all_of(grouping_column))

  # In case grouping factor is a numeric value, this will allow it to
  # show as discrete groupings instead of a gradient on the plot
  clinical[[grouping_column]] <- as.factor(clinical[[grouping_column]])

  pca_counts <- pca_counts %>%
    inner_join(clinical, by = id_column)

  pca_plot <- ggplot(
    pca_counts,
    aes(x = PC1, y = PC2, colour = !!sym(grouping_column))
  ) +
    geom_point()

} else {
  pca_plot <- ggplot(
    pca_counts,
    aes(x = PC1, y = PC2)
  ) +
    geom_point()
}

ggsave(pca_plot, filename = plot_outpath)
