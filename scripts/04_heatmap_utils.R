library(ComplexHeatmap)
library(tibble)
library(dplyr)


# FUNCTION: get_group_means
# ------------------------------------------------------------
# Computes average gene expression per group (e.g., NR vs R)
# using RNA counts from a Seurat object.
# Returns a dataframe: Gene | Group1 | Group2 | ...
# ----------------------------------------------------------
get_group_means <- function(seurat_obj, group_var) {
  groups <- unique(seurat_obj@meta.data[[group_var]])
  mean_list <- list()

  for (grp in groups) {
    cells <- colnames(seurat_obj)[seurat_obj@meta.data[[group_var]] == grp]
    avg <- rowMeans(seurat_obj@assays$RNA@counts[, cells, drop = FALSE])
    mean_list[[grp]] <- avg
  }

  df <- as.data.frame(mean_list)
  df <- rownames_to_column(df, "Gene")
  df
}

# --------------------------------------------------------------------
# FUNCTION: extract_genes
# Extracts gene symbols from enrichment outputs.
# Handles gene lists separated by "/" or ",".
# ------------------------------------------------------------
extract_genes <- function(df, column = "geneID") {
  g <- gsub("/", ",", df[[column]])
  g <- unlist(strsplit(g, ","))
  unique(trimws(g))
}

# ------------------------------------------------------------
# FUNCTION: plot_heatmap
# Plots heatmap using ComplexHeatmap.
# ------------------------------------------------------------
plot_heatmap <- function(mat, annotation = NULL, title = "") {
  Heatmap(
    as.matrix(mat),
    name = "Expression",
    col = colorRampPalette(c("blue", "white", "red"))(50),
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    top_annotation = annotation,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 8),
    column_title = title
  )
}
