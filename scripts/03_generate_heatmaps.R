library(Seurat)
library(readxl)
library(ComplexHeatmap)

source("R/heatmap_utils.R")


# Load Seurat object
# ------------------------------------------------------------
seurat <- readRDS("path/to/your/seurat_object.rds")

# ------------------------------------------------------------
# Subset specific cell type + responder categories
# ------------------------------------------------------------
sub_nr <- subset(seurat, Responder_Status == "NR" & cell_types_main == "Plasma/Dendritic")
sub_r  <- subset(seurat, Responder_Status == "R"  & cell_types_main == "Plasma/Dendritic")

# ------------------------------------------------------------
# Load enrichment results
up_df   <- read_excel("data/enrichment/up.xlsx")
down_df <- read_excel("data/enrichment/down.xlsx")

# Extract genes
up_genes   <- extract_genes(up_df)
down_genes <- extract_genes(down_df)

# ------------------------------------------------------------
# Compute mean expression for NR and R
# ------------------------------------------------------------
means_all <- get_group_means(seurat, "Responder_Status")

# ------------------------------------------------------------
# Create matrices for heatmaps
# ------------------------------------------------------------
mat_up <- means_all[means_all$Gene %in% up_genes, ]
mat_down <- means_all[means_all$Gene %in% down_genes, ]

rownames(mat_up) <- mat_up$Gene
rownames(mat_down) <- mat_down$Gene

mat_up$Gene <- NULL
mat_down$Gene <- NULL

# ------------------------------------------------------------
# Annotation for NR and R
# ------------------------------------------------------------
annotation = HeatmapAnnotation(
  Group = c("NR", "R"),
  col = list(Group = c("NR" = "brown", "R" = "green"))
)


# Generate heatmaps
heat1 <- plot_heatmap(mat_up, annotation, "Upregulated pathways")
heat2 <- plot_heatmap(mat_down, annotation, "Downregulated pathways")


# Save outputs
# ------------------------------------------------------------
pdf("data/heatmap_outputs/upregulated_heatmap.pdf")
draw(heat1)
dev.off()

pdf("data/heatmap_outputs/downregulated_heatmap.pdf")
draw(heat2)
dev.off()
