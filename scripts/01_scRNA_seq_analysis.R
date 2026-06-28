# Use the code below to read packages from requirements.txt and install if required
# Read package list
packages <- readLines("requirements.txt")

# Install missing packages

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# Apply installation function to all packages
sapply(packages, install_if_missing)

#Load necessary libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(DoubletFinder)
library(SingleR)

# Define paths to input files (adjust paths as necessary)
# Modify these paths to point to the dataset
zip_file_path <- "path_to_expression_data.txt.gz" #Replace with actual file path
patient_info_path <- "path_to_patient_info.txt.gz"  # Replace with the actual file path

# Load expression data (TPM counts) from a gzipped file
# Using gzfile to specify the file format
zip_connection <- gzfile(zip_file_path, "r")
tpm_counts <- read.table(zip_file_path, header = TRUE, row.names = 1, skip = 1, sep = "\t")

# Display the first few rows of the expression data (for verification)
head(tpm_counts)

# Load patient data
patient_info <- read.table(patient_info_path, header = TRUE, sep = "\t")

# Create Seurat object with expression data and metadata
seurat_obj <- CreateSeuratObject(counts = tpm_counts, meta.data = patient_info)

# Removing large objects to save memory
remove(tpm_counts) #Removes PM counts object
remove(patient_info) #Removes patient metadata object


# Filter only those treated by anti-PD1
# Adding cell names to metadata for easier reference
seurat_obj@meta.data$cell_names <- rownames(seurat_obj@meta.data)

# Splitting cell names to extract sample IDs easily
split_names <- strsplit((seurat_obj@meta.data$cell_names), "\\.")
#split_names <- unlist(strsplit(seurat_obj@meta.data['cell_names'], "\\."))
sample_name <- sapply(split_names, function(x) x[1])

# Add Sample IDs to metadata as a new column
seurat_obj@meta.data$Sample_id = sample_name

# (Optional) Uncomment the next line to remove the Sample_id column if necessary
#seurat_obj@meta.data$Sample_id <- NULL

#Set Sample_id as the active identity class (Idents) for subsetting and analysis
Idents(seurat_obj) <- "Sample_id"
# The `Idents()` function sets the active identity class in the Seurat object.
# This is critical for subsetting or grouping cells based on metadata features.
# NOTE: You can change the identity class (Idents) to any column in the metadata, depending on your analysis goal.

#Create a seurat obj for only anti PD-1 samples
samples_to_remove <- c('Pre_P1', 'Pre_P4', 'Post_P4', 'Pre_P6', 'Pre_P7', 'Post_P7', 'Pre_P8', 
                       'Post_P8', 'Post_P13', 'Pre_P26', 'Pre_P28', 'Post_P28', 'Post_P28_2',
                       'Post_P4_T_enriched', 'Post_P13_T_enriched', 'Post_P8_T_enriched')
seurat_obj1 <- seurat_obj[, !(seurat_obj$Sample_id %in% samples_to_remove)]
length(rownames(seurat_obj1@meta.data))
length(colnames(seurat_obj1))

#removing the old seurat object to save memory
remove(seurat_obj)


#Rearranging column order (for convenience)
seurat_obj1@meta.data <- seurat_obj1@meta.data[, c(40, 1:39)]

# Removes the first row for debugging
seurat_obj1@meta.data <- seurat_obj1@meta.data[-1, ]
seurat_obj1 <- seurat_obj1[,-1]



#Adding resopnder status column
seurat_obj1@meta.data <- seurat_obj1@meta.data %>%
  mutate(Responder_Status = ifelse(Sample_id %in% c('Post_P1', 'Post_P5_2', 'Post_P17', 'Post_P19',
                                                    'Post_P21', 'Pre_P24', 'Pre_P29', 'Pre_P33', 
                                                    'Pre_P35', 'Post_P17_myeloid_enriched',
                                                    'Post_P19_myeloid_enriched'),
                                   'Responder', 'Non-Responder'))


# Rearranging columns for convenience
seurat_obj1@meta.data <- seurat_obj1@meta.data[, c(41, 1:40)]
seurat_obj1@meta.data$Source <- 'Melanoma'
seurat_obj1@meta.data <- seurat_obj1@meta.data[, c(42, 1:41)]
---- Quality control(QC) ----

# Checking for percentage of mitochondrial and ribosomal genes
seurat_obj1[["percent.mt"]] <- PercentageFeatureSet(seurat_obj1, pattern = "^MT-")
seurat_obj1[["percent.rb"]] <- PercentageFeatureSet(seurat_obj1, pattern = "^RP[SL]")

#Changing the identity class
Idents(seurat_obj1) <- "Source"
# NOTE: Change the identity class depending on your analysis needs (e.g., clustering, sample IDs, cell types).

# Removing some unwanted variables to save memory
remove(split_names)
remove(sample_name)
remove(samples_to_remove)

# QC plots
        features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb"),
        ncol = 4,pt.size = 0.1) & 
  theme(plot.title = element_text(size=10))
FeatureScatter(seurat_obj1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')


#Filtering the outliers based on the plots---
seurat_obj1 <- subset(seurat_obj1, 
                      subset = nFeature_RNA >200 & nFeature_RNA < 4000 & percent.mt <2 & 
                        percent.rb < 10)

# Normalizes the counts using log transformations
seurat_obj1 <- NormalizeData(seurat_obj1)

#Finding 2000 highly variable features
seurat_obj1 <- FindVariableFeatures(seurat_obj1, selection.method = 'vst', nfeatures = 2000)

top10 <- head(VariableFeatures(seurat_obj1),10)

# Plot variable features
plot1 <- VariableFeaturePlot(seurat_obj1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

#Scaling
all.genes <- rownames(seurat_obj1)
seurat_obj1 <- ScaleData(seurat_obj1, features = all.genes)
# ---- Linear dimensionality reduction ----
# Perform Principal component analysis
seurat_obj1 <- RunPCA(seurat_obj1, features = VariableFeatures(object = seurat_obj1))
# Visualize the heatmap of PCA loadings for the first principal component
# dims = 1: Focus on the first principal component (this could be changed for other principal components)
DimHeatmap(seurat_obj1, dims = 1, cells = 500, balanced  = TRUE)

# Visualize PCA using scatterplot
DimPlot(seurat_obj1, reduction = "pca")

# Determine the number of significant PCs using an elbow plot
ElbowPlot(seurat_obj1)

# ---- Clustering ----
# Identify the nearest neighbors using the first 15 principal components (select PC based on variation)
seurat_obj1 <- FindNeighbors(seurat_obj1, dims = 1:15)
# Perform clustering using different resolutions
seurat_obj1 <- FindClusters(seurat_obj1, resolution = c(0.1,0.3,0.5,0.7,1))
# `FindClusters()` applies a graph-based clustering algorithm.
# The resolution parameter controls cluster granularity:
# - Lower values (e.g., 0.1) result in fewer clusters.
# - Higher values (e.g., 1) result in more clusters.
# Here, multiple resolutions are tested for comparison.
# The new clusters appear in seurat_obj1@meta.data as RNA_snn_res.


# Set the cluster identity for visualization and analysis
Idents(seurat_obj1) <- 'RNA_snn_res.0.5'
# This sets the clustering identity to the resolution of 0.5.
# Change this based on which resolution gives the best separation.

# ---- UMAP Visualization ----
#set seed for reproducibility
set.seed(123)
seurat_obj1 <- RunUMAP(seurat_obj1, dims = 1:15)

# Visualize the clusters in UMAP space
DimPlot(seurat_obj1, reduction = 'umap')

# ---- FINDING DOUBLETS ----
# finding and removing doublets is optional and not recommended as doublets may contain important mutations
#pK Identification (no ground-truth)--
sweep.res.list <- paramSweep_v3(seurat_obj1, PCs = 1:15, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

#plot pK vs BCmetric
ggplot(bcmvn, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()


#Selecting the optimal pK value 
optimal.pK <- bcmvn %>%
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
optimal.pK <- as.numeric(as.character(optimal.pK[[1]]))


#Homotypic doublet proportion estimate
annotations <- seurat_obj1@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)

nExp.poi <- round(0.076 * nrow(seurat_obj1@meta.data))
nExp.poi.adj <- round(nExp.poi * (1-homotypic.prop))

seurat_obj1 <- doubletFinder_v3(seu = seurat_obj1, PCs = 1:15,
                                pK = optimal.pK, nExp = nExp.poi.adj,
                                reuse.pANN = FALSE, sct = FALSE)


names(seurat_obj1@meta.data)
colnames(seurat_obj1@meta.data)[52] <- 'Doublets'

#Visualize doublets
DimPlot(seurat_obj1, reduction = 'umap', group.by = 'Doublets') + 
  ggtitle('UMAP - Doublets vs Singlets')

# ---- Visualizing UMAP projections to assess clustering and check batch effects if any ----
#DimPlot(seurat_obj1,label.size = 4,repel = T,label = T)
DimPlot(seurat_obj1, reduction = 'umap', group.by = 'RNA_snn_res.0.3')
#DimPlot(seurat_obj1, reduction = 'umap', group.by = 'SingleR.labels.main')
DimPlot(seurat_obj1, reduction = 'umap', group.by = 'Responder_Status')
DimPlot(seurat_obj1, reduction = 'umap', group.by = 'orig.ident')

#---Save the seurat object---
# saving the object to avoid running previous steps again (optional)
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")
install.packages("remotes")
remotes::install_github("satijalab/seurat-data")
library(SeuratDisk)
#library(SeuratData)
SaveH5Seurat(seurat_obj1, overwrite = TRUE)

#---- Cell type annotation ----
#Manual cell type identification FindConservedMarkers function
DefaultAssay(seurat_obj1)
#BiocManager::install('multtest')
#install.packages('metap')

#FindConservedMarkers
markers_cluster11 <- FindConservedMarkers(seurat_obj1, ident.1 = 11,
                                          grouping.var = 'Responder_Status')
#Repeat the above step for each clusters to find conserved markers for each clusters
#Rename Clusters
seurat_obj1 <- RenameIdents(seurat_obj1, `10`= 'B Cells' )
DimPlot(seurat_obj1, reduction = 'umap', label = TRUE)
# Automatic cell type annotation with SingleR

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("celldex")

ref <- celldex::HumanPrimaryCellAtlasData()
seurat_counts <- GetAssayData(seurat_obj1, slot = 'counts')

pred1 <- SingleR(test = seurat_counts, ref = ref, labels = ref$label.main)
pred
seurat_obj1$SingleR.labels.main <- pred1$labels[match(rownames(seurat_obj1@meta.data),
                                                      rownames(pred1))]

#Plot to observe the SingleR generated labels
DimPlot(seurat_obj1, reduction = 'umap', group.by = 'SingleR.labels.main')
remove(pred1)


#Adding new responder_status for convenience (optional, this was changed to get better labels during plotting)
colnames(seurat_obj1@meta.data)[2] <- 'Response'
seurat_obj1@meta.data <- seurat_obj1@meta.data %>%
  mutate(Responder_Status = ifelse(Sample_id %in% c('Post_P1', 'Post_P5_2', 'Post_P17', 'Post_P19',
                                                    'Post_P21', 'Pre_P24', 'Pre_P29', 'Pre_P33', 
                                                    'Pre_P35', 'Post_P17_myeloid_enriched',
                                                    'Post_P19_myeloid_enriched'),
                                   'R', 'NR'))

#Add cell type column
seurat_obj1@meta.data$cell_type <- paste(Idents(seurat_obj1), 
                                         seurat_obj1@meta.data$Responder_Status, 
                                         sep = '_')
seurat_obj1@meta.data$cell_type_main <- Idents(seurat_obj1)
# seurat_obj1@meta.data <- seurat_obj1@meta.data[, c(57, 1:56)]

# Plotting umap after labelling each cell types
DimPlot(seurat_obj1, reduction = 'umap', label = TRUE)+
  ggtitle("Cell Types")+
  theme(plot.title = element_text(hjust = 0.5))

# ---- Differential Expression Analysis Between Conditions ----

# FindMarkers between two conditions
all.cell.markers <- FindMarkers(seurat_obj1, ident.1 = 'NR', 
                                ident.2 = 'R')


# Filter results to include only genes with a p-value < 0.05
filtered_cells <- all.cell.markers[all.cell.markers$p_val < 0.05, ]
#nrow(seurat_obj1@assays$RNA@counts)

# Export the complete differential expression results to a CSV file
write.csv(filtered_cells, 'DEAnalysis.csv', row.names = TRUE) # Replace file name as needed


# ---- Subset Upregulated and Downregulated Genes ----
# Extract upregulated genes (log2 fold change > 1)
upregulated <- filtered_cells[filtered_cells$avg_log2FC > 1, ]

# Extract downregulated genes (log2 fold change < -1)
downregulated <- filtered_cells[filtered_cells$avg_log2FC < -1, ]

# Save the list of downregulated genes to a CSV file
write.csv(upregulated, 'Upregulated.csv', row.names = TRUE)

# ---- Save the Current Workspace ----

# Save the entire R workspace, including all variables, objects, and loaded data,
# to a file named 'melanoma_data_recent.RData'. This file can be loaded later to
# restore the saved state of the session.
save.image('filename.RData')
