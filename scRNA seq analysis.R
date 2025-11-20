# Load the required libraries
library(tidyverse)
library(writexl)
library(stringr)
library(org.Hs.eg.db)
library(dplyr)
library(clusterProfiler)

# Ensure that the Seurat object (seurat_obj1) is loaded into the global environment 
# before running this script. This can be done by loading it from an RDS file, e.g.,:
# seurat_obj1 <- readRDS("path/to/your/seurat_object.rds")

#View metadata of the seurat object
View(seurat_obj1@meta.data)
#Assign the ident
seurat_obj1@meta.data$cell_types_main = Idents(seurat_obj1)

# Subset the Seurat object to include only a specific cell type, 
# based on the analysis goals (e.g., "Plasma/Dendritic" cells in this case)
Plasma_i = subset(seurat_obj1,cell_types_main == 'Plasma/Dendritic')

# View the metadata of subsetted object
View(Plasma_i@meta.data)

# Set the active identities of the subsetted Seurat object to "Responder_Status"
Idents(Plasma_i)<- "Responder_Status"

# Perform differential expression analysis between two responder groups
Plasma_i.markers <- FindMarkers(Plasma_i, ident.1 = 'NR', 
                                ident.2 = 'R')

# Add the gene names as a new column in the differential expression results
Plasma_i.markers$Genes <- rownames(Plasma_i.markers)



# ---- Enrichment Analysis ----

# Annotate differentially expressed genes (diffexp) with labels:
# 'UP' for upregulated (avg_log2FC > 1 and p_val < 0.5),
# 'DOWN' for downregulated (avg_log2FC < -1 and p_val < 0.5),
# 'NO' for non-significant changes.
Plasma_i.markers <- Plasma_i.markers %>% mutate(diffexpressed = case_when(
  avg_log2FC > 1 & p_val < 0.5 ~ 'UP',         # Upregulated genes
  avg_log2FC < -1 & p_val < 0.5 ~ 'DOWN',      # Downregulated genes
  TRUE ~ 'NO'                                 # Non-significant genes
))

# Extract the gene names to use for downstream analysis
genes_in_data <- Plasma_i.markers$Genes

# Load background genes from the Reactome pathway database (.gmt file)
# this could be downloaded from https://reactome.org/
# Ensure the file path is correct, and adjust as necessary for your system.
pwl2 <- read.gmt('c2.cp.reactome.v2023.2.Hs.symbols.gmt')

# Filter the Reactome data to include only genes present in the input gene set
pwl2 <- pwl2[pwl2$gene %in% genes_in_data,]
bg_genes <- pwl2  # Save the filtered Reactome data as background genes

# Keep only the differentially expressed genes (UP and DOWN) for enrichment analysis
Plasma_i.markers <- Plasma_i.markers[Plasma_i.markers$diffexpressed != 'NO', ]

# Export the annotated differentially expressed genes to an Excel file (if required)
write_xlsx(Plasma_i.markers, "Differential_Expression_Results.xlsx")

# Split the differentially expressed genes into separate lists by their status ('UP' and 'DOWN')
deg_results <- split(Plasma_i.markers, Plasma_i.markers$diffexpressed)



# ---- Run Enrichr Analysis ----

# Perform enrichment analysis using the enricher function from clusterProfiler
# Iterate through the lists of differentially expressed genes ('UP' and 'DOWN')
res <- lapply(names(deg_results),
              function(x) enricher(gene = deg_results[[x]]$Genes,  # Input genes for enrichment
                                   TERM2GENE = bg_genes))         # Background genes
names(res) <- names(deg_results)  # Assign names ('UP' and 'DOWN') to the results


# Extract the enrichment results and combine them into a single data frame
res_df <- lapply(names(res), function(x) rbind(res[[x]]@result))  # Access the result slot
names(res_df) <- names(res)  # Assign names to the individual result data frames
res_df <- do.call(rbind, res_df)  # Combine all results into one data frame


# Clean up the enrichment results by modifying the 'diffexpressed' column (for neat labels - optional)
# Extract the differential expression label ('UP' or 'DOWN') from row names
res_df <- res_df %>% mutate(diffexpressed = gsub('\\.REACTOME.*$', '', rownames(res_df)))



# ---- Filter Enrichment Results ----

# Filter enrichment results based on the adjusted p-value (p.adjust < 0.05)
target_res <- unique(res_df$ID[res_df$p.adjust < 0.05])  # Extract significant pathways
target_res_df <- res_df[res_df$ID %in% target_res, ]  # Subset significant results



# ---- Clean Descriptions ----

# Use string manipulation to clean up pathway descriptions
# Remove the 'REACTOME_' prefix and replace underscores ('_') with spaces
target_res_df <- target_res_df %>%
  mutate(Description = str_replace_all(Description, c('REACTOME_' = '', '_' = ' ')))



# ---- Select Upregulated and Downregulated Pathways ----

# Filter the enrichment results to separate upregulated pathways
res_df_up <- target_res_df %>% filter(diffexpressed == 'UP') %>%
  dplyr::select(!c('diffexpressed'))  # Remove the 'diffexpressed' column for upregulated pathways

# Filter the enrichment results to separate downregulated pathways
res_df_down <- target_res_df %>% filter(diffexpressed == 'DOWN') %>%
  dplyr::select(!c('diffexpressed'))  # Remove the 'diffexpressed' column for downregulated pathways

# Assign pathway IDs as row names for upregulated pathways
rownames(res_df_up) <- res_df_up$ID

# Assign pathway IDs as row names for downregulated pathways
rownames(res_df_down) <- res_df_down$ID


# ---- Export Enrichment Results ----

# Save upregulated enrichment results to an Excel file
write_xlsx(res_df_up, "Upregulated pathways.xlsx")

# Save downregulated enrichment results to an Excel file
write_xlsx(res_df_down, "Downregulated pathways.xlsx")  # Replace "Nk cells" with the relevant cell type if needed




# ---- Visualize Top Enriched Pathways ----

# Create a bar plot to visualize the top 10 enriched pathways (example for downregulated pathways)
ggplot(data = res_df_up, aes(x = size, y = pathway, fill = padj)) +  # Define data and aesthetics
  geom_bar(stat = "identity") +  # Use bar plot
  labs(title = "Upregulated Pathways",  # Add title
       x = "Count",  # X-axis label
       y = "Pathway",  # Y-axis label
       fill = "p.adjust") +  # Legend label for the fill scale
  scale_fill_gradient(low = "blue", high = "red") +  # Define color gradient for adjusted p-value
  theme_minimal()  # Apply minimal theme for clean visualization
