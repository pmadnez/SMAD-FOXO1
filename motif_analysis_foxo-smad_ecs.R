library(Signac)
library(Seurat)
library(SummarizedExperiment)
library(GenomicRanges)
library(JASPAR2020)
library(BSgenome.Mmusculus.UCSC.mm10)
library(chromVAR)
library(chromVARmotifs)
library(motifmatchr)
library(TFBSTools)
setwd("Z:/data/")
#load pre-processed Seurat object
final_set <- readRDS(file = "final-set_with_daa.rds")
main.chroms <- standardChromosomes(BSgenome.Mmusculus.UCSC.mm10)
keep.peaks <- as.logical(seqnames(granges(final_set)) %in% main.chroms)
final_set[["peaks"]] <- subset(final_set[["peaks"]], features = rownames(final_set[["peaks"]])[keep.peaks])

###subset for ECs
endothelial_cells <- subset(final_set, idents = c("EC1", "EC2", "EC3", "EC4", "EC5", "EC6", "EC7"))

###get motifs

#Extract ATAC-seq count matrix from Seurat object
counts_matrix <- GetAssayData(endothelial_cells, assay = "peaks", slot = "counts")
# Filter out peaks (rows) where there are no counts in any samples
non_zero_peaks <- rowSums(counts_matrix) > 0

##filter seurat object accordingly
endothelial_cells_filtered <- subset(endothelial_cells, features = rownames(counts_matrix)[non_zero_peaks])
# Subset the SummarizedExperiment object to keep only non-zero peaks
counts_matrix_filtered <- counts_matrix[non_zero_peaks, ]
# Step 2: Extract peak information (GRanges) from Seurat object
peakSet <- granges(endothelial_cells_filtered)

#Create a SummarizedExperiment object for chromVAR
chromVarCounts <- SummarizedExperiment(assays = list(counts = counts_matrix_filtered), rowRanges = peakSet)

# Add GC bias to the object (needed for chromVAR)
chromVarCounts <- addGCBias(chromVarCounts, genome = BSgenome.Mmusculus.UCSC.mm10)

###Get JASPAR PFMs for desired motifs:
motif_ids <- c("MA0535.1","MA0480.1")  # Replace with your list of motif IDs
motif_ids
motifs <- lapply(motif_ids, function(id) getMatrixByID(JASPAR2020, ID = id))
all(sapply(motifs, function(m) class(m) == "PFMatrix"))

# Convert the list of motifs into a PFMatrixList
motifs_combined <- do.call(PFMatrixList, motifs)


# Match these motifs to the peak set
motif_ix <- matchMotifs(motifs_combined, peakSet, genome = BSgenome.Mmusculus.UCSC.mm10)


# Run chromVAR again with the filtered object
deviations <- computeDeviations(object = chromVarCounts, annotations = motif_ix)


# Get the motif activity (deviation scores) from the deviations object
motif_activities <- assay(deviations)  # Deviation scores are stored in an assay

# Get the cell type annotations from your Seurat object
cell_types <- endothelial_cells_filtered@meta.data$predicted.id

# Ensure cell types match the number of cells in motif_activities
stopifnot(length(cell_types) == ncol(motif_activities))

# Combine motif activities with cell type annotations into a data frame
library(tidyverse)

# Convert motif activities to data frame and transpose for easier plotting
motif_activity_df <- as.data.frame(t(motif_activities)) %>%
  mutate(CellType = cell_types)

# Check the structure of the motif activity data frame
head(motif_activity_df)
names(motif_activity_df) <- c("MA0535.1", "MA0480.1", "CellType")

## make violin plot to display motif activity
ggplot(motif_activity_df, aes(x = CellType, y = MA0535.1, fill = CellType)) +
  geom_violin() +
  geom_boxplot(width = 0.1, color = "black", alpha = 0.5) + 
  theme_minimal() +
  labs(title = "Foxo1 Motif Activity by Cell Type",
       x = "Cell Type", y = "Motif Activity (Deviation)") +
  scale_y_continuous(limits = c(-0.1, 0.1)) +  # Adjust limits based on your data
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Add motif activity to Seurat metadata
endothelial_cells_filtered$Smad1_Motif_Activity <- motif_activity_df$MA0535.1
endothelial_cells_filtered$Foxo1_Motif_Activity <- motif_activity_df$MA0480.1
# Extract UMAP coordinates
umap_data <- as.data.frame(Embeddings(endothelial_cells_filtered, "umap"))
umap_data$CellType <- endothelial_cells_filtered$predicted.id  # Add cell type information
umap_data$Foxo1_Motif_Activity <- endothelial_cells_filtered$Foxo1_Motif_Activity  # Add motif activity
umap_data$Smad1_Motif_Activity <- endothelial_cells_filtered$Smad1_Motif_Activity


##make umap plots of cell clusters and motif activities

#create small axis for umap plots
axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(2, "cm")
)

umap_data <- umap_data %>%
  mutate(CellType = recode(CellType, "EC8" = "EC7"))
# Plot cell cluster umap
p1 <- ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = CellType)) +
  geom_point(size = 1.25) +  # Adjust point size
  theme_minimal() +
  guides(x = axis, y = axis) +
  theme(axis.line = element_line(arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size = 1),
        axis.title = element_text(hjust = 0),
        aspect.ratio = 1)+
  labs(
    x = "UMAP1",
    y = "UMAP 2"
  ) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  NULL
p1
p2 <- ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = Foxo1_Motif_Activity)) +
  geom_point(size = 1.25, alpha = 1) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-0.08, 0.08)) +
  theme_minimal() +
  guides(x = axis, y = axis) +
  theme(axis.line = element_line(arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size = 1),
        axis.title = element_text(hjust = 0),
        aspect.ratio = 1)+
  labs(
    x = "UMAP1",
    y = "UMAP 2"
  ) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  NULL
p2
p3 <- ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = Smad1_Motif_Activity)) +
  geom_point(size = 1.25, alpha = 1) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-0.06, 0.06)) +
  theme_minimal() +
  guides(x = axis, y = axis) +
  theme(axis.line = element_line(arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size = 1),
        axis.title = element_text(hjust = 0),
        aspect.ratio = 1)+
  labs(
    x = "UMAP1",
    y = "UMAP 2"
  ) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  NULL
p3
# Save plots
ggsave("umap_EC-cell-types_final.png", plot = p1, dpi = 300)
ggsave("umap_FOXO-activity_final.png", plot = p2, dpi = 300)
ggsave("umap_SMAD-activity_final.png", plot = p3, dpi = 300)

##make plot of cell proportions per condition
# Extract cell types and conditions from metadata
cell_types <- endothelial_cells_filtered$predicted.id
conditions <- endothelial_cells_filtered$orig.ident  # Use active.ident for conditions

# Create a data frame with cell types and conditions
data_for_proportions <- data.frame(CellType = cell_types, Condition = conditions)


proportions <- data_for_proportions %>%
  group_by(Condition, CellType) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  mutate(Proportion = Count / sum(Count))
# Create the bar plot
p4 <- ggplot(proportions, aes(x = Condition, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  labs(y = "Ratio cell type/ all cells") +
  scale_y_continuous(labels = scales::percent) +  # Show proportions as percentages
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Display the plot
print(p4)

