# =============================================================================
# Figure 1: CD8 T Cell Landscape in Human Carotid Atherosclerotic Plaques
# =============================================================================
#
# Description:
#   This script generates all panels for Figure 1, which characterises the
#   CD8 T cell landscape in human carotid atherosclerotic plaques using
#   single-cell RNA sequencing (scRNA-seq). It covers:
#     - UMAP visualisation of the full plaque immune atlas and the CD8 subset
#     - Gene signature scoring with UCell (senescence, exhaustion, cytotoxicity, etc.)
#     - Hexagon-binned UMAP overlays for each module score (Spectre)
#     - Spearman correlation matrix of module scores across clusters
#     - Export of source data tables (cluster abundance, mean scores, per-cell scores)
#
# Input files required:
#   - ~/Desktop/T2D_Myeloid/data/Giannarelli_carotid.rds
#       Seurat object: full immune-cell atlas of human carotid plaques
#       (published dataset, Giannarelli lab)
#   - ~/Desktop/T2D_CD8/data/CD8_T2D_ms_version.rds
#       Seurat object: CD8 T cell subset from carotid plaques (this study)
#   - ~/Desktop/T2D_CD8/figures/Figure_1/tables/Signatures_CD8_final.xlsx
#       Excel table with gene lists for each UCell module
#   - ~/Desktop/T2D_CD8/figures/Figure_1/tables/REACTOME_SASP.1.Hs.tsv
#       REACTOME Senescence-Associated Secretory Phenotype (SASP) gene set
#
# Output files — plots:
#   - plots/UMAP_all_plaques.pdf                            (Figure 1a)
#   - figures/Figure_1/plots/annotated.pdf                  (Figure 1b)
#   - figures/Figure_1/plots/<MODULE>_Hexagon_w_legend.pdf  (Figure 1d, one per module)
#   - figures/Figure_1/plots/corr_plot_*_UPDATED.pdf        (Figure 1f, several variants)
#
# Output files — source data / tables:
#   - figures/Figure_1/tables/cluster_abundance.xlsx            (Figure 1c)
#   - figures/Figure_1/tables/ModuleGene_AverageExpression.xlsx (Figure 1e)
#   - figures/Figure_1/tables/mean_score_per_celltype.xlsx      (Figures 1e & 1g)
#   - figures/Figure_1/tables/scores_per_cell.xlsx              (Figure 1h)
#
# =============================================================================

# --- Libraries ----------------------------------------------------------------
library(Seurat)         # scRNA-seq object handling and dimensionality reduction
library(SeuratDisk)     # reading/writing Seurat objects from/to disk formats
library(SeuratWrappers) # additional Seurat utility wrappers
library(ggplot2)        # general plotting
library(patchwork)      # combining multiple ggplots
library(dplyr)          # data manipulation
library(data.table)     # fast data.frame operations
library(readxl)         # reading Excel files
library(RColorBrewer)   # colour palettes (e.g. RdBu for correlations)
library(colorspace)     # colour manipulation (desaturation of palette)

set.seed(123) # set global random seed for reproducibility


# ==============================================================================
# Figure 1a — UMAP of the full carotid plaque immune atlas
# ==============================================================================

# Load the published carotid plaque Seurat object (Giannarelli et al.)
Giannarelli_carotid <- readRDS("~/Desktop/T2D_Myeloid/data/Giannarelli_carotid.rds")
unique(Giannarelli_carotid$annotation_major)  # inspect major cell type labels
DimPlot(Giannarelli_carotid, group.by = "annotation_major", label = T, repel = T)  # quick QC plot

# Define a fixed colour map for each major cell type (for consistency across figures)
caroid_color_map <- c(
  "B cell"       = "#cab2d6",
  "GD"           = "#ff7f00",
  "Myeloid"      = "#e31a1c",
  "Other T cell" = "#1f78b4",
  "ILC"          = "#fdbf6f",
  "NK cell"      = "#fb9a99",
  "CD4 T cell"   = "#a6cee3",
  "CD8 T cell"   = "#b2df8a",
  "DP T cell"    = "#33a02c"
)

# Export final UMAP (no axes, no legend, no title — formatted for figure)
DimPlot(Giannarelli_carotid, group.by = "annotation_major", cols = caroid_color_map, shuffle = T, alpha = 0.5,
        pt.size = 2) & NoAxes() & NoLegend() & labs(title = NULL)
ggsave("plots/UMAP_all_plaques.pdf",
       units = "cm", height = 28, width = 26)


# ==============================================================================
# Figure 1b — UMAP of the CD8 T cell subset with cluster annotations
# ==============================================================================

# Load the CD8 T cell Seurat object (this study)
CD8_T2D_ms_version <- readRDS("~/Desktop/T2D_CD8/data/CD8_T2D_ms_version.rds")
Idents(CD8_T2D_ms_version) <- "final_annotation"  # set active identity to cluster labels

# Generate a colour palette for 14 CD8 clusters by interpolating 8 seed colours
# Colours are slightly desaturated for a cleaner look
sampled_colors <- c("#BC3C2999", "#0072B599", "#FFDC9199", "#20854E99", "#6F99AD99", "#EE4C9799", "#E1872799", "#6a3d9a")
sampled_colors <- desaturate(sampled_colors, amount = 0.25)
num_colors <- 14
color_palette <- colorRampPalette(sampled_colors)
created_colors <- color_palette(num_colors)
print(created_colors)  # inspect generated colours
cell_populations <- levels(CD8_T2D_ms_version)
color_mapping <- as.list(setNames(created_colors, cell_populations))  # named list: cluster -> colour

DimPlot(CD8_T2D_ms_version, label = F, group.by = "final_annotation", repel = F, cols = color_mapping, pt.size = 2)  & NoLegend() & NoAxes() & labs(title = NULL)
ggsave("~/Desktop/T2D_CD8/figures/Figure_1/plots/annotated.pdf",  units = "cm", height = 22, width = 26)


# ==============================================================================
# Figure 1c — Cluster abundance per patient (source data)
# ==============================================================================

# Extract metadata, count cells per patient per cluster, and pivot to wide format
md <- CD8_T2D_ms_version@meta.data %>% as.data.table
cluster_abundance <- md[, .N, by = c("Pat_ID", "final_annotation")] %>%
  dcast(., Pat_ID ~ final_annotation, value.var = "N")
openxlsx::write.xlsx(cluster_abundance, "~/Desktop/T2D_CD8/figures/Figure_1/tables/cluster_abundance.xlsx", asTable =F, rowNames =T)


# ==============================================================================
# Figure 1d — UCell gene signature scoring and hexagon-binned UMAP overlays
# ==============================================================================

# --- Load gene signatures -----------------------------------------------------
library(UCell)  # AUCell-based gene module scoring for single cells

# Read gene signature table (one column per module, gene names as rows)
Signatures_CD8_final <- read_excel("~/Desktop/T2D_CD8/figures/Figure_1/tables/Signatures_CD8_final.xlsx")
SENMAYO           <- Signatures_CD8_final$SENMAYO
EXHAUSTION        <- Signatures_CD8_final$EXHAUSTION
T_CELL_SENESCENCE <- Signatures_CD8_final$`T CELL SENESCENCE`
CYTOTOXIC         <- Signatures_CD8_final$CYTOTOXIC
TERM_MEMORY       <- Signatures_CD8_final$TERM_MEMORY
PROLIFERATION     <- Signatures_CD8_final$PROLIFERATION
NAIVE             <- Signatures_CD8_final$NAIVE
MEMORY            <- Signatures_CD8_final$MEMORY
RESIDENT          <- Signatures_CD8_final$RESIDENT

# Load the REACTOME SASP gene set
# The gene list is stored as a comma-separated string in row 17 of the TSV
SASP <- read.table("~/Desktop/T2D_CD8/figures/Figure_1/tables/REACTOME_SASP.1.Hs.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
SASP <- as.character(SASP[17, 2])
SASP <- unlist(strsplit(SASP, split = ","))
SASP <- unique(SASP)  # remove duplicates

# Bundle all signatures into a named list for UCell scoring
signatures <- list(
  SENMAYO          = c(SENMAYO),
  EXHAUSTION       = c(EXHAUSTION),
  T_CELL_SENESCENCE = c(T_CELL_SENESCENCE),
  CYTOTOXIC        = c(CYTOTOXIC),
  SASP             = c(SASP),
  TERM_MEMORY      = c(TERM_MEMORY),
  PROLIFERATION    = c(PROLIFERATION),
  NAIVE            = c(NAIVE),
  MEMORY           = c(MEMORY),
  RESIDENT         = c(RESIDENT))

# Score each cell for all modules using UCell (kNN-smoothed AUC scores)
DefaultAssay(CD8_T2D_ms_version) <- "RNA"
CD8_T2D_ms_version <- AddModuleScore_UCell(CD8_T2D_ms_version, features = signatures, assay = "RNA")
names(signatures)  # confirm all modules were scored

# --- Export hexagon-binned UMAP overlays with Spectre -------------------------
library(Spectre)       # flow/CyTOF analysis toolkit; used here for hex-bin UMAP plots
library(RColorBrewer)

# Extract UMAP coordinates as a data.frame
umap_coordinates <- as.data.frame(Embeddings(CD8_T2D_ms_version, reduction = "umap"))

# Extract kNN-smoothed UCell scores from metadata for each module
senmayo_scores        <- CD8_T2D_ms_version@meta.data$SENMAYO_UCell_kNN
exhaustion_scores     <- CD8_T2D_ms_version@meta.data$EXHAUSTION_UCell_kNN
t_senescence_scores   <- CD8_T2D_ms_version@meta.data$T_CELL_SENESCENCE_UCell_kNN
cytotoxic_scores      <- CD8_T2D_ms_version@meta.data$CYTOTOXIC_UCell_kNN
sasp_scores           <- CD8_T2D_ms_version@meta.data$SASP_UCell_kNN
effector_memroy_scores <- CD8_T2D_ms_version@meta.data$TERM_MEMORY_UCell_kNN
proliferation_scores  <- CD8_T2D_ms_version@meta.data$PROLIFERATION_UCell_kNN
naive_scores          <- CD8_T2D_ms_version@meta.data$NAIVE_UCell_kNN
memory_scores         <- CD8_T2D_ms_version@meta.data$MEMORY_UCell_kNN
resident_scores       <- CD8_T2D_ms_version@meta.data$RESIDENT_UCell_kNN

# Combine UMAP coordinates and all module scores into a single data.table
dat <- cbind(umap_coordinates,
             SENMAYO         = senmayo_scores,
             EXHAUSTION      = exhaustion_scores,
             T_CELL_SENESCENCE = t_senescence_scores,
             CYTOTOXIC       = cytotoxic_scores,
             SASP            = sasp_scores,
             TERM_MEMORY     = effector_memroy_scores,
             PROLFERATION    = proliferation_scores,   # note: intentional typo kept for column consistency
             NAIVE           = naive_scores,
             MEMORY          = memory_scores,
             RESIDENT        = resident_scores)
dat <- as.data.table(dat)

# Generate and save hexagon-binned UMAP plots for each module score
# Using RdBu palette: red = high score, blue = low score

# SENMAYO
p1 <- make.colour.plot(dat = dat,
                       x.axis = "UMAP_1",
                       y.axis = "UMAP_2",
                       col.axis = "SENMAYO",
                       dot.size = 1.5,
                       hex = TRUE,
                       save.to.disk = FALSE) +
  scale_fill_distiller(palette = "RdBu", direction = -1) & NoAxes()
ggsave("~/Desktop/T2D_CD8/figures/Figure_1/plots/SENMAYO_Hexagon_w_legend.pdf",
       plot = p1, units = "cm", height = 48, width = 48)

# EXHAUSTION
p2 <- make.colour.plot(dat = dat,
                       x.axis = "UMAP_1",
                       y.axis = "UMAP_2",
                       col.axis = "EXHAUSTION",
                       dot.size = 1.5,
                       hex = TRUE,
                       save.to.disk = FALSE) +
  scale_fill_distiller(palette = "RdBu", direction = -1) & NoAxes()
ggsave("~/Desktop/T2D_CD8/figures/Figure_1/plots/EXHAUSTION_Hexagon_w_legend.pdf",
       plot = p2, units = "cm", height = 48, width = 48)

# T CELL SENESCENCE
p3 <- make.colour.plot(dat = dat,
                       x.axis = "UMAP_1",
                       y.axis = "UMAP_2",
                       col.axis = "T_CELL_SENESCENCE",
                       dot.size = 1.5,
                       hex = TRUE,
                       save.to.disk = FALSE) +
  scale_fill_distiller(palette = "RdBu", direction = -1) & NoAxes()
ggsave("~/Desktop/T2D_CD8/figures/Figure_1/plots/T_CELL_SENESCENCE_Hexagon_w_legend.pdf",
       plot = p3, units = "cm", height = 48, width = 48)

# CYTOTOXIC
p4 <- make.colour.plot(dat = dat,
                       x.axis = "UMAP_1",
                       y.axis = "UMAP_2",
                       col.axis = "CYTOTOXIC",
                       dot.size = 1.5,
                       hex = TRUE,
                       save.to.disk = FALSE) +
  scale_fill_distiller(palette = "RdBu", direction = -1) & NoAxes()
ggsave("~/Desktop/T2D_CD8/figures/Figure_1/plots/CYTOTOXIC_Hexagon_w_legend.pdf",
       plot = p4, units = "cm", height = 48, width = 48)

# SASP
p5 <- make.colour.plot(dat = dat,
                       x.axis = "UMAP_1",
                       y.axis = "UMAP_2",
                       col.axis = "SASP",
                       dot.size = 1.5,
                       hex = TRUE,
                       save.to.disk = FALSE) +
  scale_fill_distiller(palette = "RdBu", direction = -1) & NoAxes()
ggsave("~/Desktop/T2D_CD8/figures/Figure_1/plots/SASP_Hexagon_w_legend.pdf",
       plot = p5, units = "cm", height = 48, width = 48)

# TERM_MEMORY (terminally differentiated effector memory)
p6 <- make.colour.plot(dat = dat,
                       x.axis = "UMAP_1",
                       y.axis = "UMAP_2",
                       col.axis = "TERM_MEMORY",
                       dot.size = 1.5,
                       hex = TRUE,
                       save.to.disk = FALSE) +
  scale_fill_distiller(palette = "RdBu", direction = -1) & NoAxes()
ggsave("~/Desktop/T2D_CD8/figures/Figure_1/plots/TERM_MEMORY_Hexagon_w_legend.pdf",
       plot = p6, units = "cm", height = 48, width = 48)

# PROLIFERATION
p7 <- make.colour.plot(dat = dat,
                       x.axis = "UMAP_1",
                       y.axis = "UMAP_2",
                       col.axis = "PROLFERATION",  # note: column name carries forward the typo from cbind above
                       dot.size = 1.5,
                       hex = TRUE,
                       save.to.disk = FALSE) +
  scale_fill_distiller(palette = "RdBu", direction = -1) & NoAxes()
ggsave("~/Desktop/T2D_CD8/figures/Figure_1/plots/PROLFERATION_Hexagon_w_legend.pdf",
       plot = p7, units = "cm", height = 48, width = 48)

# NAIVE
p8 <- make.colour.plot(dat = dat,
                       x.axis = "UMAP_1",
                       y.axis = "UMAP_2",
                       col.axis = "NAIVE",
                       dot.size = 1.5,
                       hex = TRUE,
                       save.to.disk = FALSE) +
  scale_fill_distiller(palette = "RdBu", direction = -1) & NoAxes()
ggsave("~/Desktop/T2D_CD8/figures/Figure_1/plots/NAIVE_Hexagon_w_legend.pdf",
       plot = p8, units = "cm", height = 48, width = 48)

# MEMORY
p9 <- make.colour.plot(dat = dat,
                       x.axis = "UMAP_1",
                       y.axis = "UMAP_2",
                       col.axis = "MEMORY",
                       dot.size = 1.5,
                       hex = TRUE,
                       save.to.disk = FALSE) +
  scale_fill_distiller(palette = "RdBu", direction = -1) & NoAxes()
ggsave("~/Desktop/T2D_CD8/figures/Figure_1/plots/MEMORY_Hexagon_w_legend.pdf",
       plot = p9, units = "cm", height = 48, width = 48)

# RESIDENT (tissue-resident memory)
p10 <- make.colour.plot(dat = dat,
                        x.axis = "UMAP_1",
                        y.axis = "UMAP_2",
                        col.axis = "RESIDENT",
                        dot.size = 1.5,
                        hex = TRUE,
                        save.to.disk = FALSE) +
  scale_fill_distiller(palette = "RdBu", direction = -1) & NoAxes()
ggsave("~/Desktop/T2D_CD8/figures/Figure_1/plots/RESIDENT_Hexagon_w_legend.pdf",
       plot = p10, units = "cm", height = 48, width = 48)


# ==============================================================================
# Figure 1e — Average gene expression per module (source data)
# ==============================================================================

# Calculate the mean expression of each signature gene across all cells
DefaultAssay(CD8_T2D_ms_version) <- "RNA"
all_genes    <- unique(unlist(signatures))  # all unique genes across all modules
genes_present <- all_genes[all_genes %in% rownames(CD8_T2D_ms_version)]  # keep only detected genes

# Add a dummy grouping variable to compute a single average across the whole object
CD8_T2D_ms_version$wholeobject <- "wholeobject"
avg_list <- AverageExpression(
  CD8_T2D_ms_version,
  assays   = "RNA",
  features = genes_present,
  group.by = "wholeobject"
)
avg_mat <- avg_list$RNA  # RNA assay average expression matrix

# Create a long-format table: Module | Gene | Average expression
gene_module_df <- tibble::enframe(signatures, name = "Module", value = "Gene") %>%
  tidyr::unnest(cols = Gene) %>%
  mutate(Gene = as.character(Gene)) %>%
  filter(Gene %in% genes_present)
avg_df <- as.data.frame(avg_mat) %>%
  tibble::rownames_to_column("Gene")
out_df <- gene_module_df %>%
  left_join(avg_df, by = "Gene") %>%
  arrange(Module, Gene)
openxlsx::write.xlsx(out_df,
                     file = "~/Desktop/T2D_CD8/figures/Figure_1/tables/ModuleGene_AverageExpression.xlsx", rowNames = FALSE)

# Calculate mean UCell score per cluster (used for the heatmap in Figure 1e)

#subset md
md <- CD8_T2D_ms_version@meta.data %>% as.data.table
head(md)
subset_md <- md[, .SD, .SDcols = c("Pat_ID", "conditions", "final_annotation", grep("UCell_kNN", colnames(md), value = TRUE))]
mean_scores_by_cell_type <- subset_md %>%
  group_by(final_annotation) %>%
  summarise(across(contains("UCell_kNN"), \(x) mean(x, na.rm = TRUE)))
openxlsx::write.xlsx(mean_scores_by_cell_type,
                     "~/Desktop/T2D_CD8/figures/Figure_1/tables/mean_score_per_celltype.xlsx")


# ==============================================================================
# Figure 1f — Spearman correlation matrix of module scores across clusters
# ==============================================================================

library(dplyr)
library(data.table)
library(corrplot)
library(RColorBrewer)

# Extract metadata as a data.table
md <- as.data.table(CD8_T2D_ms_version@meta.data)
md$signature_1SENEPY_SCORE_UCELL_kNN  # inspect the SENEPY combined signature column

# Define module score columns of interest
scores <- c(
  "SENMAYO_UCell_kNN",
  "EXHAUSTION_UCell_kNN",
  "T_CELL_SENESCENCE_UCell_kNN",
  "CYTOTOXIC_UCell_kNN",
  "SASP_UCell_kNN",
  "TERM_MEMORY_UCell_kNN",
  "PROLIFERATION_UCell_kNN",
  "NAIVE_UCell_kNN",
  "MEMORY_UCell_kNN",
  "RESIDENT_UCell_kNN",
  "signature_1SENEPY_SCORE_UCELL_kNN"  # composite senescence-cytotoxicity signature
)

# Guard against missing columns (e.g. if scores were added in a different run)
scores_present <- intersect(scores, colnames(md))
stopifnot(length(scores_present) >= 2)

# Compute mean score per cluster — correlation is run on cluster-level means,
# not individual cells, to reduce noise from cell-to-cell variation
mean_by_cluster <- md[, lapply(.SD, mean, na.rm = TRUE),
                      by = .(final_annotation),
                      .SDcols = scores_present] %>%
  as.data.frame()

# Keep only numeric score columns; drop any with zero variance (would produce NA in cor())
corr_data <- mean_by_cluster %>%
  select(all_of(scores_present)) %>%
  select(where(~ sd(.x, na.rm = TRUE) > 0))

# Compute Spearman correlation matrix (robust to non-normal distributions)
cor_matrix <- cor(corr_data, method = "spearman", use = "pairwise.complete.obs")

# Shared colour palette for all corrplot variants
color_palette <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(200))

# --- Pie chart corrplot -------------------------------------------------------
pdf("~/Desktop/T2D_CD8/figures/Figure_1/plots/corr_plot_pies_clusterMeans_UPDATED.pdf",
    width = 32/2.54, height = 32/2.54)
corrplot(
  cor_matrix,
  method  = "pie",
  type    = "upper",
  col     = color_palette,
  tl.col  = "black",
  tl.cex  = 1,
  tl.srt  = 45,
  cl.cex  = 0.8
)
dev.off()

# --- Circle corrplot with correlation coefficients ----------------------------
pdf("~/Desktop/T2D_CD8/figures/Figure_1/plots/corr_plot_numbers_clusterMeans_circle_UPDATED.pdf",
    width = 32/2.54, height = 32/2.54)
corrplot(
  cor_matrix,
  method       = "circle",
  type         = "upper",
  col          = color_palette,
  tl.col       = "black",
  tl.cex       = 1,
  tl.srt       = 45,
  cl.cex       = 0.8,
  addCoef.col  = "white"
)
dev.off()

# --- Square corrplot with correlation coefficients ----------------------------
pdf("~/Desktop/T2D_CD8/figures/Figure_1/plots/corr_plot_numbers_clusterMeans_square_UPDATED.pdf",
    width = 32/2.54, height = 32/2.54)
corrplot(
  cor_matrix,
  method      = "square",
  type        = "upper",
  col         = color_palette,
  tl.col      = "black",
  tl.cex      = 1,
  tl.srt      = 45,
  cl.cex      = 0.8,
  addCoef.col = "white"
)
dev.off()

# --- Minimal circle corrplot (no axis labels, for figure inset) ---------------
pdf("~/Desktop/T2D_CD8/figures/Figure_1/plots/corr_plot_numbers_clusterMeans_circle_no_labels_UPDATED.pdf",
    width = 32/2.54, height = 32/2.54)
corrplot(
  cor_matrix,
  method      = "circle",
  type        = "upper",
  col         = color_palette,
  tl.pos      = "n",      # suppress text labels
  cl.pos      = "n",      # suppress colour legend
  addCoef.col = "white",
  number.cex  = 2
)
dev.off()


# ==============================================================================
# Figure 1g — Mean UCell score per cell type (source data)
# ==============================================================================

# Subset metadata to patient ID, condition, cluster annotation, and all kNN-smoothed scores
md <- CD8_T2D_ms_version@meta.data %>% as.data.table
head(md)
subset_md <- md[, .SD, .SDcols = c("Pat_ID", "conditions", "final_annotation", grep("_kNN", colnames(md), value = TRUE))]

# Compute mean score per cluster across all cells
mean_scores_by_cell_type <- subset_md %>%
  group_by(final_annotation) %>%
  summarise(across(contains("UCell_kNN"), \(x) mean(x, na.rm = TRUE)))
openxlsx::write.xlsx(mean_scores_by_cell_type,
                     "~/Desktop/T2D_CD8/figures/Figure_1/tables/mean_score_per_celltype.xlsx")


# ==============================================================================
# Figure 1h — Per-cell UCell scores (source data)
# ==============================================================================

# Export single-cell level scores for statistical testing in downstream tools
md <- CD8_T2D_ms_version@meta.data %>% as.data.table

# Retain patient ID, condition, cluster label, and all raw UCell scores
head(md)
scores_per_cell <- md[, .SD, .SDcols = c("Pat_ID", "conditions", "final_annotation", grep("UCell", colnames(md), value = TRUE))]
head(scores_per_cell)

openxlsx::write.xlsx(scores_per_cell,
                     "~/Desktop/T2D_CD8/figures/Figure_1/tables/scores_per_cell.xlsx")
