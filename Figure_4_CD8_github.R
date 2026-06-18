# =============================================================================
# Figure 4: Pseudotime Trajectory Analysis of CD8 T Cell Differentiation
# =============================================================================
#
# Description:
#   This script infers and visualises differentiation trajectories within the
#   CD8 T cell atlas using Slingshot (pseudotime) and tradeSeq (trajectory-
#   associated gene expression). It covers:
#     - Slingshot trajectory inference from the Naïve cluster as root state,
#       producing multiple lineages across the UMAP
#     - Per-lineage pseudotime UMAP overlays
#     - Module score dynamics along pseudotime (loess smoothed)
#     - tradeSeq generalised additive model (GAM) fitting for trajectory-
#       associated transcription factor expression
#     - Heatmap of TF expression trends along the trajectory
#     - Individual TF expression curves and Nebulosa density plot of ZEB2
#
# Input files required:
#   - ~/Desktop/T2D_CD8/data/CD8_T2D_ms_version_final_modules.rds
#       Seurat object: CD8 T cells with final cluster annotations and UCell scores
#
# Output files — plots:
#   - plots/all_lineages.pdf                       (Figure 4a)
#   - plots/lineage_<name>.pdf                     (Figure 4b, per lineage, with legend)
#   - plots/lineage_<name>_no_legend.pdf           (Figure 4b, per lineage, clean)
#   - plots/<MODULE>_Module_Pseudotime.pdf         (Figure 4c, per module)
#   - plots/<MODULE>_Module_Pseudotime_clean.pdf   (Figure 4c, per module, no axes)
#   - plots/Heatmap.jpg                            (Figure 4d)
#   - plots/Gene_Expression_Pseudotime.pdf         (Figure 4e)
#   - plots/ZEB2_Nebulosa_fixed.pdf                (Figure 4f)
#
# Output files — data objects / tables:
#   - data/gam_final.rds            (fitted tradeSeq GAM model)
#   - tables/at_sig_lin2.xlsx       (trajectory-associated genes for lineage 2)
#
# =============================================================================

# --- Libraries ----------------------------------------------------------------
library(Seurat)             # scRNA-seq object handling
library(slingshot)          # pseudotime / trajectory inference
library(tradeSeq)           # trajectory-based differential expression (fitGAM)
library(TrajectoryUtils)    # trajectory data structure utilities
library(SingleCellExperiment) # Bioconductor data container
library(RColorBrewer)       # colour palettes
library(scales)             # axis scaling helpers
library(viridis)            # perceptually uniform colour scales (inferno for pseudotime)
library(ggplot2)            # plotting
library(dplyr)              # data manipulation
library(ggthemes)           # additional ggplot themes (tableau colour scale)
library(tidyr)              # data reshaping (pivot_longer / gather)
set.seed(123)               # global random seed for reproducibility

#setwd 
setwd("~/Desktop/T2D_CD8/figures/Figure_4")


# ==============================================================================
# Figure 4a — Slingshot trajectory inference and UMAP overlay
# ==============================================================================

# Generate colour palette for 14 CD8 clusters (same as Figure 1b)
sampled_colors <- c("#BC3C2999", "#0072B599", "#FFDC9199", "#20854E99", "#6F99AD99", "#EE4C9799", "#E1872799", "#6a3d9a")
num_colors     <- 14
color_palette  <- colorRampPalette(sampled_colors)
created_colors <- color_palette(num_colors)
print(created_colors)

# Load CD8 Seurat object with final cluster annotations and module scores
CD8_T2D_ms_version <- readRDS("~/Desktop/T2D_CD8/data/CD8_T2D_ms_version_final_modules.rds")
DimPlot(CD8_T2D_ms_version, group.by = "final_annotation", label = T, repel = T, cols = created_colors)

# Convert to SingleCellExperiment (required by Slingshot)
CD8_sce <- as.SingleCellExperiment(CD8_T2D_ms_version, assay = "RNA")

# Inspect available reduced dimensions
reducedDimNames(CD8_sce)

# Preview PCA layout
pca <- reducedDim(CD8_sce, "PCA")[, 1:2]
plot(pca, col = rgb(0,0,0,.5), pch=16, asp = 1)

# Preview UMAP layout
UMAP <- reducedDim(CD8_sce, "UMAP")[, 1:2]
plot(UMAP, col = rgb(0,0,0,.5), pch=16, asp = 1)

# Run Slingshot trajectory inference on the UMAP embedding
# start.clus = Naïve cluster (biologically the progenitor state)
# stretch = 1 allows curves to extend slightly beyond the data cloud
unique(CD8_sce@colData@listData[["final_annotation"]])
CD8_sce <- slingshot(CD8_sce, clusterLabels = 'final_annotation', reducedDim = 'UMAP',
                     start.clus = "C1- CCR7+ IL7R+ Naïve", stretch = 1)

# Colour cells by pseudotime along the first lineage for a quick QC plot
colors  <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(CD8_sce$slingPseudotime_1, breaks=100)]

plot(reducedDims(CD8_sce)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(CD8_sce), lwd=2, col='black')

# Plot all lineages with cluster colours (for visual annotation)
unique(CD8_sce@colData@listData[["final_annotation"]])
plot(reducedDims(CD8_sce)$UMAP, col = created_colors[CD8_sce$final_annotation], pch=16, asp = 1)
lines(SlingshotDataSet(CD8_sce), lwd=2, type = 'lineages', col = 'black')  # lineage connectors only
plot(reducedDims(CD8_sce)$UMAP, col = created_colors[CD8_sce$final_annotation], pch=16, asp = 1)
lines(SlingshotDataSet(CD8_sce), lwd=2, col='black')  # smooth curves

# Export final version for figure (no axes, no box, 12 × 12 cm)
pdf("plots/all_lineages.pdf", width = 12 / 2.54, height = 12 / 2.54)
plot(
  reducedDims(CD8_sce)$UMAP,
  col   = created_colors[CD8_sce$final_annotation],
  pch   = 16,
  cex   = 0.5,
  asp   = 1.5,
  xaxt  = 'n',  # suppress x axis
  yaxt  = 'n',  # suppress y axis
  xlab  = '',
  ylab  = '',
  bty   = 'n'   # no surrounding box
)
lines(SlingshotDataSet(CD8_sce), lwd = 2, col = 'black')
dev.off()


# ==============================================================================
# Figure 4b — Per-lineage pseudotime UMAP overlays
# ==============================================================================

# Extract smoothed lineage curves and pseudotime matrix
curves     <- slingCurves(CD8_sce)    # one curve per lineage
pt         <- slingPseudotime(CD8_sce) # cells × lineages pseudotime matrix
nms        <- colnames(pt)             # lineage names

# UMAP coordinates for coloured scatter plot
umap_coords             <- as.data.frame(reducedDim(CD8_sce, type = "UMAP"))
colnames(umap_coords)   <- c("UMAP1", "UMAP2")

plots           <- list()
plots_no_legend <- list()

# Generate one plot per lineage (with and without legend/title)
for (i in seq_along(nms)) {
  lineage <- nms[i]

  # Assign pseudotime for this lineage; cells not on this lineage get NA (shown in grey)
  umap_coords$pseudotime <- pt[, lineage]
  umap_coords$pseudotime[is.na(umap_coords$pseudotime)] <- NA

  # Full version: legend + title
  p <- ggplot(umap_coords, aes(x = UMAP1, y = UMAP2, color = pseudotime)) +
    geom_point(size = 1.5, alpha = 0.8) +
    scale_color_viridis(option = "inferno", na.value = "gray80", name = "Pseudotime") +
    theme_void() +
    theme(
      legend.position = "right",
      plot.title      = element_text(hjust = 0.5, size = 14)
    ) +
    ggtitle(paste("Lineage:", lineage)) +
    coord_fixed(ratio = 1.5)

  # Overlay the smoothed lineage curve as a black path
  lineage_curve         <- as.data.frame(curves[[i]]$s)
  colnames(lineage_curve) <- c("UMAP1", "UMAP2")
  p <- p + geom_path(data = lineage_curve, aes(x = UMAP1, y = UMAP2),
                     color = "black", size = 1)

  plots[[lineage]] <- p

  # Clean version: no legend, no title (for figure panels)
  p_no_legend <- p +
    theme(legend.position = "none", plot.title = element_blank())
  plots_no_legend[[lineage]] <- p_no_legend
}

# Save all versions to disk
for (name in names(plots)) {
  ggsave(filename = file.path("plots", paste0("lineage_", name, ".pdf")),
         plot = plots[[name]], width = 12, height = 12, units = "cm")
}
for (name in names(plots_no_legend)) {
  ggsave(filename = file.path("plots", paste0("lineage_", name, "_no_legend.pdf")),
         plot = plots_no_legend[[name]], width = 12, height = 12, units = "cm")
}


# ==============================================================================
# Figure 4c — Module scores along pseudotime (loess smoothed)
# ==============================================================================

# Build a data frame: pseudotime (lineage 2) + all module scores + cluster annotation
# Lineage 2 is selected as it captures the EMRA/senescent differentiation axis
expression_data_long <- data.frame(
  Pseudotime        = pt[, 2],
  NAIVE             = CD8_sce$NAIVE_UCell_kNN,
  MEMORY            = CD8_sce$MEMORY_UCell_kNN,
  RESIDENT          = CD8_sce$RESIDENT_UCell_kNN,
  TERM_MEM          = CD8_sce$TERM_MEMORY_UCell_kNN,
  T_CELL_SENESCENCE = CD8_sce$T_CELL_SENESCENCE_UCell_kNN,
  CYTOTOXIC         = CD8_sce$CYTOTOXIC_UCell_kNN,
  SENMAYO           = CD8_sce$SENMAYO_UCell_kNN,
  SASP              = CD8_sce$SASP_UCell_kNN,
  SENEPY            = CD8_sce$signature_1SENEPY_SCORE_UCELL_kNN,
  PROLIFERATION     = CD8_sce$PROLIFERATION_UCell_kNN,
  EXHAUSTION        = CD8_sce$EXHAUSTION_UCell_kNN,
  Annotation        = CD8_sce$final_annotation
)
head(expression_data_long)

# Pivot to long format for faceted ggplot (one row per cell per module)
expression_data_long <- expression_data_long %>%
  tidyr::pivot_longer(cols = c(NAIVE, MEMORY, RESIDENT, TERM_MEM, T_CELL_SENESCENCE, CYTOTOXIC,
                               SENMAYO, SENEPY, SASP, PROLIFERATION, EXHAUSTION),
                      names_to  = "Module",
                      values_to = "Expression")

# Recreate colour mapping for cluster annotations
Idents(CD8_T2D_ms_version) <- "final_annotation"
sampled_colors <- c("#BC3C2999", "#0072B599", "#FFDC9199", "#20854E99", "#6F99AD99", "#EE4C9799", "#E1872799", "#6a3d9a")
num_colors     <- 14
color_palette  <- colorRampPalette(sampled_colors)
created_colors <- color_palette(num_colors)
print(created_colors)
cell_populations <- levels(CD8_T2D_ms_version)
color_mapping    <- as.list(setNames(created_colors, cell_populations))

# Points coloured by cluster, loess trend line in black — all modules overlaid
p3 <- ggplot(expression_data_long, aes(x = Pseudotime, y = Expression, color = Annotation)) +
  geom_point(alpha = 0.7, size = 1.5) +
  geom_smooth(se = TRUE, method = "loess", color = "black") +  # single overall loess
  theme_classic() +
  xlab("Pseudotime") +
  ylab("Module Score") +
  ggtitle("ApoB Module Score Along Slingshot Pseudotime") +
  scale_color_manual(values = color_mapping) +
  theme(legend.position = "right")
p3

# Same as above but faceted — one panel per module (free y-axis scale)
p4 <- ggplot(expression_data_long, aes(x = Pseudotime, y = Expression, color = Annotation)) +
  geom_point(alpha = 0.7, size = 1.5) +
  geom_smooth(se = TRUE, method = "loess", color = "black") +
  theme_classic() +
  xlab("Pseudotime") +
  ylab("Module Score") +
  ggtitle("Module Score Along Slingshot Pseudotime") +
  scale_color_manual(values = color_mapping) +
  theme(legend.position = "right") +
  facet_wrap(~Module, scales = "free_y")
p4

# Export combined plots
ggsave("plots/Cuves_points_coloured_merge.pdf",     plot = p3, width = 6,  height = 4,  dpi = 300)
ggsave("plots/Cuves_points_coloured_split.pdf",     plot = p4, width = 18, height = 10, dpi = 300)

output_dir <- file.path(getwd(), "plots")

# Save one plot per module with annotation colouring + loess line
for (mod in unique(expression_data_long$Module)) {
  p <- ggplot(expression_data_long %>% filter(Module == mod),
              aes(x = Pseudotime, y = Expression, color = Annotation)) +
    geom_point(alpha = 0.7, size = 1.5) +
    geom_smooth(se = TRUE, method = "loess", color = "black") +
    theme_classic() +
    xlab("Pseudotime") +
    ylab("Module Score") +
    ggtitle(paste(mod, "Module Score Along Slingshot Pseudotime")) +
    scale_color_manual(values = color_mapping) +
    theme(legend.position = "right")

  filename_pdf <- file.path(output_dir, paste0(mod, "_Module_Pseudotime.pdf"))
  ggsave(filename_pdf, plot = p, width = 8, height = 6)
  message("Saved: ", filename_pdf)
}

# Clean version per module: no legend, no axes (for figure insets)
for (mod in unique(expression_data_long$Module)) {
  p <- ggplot(expression_data_long %>% filter(Module == mod),
              aes(x = Pseudotime, y = Expression, color = Annotation)) +
    geom_point(alpha = 0.7, size = 1.5) +
    geom_smooth(se = TRUE, method = "loess", color = "black") +
    theme_classic() +
    xlab("Pseudotime") +
    ylab("Module Score") +
    scale_color_manual(values = color_mapping) +
    theme(legend.position = "none") & NoAxes()

  filename_pdf <- file.path(output_dir, paste0(mod, "_Module_Pseudotime_clean.pdf"))
  ggsave(filename_pdf, plot = p, width = 8, height = 6)
  message("Saved: ", filename_pdf)
}


# ==============================================================================
# Figure 4d — tradeSeq: trajectory-associated transcription factor expression
# ==============================================================================
# tradeSeq fits a negative-binomial GAM per gene with spline terms for
# pseudotime, one per lineage, allowing gene expression trends to differ
# between differentiation paths.
# Reference: Van den Berge et al., Nature Communications 2020.

# Identify top 100 cluster-specific marker genes (candidates for GAM fitting)
cluister_markers <- FindAllMarkers(CD8_T2D_ms_version, logfc.threshold = 0.25, only.pos = T)
top100 <- cluister_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 100) %>%
  ungroup() %>%
  distinct(gene, .keep_all = TRUE)
top100_genes <- top100$gene

# Extract raw count matrix (required by tradeSeq; GAM is fitted on counts, not normalised data)
counts <- assay(CD8_sce, "counts")
head(rownames(counts))
head(top100_genes)
test <- intersect(top100_genes, rownames(counts))
counts <- counts[rownames(counts) %in% top100_genes, ]
missing_genes <- setdiff(top100_genes, rownames(counts))
print(missing_genes)  # inspect genes absent from the count matrix

# Pseudotime and lineage weight matrices for the first 5 lineages
pseudotime   <- slingPseudotime(CD8_sce, na = F)[,c(1,2,3,4,5)]
cellWeights  <- slingCurveWeights(CD8_sce)[,c(1,2,3,4,5)]

# Include patient ID as a random effect to account for inter-donor variability
model  <- data.frame(sample = CD8_T2D_ms_version$Pat_ID)
design <- model.matrix(~sample, model)

# Remove cells with zero total weight (not assigned to any lineage)
cellWeights <- cellWeights[which(rowSums(cellWeights) > 0), ]
pseudotime  <- pseudotime[rownames(cellWeights), ]
counts      <- counts[, c(which(colnames(counts) %in% row.names(cellWeights)))]
design      <- design[which(row.names(design) %in% row.names(cellWeights)), ]

# Fit GAM: 12 knots per lineage spline; donor-level nuisance modelled via U matrix
gam <- fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights,
               U = design, nknots = 12, verbose = T, parallel = F, genes = top100_genes)
saveRDS(gam, file = "~/Desktop/T2D_CD8/figures/Figure_3/data/gam_final.rds")
gam <- readRDS("~/Desktop/T2D_CD8/figures/Figure_3/data/gam_final.rds")

# Start-vs-end test: identifies genes with maximal change along each lineage
startRes    <- startVsEndTest(gam, lineages = T)
oStart      <- order(startRes$waldStat, decreasing = TRUE)
sigGeneStart <- names(gam)[oStart[3]]

# Association test: identifies genes that vary with pseudotime on lineage 2
# with a minimum log2 fold change threshold of 0.5
association  <- associationTest(gam, lineages = TRUE, l2fc = 0.5, global = T)

# Candidate transcription factors of interest (senescence/differentiation regulators)
genes_oi <- c("TBX21", "ZEB2", "PRDM1", "RORA", "ID2", "HOPX", "BCL6", "RUNX3",
              "EOMES", "IRF4", "TOX", "FOXO1", "ID3", "TCF7", "STAT3", "BCL11B", "GATA3")

# Filter to genes significantly associated with lineage 2 (FDR < 0.05)
at_sig_lin2   <- association[which(p.adjust(association$pvalue_2) <= 0.05),]
genes_to_export <- at_sig_lin2
genes_to_export$gene <- rownames(genes_to_export)
openxlsx::write.xlsx(genes_to_export, asTable = T, "tables/at_sig_lin2.xlsx")

# Identify the overlap between significant genes and the pre-defined TF list
overlapping_genes <- intersect(genes_oi, rownames(at_sig_lin2))
overlapping_genes

# Predict smoothed expression for each overlapping TF along lineage 1
# (nPoints = 100 interpolated pseudotime points)
yhatSmooth <- predictSmooth(gam, gene = overlapping_genes, nPoints = 100, tidy = FALSE)

# Heatmap of predicted expression (z-scored per gene, lineage 1 only)
png("plots/Heatmap.jpg", width = 8, height = 8, units = "cm", res = 300)
lin1 <- pheatmap::pheatmap(
  yhatSmooth[, 1:100],
  scale         = "row",          # z-score each gene for comparability
  cluster_cols  = FALSE,          # preserve pseudotime order on x-axis
  show_rownames = TRUE,
  show_colnames = FALSE,
  color         = inferno(256),
  legend        = FALSE,
  border_color  = NA
)
dev.off()


# ==============================================================================
# Figure 4e — Smoothed TF expression curves along pseudotime (lineage 2)
# ==============================================================================

# Extract log-normalised expression of key TFs (chosen from Figure 4d)
genes_to_plot   <- c("TOX", "EOMES", "ZEB2", "TBX21")
expression_data <- as.data.frame(t(logcounts(CD8_sce[genes_to_plot, ])))

# Add lineage 2 pseudotime and cluster annotation as covariates
expression_data$pseudotime      <- CD8_sce$slingPseudotime_2
expression_data$final_annotation <- CD8_sce$final_annotation

# Convert to long format for ggplot
expression_data_long <- tidyr::gather(expression_data, key = "gene", value = "expression",
                                       -pseudotime, -final_annotation)

# Smoothed expression curves: one line per TF, shaded confidence intervals
plot <- ggplot(expression_data_long, aes(x = pseudotime, y = expression, color = gene, fill = gene)) +
  geom_smooth(se = TRUE, alpha = 0.2) +
  theme_classic() +
  xlab("Pseudotime") +
  ylab("Expression") +
  labs(title = NULL) +
  scale_color_tableau() +  # distinct, publication-ready colour scale
  scale_fill_tableau() +
  theme(legend.position = "none")

ggsave(
  filename = "plots/Gene_Expression_Pseudotime.pdf",
  plot     = plot,
  width = 8, height = 8, units = "cm", dpi = 300
)


# ==============================================================================
# Figure 4f — Nebulosa kernel density plot of ZEB2 expression
# ==============================================================================
# Nebulosa estimates gene expression density on the UMAP using a kernel
# density approach, which reduces the visual impact of dropout noise.
# Reference: Lopez-Delisle & Bhatt, Bioinformatics 2021.

neb_plot <- Nebulosa::plot_density(CD8_T2D_ms_version, features = "ZEB2",
                                    pal = "inferno", size = 1.5) &
  NoLegend() &
  NoAxes() &
  labs(title = NULL)
neb_plot <- neb_plot + coord_fixed(ratio = 1.5)  # consistent aspect ratio with other UMAPs
ggsave("plots/ZEB2_Nebulosa_fixed.pdf", plot = neb_plot, units = "cm", height = 12, width = 12)
