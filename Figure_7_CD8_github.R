# =============================================================================
# Figure 7: ZEB2 Regulates CD8 T Cell Differentiation In Vivo —
#           CyTOF Analysis of Mouse PBMCs at 16 Weeks
# =============================================================================
#
# Description:
#   This script analyses mass cytometry (CyTOF) data from peripheral blood
#   mononuclear cells (PBMCs) of CD8-Cre+ Zeb2^wt/wt (WT) and
#   CD8-Cre+ Zeb2^flox/flox (KO) mice harvested at 16 weeks.
#   The CATALYST package is used throughout for CyTOF-specific preprocessing,
#   clustering, and visualisation. The analysis covers:
#     - FlowSOM clustering of the full PBMC panel (10×10 SOM, 15 metaclusters)
#     - Quality filtering (doublet and low-QC cluster removal)
#     - UMAP dimensionality reduction and broad immune cell type annotation
#     - CD8 T cell sub-clustering into Naive/CM, EM, EMRA, and Cytotoxic EMRA
#     - Proportional analysis of CD8 subsets per sample and per genotype
#     - Figure export: UMAP, proportion bar charts, marker heatmap, bubble plot
#
# Input files required:
#   - FCS files (16W timepoint):
#       /Volumes/research/.../PBMC_enrolled_mice/Cleaned data/*.fcs
#   - Metadata Excel file:
#       /Volumes/research/.../2026.04.10_Metadata_PBMC_CyTOF.xlsx
#       Required columns: file_name, sample_id, genotype, Timepoint, batch, mouse_id
#
# Output files — plots:
#   - ExtFig6F_immune_heatmap.pdf       (Extended Figure 6f: full PBMC marker heatmap)
#   - Fig6h_UMAP_and_proportions.pdf    (Figure 6h: PBMC UMAP + proportion bar)
#   - Fig7a_UMAP_and_CD8_proportions.pdf (Figure 7a: CD8 subset UMAP + proportion bar)
#   - ExtFig7a_CD8_bubble_plot.pdf      (Extended Figure 7a: CD8 subset marker bubble plot)
#   - ExtFig7b_CD8_stacked_bar_WTvsKO.pdf (Extended Figure 7b: per-sample CD8 proportions)
#
# Reporting standards followed:
#   Nature Portfolio reporting guidelines. FlowSOM and UMAP parameters are
#   documented inline. Cofactor = 5 (standard for PBMC CyTOF; lower than
#   cofactor 15 used for tissue samples due to higher signal intensity).
#   All source data are represented in the exported proportion tables.
#
# =============================================================================

# --- Libraries ----------------------------------------------------------------
library(CATALYST)       # CyTOF-specific: prepData, cluster, plotExprHeatmap, plotDR
library(flowCore)       # low-level FCS file I/O
library(SingleCellExperiment) # Bioconductor data container used by CATALYST
library(ggplot2)        # plotting
library(dplyr)          # data manipulation
library(tidyr)          # data reshaping
library(readxl)         # reading Excel metadata files
library(scales)         # axis formatting (percent_format)
library(cowplot)        # plot composition (plot_grid) and publication theme
library(matrixStats)    # fast row-wise operations (rowMedians)
library(RColorBrewer)   # colour palettes

# Define a clean publication theme applied globally across all plots
theme_pub <- theme_cowplot(font_size = 12) +
  theme(strip.background = element_rect(fill = "grey90", color = NA),
        legend.position  = "right")
theme_set(theme_pub)

set.seed(1234) # global random seed for reproducibility

# --- Paths --------------------------------------------------------------------
# Update these paths if running on a different machine or file server

fcs_dir   <- "/Volumes/research/giannc02lab/homes/eberhn01/PROJECTS/ZEB2 MOUSE/CyTOF/PBMC_enrolled_mice/Cleaned data"
meta_file <- "/Volumes/research/giannc02lab/homes/eberhn01/PROJECTS/ZEB2 MOUSE/CyTOF/R analysis/CyTOF analysis/2026.04.10_Metadata_PBMC_CyTOF.xlsx"
fig_dir   <- "/Volumes/research/giannc02lab/homes/eberhn01/PROJECTS/ZEB2 MOUSE/CyTOF/R analysis/16W analysis/figures"

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# --- Colour palettes ----------------------------------------------------------

# One colour per broad immune cell type (7 types → interpolated Set1 palette)
type_colors <- setNames(
  colorRampPalette(brewer.pal(9, "Set1"))(7),
  sort(c("B cells", "CD4 T cells", "CD8 T cells",
         "Classical Monocytes", "DCs", "Mixed Lymphocytes",
         "Non-classical Monocytes"))
)

# Fixed colours for the four CD8 T cell subsets (used consistently across all panels)
cd8_colors <- c(
  "Naive/CM"       = "#2196F3",
  "EM"             = "#4CAF50",
  "EMRA"           = "#9C27B0",
  "Cytotoxic EMRA" = "#F44336"
)


# =============================================================================
# STEP 1 — Load and filter metadata
# =============================================================================
# Only the 16-week timepoint is used here; baseline samples are excluded.
# Genotype is recoded to a two-level factor: WT vs KO.

md <- read_xlsx(meta_file) %>%
  mutate(
    condition = ifelse(grepl("flox/flox", genotype), "KO", "WT"),  # KO = Zeb2^flox/flox
    condition = factor(condition, levels = c("WT", "KO")),
    Timepoint = factor(Timepoint, levels = c("baseline", "16W")),
    batch     = factor(batch),
    mouse_id  = as.character(mouse_id)
  ) %>%
  dplyr::filter(Timepoint == "16W")  # restrict to 16-week endpoint

# Match each metadata row to its FCS file by filename
fcs_files <- list.files(fcs_dir, pattern = "\\.fcs$", full.names = TRUE)
md        <- md %>% mutate(fcs_path = setNames(fcs_files, basename(fcs_files))[file_name])

# Halt with an informative error if any expected FCS files are missing
missing <- md$file_name[is.na(md$fcs_path)]
if (length(missing) > 0) stop("FCS files not found:\n", paste(missing, collapse = "\n"))

cat("Samples loaded:", nrow(md), "\n")
print(md %>% count(condition))  # confirm WT/KO sample counts


# =============================================================================
# STEP 2 — Define CyTOF panel and build SingleCellExperiment
#          Arcsinh transformation with cofactor = 5
# =============================================================================
# Markers are split into two classes:
#   "type"  — lineage markers used for cell type clustering (26 markers)
#   "state" — functional/activation markers not used for clustering (8 markers)
# Cofactor = 5 is standard for PBMC CyTOF (higher signal intensity than tissue).

panel <- data.frame(
  fcs_colname = c(
    "Y89Di",    "Cd114Di",  "Cd116Di",  "Nd143Di",  "Nd146Di",  "Nd148Di",
    "Nd150Di",  "Eu151Di",  "Sm149Di",  "Sm152Di",  "Gd155Di",  "Gd158Di",
    "Dy161Di",  "Dy163Di",  "Dy164Di",  "Er168Di",  "Er170Di",  "Tm169Di",
    "Lu175Di",  "Yb172Di",  "Yb173Di",  "Pt194Di",  "Pt195Di",  "Pt196Di",
    "Pt198Di",  "Bi209Di",
    "Pr141Di",  "Er166Di",  "Er167Di",  "Ho165Di",  "Sm154Di",
    "Tb159Di",  "Yb174Di",  "Yb176Di"
  ),
  antigen = c(
    "CD45",   "B220",    "CD4",    "TCRb",    "CD43",   "CD11b",
    "CD27",   "IgM",     "CD19",   "CD49b",   "CX3CR1", "CD80",
    "CD11c",  "Siglec-F","PD-L1",  "CD138",   "NK1.1",  "CD21",
    "XCR1",   "TCRgd",   "F4/80",  "CD8a",    "Ly-6G",  "CD3",
    "Ly-6C",  "MHC-II",
    # State markers (functional, not used for type clustering):
    "CD44",   "CD62L",   "GZMB",   "CD127",   "CTLA-4",
    "PD-1",   "CD25",    "CCR7"
  ),
  marker_class = c(rep("type", 26), rep("state", 8)),
  stringsAsFactors = FALSE
)

# prepData: reads FCS files, applies arcsinh transform, stores in SCE format
sce <- prepData(
  x        = md$fcs_path,
  panel    = panel,
  md       = as.data.frame(md),
  md_cols  = list(file    = "file_name",
                  id      = "sample_id",
                  factors = c("condition", "Timepoint", "batch", "mouse_id")),
  cofactor  = 5,
  FACS      = FALSE,  # CyTOF data (not fluorescence)
  transform = TRUE    # apply arcsinh transformation
)

cat("SCE built:", ncol(sce), "cells,", nrow(sce), "markers\n")


# =============================================================================
# STEP 3 — FlowSOM clustering on type markers
#          10×10 SOM grid → up to 20 metaclusters (meta15 used downstream)
# =============================================================================
# FlowSOM first trains a self-organising map on the type marker space,
# then applies consensus hierarchical clustering to collapse SOM nodes into
# metaclusters. meta15 (15 metaclusters) is selected based on visual inspection
# of the marker heatmap.

sce <- cluster(sce, features = "type", xdim = 10, ydim = 10,
               maxK = 20, seed = 1234)


# =============================================================================
# STEP 4 — Remove doublets and low-QC cells
#
# Identified by inspecting the marker heatmap per meta15 cluster:
#   Cluster 12 — doublets: uniformly high across all lineage markers
#   Cluster 7  — low QC:   CD45 low, near-zero expression across all markers
# =============================================================================

sce <- sce[, !cluster_ids(sce, "meta15") %in% c("12", "7")]
cat("After QC removal:", ncol(sce), "cells\n")


# =============================================================================
# STEP 5 — UMAP dimensionality reduction on the full PBMC dataset
#          Downsampled to 2000 cells per sample for speed; type markers used
# =============================================================================

sce <- runDR(sce, dr = "UMAP", cells = 2000, features = "type",
             n_neighbors = 50, min_dist = 0.3, spread = 3, seed = 1234)


# =============================================================================
# STEP 6 — Broad immune cell type annotation
#
# Assignments based on canonical marker expression per meta15 cluster:
#   B cells:                 B220+  CD19+  IgM+
#   CD4 T cells:             CD3+   CD4+   TCRb+
#   CD8 T cells:             CD3+   CD8a+  TCRb+
#   Classical Monocytes:     CD11b+ Ly-6C hi
#   Non-classical Monocytes: CD11b+ CX3CR1 hi, Ly-6C lo
#   DCs:                     CD11c+ MHC-II+
#   Mixed Lymphocytes:       TCRb+  B220+  MHC-II hi
# =============================================================================

sce <- mergeClusters(
  sce,
  k     = "meta15",
  table = data.frame(
    old_cluster = c("1","2","3","4","5","6","8","9","10","11","13","14","15"),
    new_cluster = c(
      "CD4 T cells",             # 1
      "CD8 T cells",             # 2
      "B cells",                 # 3
      "B cells",                 # 4
      "B cells",                 # 5
      "Mixed Lymphocytes",       # 6
      "B cells",                 # 8
      "Non-classical Monocytes", # 9
      "CD8 T cells",             # 10
      "Classical Monocytes",     # 11
      "DCs",                     # 13
      "Classical Monocytes",     # 14
      "Non-classical Monocytes"  # 15
    ),
    stringsAsFactors = FALSE
  ),
  id = "cell_type", overwrite = TRUE
)

cat("Cell type counts:\n")
print(table(cluster_ids(sce, "cell_type")))


# =============================================================================
# STEP 7 — CD8 T cell sub-clustering
#          Subset to CD8 T cells → recluster on differentiation/functional markers
# =============================================================================
# A smaller 8×8 SOM is used for the CD8 subset (fewer cells, fewer markers).
# State markers related to memory/effector differentiation are used here instead
# of the broad lineage markers.

sce_cd8 <- sce[, cluster_ids(sce, "cell_type") == "CD8 T cells"]
cat("CD8 T cells extracted:", ncol(sce_cd8), "\n")

# CD8 differentiation and functional markers
cd8_features <- c("CD44","CD62L","CCR7","CD27","CX3CR1",
                  "GZMB","PD-1","CTLA-4","CD127","CD25","CD43")
cd8_features <- cd8_features[cd8_features %in% rownames(sce_cd8)]  # safety: keep present markers

sce_cd8 <- cluster(sce_cd8, features = cd8_features,
                   xdim = 8, ydim = 8, maxK = 12, seed = 1234)

# Initial UMAP for QC inspection before Low QC cluster removal
sce_cd8 <- runDR(sce_cd8, dr = "UMAP", cells = NULL,
                 features = cd8_features,
                 n_neighbors = 30, min_dist = 0.1, seed = 1234)


# =============================================================================
# STEP 8 — CD8 subset annotation and Low QC removal
#
# meta7 cluster assignments (from inspection of CD8 marker heatmap):
#   1 → CM (central memory: CD44hi CD62Lhi CCR7hi CD27hi)
#   2 → EM (effector memory: CD44hi CD62Llo CX3CR1+)
#   3 → EM
#   4 → EMRA (terminally differentiated: CD44hi CD62Llo CX3CR1hi)
#   5 → Low QC (near-zero expression; removed)
#   6 → Cytotoxic EMRA (GZMBhi)
#   7 → Cytotoxic EMRA
# CM is relabelled Naive/CM to reflect the mixed naive/central-memory phenotype
# in mouse peripheral blood.
# =============================================================================

sce_cd8 <- mergeClusters(
  sce_cd8,
  k     = "meta7",
  table = data.frame(
    old_cluster = c("1","2","3","4","5","6","7"),
    new_cluster = c("CM","EM","EM","EMRA","Low QC","Cytotoxic EMRA","Cytotoxic EMRA"),
    stringsAsFactors = FALSE
  ),
  id = "cd8_subset", overwrite = TRUE
)

# Remove low-quality cells
sce_cd8 <- sce_cd8[, cluster_ids(sce_cd8, "cd8_subset") != "Low QC"]

# Rename CM → Naive/CM in the final label set
sce_cd8 <- mergeClusters(
  sce_cd8,
  k     = "cd8_subset",
  table = data.frame(
    old_cluster = c("CM","EM","EMRA","Cytotoxic EMRA"),
    new_cluster = c("Naive/CM","EM","EMRA","Cytotoxic EMRA"),
    stringsAsFactors = FALSE
  ),
  id = "cd8_subset", overwrite = TRUE
)

# Re-run UMAP on the final clean CD8 cell set
sce_cd8 <- runDR(sce_cd8, dr = "UMAP", cells = NULL,
                 features = cd8_features,
                 n_neighbors = 30, min_dist = 0.1, seed = 1234)

cat("Final CD8 subset counts:\n")
print(table(cluster_ids(sce_cd8, "cd8_subset")))


# =============================================================================
# STEP 9 — Compute cell type and CD8 subset proportions
# =============================================================================

# --- Overall PBMC proportions (pooled across all samples) --------------------
df_cells  <- data.frame(cell_type = as.character(cluster_ids(sce, "cell_type")),
                         stringsAsFactors = FALSE)
ct_order   <- df_cells %>% count(cell_type) %>% arrange(desc(n)) %>% pull(cell_type)
prop_total <- df_cells %>%
  count(cell_type) %>%
  mutate(prop      = n / sum(n),
         cell_type = factor(cell_type, levels = rev(ct_order)),
         x         = "All")  # dummy x-axis variable for the stacked bar

# --- CD8 subset proportions per sample and per condition ---------------------
cd8_prop_df <- data.frame(
  sample_id  = sce_cd8$sample_id,
  cd8_subset = cluster_ids(sce_cd8, "cd8_subset"),
  condition  = sce_cd8$condition,
  stringsAsFactors = FALSE
) %>%
  count(sample_id, condition, cd8_subset) %>%
  group_by(sample_id) %>%
  mutate(prop = n / sum(n)) %>%  # normalise to total CD8 T cells per sample
  ungroup() %>%
  mutate(cd8_subset = factor(cd8_subset, levels = names(cd8_colors)),
         condition  = factor(condition,  levels = c("WT", "KO")))

# Mean proportion per condition (used in the stacked bar in Figure 7a)
cd8_prop_mean <- cd8_prop_df %>%
  group_by(condition, cd8_subset) %>%
  summarise(prop = mean(prop), .groups = "drop") %>%
  mutate(cd8_subset = factor(cd8_subset, levels = names(cd8_colors)),
         condition  = factor(condition,  levels = c("WT", "KO")))


# =============================================================================
# FIGURES
# =============================================================================

# ── Extended Figure 6f — Full PBMC immune landscape heatmap -------------------
# Scaled median arcsinh expression per cell type; rows = markers, columns = types.
# bar = TRUE shows cell counts per type; perc = TRUE shows relative proportions.

all_markers_ordered <- c(
  "CD45","B220","CD4","TCRb","CD43","CD11b","CD27",
  "IgM","CD19","CD49b","CX3CR1","CD80","CD11c","Siglec-F",
  "PD-L1","CD138","NK1.1","CD21","XCR1","TCRgd","F4/80",
  "CD8a","Ly-6G","CD3","Ly-6C","MHC-II",
  "CD44","CD62L","GZMB","CD127","CTLA-4","PD-1","CD25","CCR7"
)
all_markers_ordered <- all_markers_ordered[all_markers_ordered %in% rownames(sce)]  # safety check

pdf(file.path(fig_dir, "ExtFig6F_immune_heatmap.pdf"), width = 14, height = 9)
print(plotExprHeatmap(sce, features = all_markers_ordered,
                      by = "cluster_id", k = "cell_type",
                      scale = "last",       # scale across clusters (last = z-score columns)
                      bars = TRUE, perc = TRUE,
                      row_clust = FALSE,    # preserve manual marker order
                      col_clust = TRUE))    # cluster cell types by expression similarity
dev.off()


# ── Figure 6h — PBMC UMAP + overall proportion bar chart ----------------------
# Left panel: UMAP coloured by cell type.
# Right panel: single stacked bar showing proportions of all immune cell types.

p_umap_6h <- plotDR(sce, dr = "UMAP", color_by = "cell_type") +
  scale_color_manual(values = type_colors, name = NULL) +
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 1)) +
  theme_void() +
  theme(legend.position = "right", legend.text = element_text(size = 9))

p_bar_6h <- ggplot(prop_total, aes(x = x, y = prop, fill = cell_type)) +
  geom_col(width = 0.6, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = type_colors, guide = "none") +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0),
                     name = "Proportion (% of Live cells)") +
  labs(x = NULL) + theme_pub +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggsave(file.path(fig_dir, "Fig6h_UMAP_and_proportions.pdf"),
       plot_grid(p_umap_6h, p_bar_6h, nrow = 1, rel_widths = c(3, 1),
                 align = "h", axis = "tb"),
       width = 11, height = 6)


# ── Figure 7a — CD8 T cell subset UMAP + mean proportion bar (WT vs KO) ------
# Left panel: UMAP of CD8 T cells coloured by subset.
# Right panel: mean-proportion stacked bar per genotype.
# Genotype labels use bquote() for superscript formatting.

p_umap_7a <- plotDR(sce_cd8, dr = "UMAP", color_by = "cd8_subset") +
  scale_color_manual(values = cd8_colors, name = "Subset") +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme_void() +
  theme(legend.position = "right", legend.text = element_text(size = 9))

p_bar_7a <- ggplot(cd8_prop_mean, aes(x = condition, y = prop, fill = cd8_subset)) +
  geom_col(width = 0.6, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = cd8_colors, guide = "none") +
  scale_x_discrete(labels = c(
    "WT" = bquote("CD8-Cre"^"+" ~ "Zeb2"^"wt/wt"),
    "KO" = bquote("CD8-Cre"^"+" ~ "Zeb2"^"flox/flox")
  )) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  labs(x = NULL, y = "Proportion (% of CD8 T cells)") + theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

ggsave(file.path(fig_dir, "Fig7a_UMAP_and_CD8_proportions.pdf"),
       plot_grid(p_umap_7a, p_bar_7a, nrow = 1, rel_widths = c(3, 1),
                 align = "h", axis = "tb"),
       width = 11, height = 6)


# ── Extended Figure 7a — CD8 subset marker bubble plot ------------------------
# Dot size   = % cells with arcsinh expression > 0.5 (threshold for "positive")
# Dot colour = min-max scaled median expression (comparable across markers)

EXPR_THRESHOLD <- 0.5  # arcsinh expression cutoff for "expressing" classification
expr_mat <- assay(sce_cd8, "exprs")
subsets  <- cluster_ids(sce_cd8, "cd8_subset")

# Compute per-subset median expression and % positive for each CD8 feature
bubble_dat <- do.call(rbind, lapply(levels(subsets), function(ss) {
  sub_mat <- expr_mat[cd8_features, subsets == ss, drop = FALSE]
  data.frame(subset   = ss,
             marker   = cd8_features,
             med_expr = rowMedians(sub_mat),
             pct_expr = rowMeans(sub_mat > EXPR_THRESHOLD) * 100,
             stringsAsFactors = FALSE)
})) %>%
  group_by(marker) %>%
  # Min-max scale per marker so colours are comparable across markers
  mutate(scaled_med = (med_expr - min(med_expr)) /
                      (max(med_expr) - min(med_expr) + 1e-9)) %>%  # 1e-9 avoids division by zero
  ungroup() %>%
  mutate(subset = factor(subset, levels = c("EMRA","Cytotoxic EMRA","EM","Naive/CM")),
         marker = factor(marker, levels = cd8_features))

ggsave(file.path(fig_dir, "ExtFig7a_CD8_bubble_plot.pdf"),
       ggplot(bubble_dat, aes(x = marker, y = subset,
                              size = pct_expr, color = scaled_med)) +
         geom_point() +
         scale_size_continuous(name = "% Expressing", range = c(1, 10),
                               breaks = c(25, 50, 75)) +
         scale_color_gradientn(name = "Scaled Median\nExpression",
                               colors = c("#2166AC","white","#D6604D"),
                               limits = c(0, 1)) +
         labs(x = NULL, y = NULL) + theme_pub +
         theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
               axis.text.y = element_text(size = 10)),
       width = 8, height = 4)


# ── Extended Figure 7b — Per-sample CD8 subset proportions (WT vs KO) --------
# Individual stacked bar per mouse, faceted by genotype.
# Shows within-genotype variability in CD8 subset composition.

ggsave(file.path(fig_dir, "ExtFig7b_CD8_stacked_bar_WTvsKO.pdf"),
       ggplot(cd8_prop_df, aes(x = sample_id, y = prop, fill = cd8_subset)) +
         geom_col(width = 0.8, color = "white", linewidth = 0.2) +
         facet_wrap(~ condition, scales = "free_x",
                    labeller = labeller(condition = c(
                      "WT" = "CD8-Cre+ Zeb2^WT/WT",
                      "KO" = "CD8-Cre+ Zeb2^flox/flox"
                    ))) +
         scale_fill_manual(values = cd8_colors, name = "Subset") +
         scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
         labs(x = NULL, y = "Proportion (% of CD8 T cells)") + theme_pub +
         theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 9),
               strip.text   = element_text(face = "bold", size = 10)),
       width = 10, height = 5)


cat("\nAll figures saved to:", fig_dir, "\n")
