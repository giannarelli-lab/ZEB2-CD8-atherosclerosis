# =============================================================================
# Figure 3: Differential Abundance and Module Associations in TEMRA CD8 T Cells
# =============================================================================
#
# Description:
#   This script generates panels for Figure 3, which examines how CD8 T cell
#   subset composition and functional module scores associate with clinical
#   variables. It covers:
#     - Export of mean UCell module scores per patient (all cells and TEMRA only)
#     - Radar charts comparing module correlations with clinical metadata
#       (age, BMI, pathology, diabetes diagnosis, etc.)
#     - Differential abundance testing with MiloR:
#         * Plaque pathology (Fibroatheroma vs Fibrocalcific)
#         * Diabetes status (Type 2 vs No Diabetes)
#     - Export of cluster abundance tables stratified by pathology and condition
#
# Input files required:
#   - ~/Desktop/T2D_CD8/data/CD8_T2D_ms_version_final_modules.rds
#       Seurat object: CD8 T cells with final cluster annotations and UCell scores
#   - ~/Desktop/T2D_CD8/data/CD8_T2D_ms_version.rds
#       Seurat object: CD8 T cells (base version, used for T2D MiloR analysis)
#   - ~/Desktop/T2D_CD8/figures/Figure_2/tables/TEMRA_modules_radarchart.xlsx
#       Pre-computed Spearman correlation coefficients for TEMRA module radar chart
#
# Output files — plots:
#   - figures/Figure_2/TEMRA_radar_Functional.pdf     (Figure 3a)
#   - figures/Figure_2/TEMRA_radar_Dysfunctional.pdf  (Figure 3b)
#   - figures/Figure_2/plots/beeswarm_pathology.pdf   (Figure 3d)
#   - figures/Figure_2/plots/nh_graph_pathology.pdf   (Figure 3d)
#   - plots/beeswarm_T2D.pdf                          (Figure 3e)
#   - plots/nh_graph_T2D.pdf                          (Figure 3e)
#
# Output files — source data / tables:
#   - figures/Figure_2/tables/patient_ucell_means.xlsx         (Figures 3a/b, all cells)
#   - figures/Figure_2/tables/patient_ucell_means_TEMRA.xlsx   (Figures 3a/b, TEMRA)
#   - figures/Figure_2/tables/MiloR_pathology/da_results.csv   (Figure 3d)
#   - tables/MiloR_T2D/da_results.csv                         (Figure 3e)
#   - figures/Figure_2/tables/cluster_abundance_with_pathology.xlsx (Figures 3f/g)
#
# =============================================================================

# --- Libraries ----------------------------------------------------------------
library(Seurat)           # scRNA-seq object handling
library(SeuratDisk)       # object serialisation
library(SeuratWrappers)   # Seurat utility wrappers
library(ggplot2)          # plotting
library(patchwork)        # combining ggplots
library(dplyr)            # data manipulation
library(data.table)       # fast tabular operations
library(readxl)           # reading Excel files
library(ggalluvial)       # alluvial / Sankey plots
library(miloR)            # differential abundance testing on KNN graphs
library(SummarizedExperiment) # Bioconductor data container (used by miloR)

set.seed(123) # global random seed for reproducibility


# ==============================================================================
# Figures 3a & 3b — Mean UCell module scores per patient (source data)
# ==============================================================================

# Load the CD8 Seurat object that includes final cluster annotations and UCell scores
CD8_T2D_ms_version_final_modules <- readRDS("~/Desktop/T2D_CD8/data/CD8_T2D_ms_version_final_modules.rds")
head(CD8_T2D_ms_version_final_modules@meta.data)
unique(CD8_T2D_ms_version_final_modules$Pat_ID)  # inspect donor IDs

# Define the UCell module columns to summarise
ucell_modules <- c(
  "SENMAYO_UCell_kNN",
  "EXHAUSTION_UCell_kNN",
  "T_CELL_SENESCENCE_UCell_kNN",
  "CYTOTOXIC_UCell_kNN",
  "SENMAYO_CYTOTOXIC_UCell_kNN",      # composite: SENMAYO × CYTOTOXIC
  "TSEN_CYTOTOXIC_UCell_kNN",         # composite: T_CELL_SENESCENCE × CYTOTOXIC
  "SASP_UCell_kNN",
  "TERM_MEMORY_UCell_kNN",
  "PROLIFERATION_UCell_kNN",
  "NAIVE_UCell_kNN",
  "MEMORY_UCell_kNN",
  "RESIDENT_UCell_kNN",
  "signature_1SENEPY_SCORE_UCELL_kNN" # combined senescence-cytotoxicity signature
)

# Extract metadata for all cells and compute per-patient means across all CD8 clusters
md <- CD8_T2D_ms_version_final_modules@meta.data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell_id")

patient_ucell_means <- md %>%
  group_by(Pat_ID) %>%
  summarise(
    across(all_of(ucell_modules), ~ mean(.x, na.rm = TRUE)),
    n_cells = n()
  )

openxlsx::write.xlsx(patient_ucell_means, asTable = T,
                     "~/Desktop/T2D_CD8/figures/Figure_2/tables/patient_ucell_means.xlsx")

# --- TEMRA subset: compute per-patient module means for TEMRA clusters only ---
# TEMRA = terminally differentiated effector memory RA cells (three clusters)
Idents(CD8_T2D_ms_version_final_modules) <- "final_annotation"
levels(CD8_T2D_ms_version_final_modules)

TEMRA <- subset(CD8_T2D_ms_version_final_modules,
                idents = c("C2- CX3CR1+ ADGRG1+ EMRA", "C8- HAVCR+ HOPX+ EMRA", "C10 - FCGR3B+ B3GAT1+ EMRA"))

# Reuse the same module list for TEMRA
ucell_modules <- c(
  "SENMAYO_UCell_kNN",
  "EXHAUSTION_UCell_kNN",
  "T_CELL_SENESCENCE_UCell_kNN",
  "CYTOTOXIC_UCell_kNN",
  "SENMAYO_CYTOTOXIC_UCell_kNN",
  "TSEN_CYTOTOXIC_UCell_kNN",
  "SASP_UCell_kNN",
  "TERM_MEMORY_UCell_kNN",
  "PROLIFERATION_UCell_kNN",
  "NAIVE_UCell_kNN",
  "MEMORY_UCell_kNN",
  "RESIDENT_UCell_kNN",
  "signature_1SENEPY_SCORE_UCELL_kNN"
)

md <- TEMRA@meta.data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell_id")

patient_ucell_means_TEMRA <- md %>%
  group_by(Pat_ID) %>%
  summarise(
    across(all_of(ucell_modules), ~ mean(.x, na.rm = TRUE)),
    n_cells = n()
  )

openxlsx::write.xlsx(patient_ucell_means_TEMRA, asTable = T,
                     "~/Desktop/T2D_CD8/figures/Figure_2/tables/patient_ucell_means_TEMRA.xlsx")

# --- Radar charts for TEMRA module associations with clinical variables --------
library(readxl)
library(dplyr)
library(tibble)
library(fmsb)  # radar / spider chart plotting

# 1) Read pre-computed correlation table (module vs clinical variable)
corr_raw <- read_excel("~/Desktop/T2D_CD8/figures/Figure_2/tables/TEMRA_modules_radarchart.xlsx")

# Identify the baseline (row label) column — may be named "Baseline" or "...1"
baseline_col <- if ("Baseline" %in% names(corr_raw)) {
  "Baseline"
} else if ("...1" %in% names(corr_raw)) {
  "...1"
} else {
  names(corr_raw)[1]
}

# 2) Convert to numeric matrix with clinical variables as row names
corr_df <- corr_raw %>%
  rename(Baseline = all_of(baseline_col)) %>%
  mutate(Baseline = trimws(as.character(Baseline))) %>%
  mutate(across(-Baseline, ~ as.numeric(gsub(",", ".", as.character(.x))))) %>%
  column_to_rownames("Baseline")

stopifnot(nrow(corr_df) >= 3)  # radar chart requires at least 3 axes
stopifnot(ncol(corr_df) >= 3)

# 3) Helper: resolve module column names with partial prefix matching
#    (handles minor naming differences between Excel and Seurat metadata)
resolve_cols <- function(target_names, available_names, prefix_len = 25) {
  resolved <- vapply(target_names, function(tn) {
    if (tn %in% available_names) return(tn)
    pref <- substr(tn, 1, min(prefix_len, nchar(tn)))
    hits <- which(startsWith(available_names, pref))
    if (length(hits) == 1) return(available_names[hits])
    NA_character_
  }, character(1))

  missing <- target_names[is.na(resolved)]
  if (length(missing) > 0) {
    message("These requested module columns were NOT found (even with prefix matching):\n  ",
            paste(missing, collapse = ", "))
    message("\nAvailable columns are:\n  ",
            paste(available_names, collapse = ", "))
    stop("Fix column names or adjust module lists.")
  }
  unname(resolved)
}

# 4) Define the two module groups to show in separate radar charts
#    Functional: memory, proliferation, and residency signatures
functional_modules_wanted <- c(
  "NAIVE_UCell_kNN",
  "PROLIFERATION_UCell_kNN",
  "MEMORY_UCell_kNN",
  "TERM_MEMORY_UCell_kNN",
  "RESIDENT_UCell_kNN"
)

#    Dysfunctional: senescence, exhaustion, cytotoxicity, and SASP signatures
dysfunctional_modules_wanted <- c(
  "EXHAUSTION_UCell_kNN",
  "SENMAYO_UCell_kNN",
  "CYTOTOXIC_UCell_kNN",
  "SASP_UCell_kNN",
  "T_CELL_SENESCENCE_UCell_kNN",
  "signature_1SENEPY_SCORE_UCELL_kNN"
)

functional_modules    <- resolve_cols(functional_modules_wanted, colnames(corr_df))
dysfunctional_modules <- resolve_cols(dysfunctional_modules_wanted, colnames(corr_df))

# 5) Radar plot function
#    Each line = one module; axes = clinical variables; values = Spearman rho
plot_module_radar <- function(corr_df, modules, title,
                              axis_min = -1, axis_max = 1,
                              baseline_order = NULL,
                              legend_labels  = NULL,
                              cols_border,
                              alpha_fill = 0.18) {

  stopifnot(all(modules %in% colnames(corr_df)))

  df <- corr_df[, modules, drop = FALSE]

  # Optionally reorder axes (clinical variables) to a specified sequence
  if (!is.null(baseline_order)) {
    rn <- gsub(" ", " ", rownames(df))  # replace non-breaking spaces
    rn <- trimws(rn)
    rownames(df) <- rn

    bo <- gsub(" ", " ", baseline_order)
    bo <- trimws(bo)

    keep <- bo[bo %in% rownames(df)]
    df   <- df[keep, , drop = FALSE]
  }

  if (nrow(df) < 3)
    stop("Radar plot needs >= 3 baseline axes; only found: ", nrow(df))

  # fmsb expects rows: max, min, then data rows
  df_t <- as.data.frame(t(df))

  radar_df <- rbind(
    rep(axis_max, ncol(df_t)),   # row 1: axis maximum
    rep(axis_min, ncol(df_t)),   # row 2: axis minimum
    df_t                         # rows 3+: one per module
  )
  colnames(radar_df) <- rownames(df)

  n <- nrow(df_t)

  cols_border <- cols_border[seq_len(n)]
  cols_fill   <- grDevices::adjustcolor(cols_border, alpha.f = alpha_fill)

  if (is.null(legend_labels)) {
    legend_labels <- rownames(df_t)
    legend_labels <- gsub("_UCell_kNN$", "", legend_labels)  # strip suffix for cleaner legend
    legend_labels <- gsub("_", " ", legend_labels)
  }

  op <- par(mar = c(1.2, 1.2, 3, 12), xpd = NA)  # wide right margin for legend
  on.exit(par(op), add = TRUE)

  fmsb::radarchart(
    radar_df,
    axistype     = 1,
    title        = title,
    seg          = 4,         # 4 grid segments between min and max
    caxislabels  = c(axis_min, axis_min/2, 0, axis_max/2, axis_max),
    pcol         = cols_border,
    pfcol        = NA,        # no fill (lines only)
    plwd         = 2,
    plty         = 1,
    cglcol       = "black",
    cglty        = 1,
    axislabcol   = NA,
    vlcex        = 0.85
  )

  legend(
    x        = "right",
    inset    = c(-0.42, 0),
    legend   = legend_labels,
    bty      = "n",
    pch      = 15,
    pt.cex   = 1.2,
    pt.bg    = cols_fill,
    col      = cols_border,
    cex      = 0.8
  )
}

# 6) Define module line colours (one colour per module, consistent across charts)
cols_dysfunctional <- c(
  "#31A354",  # EXHAUSTION — green
  "#8E0152",  # SENMAYO — magenta-red
  "#CC4C02",  # CYTOTOXIC — burnt orange
  "#B2182B",  # SASP — strong red
  "#67001F",  # T CELL SENESCENCE — dark crimson
  "black"     # SENEPY composite
)

cols_functional <- c(
  "#08519C",  # NAIVE — deep blue
  "#3182BD",  # PROLIFERATION — medium blue
  "#5E3C99",  # MEMORY — purple
  "#006D2C",  # TERM_MEMORY — dark green
  "#41B6C4"   # RESIDENT — cyan
)

# 7) Clinical variable axis order (controls which metadata appear on the radar axes)
baseline_order <- c(
  "Age (years)",
  "Gender",
  "BMI classification",
  "Tobacco Use",
  "Hypertension",
  "High Lipids",
  "Carotid Stenosis",
  "Pathology",
  "Diagnosis group for CD8 project"
)

# 8) Export radar charts
pdf(file = "~/Desktop/T2D_CD8/figures/Figure_2/TEMRA_radar_Functional.pdf",
    width = 7.5, height = 7.5, useDingbats = FALSE)
plot_module_radar(
  corr_df,
  functional_modules,
  title          = NULL,
  baseline_order = baseline_order,
  cols_border    = cols_functional,
  alpha_fill     = 0.20
)
dev.off()

pdf(file = "~/Desktop/T2D_CD8/figures/Figure_2/TEMRA_radar_Dysfunctional.pdf",
    width = 7.5, height = 7.5, useDingbats = FALSE)
plot_module_radar(
  corr_df,
  dysfunctional_modules,
  title          = NULL,
  baseline_order = baseline_order,
  cols_border    = cols_dysfunctional,
  alpha_fill     = 0.18,
)
dev.off()


# ==============================================================================
# Figures 3d & 3e — MiloR differential abundance testing
# ==============================================================================
# MiloR identifies cell neighbourhoods (KNN graph-based) that are differentially
# abundant between biological conditions, without requiring hard cluster boundaries.
# Reference: Dann et al., Nature Biotechnology 2022.

library(scater)  # provides as.SingleCellExperiment and plotReducedDim

# --- MiloR: plaque pathology (Fibroatheroma vs Fibrocalcific) -----------------

Idents(CD8_T2D_ms_version_final_modules) <- 'pathology'
levels(CD8_T2D_ms_version_final_modules)

# Subset to the two pathology groups and convert to SingleCellExperiment
CD8_sce <- as.SingleCellExperiment(
  subset(CD8_T2D_ms_version_final_modules, idents = c('Fibroatheroma', 'Fibrocalcific')),
  assay = "integrated"
)
milo <- Milo(CD8_sce)

# Build KNN graph (k = 20 neighbours, d = 10 PCA dimensions)
milo <- buildGraph(milo, k = 20, d = 10, reduced.dim = "PCA")
# Define neighbourhoods: one per cell sampled at prop = 1 (all cells as index cells)
milo <- makeNhoods(milo, prop = 1, k = 20, d=10, refined = TRUE, reduced_dims = "PCA")
plotNhoodSizeHist(milo)  # QC: distribution of neighbourhood sizes
# Count cells per sample in each neighbourhood
milo <- countCells(milo, meta.data = as.data.frame(colData(milo)), sample="Pat_ID")
head(nhoodCounts(milo))

unique(CD8_T2D_ms_version_final_modules$Pat_ID)
unique((CD8_T2D_ms_version_final_modules$pathology))

# Build experimental design data frame (one row per sample)
milo_design <- data.frame(colData(milo))[,c("Pat_ID", "pathology")]
milo_design$Pat_ID   <- as.factor(milo_design$Pat_ID)
milo_design$pathology <- as.factor(milo_design$pathology)
milo_design <- distinct(milo_design)
rownames(milo_design) <- milo_design$Pat_ID
milo_design

# Compute neighbourhood distance (required for testNhoods)
milo <- calcNhoodDistance(milo, d=10, reduced.dim = "PCA")
# Test for differential abundance using a negative binomial GLM
da_results <- testNhoods(milo, design = ~ pathology, design.df = milo_design)
head(da_results)
ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)  # QC: p-value distribution

milo <- buildNhoodGraph(milo)  # build graph for plotting
umap_pl    <- plotReducedDim(milo, dimred = "UMAP", colour_by="pathology", text_by = "celltype", text_size = 3, point_size=0.5) +guides(fill="none")
nh_graph_pl <- plotNhoodGraphDA(milo, da_results, layout="UMAP", alpha=0.1)
nh_graph_pl$layers[[1]] <- NULL  # remove default first layer (replaced below with custom points)
nh_graph_pl

# Annotate neighbourhoods with the most abundant cell type in each neighbourhood
da_results <- annotateNhoods(milo, da_results, coldata_col = "final_annotation")
head(da_results)
ggplot(da_results, aes(final_annotation_fraction)) + geom_histogram(bins=50)  # QC: annotation purity

# Beeswarm plot: logFC per neighbourhood, grouped by annotated cell type
plotDAbeeswarm(da_results, group.by = "final_annotation")
ggsave("~/Desktop/T2D_CD8/figures/Figure_2/plots/beeswarm_pathology.pdf",
       units = "cm", height = 22, width = 26)

# Custom neighbourhood graph plot with manual point styling
nh_graph_pl <- plotNhoodGraphDA(milo, da_results, layout="UMAP", alpha = 0.1)
plot_data       <- nh_graph_pl$data
white_points    <- plot_data
colored_points  <- plot_data[plot_data$colour_by != 0, ]  # highlight DA neighbourhoods
nh_graph_pl <- nh_graph_pl +
  geom_point(data = white_points, aes(x = x, y = y), color = "black", size = 5, shape = 21, stroke = 0.2, fill = "lightgrey", show.legend = FALSE)
nh_graph_pl <- nh_graph_pl +
  geom_point(data = colored_points, aes(x = x, y = y, fill = colour_by), size = 8, shape = 21, stroke = 0.2, show.legend = FALSE)
# Diverging scale: green = enriched in Fibrocalcific, purple = enriched in Fibroatheroma
nh_graph_pl <- nh_graph_pl +
  scale_fill_gradient2(low = "#1b9e77", mid = "white", high = "#7570b3", midpoint = 0)
nh_graph_pl <- nh_graph_pl +
  theme_minimal() +
  theme(legend.position = "right") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
nh_graph_pl$layers[[1]] <- NULL  # remove the default Milo layer (already re-added above)
print(nh_graph_pl)
nh_graph_pl & NoLegend() & NoAxes()
ggsave("~/Desktop/T2D_CD8/figures/Figure_2/plots/nh_graph_pathology.pdf",
       units = "cm", height = 22, width = 26)
nh_graph_pl
ggsave("~/Desktop/T2D_CD8/figures/Figure_2/plots/nh_graph_pathology_w_legend.pdf",
       units = "cm", height = 22, width = 26)

# Export DA results per cell type (one CSV per annotation, sorted by logFC)
head(da_results)
dir.create("~/Desktop/T2D_CD8/Figures/Figure_2/tables/MiloR_pathology")
setwd("~/Desktop/T2D_CD8/figures/Figure_2/tables/MiloR_pathology")
getwd()

sorted_df          <- da_results[order(da_results$logFC), ]
unique_annotations <- unique(sorted_df$final_annotation)
for (final_annotation in unique_annotations) {
  annotation_df <- sorted_df[sorted_df$final_annotation == final_annotation, ]
  file_name     <- paste0(final_annotation, "_sorted_by_FC.csv")
  write.table(annotation_df, file = file_name, sep = ",", row.names = FALSE)
}
write.csv(da_results, "da_results.csv")


# --- MiloR: diabetes status (Type 2 vs No Diabetes) --------------------------

library(miloR)
library(SummarizedExperiment)
library(scater)

# Load base CD8 object (without final_modules, for T2D comparison)
CD8_T2D_ms_version <- readRDS("~/Desktop/T2D_CD8/data/CD8_T2D_ms_version.rds")

Idents(CD8_T2D_ms_version) <- 'conditions'
levels(CD8_T2D_ms_version)

# Subset to diabetic (Type 2) and non-diabetic donors
CD8_sce <- as.SingleCellExperiment(
  subset(CD8_T2D_ms_version, idents = c('Type 2', 'No Diabetes')),
  assay = "integrated"
)
milo <- Milo(CD8_sce)
milo <- buildGraph(milo, k = 20, d = 10, reduced.dim = "PCA")
milo <- makeNhoods(milo, prop = 1, k = 20, d=10, refined = TRUE, reduced_dims = "PCA")
plotNhoodSizeHist(milo)
milo <- countCells(milo, meta.data = as.data.frame(colData(milo)), sample="Pat_ID")
head(nhoodCounts(milo))

unique(CD8_T2D_ms_version$Pat_ID)
unique((CD8_T2D_ms_version$conditions))

# Build design data frame
milo_design <- data.frame(colData(milo))[,c("Pat_ID", "conditions")]
milo_design$Pat_ID <- as.factor(milo_design$Pat_ID)
milo_design$T2D    <- as.factor(milo_design$conditions)
milo_design <- distinct(milo_design)
rownames(milo_design) <- milo_design$Pat_ID

milo <- calcNhoodDistance(milo, d=10, reduced.dim = "PCA")
da_results <- testNhoods(milo, design = ~ conditions, design.df = milo_design)
head(da_results)
ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)

milo <- buildNhoodGraph(milo)
umap_pl     <- plotReducedDim(milo, dimred = "UMAP", colour_by="conditions", text_by = "celltype", text_size = 3, point_size=0.5) +guides(fill="none")
nh_graph_pl <- plotNhoodGraphDA(milo, da_results, layout="UMAP", alpha=0.1)
nh_graph_pl$layers[[1]] <- NULL
nh_graph_pl

da_results <- annotateNhoods(milo, da_results, coldata_col = "final_annotation")
head(da_results)
ggplot(da_results, aes(final_annotation_fraction)) + geom_histogram(bins=50)

plotDAbeeswarm(da_results, group.by = "final_annotation")
ggsave("plots/beeswarm_T2D.pdf",
       units = "cm", height = 22, width = 26)

nh_graph_pl <- plotNhoodGraphDA(milo, da_results, layout="UMAP", alpha = 0.1)
plot_data      <- nh_graph_pl$data
white_points   <- plot_data
colored_points <- plot_data[plot_data$colour_by != 0, ]
nh_graph_pl <- nh_graph_pl +
  geom_point(data = white_points, aes(x = x, y = y), color = "black", size = 5, shape = 21, stroke = 0.2, fill = "lightgrey", show.legend = FALSE)
nh_graph_pl <- nh_graph_pl +
  geom_point(data = colored_points, aes(x = x, y = y, fill = colour_by), size = 8, shape = 21, stroke = 0.2, show.legend = FALSE)
# Diverging scale: blue = enriched in No Diabetes, red = enriched in Type 2
nh_graph_pl <- nh_graph_pl +
  scale_fill_gradient2(low = "#2c7bb6", mid = "white", high = "#d7191c", midpoint = 0)
nh_graph_pl <- nh_graph_pl +
  theme_minimal() +
  theme(legend.position = "right") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
nh_graph_pl$layers[[1]] <- NULL
print(nh_graph_pl)
nh_graph_pl & NoLegend() & NoAxes()
ggsave("plots/nh_graph_T2D.pdf",
       units = "cm", height = 22, width = 26)
nh_graph_pl
ggsave("plots/nh_graph_T2D_w_legend.pdf",
       units = "cm", height = 22, width = 26)

# Export DA results per cluster (sorted by logFC, one file per annotation)
head(da_results)
dir.create("tables/MiloR_T2D")
setwd("tables/MiloR_T2D")
getwd()

sorted_df          <- da_results[order(da_results$logFC), ]
unique_annotations <- unique(sorted_df$final_annotation)
for (final_annotation in unique_annotations) {
  annotation_df <- sorted_df[sorted_df$final_annotation == final_annotation, ]
  file_name     <- paste0(final_annotation, "_sorted_by_FC.csv")
  write.table(annotation_df, file = file_name, sep = ",", row.names = FALSE)
}
write.csv(da_results, "da_results.csv")


# ==============================================================================
# Figures 3f & 3g — Cluster abundance per patient stratified by pathology and
#                   diabetes condition (source data)
# ==============================================================================

# Extract metadata for the full annotated object
md <- CD8_T2D_ms_version_final_modules@meta.data %>% as.data.table()

# Count cells per donor per cluster, retaining pathology and condition labels
cluster_abundance <- md[, .N, by = c("Pat_ID", "final_annotation", "pathology", "conditions")]

# Pivot to wide format: one row per donor, one column per cluster (fill = 0 for absent clusters)
cluster_abundance_wide <- dcast(cluster_abundance, Pat_ID + pathology + conditions ~ final_annotation, value.var = "N", fill = 0)

openxlsx::write.xlsx(cluster_abundance_wide,
                     "~/Desktop/T2D_CD8/figures/Figure_2/tables/cluster_abundance_with_pathology.xlsx",
                     asTable = FALSE, rowNames = TRUE)
