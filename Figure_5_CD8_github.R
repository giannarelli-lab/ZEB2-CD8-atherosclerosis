# =============================================================================
# Figure 5: ZEB2 as a Regulator of CD8 T Cell Senescence and Cytotoxicity
# =============================================================================
#
# Description:
#   This script characterises the functional role of ZEB2 in CD8 T cells by:
#     - Differential gene expression between ZEB2-positive and ZEB2-negative
#       CD8 T cells across all clusters (MAST test)
#     - Spearman correlation of ZEB2 (and related TF) expression with
#       functional module scores; corrplot and individual scatter plots
#     - Export of correlation coefficients and p-values as source data
#     - CyTOF-based heatmap from a ZEB2 knock-out (KO) mouse experiment,
#       showing protein-level changes across stimulation conditions and days
#
# Input files required:
#   - CD8_T2D_ms_version (Seurat object in memory from earlier sections,
#     or reload from ~/Desktop/T2D_CD8/data/CD8_T2D_ms_version.rds)
#   - ~/Desktop/T2D_CD8/figures/Figure_3/Heatmap_data/data/*.xlsx
#       Per-sample CyTOF median expression tables (ZEB2 KO experiment)
#       Expected filename format: Day<N>_<Genotype>_<Stim>.xlsx
#       (e.g. Day1_Control_Unstim.xlsx, Day3_ZEB2KO_Stim.xlsx)
#
# Output files — plots:
#   - figures/Figure_3/DEG_ZEB2/ZEB2_pos_vs_neg/[volcano + tables] (Figure 5a/c)
#   - plots/corr_plot_pies.pdf                                      (Figure 5b)
#   - plots/corr_plot_numbers.pdf                                    (Figure 5b)
#   - plots/ZEB2_T_Scen.pdf                                         (Figure 5b)
#   - plots/ZEB2_Cytotoxicity.pdf                                   (Figure 5b)
#   - plots/ZEB2_SENMAYO.pdf                                        (Figure 5b)
#   - plots/ZEB2_SASP.pdf                                           (Figure 5b)
#   - figures/Figure_3/Heatmap_data/plots/heatmap_median_log.pdf    (Figure 5e)
#
# Output files — source data / tables:
#   - figures/Figure_3/DEG_ZEB2/ZEB2_pos_vs_neg/sigs_p_adj_ZEB2_DEG.xlsx  (Figure 5a/c)
#   - figures/Figure_3/DEG_ZEB2/ZEB2_pos_vs_neg/ZEB2_DEG.xlsx             (Figure 5a/c)
#   - figures/Figure_3/tables/zeb2_modules_for_correlation.xlsx            (Figure 5b)
#   - tables/correlation_with_pvalues.xlsx                                  (Figure 5b)
#
# =============================================================================

# --- Libraries ----------------------------------------------------------------
library(Seurat)      # scRNA-seq object handling
library(dplyr)       # data manipulation
library(data.table)  # fast tabular operations
library(corrplot)    # correlation matrix visualisation
set.seed(123)        # global random seed for reproducibility


# ==============================================================================
# Figures 5a & 5c — DEG analysis: ZEB2-positive vs ZEB2-negative CD8 T cells
# ==============================================================================
# Cells are split into ZEB2+ (any detectable expression) and ZEB2- (dropout).
# MAST (Mixed-effects model of Association of Single-cell Transcriptomics) is
# used because it models the bimodal distribution (expressed vs not expressed)
# characteristic of scRNA-seq data.
# Reference: Finak et al., Genome Biology 2015.

Idents(CD8_T2D_ms_version) <- "final_annotation"
levels(CD8_T2D_ms_version)
DimPlot(CD8_T2D_ms_version, label = T, repel = T, group.by = "final_annotation")

# Split cells based on ZEB2 expression: positive = any expression > 0
ZEB_pos <- subset(x = CD8_T2D_ms_version, subset = ZEB2 > 0)
ZEB_neg <- subset(x = CD8_T2D_ms_version, subset = ZEB2 > 0, invert = T)  # invert = TRUE selects ZEB2 == 0
ZEB_pos$ZEB_status <- "positive"
ZEB_neg$ZEB_status <- "negative"

# Slim objects to RNA assay only before merging (removes large unnecessary assays)
ZEB_pos <- DietSeurat(ZEB_pos, assays = "RNA")
ZEB_neg <- DietSeurat(ZEB_neg, assays = "RNA")
ZEB2    <- merge(x = ZEB_pos, y = ZEB_neg)

# Normalise and scale the merged object for DEG testing
ZEB2 <- NormalizeData(ZEB2)
ZEB2 <- ScaleData(ZEB2)
scCustomize::VlnPlot_scCustom(ZEB2, features = "ZEB2", plot_boxplot = T, group.by = "ZEB_status")  # QC: confirm split

# Run MAST differential expression: ZEB2+ vs ZEB2-
Idents(ZEB2) <- "ZEB_status"
ZEB2_DEG    <- FindMarkers(ZEB2, test.use = "MAST", ident.1 = "positive", ident.2 = "negative")

# Quick volcano plot for inspection
EnhancedVolcano::EnhancedVolcano(ZEB2_DEG, x = "avg_log2FC", y = "p_val_adj", lab = rownames(ZEB2_DEG))

# Compute additional BH-adjusted p-values from the raw p-value column
ZEB2_DEG$p_val_adj_BH <- p.adjust(ZEB2_DEG$p_val, method = "BH")

# Export: significant genes only (p_val_adj < 0.05)
sigs_p_adj_ZEB2_DEG        <- ZEB2_DEG[ZEB2_DEG$p_val_adj < 0.05, ]
sigs_p_adj_ZEB2_DEG$genes  <- rownames(sigs_p_adj_ZEB2_DEG)
openxlsx::write.xlsx(sigs_p_adj_ZEB2_DEG, asTable = T,
                     "~/Desktop/T2D_CD8/figures/Figure_3/DEG_ZEB2/ZEB2_pos_vs_neg/sigs_p_adj_ZEB2_DEG.xlsx")

# Export: full DEG table (all genes tested)
ZEB2_DEG$genes <- rownames(ZEB2_DEG)
openxlsx::write.xlsx(ZEB2_DEG, asTable = T,
                     "~/Desktop/T2D_CD8/figures/Figure_3/DEG_ZEB2/ZEB2_pos_vs_neg/ZEB2_DEG.xlsx")


# ==============================================================================
# Figure 5b — Correlation of ZEB2 and related TFs with functional module scores
# ==============================================================================

# Reload base CD8 object (ensures clean state for correlation analysis)
CD8_T2D_ms_version <- readRDS("~/Desktop/T2D_CD8/data/CD8_T2D_ms_version.rds")

# Fetch expression of ZEB2 and four co-regulated transcription factors
zeb2_expr <- FetchData(CD8_T2D_ms_version, vars = c("ZEB2", "TBX21", "EOMES", "TOX", "TCF7", "GATA3"))

# Extract UCell module scores and cluster annotation from metadata
md <- as.data.frame(CD8_T2D_ms_version@meta.data)
md <- md %>%
  select(matches("UCELL_kNN"), final_annotation)
head(zeb2_expr)
head(md)

# Add row names as explicit column for merging
zeb2_expr <- zeb2_expr %>%
  tibble::rownames_to_column(var = "cell_barcode")
md <- md %>%
  tibble::rownames_to_column(var = "cell_barcode")

# Merge module scores and TF expression per cell
merged_df <- left_join(md, zeb2_expr, by = "cell_barcode")
head(merged_df)

openxlsx::write.xlsx(merged_df, asTable = T,
                     "~/Desktop/T2D_CD8/figures/Figure_3/tables/zeb2_modules_for_correlation.xlsx")

# --- Corrplot -----------------------------------------------------------------
library(dplyr)
library(corrplot)

# Remove composite/derived scores to avoid redundancy in the correlation matrix
merged_df$SENMAYO_CYTOTOXIC_UCell_kNN <- NULL
merged_df$TSEN_CYTOTOXIC_UCell_kNN    <- NULL

# Select numeric columns: all UCell scores + TF expression values
corr_data <- merged_df %>%
  select(matches("UCELL"), ZEB2, TBX21, EOMES, TOX, TCF7, GATA3)

# Spearman correlation matrix (cell-level; n can be large, robust method preferred)
cor_matrix    <- cor(corr_data, method = "spearman")
color_palette <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(200))

pdf("plots/corr_plot_pies.pdf", width = 32 / 2.54, height = 32 / 2.54)
corrplot(cor_matrix, method = "pie", type = "upper",
         col = color_palette, tl.col = "black", tl.cex = 1, cl.cex = 0.8)
dev.off()

pdf("plots/corr_plot_numbers.pdf", width = 32 / 2.54, height = 32 / 2.54)
corrplot(cor_matrix, method = "circle", type = "upper",
         col = color_palette, tl.col = "black", tl.cex = 1,
         cl.cex = 0.8, addCoef.col = "black")
dev.off()

# --- Individual scatter plots: ZEB2 expression vs each module score -----------
# Linear regression line with 95% CI shading (red fill)

# ZEB2 vs Cytotoxicity score
ggplot(merged_df, aes(x = ZEB2, y = CYTOTOXIC_UCell_kNN)) +
  geom_smooth(method = "lm", se = TRUE, color = "red", fill = "#fca5a5") +
  labs(x = "ZEB2 Expression", y = "Cytotoxicity Score") +
  theme_minimal() +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  theme(
    panel.grid  = element_blank(),
    axis.line   = element_line(color = "black"),
    axis.ticks  = element_line(color = "black")
  )
ggsave("plots/ZEB2_T_Scen.pdf", units = "cm", height = 12, width = 12)

# ZEB2 vs T cell senescence score
ggplot(merged_df, aes(x = ZEB2, y = T_CELL_SENESCENCE_UCell_kNN)) +
  geom_smooth(method = "lm", se = TRUE, color = "red", fill = "#fca5a5") +
  labs(x = "ZEB2 Expression", y = "T cell senescence Score") +
  theme_minimal() +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  theme(
    panel.grid  = element_blank(),
    axis.line   = element_line(color = "black"),
    axis.ticks  = element_line(color = "black")
  )
ggsave("plots/ZEB2_Cytotoxicity.pdf", units = "cm", height = 12, width = 12)

# ZEB2 vs Cytotoxicity (duplicate saved with alternative file name)
ggplot(merged_df, aes(x = ZEB2, y = CYTOTOXIC_UCell_kNN)) +
  geom_smooth(method = "lm", se = TRUE, color = "red", fill = "#fca5a5") +
  labs(x = "ZEB2 Expression", y = "Cytotoxicity Score") +
  theme_minimal() +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  theme(
    panel.grid  = element_blank(),
    axis.line   = element_line(color = "black"),
    axis.ticks  = element_line(color = "black")
  )
ggsave("plots/ZEB2_T_Scen.pdf", units = "cm", height = 12, width = 12)

# ZEB2 vs T cell senescence (duplicate)
ggplot(merged_df, aes(x = ZEB2, y = T_CELL_SENESCENCE_UCell_kNN)) +
  geom_smooth(method = "lm", se = TRUE, color = "red", fill = "#fca5a5") +
  labs(x = "ZEB2 Expression", y = "T cell senescence Score") +
  theme_minimal() +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  theme(
    panel.grid  = element_blank(),
    axis.line   = element_line(color = "black"),
    axis.ticks  = element_line(color = "black")
  )
ggsave("plots/ZEB2_Cytotoxicity.pdf", units = "cm", height = 12, width = 12)

# ZEB2 vs SENMAYO score
ggplot(merged_df, aes(x = ZEB2, y = SENMAYO_UCell_kNN)) +
  geom_smooth(method = "lm", se = TRUE, color = "red", fill = "#fca5a5") +
  labs(x = "ZEB2 Expression", y = "SENMAYO Score") +
  theme_minimal() +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  theme(
    panel.grid  = element_blank(),
    axis.line   = element_line(color = "black"),
    axis.ticks  = element_line(color = "black")
  )
ggsave("plots/ZEB2_SENMAYO.pdf", units = "cm", height = 12, width = 12)

# ZEB2 vs SASP score
ggplot(merged_df, aes(x = ZEB2, y = SASP_UCell_kNN)) +
  geom_smooth(method = "lm", se = TRUE, color = "red", fill = "#fca5a5") +
  labs(x = "ZEB2 Expression", y = "SASP Score") +
  theme_minimal() +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  theme(
    panel.grid  = element_blank(),
    axis.line   = element_line(color = "black"),
    axis.ticks  = element_line(color = "black")
  )
ggsave("plots/ZEB2_SASP.pdf", units = "cm", height = 12, width = 12)

# --- Export full correlation matrix with p-values (source data) ---------------
# Pairwise Spearman tests are computed manually to retain p-values
cor_mat <- matrix(NA, ncol = ncol(corr_data), nrow = ncol(corr_data))
p_mat   <- matrix(NA, ncol = ncol(corr_data), nrow = ncol(corr_data))
colnames(cor_mat) <- rownames(cor_mat) <- colnames(corr_data)
colnames(p_mat)   <- rownames(p_mat)   <- colnames(corr_data)

for (i in 1:ncol(corr_data)) {
  for (j in 1:ncol(corr_data)) {
    test         <- cor.test(corr_data[[i]], corr_data[[j]], method = "spearman")
    cor_mat[i, j] <- test$estimate
    p_mat[i, j]   <- test$p.value
  }
}

cor_mat <- as.data.frame(cor_mat)
p_mat   <- as.data.frame(p_mat)

# Move row names to a column for Excel export
cor_mat <- cor_mat %>% tibble::rownames_to_column(var = "Feature")
p_mat   <- p_mat   %>% tibble::rownames_to_column(var = "Feature")

openxlsx::write.xlsx(list("Correlation" = cor_mat, "P_values" = p_mat),
                     file = "tables/correlation_with_pvalues.xlsx", asTable = TRUE)


# ==============================================================================
# Figure 5e — CyTOF heatmap: ZEB2 KO vs Control across stimulation conditions
# ==============================================================================
# This section processes per-sample median CyTOF expression files from a
# mouse ZEB2 KO experiment and generates a clustered heatmap.
# Filename convention: Day<N>_<Genotype>_<Stim>.xlsx
#   Genotype: Control or ZEB2KO
#   Stim: Stim (anti-CD3/CD28) or Unstim
#   Day: 1, 3, or 5

library(readxl)
library(dplyr)
library(stringr)
library(tibble)
library(purrr)
library(ComplexHeatmap)  # publication-quality heatmaps with annotations
library(circlize)        # colour scale definition for ComplexHeatmap
library(ggplot2)

# --- Settings -----------------------------------------------------------------
data_dir <- "~/Desktop/T2D_CD8/figures/Figure_3/Heatmap_data/data"
pattern  <- "\\.xlsx$"

# Log10(x + 1) transform: compresses dynamic range of median CyTOF intensities
log_fun <- function(x) log10(x + 1)

# --- Helper functions ---------------------------------------------------------

# Reads one Excel file, log transforms, and returns per-marker medians
read_one_file_medians <- function(path) {
  df <- read_excel(path)

  # Drop event number column if present (not a marker)
  df <- df %>% select(-matches("^Event\\s*#"))

  # Handle European decimal notation (comma → period)
  df <- df %>%
    mutate(across(everything(), ~ as.numeric(str_replace(as.character(.x), ",", "."))))

  # Log transform and compute median per marker across all cells in the sample
  meds <- df %>%
    mutate(across(everything(), ~ log_fun(.x))) %>%
    summarise(across(everything(), ~ median(.x, na.rm = TRUE)))

  meds
}

# Parses experimental metadata from the filename
# Expected format: Day<N>_<Genotype>_<Stim>  (e.g. Day1_ZEB2KO_Stim)
parse_meta <- function(path) {
  nm <- tools::file_path_sans_ext(basename(path))

  m <- str_match(nm, "^(Day\\d+)_?(WT|KO)_?(Stim|Unstim)$")

  if (any(is.na(m))) {
    return(tibble(Sample = nm, Day = NA, Genotype = NA, Stim = NA))
  }

  tibble(
    Sample   = nm,
    Day      = m[,2],
    Genotype = m[,3],
    Stim     = m[,4]
  )
}

# --- Main: read and assemble matrix -------------------------------------------
files <- list.files(data_dir, pattern = pattern, full.names = TRUE)

meta     <- map_dfr(files, parse_meta)
med_list <- map(files, read_one_file_medians)

# Combine into a marker × sample matrix
mat <- bind_rows(med_list) %>%
  mutate(Sample = meta$Sample) %>%
  column_to_rownames("Sample") %>%
  as.matrix()

mat <- t(mat)  # orient: markers as rows, samples as columns

# Z-score across samples per marker (centred on 0, scaled by SD)
mat_z          <- t(scale(t(mat)))
mat_z[is.na(mat_z)] <- 0  # replace NA (from zero-variance markers) with 0

# Column metadata aligned to the matrix column order
col_meta <- meta %>%
  filter(Sample %in% colnames(mat)) %>%
  slice(match(colnames(mat), Sample))

# Top annotation bar: Day, Genotype, Stimulation condition
ha <- HeatmapAnnotation(
  Day      = col_meta$Day,
  Genotype = col_meta$Genotype,
  Stim     = col_meta$Stim
)

plot_mat <- mat_z  # use Z-scored matrix for the final figure

# Quick QC heatmap with hierarchical clustering (for inspection only)
Heatmap(
  plot_mat,
  name             = ifelse(identical(plot_mat, mat_z), "Z", "Median log"),
  cluster_rows     = TRUE,
  cluster_columns  = F,
  show_row_names   = TRUE,
  show_column_names = TRUE,
  column_names_rot = 90,
  border           = T
)

# --- Final heatmap for figure -------------------------------------------------
library(circlize)

# Diverging colour scale: black → white → red (low → mid → high Z-score)
col_fun <- colorRamp2(
  c(min(plot_mat), 0, max(plot_mat)),
  c("black", "white", "red")
)

# Cap values at ±3 to reduce influence of extreme outliers
plot_mat <- pmax(pmin(plot_mat, 3), -3)

# Manually define column order: group by genotype × stimulation, ordered by day
col_order <- c(
  "Day1_Control_Unstim",
  "Day3_Control_Unstim",
  "Day5_Control_Unstim",
  "Day1_ZEB2KO_Unstim",
  "Day3_ZEB2KO_Unstim",
  "Day5_ZEB2KO_Unstim",
  "Day1_Control_Stim",
  "Day3_Control_Stim",
  "Day5_Control_Stim",
  "Day1_ZEB2KO_Stim",
  "Day3_ZEB2KO_Stim",
  "Day5_ZEB2KO_Stim"
)
plot_mat <- plot_mat[, col_order, drop = FALSE]

rownames(plot_mat)

# Manually define row order to group markers by functional category
row_order <- c(
  "164Dy_CD69_(Dy164Di)",          # activation / residency
  "166Er_CD25_(Er166Di)",           # activation (IL-2 receptor alpha)
  "170Er_CD38_(Er170Di)",           # activation / senescence
  "144Nd_CD103_(Nd144Di)",          # tissue residency
  "143Nd_CD45RA_(Nd143Di)",         # differentiation state (EMRA marker)
  "148Nd_CD16_(Nd148Di)",           # NK-like / EMRA
  "165Ho_LAG3_(Ho165Di)",           # exhaustion
  "153Eu_PD-1_(Eu153Di)",           # exhaustion
  "156Gd_Granzyme_A_(Gd156Di)",     # cytotoxicity
  "171Yb_Granzyme_B_(Yb171Di)",     # cytotoxicity
  "147Sm_Perforin_(Sm147Di)",       # cytotoxicity
  "113In_CD57_(In113Di)",           # senescence / EMRA
  "158Gd_IL-2_(Gd158Di)",           # effector cytokine
  "141Pr_IFNG_(Pr141Di)",           # effector cytokine
  "152Sm_TNFa_(Sm152Di)",           # effector cytokine
  "173Yb_Tbet_(Yb173Di)"            # transcription factor
)
plot_mat <- plot_mat[row_order, , drop = FALSE]

pdf(
  file      = "~/Desktop/T2D_CD8/figures/Figure_3/Heatmap_data/plots/heatmap_median_log.pdf",
  width     = 12,
  height    = 10,
  useDingbats = FALSE
)

Heatmap(
  plot_mat,
  name              = ifelse(identical(plot_mat, mat_z), "Z", "Median log"),
  col               = col_fun,
  cluster_rows      = FALSE,     # row order is fixed manually above
  cluster_columns   = FALSE,     # column order is fixed manually above
  show_row_names    = TRUE,
  show_column_names = TRUE,
  column_names_rot  = 90,
  border            = FALSE,
  rect_gp           = grid::gpar(col = "white", lwd = 8)  # white grid lines between cells
)

dev.off()
