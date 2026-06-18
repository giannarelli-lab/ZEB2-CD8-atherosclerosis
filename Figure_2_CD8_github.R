# =============================================================================
# Figure 2: CyTOF Validation of the CD8 T Cell Landscape
# =============================================================================
#
# Description:
#   This script processes mass cytometry (CyTOF) data from carotid plaque
#   specimens and projects it onto the scRNA-seq CD8 T cell reference map.
#   It covers:
#     - Import and preprocessing of FCS files (marker name cleaning, arcsinh
#       transformation with cofactor = 15, metadata annotation)
#     - Conversion to Seurat, RPCA-based batch integration across donors
#     - Unsupervised clustering and broad lineage annotation (CD4/CD8/NK)
#     - Label transfer of scRNA-seq cluster identities onto CyTOF CD8 cells
#       (Seurat MapQuery using scRNA-seq as reference)
#     - DotPlot of CyTOF marker expression grouped by functional category
#     - Feature plots of cytotoxicity markers on the projected UMAP
#
# Input files required:
#   - ~/Desktop/T2D_CD8/figures/Figure_2/CyTOF_analysis/data/*.fcs
#       Raw FCS files from CyTOF experiments (NKT panel, 18 files)
#   - ~/Desktop/T2D_CD8/figures/Figure_2/CyTOF_analysis/metadata/metadata.csv
#       Sample metadata table (Panel, Sample_ID, Tissue, etc.)
#   - ~/Desktop/T2D_CD8/data/CD8_T2D_ms_version.rds
#       scRNA-seq CD8 Seurat object (used as reference for label transfer)
#
# Output files — plots:
#   - plots/Markers_panel2.jpg                             (QC DotPlot)
#   - plots/UMAP_panel2.jpg                                (QC UMAP)
#   - plots/UMAP_all_cells_by_annotation.pdf               (Figure 2b)
#   - plots/UMAP_RNA_grey_CyTOF_colored.pdf                (Figure 2b)
#   - plots/Dotplot_grouped_by_functional_category.pdf     (Figure 2c)
#   - plots/Dotplot_grouped_by_functional_category_square.pdf (Figure 2c)
#   - plots/Feature_<marker>.pdf                           (Figure 2d, per marker)
#
# Output files — objects:
#   - object/seurat.rds   (integrated CyTOF Seurat object, saved mid-analysis)
#
# =============================================================================

# --- Libraries ----------------------------------------------------------------
library('Spectre')   # CyTOF/flow cytometry analysis (FCS import, arcsinh transform, colour plots)
library('flowCore')  # low-level FCS file handling
library(dplyr)       # data manipulation
library(reticulate)  # Python interface (used for scVI environment)
set.seed(123)        # global random seed for reproducibility

# Point reticulate to the conda environment that contains scVI (used if scVI-based steps are run)
use_condaenv("/Users/hauke/opt/anaconda3/envs/scvi-env/")

# Increase memory limit for large single-cell objects passed between R and future workers
options(future.globals.maxSize = 4 * 1024^3)

# --- Directory structure setup ------------------------------------------------
setwd("~/Desktop/T2D_CD8/figures/Figure_2/CyTOF_analysis")
getwd()
PrimaryDirectory <- getwd()

# Input: directory containing raw FCS files
setwd("~/Desktop/T2D_CD8/figures/Figure_2/CyTOF_analysis/data")
InputDirectory <- getwd()
setwd(PrimaryDirectory)

# Metadata: CSV with sample annotations
setwd("~/Desktop/T2D_CD8/figures/Figure_2/CyTOF_analysis/metadata")
MetaDirectory <- getwd()
setwd(PrimaryDirectory)

# Output: create and save to SpectreExports directory
dir.create("SpectreExports")
setwd("SpectreExports")
getwd()
OutputDirectory <- getwd()
setwd(PrimaryDirectory)


# ==============================================================================
# Figure 2b — CyTOF data preprocessing, integration, and label transfer
# ==============================================================================

# --- Import FCS files ---------------------------------------------------------

setwd(InputDirectory)
list.files(InputDirectory, ".fcs")  # confirm FCS files are present

# Read all FCS files into a named list (one data.table per file)
data.list <- Spectre::read.files(file.loc = InputDirectory,
                                 file.type = ".fcs",
                                 do.embed.file.names = TRUE)

# Quality check: verify column names and row counts match across files
check <- do.list.summary(data.list)
check$name.table   # channel name overview
check$ncol.check   # column count per file
check$nrow.check   # cell count per file
data.list[[1]]     # inspect first file

# --- Clean marker names -------------------------------------------------------
# CyTOF channels are named with metal isotope prefixes (e.g. "Nd144Di_CD103")
# This helper strips those prefixes and standardises specific marker names
# for consistent merging across FCS files from different acquisition batches

megalist_cleaned_names <- list()

clean_marker_names <- function(cols) {
  gsub("\\bPerf\\b", "Perforin", gsub("Granzyme_B","GranzymeB",gsub("PD-1","PD1",gsub("^.*?_","",gsub("^.*?Di_","", cols)))))
}

for (i in 1:length(data.list)){
  current_marker_cols <- NULL
  current_marker_cols <- grep("Di_", colnames(data.list[[i]]), value = TRUE)  # select mass channels
  current_names_clean <- clean_marker_names(current_marker_cols)
  megalist_cleaned_names[[i]] <- current_names_clean
  # Replace raw channel names (columns 3 onwards) with cleaned marker names
  colnames(data.list[[i]])[3:(length(current_names_clean) + 2)] <- current_names_clean
}

# --- Identify markers common to all FCS files ---------------------------------
# Only keep markers present in every file to ensure a consistent feature space
common_markers <- Reduce(intersect, megalist_cleaned_names)
common_markers
# Append non-marker bookkeeping columns
common_markers <- append(common_markers,c("FileName", "FileNo","Time_Time","Event_length_Event_length"))

# Remove DNA intercalator channels, bead channels, and non-target masses
common_markers <- common_markers[!common_markers %in% c("DNA","93Nb","102Pd","104Pd","105Pd","106Pd","108Pd","110Pd", "127I","140Ce","157Gd","181Ta","192Os","194Pt","B2M_CD298_Plaque","209Bi")]
common_markers

# Subset each file to the common marker set
for (j in 1:length(data.list)){
  bool <- NULL
  bool <- colnames(data.list[[j]]) %in% common_markers
  data.list[[j]] <- data.list[[j]][, bool, with = FALSE]
}

# Verify column counts match
length(common_markers)
size_list <- lapply(data.list, dim)
size_list

# --- Merge files and apply arcsinh transformation -----------------------------

cell.dat <- Spectre::do.merge.files(dat = data.list)  # row-bind all files
cell.dat

# Separate marker columns from bookkeeping metadata columns
meta_cols   <- c("FileName", "FileNo", "Time_Time", "Event_length_Event_length")
marker_cols <- setdiff(colnames(cell.dat), meta_cols)
markers <- marker_cols

# Arcsinh transformation with cofactor = 15 (standard for CyTOF data)
# Compresses high-expression values while preserving low-expression resolution
to.asinh <- markers
cofactor  <- 15
cell.dat  <- do.asinh(cell.dat, to.asinh, cofactor = cofactor)

# --- Add sample metadata ------------------------------------------------------
setwd(MetaDirectory)
meta.dat <- fread("metadata.csv")
meta.dat

# Restrict to the NKT panel used for this analysis
unique(meta.dat$Panel)
meta.dat <- meta.dat[meta.dat$Panel == "NKT panel with 18files and gating", ]
unique(meta.dat$Panel)

# Merge selected metadata columns into the cell data table
colnames(meta.dat)
sample.info <- meta.dat[,c(1:5,16,22)]
cell.dat    <- do.add.cols(cell.dat, "FileName", sample.info, "FileName", rmv.ext = TRUE)

# Redefine metadata columns after merging (additional columns were added)
meta_cols <- setdiff(colnames(cell.dat), marker_cols)

# --- Convert to SingleCellExperiment then Seurat ------------------------------
library(data.table)
library(S4Vectors)
library(SingleCellExperiment)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(ggplot2)

dt <- as.data.table(cell.dat)
unique(colnames(dt))

# Identify raw (untransformed) and arcsinh-transformed marker columns
raw_marker_cols   <- markers
asinh_marker_cols <- grep("_asinh$", colnames(dt), value = TRUE)

# Helper to strip metal prefixes and the _asinh suffix from column names
clean_marker_names <- function(cols) {
  gsub("^.*?Di_", "", gsub("_asinh$", "", cols))
}

raw_names_clean   <- markers
asinh_names_clean <- clean_marker_names(asinh_marker_cols)

# Transpose to marker x cell matrices (required by SingleCellExperiment)
raw_marker_mat   <- t(as.matrix(dt[, ..raw_marker_cols]))
asinh_marker_mat <- t(as.matrix(dt[, ..asinh_marker_cols]))

rownames(raw_marker_mat)   <- raw_names_clean
rownames(asinh_marker_mat) <- asinh_names_clean

# Assign unique cell barcodes
cell_ids <- paste0("cell_", seq_len(nrow(dt)))
colnames(raw_marker_mat)   <- cell_ids
colnames(asinh_marker_mat) <- cell_ids

# Package non-marker columns as cell-level metadata
meta_cols  <- setdiff(colnames(dt), c(raw_marker_cols, asinh_marker_cols))
cell_meta  <- DataFrame(dt[, ..meta_cols], row.names = cell_ids)

# Build SingleCellExperiment: raw counts + arcsinh-transformed logcounts
sce <- SingleCellExperiment(
  assays = list(
    counts    = raw_marker_mat,
    logcounts = asinh_marker_mat
  ),
  colData = cell_meta
)

# Convert SCE to Seurat object for downstream dimensionality reduction and integration
seurat <- as.Seurat(sce, counts = "counts", data = "logcounts")
seurat

DefaultAssay(seurat)           # confirm active assay
head(seurat@meta.data)         # confirm metadata transferred
all(colnames(seurat) == colnames(sce))  # sanity check: cell names must match

# --- Batch integration with RPCA ----------------------------------------------

# Rename the default assay for clarity and remove the old placeholder
seurat[['CyTOF_assay']] <- seurat[['originalexp']]
DefaultAssay(seurat)    <- "CyTOF_assay"
seurat[['originalexp']] <- NULL

# Convert to Seurat V5 assay format, required by IntegrateLayers
seurat[["CyTOF_assay5"]] <- as(object = seurat[["CyTOF_assay"]], Class = "Assay5")
DefaultAssay(seurat) <- "CyTOF_assay5"
table(seurat$Sample_ID, seurat$Tissue)  # inspect samples per tissue
Idents(seurat) <- "Tissue"
table(seurat$Tissue)

# Split the assay by donor (Sample_ID) to enable per-donor integration
seurat[["CyTOF_assay5"]] <- split(seurat[["CyTOF_assay5"]], f = seurat$Sample_ID)

levels(seurat)
rownames(seurat)  # confirm marker names are intact

table(seurat$Sample_ID, seurat$Tissue)
seurat <- FindVariableFeatures(seurat, selection.method = "vst")
all.markers <- rownames(seurat)
seurat <- ScaleData(seurat, features = all.markers)

# PCA for initial dimensionality reduction (15 PCs capture the main CyTOF variance)
seurat <- RunPCA(seurat, approx = F, npcs = 15)
ElbowPlot(seurat)  # inspect variance explained per PC

# RPCA integration: aligns donor-specific PCA embeddings to correct batch effects
seurat <- IntegrateLayers(
  object = seurat, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "RPCA",
  verbose = T, k.anchor = 50, dims = 1:15
)
seurat <- FindNeighbors(seurat, reduction = "RPCA", dims = 1:15)
# Leiden clustering (algorithm = 4) at resolution 0.5
seurat <- FindClusters(seurat, resolution = 0.5, cluster.name = "RPCA_clusters", algorithm = 4, method = "igraph")
seurat <- RunUMAP(seurat, dims = 1:15, reduction = "RPCA", reduction.name = "umap.RPCA", a = 2, b = 0.7)
DimPlot(seurat, reduction = "umap.RPCA", group.by = c("Tissue", "Sample_ID", "RPCA_clusters"))

# Identify cluster marker proteins to guide annotation
seurat <- JoinLayers(seurat)  # re-join split layers after integration
Idents(seurat) <- "RPCA_clusters"
RPCA_markers <- FindAllMarkers(seurat, min.pct = 0.25, only.pos = T, logfc.threshold = 0.2)

# Export significant markers to help with manual cluster annotation
significant_RPCA_markers <- subset(RPCA_markers, p_val_adj < 0.05)
openxlsx::write.xlsx(significant_RPCA_markers, asTable = T,
                     "tables/significant_RPCA_markers.xlsx")
DotPlot(seurat, features = rownames(seurat), cols = "RdYlBu", group.by = "RPCA_clusters") &
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("plots/Markers_panel2.jpg",
       units = "cm", height = 18, width = 22)
DimPlot(seurat, reduction = "umap.RPCA", group.by = "RPCA_clusters", pt.size = 2, label = T)
ggsave("plots/UMAP_panel2.jpg",
       units = "cm", height = 20, width = 22)

# Save integrated object to disk before annotation (checkpoint)
saveRDS(seurat,
        "object/seurat.rds")

# --- Broad lineage annotation -------------------------------------------------

# Reload from checkpoint
seurat <- readRDS("~/Desktop/T2D_CD8/figures/Figure_2/CyTOF_analysis/object/seurat.rds")

# Manually assign each RPCA cluster to a broad lineage (CD4, CD8, or NK)
# based on marker expression from the DotPlot above
Idents(seurat) <- "RPCA_clusters"
levels(seurat)
seurat <- RenameIdents(seurat,
                       `1`  = "CD4",
                       `2`  = "CD4",
                       `3`  = "CD8",
                       `4`  = "NK",
                       `5`  = "CD8",
                       `6`  = "CD8",
                       `7`  = "CD4",
                       `8`  = "CD8",
                       `9`  = "CD8",
                       `10` = "CD8",
                       `11` = "CD4",
                       `12` = "CD8",
                       `13` = "CD8",
                       `14` = "CD8",
                       `15` = "CD4",
                       `16` = "CD4",
                       `17` = "CD4"
)
seurat[["lineage"]] <- Idents(object = seurat)
DimPlot(seurat, reduction = "umap.RPCA", group.by = "lineage", pt.size = 2, label = T)

# --- Label transfer: project CyTOF CD8 cells onto scRNA-seq reference ---------
library(ggplot2)

# Subset CyTOF data to CD8 T cells only
CD8_CytOF <- subset(seurat, idents = "CD8")
CD8_CytOF
rownames(CD8_CytOF)  # inspect available CyTOF markers

# Load scRNA-seq CD8 reference object
CD8_T2D_ms_version <- readRDS("~/Desktop/T2D_CD8/data/CD8_T2D_ms_version.rds")

# Rebuild UMAP model on the reference (return.model = TRUE required for MapQuery)
CD8_T2D_ms_version <- RunUMAP(CD8_T2D_ms_version, dims = 1:30, verbose = TRUE, return.model = T, assay = "integrated")
DimPlot(CD8_T2D_ms_version, reduction = "umap", group.by = "final_annotation", pt.size = 2, label = T) +
  scale_y_reverse() + scale_x_reverse()

# Identify overlapping features between CyTOF panel and scRNA-seq genes
matching_genes <- intersect(rownames(CD8_CytOF), rownames(CD8_T2D_ms_version))

# Rename CyTOF protein markers to their HGNC gene symbol equivalents
# so they can be used for anchor-based integration with the scRNA-seq reference
original_matrix <- GetAssayData(CD8_CytOF, assay = "CyTOF_assay", layer = "data")
rownames(original_matrix)

gene_renaming <- c(
  "CD45"    = "PTPRC",
  "CD57"    = "B3GAT1",
  "CD103"   = "ITGAE",
  "CD4"     = "CD4",
  "CD8"     = "CD8A",
  "Perforin" = "PRF1",
  "CD127"   = "IL7R",
  "CD123"   = "IL3RA",
  "TIGIT"   = "TIGIT",
  "2B4"     = "CD244",
  "CD27"    = "CD27",
  "TIM3"    = "HAVCR2",
  "NKG2C"   = "KLRC2",
  "NKp30"   = "NCR3",
  "CD56"    = "NCAM1",
  "NKG2A"   = "KLRC1",
  "KIR3DL1" = "KIR3DL1",
  "CD69"    = "CD69",
  "LAG3"    = "LAG3",
  "NKG2D"   = "KLRK1",
  "CCR7"    = "CCR7",
  "CD3"     = "CD3D",
  "NKp44"   = "NCR2",
  "DNAM1"   = "CD226",
  "HLADR"   = "HLA-DRA",
  "PD1"     = "PDCD1",
  "GranzymeB" = "GZMB"
)

# Apply renaming to the expression matrix
new_matrix              <- original_matrix
rownames(new_matrix)    <- plyr::mapvalues(
  rownames(new_matrix),
  from = names(gene_renaming),
  to   = gene_renaming,
  warn_missing = FALSE
)

# Store the renamed matrix as a new assay for integration
new_assay                        <- CreateAssayObject(data = new_matrix)
CD8_CytOF[["assay_for_mapping"]] <- new_assay
DefaultAssay(CD8_CytOF)          <- "assay_for_mapping"
DefaultAssay(CD8_CytOF)

# Re-check shared features after renaming
matching_genes <- intersect(rownames(CD8_CytOF), rownames(CD8_T2D_ms_version))

# Find transfer anchors between the scRNA-seq reference and CyTOF query
transfer.anchors <- FindTransferAnchors(reference = CD8_T2D_ms_version, query = CD8_CytOF, dims = 1:15,
                                        reference.reduction = "pca", features = matching_genes)

# Project CyTOF cells onto the scRNA-seq UMAP and transfer cluster labels
CD8_CytOF.query <- MapQuery(anchorset = transfer.anchors, reference = CD8_T2D_ms_version, query = CD8_CytOF,
                            refdata = list(final_annotation = "final_annotation"), reference.reduction = "pca",
                            reduction.model = "umap")

DimPlot(CD8_CytOF.query, group.by = "predicted.final_annotation", label = TRUE, repel = T) +
  scale_y_reverse() + scale_x_reverse()

# Define consistent colour map for CD8 cluster annotations (shared with scRNA-seq figures)
cd8_color_map <- c(
  "C0- GZMK+ EM"                  = "#BC3C29",
  "C1- CCR7+ IL7R+ Naïve"         = "#565974",
  "C2- CX3CR1+ ADGRG1+ EMRA"      = "#137AB2",
  "C3- IL-7R+ LYAR+ EM"           = "#9CB39E",
  "C4- CCL4+ EM"                  = "#DCCE86",
  "C5- IL-2RA+ CTLA4+ RM"         = "#649F62",
  "C6- KLRB1+ MAIT"               = "#328963",
  "C7 -IFN-g+ CCL4L2+ EM"         = "#5C9497",
  "C8- HAVCR+ HOPX+ EMRA"         = "#9681A6",
  "C09- ITGA1+ NEAT1+ RM"         = "#DA579A",
  "C10 - FCGR3B+ B3GAT1+ EMRA"   = "#E9626B",
  "C11- IFIT1+ ISG20+ IFN-response" = "#E2822F",
  "C12- GATA3+ Naïve mixed"        = "#AA645C",
  "C13- LMNA+ EM"                  = "#6A3D9A"
)

DimPlot(CD8_T2D_ms_version, group.by = "final_annotation", label = TRUE, repel = T,
        cols = cd8_color_map) +
  scale_y_reverse() + scale_x_reverse()

DimPlot(CD8_CytOF.query, group.by = "predicted.final_annotation", label = TRUE, repel = T,
        cols = cd8_color_map) +
  scale_y_reverse() + scale_x_reverse()

# --- Export UMAP overlays for Figure 2b ---------------------------------------
library(ggplot2)
library(patchwork)

# Plot 1: all cells coloured by predicted cluster annotation
umap1_plot <- ggplot(merged_umap, aes(x = UMAP_1, y = UMAP_2, color = renamecelltype)) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_color_manual(values = cd8_color_map) +
  scale_y_reverse() +
  scale_x_reverse() +
  labs(title = NULL) +
  theme_minimal() +
  theme(
    panel.grid  = element_blank(),
    axis.text   = element_blank(),
    axis.ticks  = element_blank(),
    axis.title  = element_blank(),
    plot.title  = element_text(hjust = 0.5)
  ) & NoLegend()

ggsave(
  filename = "plots//UMAP_all_cells_by_annotation.pdf",
  plot     = umap1_plot,
  width = 9, height = 10, units = "cm", dpi = 300
)

# Plot 2: scRNA-seq cells shown in grey, CyTOF cells coloured by predicted annotation
# This visualisation highlights where CyTOF cells land within the scRNA-seq landscape
umap2_plot <- ggplot() +
  geom_point(
    data  = subset(merged_umap, source == "CD8_T2D_ms_version"),
    aes(x = UMAP_1, y = UMAP_2),
    color = "grey80", size = 1.5, alpha = 0.4
  ) +
  geom_point(
    data  = subset(merged_umap, source == "CD8_CytOF.query"),
    aes(x = UMAP_1, y = UMAP_2, color = renamecelltype),
    size = 1.5, alpha = 0.9
  ) +
  scale_color_manual(values = cd8_color_map) +
  scale_y_reverse() +
  scale_x_reverse() +
  labs(title = NULL) +
  theme_minimal() +
  theme(
    panel.grid  = element_blank(),
    axis.text   = element_blank(),
    axis.ticks  = element_blank(),
    axis.title  = element_blank(),
    plot.title  = element_text(hjust = 0.5)
  ) & NoLegend()

ggsave(
  filename = "plots/UMAP_RNA_grey_CyTOF_colored.pdf",
  plot     = umap2_plot,
  width = 9, height = 10, units = "cm", dpi = 300
)


# ==============================================================================
# Figure 2c — DotPlot of CyTOF markers grouped by functional category
# ==============================================================================

library(dplyr)
library(ggplot2)
library(scales)
library(grid)

# Define markers grouped by functional category (determines facet order)
marker_categories <- list(
  "Memory / Differentiation state / Residency" = c(
    "CD27", "CCR7", "CD69", "CD127", "CD123", "CD103"
  ),
  "EMRA" = c(
    "CD57", "CD45RA", "KIR3DL1"
  ),
  "Cytotoxicity" = c(
    "Perforin", "GranzymeB"
  ),
  "Exhaustion / Inhibitory" = c(
    "TIGIT", "PD1", "TIM3", "LAG3"
  ),
  "Innate-like receptors" = c(
    "NKG2D", "2B4", "DNAM1", "NKG2A", "CD56", "NKp44", "NKp30", "NKG2C"
  )
)

# Flatten to ordered vector for x-axis placement
all_markers <- unname(unlist(marker_categories))

# Build a lookup: marker -> category (used for facet assignment)
marker_group_lut <- rep(names(marker_categories), lengths(marker_categories))
names(marker_group_lut) <- all_markers

# Safety check: remove markers not found in the CyTOF assay
present_markers <- all_markers[all_markers %in% rownames(CD8_CytOF.query[["CyTOF_assay"]])]
missing_markers <- setdiff(all_markers, present_markers)
if (length(missing_markers) > 0) {
  message("Markers not found and removed: ", paste(missing_markers, collapse = ", "))
}
all_markers <- present_markers

# Recreate lookup table after removing absent markers
marker_group_lut <- marker_group_lut[all_markers]

# Generate raw DotPlot object from Seurat (data extracted below for custom ggplot)
dp_raw <- DotPlot(
  CD8_CytOF.query,
  features  = all_markers,
  group.by  = "predicted.final_annotation",
  cols      = "RdBu",
  dot.scale = 6
)

# Extract underlying data and annotate with marker group and ordered factors
dp_data <- dp_raw$data %>%
  mutate(
    marker_group  = factor(
      marker_group_lut[as.character(features.plot)],
      levels = names(marker_categories)
    ),
    features.plot = factor(features.plot, levels = all_markers),
    id            = factor(id, levels = cluster_lvls)
  )

# Colour y-axis labels to match the CD8 cluster colour scheme
y_label_colours <- cd8_color_map[cluster_lvls]
y_label_colours[is.na(y_label_colours)] <- "grey40"

# Custom ggplot DotPlot: dot size = % cells expressing, colour = scaled expression
dotplot <- ggplot(
  dp_data,
  aes(x = features.plot, y = id, size = pct.exp, color = avg.exp.scaled)
) +
  geom_point() +

  # Diverging colour scale centred at 0, capped at ±2.5
  scale_color_distiller(
    palette   = "RdBu",
    direction = -1,
    limits    = c(-2.5, 2.5),
    oob       = squish,
    name      = "Scaled\nexpression",
    guide     = guide_colorbar(barwidth = 0.5, barheight = 4, title.position = "top")
  ) +

  scale_size_continuous(
    name   = "% cells\nexpressing",
    range  = c(0.3, 10),
    breaks = c(10, 25, 50, 75, 100),
    guide  = guide_legend(override.aes = list(colour = "grey50"))
  ) +

  # Facet columns by marker functional category; each panel scales to its own width
  facet_grid(. ~ marker_group, scales = "free_x", space = "free_x") +

  labs(x = NULL, y = NULL) +

  theme_minimal(base_size = 9) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y      = element_text(size = 7.5, colour = y_label_colours, face = "bold"),
    panel.grid.major = element_line(colour = "grey93", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(colour = "grey35", fill = NA, linewidth = 0.5),
    strip.background = element_rect(fill = "grey90", colour = "grey35", linewidth = 0.5),
    strip.text       = element_text(face = "bold", size = 8),
    panel.spacing    = unit(0.25, "cm"),
    legend.position  = "right",
    legend.title     = element_text(size = 7.5, hjust = 0.5),
    legend.text      = element_text(size = 7),
    legend.key.size  = unit(0.35, "cm"),
    legend.box       = "vertical",
    legend.spacing.y = unit(0.3, "cm"),
    plot.margin      = margin(6, 6, 6, 6)
  )

ggsave("plots/Dotplot_grouped_by_functional_category.pdf",
       plot = dotplot, width = 48, height = 10, units = "cm", dpi = 300)

ggsave("plots/Dotplot_grouped_by_functional_category_square.pdf",
       plot = dotplot, width = 32, height = 18, units = "cm", dpi = 300)


# ==============================================================================
# Figure 2d — Feature plots of cytotoxicity markers on projected UMAP
# ==============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(viridis)

# Scale the CyTOF assay before feature extraction
DefaultAssay(CD8_CytOF.query) <- "CyTOF_assay"
CD8_CytOF.query <- ScaleData(CD8_CytOF.query)

# Define markers to visualise (cytotoxicity panel)
marker_list <- c("Perforin", "GranzymeB")
marker_list <- intersect(marker_list, rownames(CD8_CytOF.query))  # keep only present markers

# Attach cell barcodes to UMAP coordinate data frames for joining
umap_coords_CYTOF$cell <- rownames(umap_coords_CYTOF)
umap_coords_RNA$cell   <- rownames(umap_coords_RNA)

# Diagnostic: confirm CyTOF UMAP rows match the Seurat object
cat("CyTOF UMAP rows:", nrow(umap_coords_CYTOF), "\n")
cat("CyTOF Seurat cells:", ncol(CD8_CytOF.query), "\n")
cat("Shared cells:", length(intersect(umap_coords_CYTOF$cell, colnames(CD8_CytOF.query))), "\n")

dir.create("plots", showWarnings = FALSE)

# Loop through each marker and create a UMAP feature plot
# scRNA-seq cells are shown in grey as a background reference
for (marker in marker_list) {

  # Fetch scaled expression values for this marker
  expr_df <- FetchData(CD8_CytOF.query, vars = marker, layer = "scale.data")
  expr_df$cell <- rownames(expr_df)
  colnames(expr_df)[1] <- "marker_expr"

  # Join expression values to the CyTOF UMAP coordinate table
  umap_coords_CYTOF_marker <- umap_coords_CYTOF %>%
    dplyr::select(-dplyr::any_of("marker_expr")) %>%
    dplyr::left_join(expr_df, by = "cell")

  # Prepare the RNA background (expression set to 0 for grey colouring)
  umap_coords_RNA_marker <- umap_coords_RNA %>%
    dplyr::select(-dplyr::any_of("marker_expr"))
  umap_coords_RNA_marker$marker_expr <- 0

  # Align columns before row-binding RNA and CyTOF data
  common_cols <- intersect(colnames(umap_coords_RNA_marker), colnames(umap_coords_CYTOF_marker))
  merged_marker_df <- rbind(
    umap_coords_RNA_marker[, common_cols],
    umap_coords_CYTOF_marker[, common_cols]
  )

  # Separate into RNA background and CyTOF foreground for layered plotting
  cytof_cells <- merged_marker_df %>%
    dplyr::filter(source == "CD8_CytOF.query") %>%
    dplyr::arrange(marker_expr)  # plot low-expression cells first (avoid overplotting)
  rna_cells <- merged_marker_df %>%
    dplyr::filter(source == "CD8_T2D_ms_version")

  p <- ggplot() +
    geom_point(
      data  = rna_cells,
      aes(x = umap_1, y = umap_2),
      color = "grey80", size = 1, alpha = 0.4
    ) +
    geom_point(
      data  = cytof_cells,
      aes(x = umap_1, y = umap_2, color = marker_expr),
      size = 1.5
    ) +
    scale_color_gradientn(colors = viridis::inferno(256), na.value = "grey90") +
    scale_x_reverse() +
    scale_y_reverse() +
    labs(title = marker, color = marker) +
    theme_minimal() +
    theme(
      panel.grid  = element_blank(),
      axis.title  = element_blank(),
      axis.text   = element_blank(),
      axis.ticks  = element_blank()
    )

  ggsave(
    filename = paste0("plots/Feature_", marker, ".pdf"),
    plot = p, width = 12, height = 10, units = "cm", dpi = 300
  )
}
