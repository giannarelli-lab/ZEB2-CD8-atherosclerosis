# ZEB2 Directs Senescent and Cytotoxic Terminal Differentiation of CD8⁺ T Cells in Atherosclerosis

---

## Abstract

Senescent cytotoxic CD8⁺ T cells accumulate in atherosclerotic plaques and are associated with adverse cardiovascular outcomes, yet the transcriptional mechanisms underlying this pathogenic state remain poorly understood. Using single-cell transcriptomic profiling of human carotid atherosclerotic plaques, we identified a population of CD8⁺ TEMRA cells exhibiting a coordinated senescence and cytotoxic gene program that preferentially localise to plaque shoulder and luminal regions associated with inflammation and plaque vulnerability. Trajectory and gene-regulatory network analyses identified ZEB2 (zinc finger E-box binding homeobox 2), encoded within a coronary artery disease susceptibility locus, as a master regulator of this state. CRISPR-mediated deletion of ZEB2 in primary human CD8⁺ T cells impaired Granzyme B-mediated cytotoxicity and senescence. Accordingly, CD8⁺ T cell–specific or Granzyme B-restricted *Zeb2* deletion in hypercholesterolaemic mice reduced atherosclerosis burden, macrophage accumulation, plaque necrosis, and inflammation. These findings establish ZEB2 as a transcriptional regulator that couples cytotoxic differentiation and senescence in CD8⁺ T cells and establishes this program as a driver of atherosclerosis progression.

---

## Repository overview

This repository contains the R code used to generate all main and extended figures in the manuscript. Each script corresponds to one figure and is self-contained. Scripts are numbered to match figure numbers in the paper.

```
.
├── Figure_1_CD8_github.R   # scRNA-seq CD8 landscape, UCell module scoring, correlations
├── Figure_2_CD8_github.R   # CyTOF processing, RPCA integration, label transfer
├── Figure_3_CD8_github.R   # TEMRA module associations, MiloR differential abundance
├── Figure_4_CD8_github.R   # Slingshot pseudotime, tradeSeq GAM, ZEB2 Nebulosa
├── Figure_5_CD8_github.R   # ZEB2 DEG analysis, TF-module correlations, KO heatmap
├── Figure_7_CD8_github.R   # Mouse PBMC CyTOF (ZEB2 KO), CATALYST pipeline
└── README.md
```

Each script contains:
- A header block listing input files, output files (by panel), and reporting standards
- Section headers aligned to individual figure panels
- Inline comments explaining analytical choices

---

## Data availability

| Dataset | Description | Accession |
|---------|-------------|-----------|
| Human carotid scRNA-seq (Giannarelli et al.) | Full PBMC/plaque immune atlas (reference dataset) | [accession] |
| Human carotid CD8 scRNA-seq (this study) | CD8 T cell subset, all donors | [accession] |
| Human CyTOF — NKT panel | FCS files, carotid plaque specimens | [accession] |
| Mouse PBMC CyTOF — 16 W | FCS files, *Zeb2* KO/WT mice | [accession] |

---

## System requirements

### Software

| Software | Version tested | Notes |
|----------|---------------|-------|
| R | ≥ 4.3 | |
| Bioconductor | ≥ 3.18 | |

### R packages

#### Single-cell RNA-seq (Figures 1–5)

| Package | Version tested | Source |
|---------|---------------|--------|
| Seurat | ≥ 5.0 | CRAN |
| SeuratDisk | ≥ 0.0.0.9021 | GitHub (mojaveazure) |
| SeuratWrappers | ≥ 0.3.5 | GitHub (satijalab) |
| UCell | ≥ 2.6 | Bioconductor |
| slingshot | ≥ 2.10 | Bioconductor |
| tradeSeq | ≥ 1.16 | Bioconductor |
| TrajectoryUtils | ≥ 1.10 | Bioconductor |
| SingleCellExperiment | ≥ 1.24 | Bioconductor |
| miloR | ≥ 2.0 | Bioconductor |
| scater | ≥ 1.30 | Bioconductor |
| Nebulosa | ≥ 1.12 | Bioconductor |

#### CyTOF (Figures 2 & 7)

| Package | Version tested | Source |
|---------|---------------|--------|
| CATALYST | ≥ 1.26 | Bioconductor |
| flowCore | ≥ 2.14 | Bioconductor |
| Spectre | ≥ 1.1.0 | GitHub (ImmuneDynamics) |

#### Visualisation and utilities

| Package | Version tested | Source |
|---------|---------------|--------|
| ggplot2 | ≥ 3.5 | CRAN |
| patchwork | ≥ 1.2 | CRAN |
| cowplot | ≥ 1.1 | CRAN |
| ComplexHeatmap | ≥ 2.18 | Bioconductor |
| corrplot | ≥ 0.92 | CRAN |
| EnhancedVolcano | ≥ 1.20 | Bioconductor |
| fmsb | ≥ 0.7 | CRAN |
| ggalluvial | ≥ 0.12 | CRAN |
| pheatmap | ≥ 1.0.12 | CRAN |
| viridis | ≥ 0.6 | CRAN |
| RColorBrewer | ≥ 1.1 | CRAN |
| colorspace | ≥ 2.1 | CRAN |
| dplyr | ≥ 1.1 | CRAN |
| data.table | ≥ 1.15 | CRAN |
| tidyr | ≥ 1.3 | CRAN |
| readxl | ≥ 1.4 | CRAN |
| openxlsx | ≥ 4.2 | CRAN |
| matrixStats | ≥ 1.3 | CRAN |
| reticulate | ≥ 1.36 | CRAN |
| scCustomize | ≥ 2.1 | CRAN |

### Hardware

All analyses were run on a macOS system with ≥ 32 GB RAM. 

---

## Reproducibility

- All scripts call `set.seed(123)` (or `set.seed(1234)` for Figure 7) at the top.
- Random seeds for UMAP and FlowSOM are passed explicitly to each function call.
- Package versions used for the final figures are listed above. Minor version differences may produce slightly different UMAP layouts but will not affect statistical conclusions.

---

## License

This code is released under the **MIT License** — see [LICENSE](LICENSE) for details.

---

## Contact

For questions about the code or data, please open a [GitHub Issue](https://github.com/giannarelli-lab/ZEB2-CD8-atherosclerosis/issues)
