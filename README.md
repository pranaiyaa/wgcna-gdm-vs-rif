# WGCNA Analysis: GDM vs RIF

This project analyzes gene co-expression networks in patients with Gestational Diabetes Mellitus (GDM) and Recurrent Implantation Failure (RIF) using microarray data. The goal is to identify gene modules correlated with disease traits and extract hub genes relevant to each condition.


# Pipeline Overview

# 1. Normalization
- Platform: Affymetrix HTA 2.0 Array
- Method: RMA normalization using `oligo`
- Gene symbol mapping via `hta20transcriptcluster.db`
- Aggregated probes to gene-level expression

# 2. Differential Expression Analysis
- Tool: `limma`
- Groups: GDM vs Control
- Outputs: 
  - Full DEGs
  - Filtered DEGs (logFC > 0.1 & adj.p < 0.01)
  - Volcano plot
  - Heatmap of top 50 genes

# 3. WGCNA (RIF condition)
- Sample clustering & soft-thresholding
- Network construction with `blockwiseModules`
- Module eigengene analysis
- Trait correlation heatmaps
- Extraction of condition-specific modules and genes


# Tools Used

| Purpose | Tools |
|--------|--------|
| Normalization | `oligo`, `hta20transcriptcluster.db` |
| DEG Analysis | `limma`, `ggplot2`, `pheatmap` |
| Network Analysis | `WGCNA`, `genefilter` |
| Metadata handling | `tidyverse`, `readxl`, `AnnotationDbi` |


# Repository Structure
wgcna-gdm-vs-rif/
├── data/ # Input data files (CELs, metadata, processed expression)
│ ├── raw/ # Raw CEL files or GEO data
│ ├── normalized_expression.csv
│ └── traits_metadata.csv
│
├── scripts/ # R scripts for each stage of the pipeline
│ ├── 01_normalization.R
│ ├── 02_differential_expression.R
│ └── 03_wgcna_rif_pipeline.R
│
├── results/ # Output results, plots, and module data
│ ├── volcano_plot.png
│ ├── heatmap_top50_genes.png
│ ├── module_trait_heatmap.png
│ └── genes_by_module/ # CSVs of genes in key modules
│
├── README.md # Project overview and usage guide


