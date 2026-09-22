# EdgeR Differential Expression Analysis

## Overview

This repository contains an end-to-end differential expression analysis of FFPE (Formalin-Fixed Paraffin-Embedded) raw count data using the EdgeR package in R. The analysis compares gene expression between **resistance** and **sensitive** groups from the GSE336901 dataset.

## Dataset

- **Source:** GSE336901_Fig3_FFPE_raw_counts_matrix_GEO.csv
- **Samples:** 10 total (5 resistance group, 5 sensitive group)
- **Data Type:** Raw gene expression counts from FFPE tissue samples

## Analysis Summary

| Metric | Value |
|--------|-------|
| Total genes analyzed | 13,213 (after filtering) |
| Significant genes (FDR < 0.05, \|logFC\| > 1) | 81 |
| Upregulated in Resistance | 31 genes |
| Upregulated in Sensitive | 50 genes |

## Files

### Input
- `GSE336901_Fig3_FFPE_raw_counts_matrix_GEO.csv` - Raw count matrix with gene annotations

### Analysis Script
- `deg_analysis.R` - R script performing the complete EdgeR analysis

### Output
- `deg_results.csv` - Full results table with all 13,213 genes (logFC, logCPM, PValue, FDR)
- `significant_genes.csv` - Filtered significant genes (81 genes meeting criteria)
- `ma_plot.png` - MA plot visualization showing expression vs fold change
- `volcano_plot.png` - Volcano plot showing statistical significance vs fold change
- `heatmap_top50.png` - Heatmap of top 50 differentially expressed genes

## Top 10 Most Significant Differentially Expressed Genes

| Gene | logFC | logCPM | PValue | FDR |
|------|-------|--------|--------|-----|
| MSTRG.6117\|SFTPA2 | 6.36 | 6.98 | 8.07e-08 | 0.00054 |
| MSTRG.35943\|TPRG1 | -6.45 | 6.70 | 8.24e-08 | 0.00054 |
| MSTRG.40337\|SPOCK1 | -9.63 | 4.32 | 2.68e-06 | 0.0118 |
| MSTRG.35768\|RNA5SP149 | 4.59 | 8.42 | 9.67e-06 | 0.0118 |
| MSTRG.47827\|MMP16 | -5.22 | 4.43 | 1.07e-05 | 0.0118 |
| MSTRG.11999\|SLC41A2 | -9.18 | 3.91 | 1.26e-05 | 0.0118 |
| MSTRG.4483\|RNA5S2 | 4.47 | 11.73 | 1.70e-05 | 0.0118 |
| MSTRG.4486\|RNA5S5 | 4.47 | 11.79 | 1.72e-05 | 0.0118 |
| MSTRG.4491\|RNA5S10 | 4.47 | 11.82 | 1.72e-05 | 0.0118 |
| MSTRG.4496\|RNA5S15 | 4.46 | 11.84 | 1.80e-05 | 0.0118 |

## Requirements

- R (version 4.5.2 or higher)
- R packages:
  - `edgeR` - For differential expression analysis
  - `ggplot2` - For data visualization
  - `pheatmap` - For heatmap generation
  - `RColorBrewer` - For color palettes

## Installation

```r
# Install required packages
install.packages(c("edgeR", "ggplot2", "pheatmap", "RColorBrewer"))
```

## Usage

```bash
# Run the analysis
Rscript deg_analysis.R
```

Or within R:
```r
source("deg_analysis.R")
```

## Analysis Pipeline

1. **Data Loading** - Read CSV file and extract count matrix
2. **Filtering** - Remove low-expression genes using `filterByExpr()`
3. **Normalization** - TMM (Trimmed Mean of M-values) normalization
4. **Dispersion Estimation** - Calculate biological variability
5. **Statistical Testing** - Exact test for two-group comparison
6. **Results Extraction** - Get logFC, logCPM, PValue, FDR for all genes
7. **Significance Filtering** - Identify DEGs (FDR < 0.05, |logFC| > 1)
8. **Visualization** - Generate MA plot, volcano plot, and heatmap

## Methodology

### EdgeR Approach
- Uses negative binomial distribution to model count data
- Accounts for biological variability through dispersion estimation
- Implements exact test for comparing two groups
- Applies FDR correction for multiple testing

### Significance Criteria
- **FDR < 0.05**: False Discovery Rate below 5% (statistical significance)
- **|logFC| > 1**: Absolute log2 fold change greater than 1 (biological relevance)

## License

This project is available as open source under the MIT License.

## Author

Generated using EdgeR (R) for differential expression analysis of FFPE gene expression data.