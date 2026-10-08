# RNA-seq Analysis of GSE116583 — Alveolar Macrophages

This project contains an R-based differential expression analysis of the GSE116583 dataset from NCBI GEO. The data consists of RPKM values from alveolar macrophages across three conditions: **Naive**, **Allo_02H** (2 hours post-allograft), and **Allo_24H** (24 hours post-allograft).

The dataset was obtained from [this paper](https://doi.org/10.1165/rcmb.2017-0430TR) titled *A Beginner's Guide to Analysis of RNA Sequencing Data*. More details can be found in the [References](#references) section below.

## Prerequisites

This project uses [`renv`](https://rstudio.github.io/renv/) for reproducible package management. To set up the environment:

```r
install.packages("renv")  # if you don't have renv installed
renv::restore()
```

The `.Rprofile` in the project root automatically activates the renv environment when R is started from this directory.

## Running the Analysis

From the project root, run:

```r
source("analysis.R")
```

Or from the command line:

```bash
Rscript analysis.R
```

## Analysis Workflow

The analysis is performed entirely in R using the following approach:

1. **Data Loading** — Reads the RPKM matrix (genes × samples) from `GSE116583_transplant.am.htseq.all.rpkm.txt`.
2. **Transformation** — Applies `log2(RPKM + 1)` transformation.
3. **Filtering** — Retains genes with `log2(RPKM+1) > 1` in at least 3 of 4 samples within at least one biological group.
4. **Differential Expression** — Uses **limma** with `eBayes(trend = TRUE)` on the log-transformed RPKM values. Three pairwise contrasts are tested:
   - `Allo_02H vs Naive`
   - `Allo_24H vs Naive`
   - `Allo_24H vs Allo_02H`
5. **Significance Thresholds** — Genes are classified as differentially expressed if `|log2FC| > 1` and `adjusted p-value (BH) < 0.05`.
6. **PCA** — Performed on the top 1000 most variable genes.
7. **Heatmap** — Top 50 DEGs from the `Allo_24H vs Naive` contrast.

> **Note:** The input data is RPKM (not raw counts), so results are exploratory in nature. If raw counts were available, DESeq2/edgeR/limma-voom would be preferred.

## Output Files

All results are written to the `results/` directory:

| File | Description |
|------|-------------|
| `limma_Allo_02H_vs_Naive.csv` | Full limma results for Allo_02H vs Naive |
| `limma_Allo_24H_vs_Naive.csv` | Full limma results for Allo_24H vs Naive |
| `limma_Allo_24H_vs_Allo_02H.csv` | Full limma results for Allo_24H vs Allo_02H |
| `Volcano_Allo_02H_vs_Naive.png` | Volcano plot for Allo_02H vs Naive |
| `Volcano_Allo_24H_vs_Naive.png` | Volcano plot for Allo_24H vs Naive |
| `Volcano_Allo_24H_vs_Allo_02H.png` | Volcano plot for Allo_24H vs Allo_02H |
| `PCA_plot.png` | PCA of log2(RPKM+1) on top 1000 variable genes |
| `Heatmap_top50_DEGs.png` | Heatmap of top 50 DEGs (Allo_24H vs Naive) |
| `sessionInfo.txt` | R session info for reproducibility |

## Project Structure

```
.
├── .Rprofile                         # Activates renv on R startup
├── analysis.R                        # Main analysis script
├── GSE116583_transplant.am.htseq.all.rpkm.txt  # Input RPKM matrix
├── renv/                             # renv environment
│   ├── activate.R
│   └── library/
├── renv.lock                         # Locked package versions
├── results/                          # Output directory (created by script)
└── README.md
```

## References

1. Koch CM, Chiu SF, Akbarpour M, Bharat A, Ridge KM, Bartom ET, Winter DR. A Beginner's Guide to Analysis of RNA Sequencing Data. *Am J Respir Cell Mol Biol*. 2018 Aug;59(2):145-157. doi: 10.1165/rcmb.2017-0430TR. PMID: 29624415; PMCID: PMC6096346.
