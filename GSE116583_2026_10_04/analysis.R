library(limma)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)

INPUT <- "GSE116583_transplant.am.htseq.all.rpkm.txt"
OUTDIR <- "results"
FDR <- 0.05
LFC <- 1
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

rpkm <- read.csv(INPUT, sep="\t", row.names = 1)
log_rpkm <- log2(rpkm + 1)

colnames(rpkm)
sample_info <- data.frame(
  sample = colnames(rpkm), 
  group = ifelse(grepl("Naive", colnames(rpkm)), "Naive", ifelse(grepl("02H", colnames(rpkm)), "Allo_02H", "Allo_24H")), 
  stringsAsFactors = FALSE
)
sample_info$group <- factor(sample_info$group, levels = c("Naive", "Allo_02H", "Allo_24H"))
