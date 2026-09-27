library(GEOquery)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(tidyverse)

counts <- read.csv("GSE344886_GEO_raw_count_6samples.tsv", sep="\t", row.names = 1)

gse <- getGEO(filename = "GSE344886_series_matrix.txt", GSEMatrix = TRUE)

metadata <- pData(gse)

metadata[, c("geo_accession","title","cell line:ch1","cell type:ch1","tissue:ch1","treatment:ch1")]

table(metadata$`treatment:ch1`)

# mapping each counts column name to the appropriate GEO accession ID
sample_map <- c(
  "BEV_1" = "GSM9986991",
  "BEV_2" = "GSM9986992",
  "BEV_3" = "GSM9986993",
  "CTRL_1" = "GSM9986994",
  "CTRL_2" = "GSM9986995",
  "CTRL_3" = "GSM9986996"
)

colnames(counts) <- sample_map[colnames(counts)]

metadata <- metadata[colnames(counts), ]

# Creating condition factors
metadata$condition <- factor(metadata$`treatment:ch1`)
table(metadata$condition)

metadata$condition <- relevel(metadata$condition,ref = "untreated control")

keep <- rowSums(counts >= 10) >= 3
counts_filtered <- counts[keep, ]

# Checking the difference between the two data
dim(counts)
dim(counts_filtered)

dds <- DESeqDataSetFromMatrix(countData = counts_filtered, colData = metadata, design = ~ condition)

dds <- DESeq(dds)

res <- results(dds, contrast = c("condition","bevacizumab-treated","untreated control"))

res <- res[order(res$padj), ]

write.csv(as.data.frame(res),"GSE344886_DESeq2_bevacizumab_vs_control.csv")

# PCA
vsd <- vst(dds, blind = FALSE)

plotPCA(vsd, intgroup = "condition")

# Sample-sample distances
sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)

rownames(sampleDistMatrix) <- colnames(vsd)
colnames(sampleDistMatrix) <- colnames(vsd)
annotation <- data.frame(Treatment = colData(vsd)$condition)
rownames(annotation) <- colnames(vsd)
pheatmap(sampleDistMatrix, annotation_col = annotation, annotation_row = annotation, main = "sample-sample distance")

# MA plot
plotMA(res, ylim = c(-5, 5), alpha = 0.05)
