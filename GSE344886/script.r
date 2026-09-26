library(GEOquery)
library(DESeq2)

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

dim(counts)
dim(counts_filtered)

dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = metadata,
  design = ~ condition
)

dds <- DESeq(dds)

res <- results(
  dds,
  contrast = c(
    "condition",
    "bevacizumab-treated",
    "untreated control"
  )
)

res <- res[order(res$padj), ]

head(res)

write.csv(
  as.data.frame(res),
  "GSE344886_DESeq2_bevacizumab_vs_control.csv"
)

# PCA
vsd <- vst(dds, blind = FALSE)

plotPCA(vsd, intgroup = "condition")

# Sample-sample distances
sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)

heatmap(
  sampleDistMatrix,
  main = "Sample-to-Sample Distances"
)

# MA plot
plotMA(
  res,
  ylim = c(-5, 5),
  alpha = 0.05
)

# Volcano plot
res_df <- as.data.frame(res)

res_df$significant <- with(
  res_df,
  !is.na(padj) &
    padj < 0.05 &
    abs(log2FoldChange) >= 1
)

plot(
  res_df$log2FoldChange,
  -log10(res_df$pvalue),
  pch = 16,
  xlab = "log2 Fold Change",
  ylab = "-log10(p-value)",
  main = "Bevacizumab vs Control"
)
