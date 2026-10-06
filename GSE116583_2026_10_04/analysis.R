## RNA-seq analysis of GSE116583 alveolar macrophage RPKM data
## Input: RPKM matrix (genes x samples), not raw counts, so use limma on log2(RPKM+1).
## Run from project root; .Rprofile activates the project renv. No new packages.
## NOTE: RPKM matrix (not raw counts) -> results are exploratory; raw counts +
## DESeq2/edgeR/limma-voom would be preferred if available.

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

## ---- Load and inspect ------------------------------------------------------
rpkm <- read.csv(INPUT, sep = "\t", row.names = 1, check.names = FALSE)
cat("Dimensions (genes x samples):", nrow(rpkm), "x", ncol(rpkm), "\n")
cat("Total missing values:", sum(is.na(rpkm)), "\n")
stopifnot(!any(duplicated(rownames(rpkm))))

sample_info <- data.frame(
  sample = colnames(rpkm),
  group = ifelse(grepl("Naive", colnames(rpkm)), "Naive",
                 ifelse(grepl("02H", colnames(rpkm)), "Allo_02H",
                        ifelse(grepl("24H", colnames(rpkm)), "Allo_24H", NA_character_))),
  stringsAsFactors = FALSE
)
stopifnot(!any(is.na(sample_info$group)))
sample_info$group <- factor(sample_info$group, levels = c("Naive", "Allo_02H", "Allo_24H"))
print(table(sample_info$group))
stopifnot(all(colnames(rpkm) == sample_info$sample))

## ---- Transform and filter ---------------------------------------------------
log_rpkm <- log2(rpkm + 1)

# Group-aware filter: gene must pass RPKM > 1 (log2(RPKM+1) > 1) in at least
# 3 of 4 samples within at least one biological group.
pass <- log_rpkm > 1
grp_pass <- sapply(levels(sample_info$group), function(g) {
  idx <- which(sample_info$group == g)
  rowSums(pass[, idx, drop = FALSE], na.rm = TRUE) >= 3
})
keep <- apply(grp_pass, 1, any)
filtered_log_rpkm <- log_rpkm[keep, ]
cat("Genes passing filter:", nrow(filtered_log_rpkm), "of", nrow(log_rpkm), "\n")

## ---- limma model ------------------------------------------------------------
design <- model.matrix(~ 0 + group, data = sample_info)
colnames(design) <- levels(sample_info$group)

fit <- lmFit(filtered_log_rpkm, design)

contrast_matrix <- makeContrasts(
  Allo_02H_vs_Naive = Allo_02H - Naive,
  Allo_24H_vs_Naive = Allo_24H - Naive,
  Allo_24H_vs_Allo_02H = Allo_24H - Allo_02H,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
# trend = TRUE: mean-variance trend appropriate for RPKM-like continuous data
fit2 <- eBayes(fit2, trend = TRUE)

## ---- Results per contrast ---------------------------------------------------
ctr_names <- colnames(contrast_matrix)
for (ctr in ctr_names) {
  res <- topTable(fit2, coef = ctr, number = Inf, adjust.method = "BH")
  res$diffexpressed <- "NO"
  res$diffexpressed[res$logFC > LFC & res$adj.P.Val < FDR] <- "UP"
  res$diffexpressed[res$logFC < -LFC & res$adj.P.Val < FDR] <- "DOWN"
  res$gene <- rownames(res)
  write.csv(res, file = file.path(OUTDIR, paste0("limma_", ctr, ".csv")), row.names = FALSE)
  cat(ctr, "- UP:", sum(res$diffexpressed == "UP"),
      " DOWN:", sum(res$diffexpressed == "DOWN"), "\n")
  
  # Top 10 by adjusted p-value among significant genes
  sig <- res[res$diffexpressed != "NO", ]
  top10 <- sig[order(sig$adj.P.Val), ][seq_len(min(10, nrow(sig))), ]
  p_volcano <- ggplot(res, aes(x = logFC, y = -log10(adj.P.Val), color = diffexpressed)) +
    geom_point(alpha = 0.4, size = 1.8) +
    scale_color_manual(values = c("DOWN" = "dodgerblue3", "NO" = "grey70", "UP" = "firebrick3")) +
    geom_vline(xintercept = c(-LFC, LFC), linetype = "dashed", color = "black", linewidth = 0.4) +
    geom_hline(yintercept = -log10(FDR), linetype = "dashed", color = "black", linewidth = 0.4) +
    geom_text_repel(data = top10, aes(label = gene), max.overlaps = 15, size = 3.5, color = "black") +
    theme_classic(base_size = 12) +
    labs(title = gsub("_", " ", ctr), x = "log2 Fold Change",
         y = "-log10 Adjusted P-value", color = "Status")
  ggsave(file.path(OUTDIR, paste0("Volcano_", ctr, ".png")), plot = p_volcano, width = 6.5, height = 5.5)
}

## ---- PCA on top 1000 variable genes ----------------------------------------
n_top <- min(1000, nrow(filtered_log_rpkm))
var_genes <- head(order(apply(filtered_log_rpkm, 1, var), decreasing = TRUE), n_top)
pca <- prcomp(t(filtered_log_rpkm[var_genes, ]), center = TRUE, scale. = FALSE)

pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  group = sample_info$group,
  sample = sample_info$sample
)
percent_var <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 4, alpha = 0.8) +
  theme_bw(base_size = 12) +
  labs(title = "PCA of Log2(RPKM+1)",
       x = paste0("PC1 (", percent_var[1], "% variance)"),
       y = paste0("PC2 (", percent_var[2], "% variance)")) +
  scale_color_brewer(palette = "Set1")

ggsave(file.path(OUTDIR, "PCA_plot.png"), plot = p_pca, width = 6, height = 5)

## ---- Heatmap of top 50 DEGs (Allo_24H vs Naive) ----------------------------
res_24h <- read.csv(file.path(OUTDIR, "limma_Allo_24H_vs_Naive.csv"))
sig <- res_24h[res_24h$diffexpressed != "NO", ]
top50_genes <- head(sig$gene[order(sig$adj.P.Val)], 50)

if (length(top50_genes) > 0) {
  mat <- filtered_log_rpkm[top50_genes, ]
  annotation_col <- data.frame(Group = sample_info$group)
  rownames(annotation_col) <- sample_info$sample
  pheatmap(
    mat,
    scale = "row",
    annotation_col = annotation_col,
    show_rownames = TRUE,
    show_colnames = FALSE,
    color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100),
    fontsize_row = 7,
    main = "Top 50 DEGs: Allo_24H vs Naive",
    filename = file.path(OUTDIR, "Heatmap_top50_DEGs.png"),
    width = 7,
    height = 8
  )
} else {
  message("No significant genes for heatmap.")
}

## ---- Reproducibility --------------------------------------------------------
sink(file.path(OUTDIR, "sessionInfo.txt"))
print(sessionInfo())
sink()
cat("--- Session info written to results/sessionInfo.txt ---\n")
