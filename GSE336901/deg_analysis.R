# EdgeR Differential Expression Analysis
# GSE336901_Fig3_FFPE_raw_counts_matrix_GEO.csv

# Load required libraries
library(edgeR)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# Set working directory
setwd("d:/DEG")

# 1. Load data
cat("Loading data...\n")
counts_data <- read.csv("GSE336901_Fig3_FFPE_raw_counts_matrix_GEO.csv", row.names = 1)

# Extract gene information
# gene_id is row names, gene_symbol is column 1
gene_symbols <- counts_data[, 1]
counts <- counts_data[, 2:ncol(counts_data)]

# Clean column names
colnames(counts) <- gsub("figure3_in_house_FFPE_", "", colnames(counts))
colnames(counts) <- gsub("_Sample_", "_", colnames(counts))
colnames(counts) <- gsub("_combined", "", colnames(counts))

# 2. Create group labels based on column names
# Resistance samples have "resistance_group" in name, sensitive have "sensitive_group"
group_labels <- ifelse(grepl("resistance", colnames(counts)), "resistance", "sensitive")
group <- factor(group_labels)
cat("Group distribution:\n")
print(table(group))

# 3. Create DGEList object
cat("\nCreating DGEList object...\n")
dge <- DGEList(counts = counts, group = group)
dge$genes <- data.frame(gene_symbol = gene_symbols)

# 4. Filter low-expression genes
cat("\nFiltering low-expression genes...\n")
keep <- filterByExpr(dge)
cat("Genes before filtering:", nrow(dge$counts), "\n")
cat("Genes after filtering:", sum(keep), "\n")
dge <- dge[keep, , keep.lib.sizes = FALSE]

# 5. Normalization (TMM)
cat("\nCalculating normalization factors...\n")
dge <- calcNormFactors(dge)
cat("Normalization factors:\n")
print(dge$samples$norm.factors)

# 6. Dispersion estimation
cat("\nEstimating dispersions...\n")
dge <- estimateDisp(dge)

# 7. Differential expression analysis
cat("\nPerforming exact test...\n")
et <- exactTest(dge)

# 8. Extract results
cat("\nExtracting results...\n")
results <- topTags(et, n = Inf)
results_df <- as.data.frame(results)

# Add gene symbol column
results_df$gene_symbol <- dge$genes$gene_symbol

# 9. Save full results
cat("\nSaving results...\n")
write.csv(results_df, "deg_results.csv", row.names = TRUE)

# 10. Identify significant genes
sig_genes <- results_df[results_df$FDR < 0.05 & abs(results_df$logFC) > 1, ]
cat("Number of significant genes (FDR < 0.05, |logFC| > 1):", nrow(sig_genes), "\n")
write.csv(sig_genes, "significant_genes.csv", row.names = TRUE)

# 11. Create MA plot using ggplot2
cat("\nCreating MA plot...\n")
ma_data <- data.frame(
  logCPM = results_df$logCPM,
  logFC = results_df$logFC,
  FDR = results_df$FDR
)
ma_data$sig <- "Not Significant"
ma_data$sig[ma_data$FDR < 0.05 & ma_data$logFC > 1] <- "Up in Resistance"
ma_data$sig[ma_data$FDR < 0.05 & ma_data$logFC < -1] <- "Up in Sensitive"

png("ma_plot.png", width = 800, height = 600)
ggplot(ma_data, aes(x = logCPM, y = logFC, color = sig)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = c("Up in Resistance" = "red", 
                               "Up in Sensitive" = "blue", 
                               "Not Significant" = "gray")) +
  theme_minimal() +
  labs(title = "MA Plot - Resistance vs Sensitive",
       x = "Log2 CPM",
       y = "Log2 Fold Change",
       color = "Status") +
  geom_hline(yintercept = 0, linetype = "dashed")
dev.off()

# 12. Create Volcano plot
cat("Creating volcano plot...\n")
volcano_data <- results_df
volcano_data$sig <- "Not Significant"
volcano_data$sig[volcano_data$FDR < 0.05 & volcano_data$logFC > 1] <- "Up in Resistance"
volcano_data$sig[volcano_data$FDR < 0.05 & volcano_data$logFC < -1] <- "Up in Sensitive"

png("volcano_plot.png", width = 800, height = 600)
ggplot(volcano_data, aes(x = logFC, y = -log10(FDR), color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up in Resistance" = "red", 
                               "Up in Sensitive" = "blue", 
                               "Not Significant" = "gray")) +
  theme_minimal() +
  labs(title = "Volcano Plot - Resistance vs Sensitive",
       x = "Log2 Fold Change",
       y = "-log10(FDR)",
       color = "Status") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed")
dev.off()

# 13. Create heatmap of top 50 DEGs
cat("Creating heatmap...\n")
top_genes <- head(order(et$table$PValue), 50)
top_counts <- cpm(dge)[top_genes, ]
top_counts_scaled <- t(scale(t(top_counts)))

annotation_col <- data.frame(Group = group)
rownames(annotation_col) <- colnames(counts)

png("heatmap_top50.png", width = 1000, height = 1200)
pheatmap(top_counts_scaled,
         annotation_col = annotation_col,
         show_rownames = FALSE,
         show_colnames = TRUE,
         fontsize_col = 10,
         main = "Top 50 Differentially Expressed Genes")
dev.off()

# 14. Summary statistics
cat("\n=== SUMMARY ===\n")
cat("Total genes analyzed:", nrow(results_df), "\n")
cat("Significant genes (FDR < 0.05):", sum(results_df$FDR < 0.05), "\n")
cat("Significant genes with |logFC| > 1:", nrow(sig_genes), "\n")
cat("Upregulated in resistance:", sum(results_df$FDR < 0.05 & results_df$logFC > 1), "\n")
cat("Upregulated in sensitive:", sum(results_df$FDR < 0.05 & results_df$logFC < -1), "\n")

# 15. Top 10 most significant genes
cat("\nTop 10 most significant genes:\n")
print(head(results_df[order(results_df$FDR), ], 10))

cat("\nAnalysis complete! Output files generated:\n")
cat("- deg_results.csv\n")
cat("- significant_genes.csv\n")
cat("- ma_plot.png\n")
cat("- volcano_plot.png\n")
cat("- heatmap_top50.png\n")