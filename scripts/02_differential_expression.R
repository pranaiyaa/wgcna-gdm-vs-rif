# ───────────────────────────────────────────────
# Differential Expression Analysis using limma
# ───────────────────────────────────────────────

library(limma)
library(ggplot2)
library(pheatmap)

# Load normalized expression
expr <- read.csv("data/normalized_expression.csv", row.names = 1)
expr_matrix <- as.matrix(expr)

# Extract condition labels from column names
sample_names <- colnames(expr_matrix)
condition <- factor(ifelse(grepl("GDM", sample_names, ignore.case = TRUE), "GDM", "Control"))

# Create design matrix and contrast
design <- model.matrix(~ 0 + condition)
colnames(design) <- levels(condition)
contrast <- makeContrasts(GDM_vs_Control = GDM - Control, levels = design)

# Fit linear model and apply empirical Bayes
fit <- lmFit(expr_matrix, design)
fit <- contrasts.fit(fit, contrast)
fit <- eBayes(fit)

# Extract DEGs
deg_results <- topTable(fit, coef = "GDM_vs_Control", adjust = "BH", number = Inf)
write.csv(deg_results, "results/deg_results_full.csv")

# Label genes
deg_results$significance <- "Not Significant"
deg_results$significance[deg_results$logFC > 0.1 & deg_results$adj.P.Val < 0.01] <- "Upregulated"
deg_results$significance[deg_results$logFC < -0.1 & deg_results$adj.P.Val < 0.01] <- "Downregulated"

# Save filtered DEGs
filtered <- deg_results[deg_results$adj.P.Val < 0.01 & abs(deg_results$logFC) > 0.1, ]
write.csv(filtered, "results/deg_results_filtered.csv")

# Plot volcano
ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val), color = significance)) +
  geom_point() +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "log2 Fold Change", y = "-log10 Adjusted P-value")
ggsave("results/volcano_plot.png")

# Heatmap of top 50 DEGs
top50_genes <- rownames(filtered[order(filtered$adj.P.Val),])[1:50]
heatmap_data <- expr_matrix[top50_genes, ]
heatmap_scaled <- t(scale(t(heatmap_data)))
annotation_col <- data.frame(Condition = condition)
rownames(annotation_col) <- colnames(heatmap_scaled)

pheatmap(heatmap_scaled, annotation_col = annotation_col, main = "Top 50 DEGs")
ggsave("results/heatmap_top50_genes.png")
