# ──────────────────────────────────────────────────────────────
# 03_wgcna_rif_pipeline.R
# Co-expression network analysis (WGCNA) for RIF vs Control
# Author: pranaiyaa
# ──────────────────────────────────────────────────────────────

# ── 1. Load Required Libraries ────────────────────────────────
if (!require("WGCNA")) install.packages("WGCNA"); library(WGCNA)
if (!require("genefilter")) install.packages("genefilter"); library(genefilter)
if (!require("readxl")) install.packages("readxl"); library(readxl)
if (!require("DESeq2")) BiocManager::install("DESeq2"); library(DESeq2)

options(stringsAsFactors = FALSE)
allowWGCNAThreads()  # Enable multithreading

# ── 2. Load and Preprocess Expression Data ─────────────────────
data_path <- "data/RIF_input_final_cleaned.xlsx"
expr_raw <- read_excel(data_path)
expr_raw <- as.data.frame(expr_raw)
rownames(expr_raw) <- expr_raw[[1]]
expr_raw[[1]] <- NULL
expr_raw[] <- lapply(expr_raw, function(x) as.numeric(as.character(x)))
datExpr <- t(as.matrix(expr_raw))

# ── 3. Optional DEG Calculation with DESeq2 ────────────────────
count_data <- t(datExpr)
col_data <- data.frame(
  row.names = colnames(count_data),
  condition = ifelse(grepl("^RIF", colnames(count_data)), "RIF", "Control")
)
dds <- DESeqDataSetFromMatrix(countData = round(count_data),
                              colData = col_data,
                              design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
deg_df <- as.data.frame(res)
deg_df$Gene <- rownames(deg_df)
write.csv(deg_df, "results/DEG_results_r.csv", row.names = FALSE)

# ── 4. Filter Bad Samples and Genes ────────────────────────────
gsg <- goodSamplesGenes(datExpr)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

# ── 5. Sample Clustering ───────────────────────────────────────
sampleTree <- hclust(dist(datExpr), method = "average")
plot(sampleTree, main = "Sample Clustering", cex = 0.6)

# Optional removal of outliers
clust <- cutreeStatic(sampleTree, cutHeight = 100, minSize = 10)
datExpr_clean <- datExpr[clust != 0, ]

# ── 6. Soft-Threshold Power Selection ──────────────────────────
powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(datExpr_clean, powerVector = powers, verbose = 5)

# Visualize scale-free topology
par(mfrow = c(1, 2))
plot(sft$fitIndices$Power, -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq,
     xlab = "Soft Threshold", ylab = "Scale Free Topology Fit", type = "n",
     main = "Scale Independence")
text(sft$fitIndices$Power, -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq,
     labels = powers, col = "red")
abline(h = 0.8, col = "red")

plot(sft$fitIndices$Power, sft$fitIndices$mean.k.,
     xlab = "Soft Threshold", ylab = "Mean Connectivity", type = "n",
     main = "Mean Connectivity")
text(sft$fitIndices$Power, sft$fitIndices$mean.k., labels = powers, col = "red")

softPower <- 10  # Set manually based on previous plot

# ── 7. Network Construction and Module Detection ───────────────
net <- blockwiseModules(datExpr_clean,
                        power = softPower,
                        TOMType = "signed",
                        networkType = "signed",
                        minModuleSize = 30,
                        mergeCutHeight = 0.25,
                        numericLabels = FALSE,
                        pamRespectsDendro = FALSE,
                        verbose = 3)

# ── 8. Plot Dendrogram with Modules ────────────────────────────
for (i in seq_along(net$dendrograms)) {
  plotDendroAndColors(net$dendrograms[[i]],
                      colors = cbind(labels2colors(net$unmergedColors[net$blockGenes[[i]]]),
                                     labels2colors(net$colors[net$blockGenes[[i]]])),
                      groupLabels = c("Unmerged", "Merged"),
                      main = paste("Dendrogram - Block", i),
                      addGuide = TRUE)
}

# ── 9. Eigengene Clustering ────────────────────────────────────
moduleColors <- labels2colors(net$colors)
MEList <- moduleEigengenes(datExpr_clean, colors = moduleColors)
MEs <- MEList$eigengenes
MEDiss <- 1 - cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")
plot(METree, main = "Clustering of Module Eigengenes")

# ── 10. Module–Trait Correlation ───────────────────────────────
sample_names <- rownames(datExpr_clean)
traitData <- data.frame(
  row.names = sample_names,
  RIF = ifelse(grepl("^RIF", sample_names), 1, 0),
  Control = ifelse(grepl("^Con", sample_names), 1, 0)
)

commonSamples <- intersect(rownames(traitData), rownames(MEs))
traitData <- traitData[commonSamples, , drop = FALSE]
MEs <- MEs[commonSamples, , drop = FALSE]

moduleTraitCor <- cor(MEs, traitData, use = "pairwise.complete.obs")
moduleTraitPval <- corPvalueStudent(moduleTraitCor, nSamples = nrow(MEs))

# Text matrix for heatmap
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPval, 1), ")", sep = "")
par(mar = c(5, 9, 4, 2) + 0.1)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(traitData),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               zlim = c(-1, 1),
               main = "Module–Trait Relationships")

# ── 11. Extract Significant Modules and Genes ──────────────────
sigModules <- rownames(moduleTraitPval)[apply(moduleTraitPval, 1, function(p) any(p < 0.05))]
rif_corr <- moduleTraitCor[, "RIF"]
control_corr <- moduleTraitCor[, "Control"]
rif_p <- moduleTraitPval[, "RIF"]
control_p <- moduleTraitPval[, "Control"]

rif_upregulated <- rownames(moduleTraitCor)[rif_corr > control_corr & rif_p < 0.05]
control_upregulated <- rownames(moduleTraitCor)[control_corr > rif_corr & control_p < 0.05]

# Clean names
rif_upregulated <- gsub("^ME", "", rif_upregulated)
control_upregulated <- gsub("^ME", "", control_upregulated)

cat("RIF Upregulated Modules:\n"); print(rif_upregulated)
cat("Control Upregulated Modules:\n"); print(control_upregulated)

# ── 12. Save Gene Lists for Key Modules ────────────────────────
modules_of_interest <- c("brown", "darkgreen", "darkgrey", "darkred", "salmon", "tan", "turquoise")
output_dir <- "results/genes_by_module"
dir.create(output_dir, showWarnings = FALSE)

# Optional: join with DEG results
for (mod in modules_of_interest) {
  genes_in_mod <- colnames(datExpr_clean)[which(moduleColors == mod)]
  gene_df <- data.frame(Gene = genes_in_mod)
  
  # Merge with DEGs if available
  merged <- merge(gene_df, deg_df, by = "Gene", all.x = TRUE)
  write.csv(merged, file.path(output_dir, paste0("Genes_RIF_", mod, ".csv")), row.names = FALSE)
}

