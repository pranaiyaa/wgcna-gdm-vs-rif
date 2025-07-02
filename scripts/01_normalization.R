# ───────────────────────────────────────────────
# Normalization of CEL files using RMA method
# ───────────────────────────────────────────────

library(oligo)
library(GEOquery)
library(tidyverse)
library(pd.hta.2.0)
library(AnnotationDbi)
library(hta20transcriptcluster.db) 

# Path to CEL files
cel_path <- "data/CEL_files/"
output_csv <- "data/normalized_expression.csv"

# Read all CEL files
cel_files <- list.celfiles(cel_path, full.names = TRUE)
raw_data <- read.celfiles(cel_files)

# Normalize using RMA
norm_data <- rma(raw_data)

# Convert to data frame
exprs_df <- as.data.frame(exprs(norm_data))
exprs_df$PROBEID <- rownames(exprs_df)

# Map probe IDs to gene symbols
gene_map <- AnnotationDbi::select(
  hta20transcriptcluster.db,
  keys = exprs_df$PROBEID,
  columns = c("SYMBOL"),
  keytype = "PROBEID"
)

# Merge and aggregate by gene symbol
merged_data <- exprs_df %>%
  left_join(gene_map, by = "PROBEID") %>%
  filter(!is.na(SYMBOL)) %>%
  select(-PROBEID) %>%
  group_by(SYMBOL) %>%
  summarise(across(everything(), mean), .groups = "drop") %>%
  rename(GeneSymbol = SYMBOL)

# Save normalized data
write.csv(merged_data, file = output_csv, row.names = FALSE)
