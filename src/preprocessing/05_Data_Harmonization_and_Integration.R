##############################################################
# Data Harmonization of TCGA-HNSC and GTEx Esophagus Mucosa
##############################################################
# Harmonize metadata and raw RNA-seq count matrices for
# TCGA Head and Neck Squamous Cell Carcinoma (HNSC) and
# GTEx Esophagus Mucosa normal tissue samples.
#
# Outputs:
# - TCGA_Harmonized_Metadata.csv
# - GTEx_Harmonized_Metadata.csv
# - TCGA_Counts_Harmonized.csv
# - GTEx_Esophagus_Mucosa_Counts_Harmonized.csv
# - HNSC_Common_Genes_Counts.csv
# - GTEx_Common_Genes_Counts.csv
# - Combined_Count_Matrix.csv
# - Combined_Metadata.csv
#
##############################################################
# Load HNSC Metadata
##############################################################

metadata <- read.csv(
  "HNSC_metadata_521.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

dim(metadata)
colnames(metadata)
head(metadata)

##############################################################
# Extract TCGA Samples
# (Primary Tumor + Solid Tissue Normal)
##############################################################

hnsc_tcga <- subset(
  metadata,
  sample_type %in% c(
    "Primary Tumor",
    "Solid Tissue Normal"
  )
)

dim(hnsc_tcga)
table(hnsc_tcga$sample_type)

##############################################################
# Generate TCGA Sample IDs
##############################################################

hnsc_tcga$Sample_ID <- substr(
  hnsc_tcga$cases,
  1,
  15
)

head(hnsc_tcga$Sample_ID)
length(unique(hnsc_tcga$Sample_ID))

##############################################################
# Load GTEx Esophagus Mucosa Metadata
##############################################################

gtex_metadata <- read.csv(
  "GTEx_Esophagus_Mucosa_Metadata.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

dim(gtex_metadata)
colnames(gtex_metadata)
head(gtex_metadata)

##############################################################
# Create Harmonized TCGA Metadata
##############################################################

hnsc_harmonized <- data.frame(
  Sample_ID = hnsc_tcga$Sample_ID,
  
  Condition = ifelse(
    hnsc_tcga$sample_type == "Primary Tumor",
    "Tumor",
    "Normal"
  ),
  
  Batch = "TCGA",
  
  Tissue = ifelse(
    hnsc_tcga$sample_type == "Primary Tumor",
    "HNSC",
    "Adjacent_Normal"
  ),
  
  Study = "TCGA",
  
  stringsAsFactors = FALSE
)

dim(hnsc_harmonized)
table(hnsc_harmonized$Condition)
head(hnsc_harmonized)

##############################################################
# Create Harmonized GTEx Metadata
##############################################################

gtex_harmonized <- data.frame(
  Sample_ID = gtex_metadata$sample,
  Condition = "Normal",
  Batch = "GTEX",
  Tissue = "Esophagus_Mucosa",
  Study = "GTEX",
  stringsAsFactors = FALSE
)

dim(gtex_harmonized)
head(gtex_harmonized)

##############################################################
# Save Harmonized Metadata
##############################################################

write.csv(hnsc_harmonized, "TCGA_Harmonized_Metadata.csv", row.names = FALSE)
write.csv(gtex_harmonized, "GTEx_Harmonized_Metadata.csv", row.names = FALSE)

cat("Harmonized metadata files saved successfully.\n")

##############################################################
# Load HNSC Raw Count Matrix
##############################################################

counts <- read.csv(
  "hnsc_counts_HUGO_28_charac_TCGA_Barcodes.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

dim(counts)
head(counts[, 1:6])
ncol(counts) - 2

##############################################################
# Standardize TCGA Sample IDs
##############################################################

tcga_sample_ids <- colnames(counts)[-c(1, 2)]
tcga_sample_ids <- substr(tcga_sample_ids, 1, 15)
colnames(counts)[-c(1, 2)] <- tcga_sample_ids

head(colnames(counts)[3:10])

##############################################################
# Verify Metadata and Count Matrix Alignment
##############################################################

sum(hnsc_harmonized$Sample_ID %in% colnames(counts))
length(hnsc_harmonized$Sample_ID)
setdiff(hnsc_harmonized$Sample_ID, colnames(counts))

##############################################################
# Harmonize HNSC Count Matrix
##############################################################

counts_harmonized <- counts[
  ,
  c(
    1,
    2,
    match(hnsc_harmonized$Sample_ID, colnames(counts))
  )
]

dim(counts_harmonized)

##############################################################
# Verify Sample Order
##############################################################

stopifnot(
  all(hnsc_harmonized$Sample_ID == colnames(counts_harmonized)[-c(1, 2)])
)

##############################################################
# Save Harmonized HNSC Count Matrix
##############################################################

write.csv(counts_harmonized, "TCGA_Counts_Harmonized.csv", row.names = FALSE)

cat("TCGA count matrix harmonized successfully.\n")

##############################################################
# GTEx side Harmonization
##############################################################

library(data.table)

header <- names(
  fread("TcgaTargetGtex_gene_expected_count", nrows = 0)
)

length(header)
head(header)

##############################################################
# Identify GTEx Esophagus Mucosa Columns
##############################################################

gtex_columns <- which(header %in% gtex_harmonized$Sample_ID)
length(gtex_columns)

##############################################################
# Check Missing GTEx Samples
##############################################################

missing_samples <- setdiff(gtex_harmonized$Sample_ID, header)
length(missing_samples)
missing_samples

##############################################################
# Extract GTEx Esophagus Mucosa Count Matrix
##############################################################

cols_to_read <- c(1, gtex_columns)

gtex_counts <- fread(
  "TcgaTargetGtex_gene_expected_count",
  select = cols_to_read
)

dim(gtex_counts)
head(gtex_counts[, 1:6])

##############################################################
# Harmonize GTEx Count Matrix
##############################################################

gtex_counts_harmonized <- gtex_counts[
  ,
  c("sample", gtex_harmonized$Sample_ID),
  with = FALSE
]

dim(gtex_counts_harmonized)

stopifnot(
  all(gtex_harmonized$Sample_ID == colnames(gtex_counts_harmonized)[-1])
)

##############################################################
# Save Harmonized GTEx Count Matrix
##############################################################

write.csv(
  gtex_counts_harmonized,
  "GTEx_Esophagus_Mucosa_Counts_Harmonized.csv",
  row.names = FALSE
)

cat("GTEx raw count matrix harmonized successfully.\n")

##############################################################
# Gene Harmonization
##############################################################

# ---- Remove Ensembl version numbers - TCGA ----
counts_harmonized$Ensembl <- sub("\\..*", "", counts_harmonized$Ensembl)
head(counts_harmonized$Ensembl)

# ---- Remove Ensembl version numbers - GTEx ----
gtex_counts_harmonized$sample <- sub("\\..*", "", gtex_counts_harmonized$sample)
head(gtex_counts_harmonized$sample)

##############################################################
# Remove Duplicate Ensembl IDs (TCGA side)
##############################################################

# This step previously failed silently due to a variable-name typo
# ('ounts_harmonized' instead of 'counts_harmonized'), which meant
# duplicates were never actually dropped even though the code below it
# looked like it had run successfully. stopifnot() now makes that
# failure impossible to miss - the script halts instead of silently
# carrying 44 duplicate Ensembl IDs into every downstream step.
counts_harmonized <- counts_harmonized[
  !duplicated(counts_harmonized$Ensembl),
]

stopifnot(sum(duplicated(counts_harmonized$Ensembl)) == 0)

dim(counts_harmonized)

##############################################################
# Identify Common Ensembl Genes
##############################################################

common_genes <- intersect(
  counts_harmonized$Ensembl,
  gtex_counts_harmonized$sample
)

length(common_genes)

stopifnot(length(common_genes) > 0)

##############################################################
# Filter Common Genes
##############################################################

counts_common <- counts_harmonized[
  counts_harmonized$Ensembl %in% common_genes,
]

gtex_counts_common <- gtex_counts_harmonized[
  gtex_counts_harmonized$sample %in% common_genes,
]

dim(counts_common)
dim(gtex_counts_common)

##############################################################
# Reorder Genes
##############################################################

counts_common <- counts_common[
  match(common_genes, counts_common$Ensembl),
]

gtex_counts_common <- gtex_counts_common[
  match(common_genes, gtex_counts_common$sample),
]

stopifnot(
  all(counts_common$Ensembl == gtex_counts_common$sample)
)

##############################################################
# Save Common Gene Count Matrices
##############################################################

write.csv(counts_common, "HNSC_Common_Genes_Counts.csv", row.names = FALSE)
write.csv(gtex_counts_common, "GTEx_Common_Genes_Counts.csv", row.names = FALSE)

cat("Common gene count matrices saved successfully.\n")

##############################################################
# Build Final Combined Count Matrix
##############################################################

# IMPORTANT: tcga_expression and gtex_expression are (re)created here,
# immediately before use, from the CURRENT counts_common /
# gtex_counts_common objects. Extracting them earlier in the script (or
# reusing them from a previous session/run) is exactly what caused the
# 755-vs-793-column mismatch before: stale objects from an old 482-sample
# run silently persisted in the R environment and got cbind'ed in here
# instead of the current 519-sample data.
tcga_expression <- counts_common[, -(1:2)]
gtex_expression <- gtex_counts_common[, -1]

# Sanity check: expression matrix sample counts must match the metadata
# they're about to be paired with, BEFORE building combined_counts -
# catches a stale/mismatched object immediately instead of after the
# fact via a failed alignment check further down.
stopifnot(
  ncol(tcga_expression) == nrow(hnsc_harmonized),
  ncol(gtex_expression) == nrow(gtex_harmonized)
)

dim(tcga_expression)
dim(gtex_expression)

##############################################################
# Combine TCGA and GTEx Count Matrices
##############################################################

combined_counts <- cbind(
  Ensembl = counts_common$Ensembl,
  tcga_expression,
  gtex_expression
)

dim(combined_counts)

##############################################################
# Build Combined Metadata
##############################################################

combined_metadata <- rbind(
  hnsc_harmonized,
  gtex_harmonized
)

dim(combined_metadata)
head(combined_metadata)

##############################################################
# Verify Final Alignment
##############################################################

stopifnot(
  nrow(combined_metadata) == (ncol(combined_counts) - 1),
  all(combined_metadata$Sample_ID == colnames(combined_counts)[-1])
)

##############################################################
# Save Final Combined Dataset
##############################################################

write.csv(combined_counts, "Combined_Count_Matrix.csv", row.names = FALSE)
write.csv(combined_metadata, "Combined_Metadata.csv", row.names = FALSE)

cat("\n=============================================\n")
cat(" Data Harmonization Completed Successfully\n")
cat("=============================================\n")
cat("Genes   :", nrow(combined_counts), "\n")
cat("Samples :", ncol(combined_counts) - 1, "\n")
cat("Tumor   :", sum(combined_metadata$Condition == "Tumor"), "\n")
cat("Normal  :", sum(combined_metadata$Condition == "Normal"), "\n")
cat("=============================================\n")