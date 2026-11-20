##############################################################
# Preprocessing, Batch Effect Assessment and ComBat-Seq Correction
##############################################################
#
# Objective:
# Perform quality assessment and batch correction for the
# integrated TCGA-HNSC and GTEx Esophagus Mucosa RNA-seq dataset
# prior to differential expression analysis.
#
# Workflow:
# 1. Load harmonized count matrix and metadata
# 2. Construct DESeq2 dataset
# 3. Filter lowly expressed genes
# 4. Perform Variance Stabilizing Transformation (VST)
# 5. Assess batch effects using PCA before correction
# 6. Apply ComBat-Seq batch correction
# 7. Assess batch correction using PCA (PC1-PC2 and PC3-PC4)
# 8. Investigate potential sample outliers
# 9. Save corrected count matrix and PCA figures
#
# Inputs:
# - Combined_Count_Matrix.csv
# - Combined_Metadata.csv
#
# Outputs:
# - PCA_Before_BatchCorrection.png
# - PCA_Before_BatchCorrection.pdf
# - Combined_Count_Matrix_ComBatSeq.csv
# - PCA_After_BatchCorrection.png
# - PCA_After_BatchCorrection.pdf
# - PCA_PC3_PC4_AfterBatchCorrection.png
# - PCA_PC3_PC4_AfterBatchCorrection.pdf
##############################################################


##############################################################
# Load Required Packages
##############################################################

library(DESeq2)
library(ggplot2)
library(ggpubr)

##############################################################
# Load Harmonized Dataset
##############################################################

counts <- read.csv(
  "Combined_Count_Matrix.csv",
  row.names = 1,
  check.names = FALSE
)

metadata <- read.csv(
  "Combined_Metadata.csv",
  stringsAsFactors = FALSE
)

##############################################################
# Prepare Metadata
##############################################################

rownames(metadata) <- metadata$Sample_ID

metadata$Sample_ID <- NULL

# Verify
head(metadata)

dim(metadata)

##############################################################
# Verify Sample Alignment
##############################################################

all(
  colnames(counts) ==
    rownames(metadata)
)

##############################################################
# Prepare Metadata Variables
##############################################################

metadata$Condition <- factor(
  metadata$Condition,
  levels = c("Normal","Tumor")
)

metadata$Batch <- factor(metadata$Batch)

##############################################################
# Create DESeq2 Dataset
##############################################################

dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts)),
  colData = metadata,
  design = ~1
)

##############################################################
# Filter Lowly Expressed Genes
##############################################################

dds <- dds[
  rowSums(counts(dds)) >= 10,
]

dim(dds)

##############################################################
# Variance Stabilizing Transformation
##############################################################

vsd <- vst(
  dds,
  blind = TRUE
)

##############################################################
# Principal Component Analysis
##############################################################

pca_data <- plotPCA(
  vsd,
  intgroup = c(
    "Batch",
    "Condition"
  ),
  returnData = TRUE
)

percentVar <- round(
  100 *
    attr(
      pca_data,
      "percentVar"
    )
)

##############################################################
# PCA Plot
##############################################################
pca_plot <- ggplot(
  pca_data,
  aes(
    PC1,
    PC2,
    color = Batch,
    shape = Condition
  )
) +
  geom_point(size = 3) +
  theme_pubr() +
  labs(
    title = "PCA Before Batch Correction",
    
    subtitle = paste(
      "Samples:", nrow(metadata),
      "| Genes:", nrow(counts)
    ),
    
    x = paste0(
      "PC1: ",
      percentVar[1],
      "% variance"
    ),
    
    y = paste0(
      "PC2: ",
      percentVar[2],
      "% variance"
    ),
    
    caption = "Variance Stabilizing Transformation (VST) before batch correction"
  )

pca_plot <- pca_plot +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 11
    ),
    plot.caption = element_text(
      hjust = 1,
      size = 9,
      face = "italic"
    )
  )

print(pca_plot)

##############################################################
# Save PCA Plot
##############################################################

ggsave(
  "PCA_Before_BatchCorrection.png",
  pca_plot,
  width = 8,
  height = 6,
  dpi = 600
)

ggsave(
  "PCA_Before_BatchCorrection.pdf",
  pca_plot,
  width = 8,
  height = 6
)

cat("\n=====================================\n")
cat(" PCA Before Batch Correction Complete\n")
cat("=====================================\n")


table(metadata$Batch, metadata$Condition)

##############################################################
# ComBat-Seq Batch Correction
##############################################################
#
# Objective:
# Correct batch effects between TCGA-HNSC and GTEx
# while preserving biological differences between
# Tumor and Normal samples.
##############################################################

##############################################################
# Load Required Package
##############################################################

library(sva)

##############################################################
# Convert Counts to Integers
##############################################################

raw_counts <- round(
  as.matrix(counts)
)

# Verify
mode(raw_counts)

dim(raw_counts)

all(raw_counts == round(raw_counts))

##############################################################
# Perform ComBat-Seq Batch Correction
##############################################################

combat_counts <- ComBat_seq(
  counts = raw_counts,
  batch = metadata$Batch,
  group = metadata$Condition
)

##############################################################
# Verify Corrected Count Matrix
##############################################################

dim(combat_counts)

head(combat_counts[,1:5])

##############################################################
# Save Corrected Count Matrix
##############################################################

write.csv(
  combat_counts,
  "Combined_Count_Matrix_ComBatSeq.csv",
  row.names = TRUE
)

cat("ComBat-Seq corrected count matrix saved successfully.\n")

#############################
#PCA After Batch Correction
###############################

##############################################################
# Create DESeq2 Object
##############################################################

dds_corrected <- DESeqDataSetFromMatrix(
  countData = combat_counts,
  colData = metadata,
  design = ~ Condition
)

##############################################################
# Variance Stabilizing Transformation
##############################################################

vsd_corrected <- vst(
  dds_corrected,
  blind = TRUE
)

##############################################################
# Perform PCA
##############################################################

pca_corrected <- plotPCA(
  vsd_corrected,
  intgroup = c("Batch", "Condition"),
  returnData = TRUE
)

percentVar <- round(
  100 * attr(
    pca_corrected,
    "percentVar"
  )
)

##############################################################
# PCA Visualization
##############################################################

pca_after_plot <- ggplot(
  pca_corrected,
  aes(
    PC1,
    PC2,
    color = Batch,
    shape = Condition
  )
) +
  geom_point(size = 3) +
  theme_pubr() +
  labs(
    title = "PCA After Batch Correction",
    
    subtitle = paste(
      "Samples:", nrow(metadata),
      "| Genes:", nrow(combat_counts)
    ),
    
    x = paste0(
      "PC1: ",
      percentVar[1],
      "% variance"
    ),
    
    y = paste0(
      "PC2: ",
      percentVar[2],
      "% variance"
    ),
    
    caption = "Variance Stabilizing Transformation after ComBat-Seq"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    ),
    
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 11
    ),
    
    plot.caption = element_text(
      hjust = 1,
      size = 9,
      face = "italic"
    )
  )

print(pca_after_plot)

##############################################################
# Save PCA Plot
##############################################################

ggsave(
  "PCA_After_BatchCorrection.png",
  pca_after_plot,
  width = 8,
  height = 6,
  dpi = 600
)

ggsave(
  "PCA_After_BatchCorrection.pdf",
  pca_after_plot,
  width = 8,
  height = 6
)

##############################################################
# Summary
##############################################################

cat("\n=====================================\n")
cat(" PCA After Batch Correction Complete\n")
cat("=====================================\n")


#####################################################
#validation analyses
#####################################################

##############################################################
# PCA Using Principal Components 3 and 4
##############################################################

# Extract VST expression matrix
vsd_matrix <- assay(vsd_corrected)

# Perform PCA
pca <- prcomp(
  t(vsd_matrix),
  center = TRUE,
  scale. = FALSE
)

# Percentage variance explained
percentVar <- round(
  100 * (pca$sdev^2 / sum(pca$sdev^2)),
  1
)

# Create dataframe
pca34_data <- data.frame(
  PC3 = pca$x[,3],
  PC4 = pca$x[,4],
  metadata
)

##############################################################
# PCA Plot (PC3 vs PC4)
##############################################################

pca34_plot <- ggplot(
  pca34_data,
  aes(
    PC3,
    PC4,
    color = Batch,
    shape = Condition
  )
) +
  geom_point(size = 3) +
  theme_pubr() +
  labs(
    title = "PCA After Batch Correction (PC3 vs PC4)",
    
    subtitle = paste(
      "Samples:", nrow(metadata),
      "| Genes:", nrow(combat_counts)
    ),
    
    x = paste0(
      "PC3: ",
      percentVar[3],
      "% variance"
    ),
    
    y = paste0(
      "PC4: ",
      percentVar[4],
      "% variance"
    ),
    
    caption = "Variance Stabilizing Transformation after ComBat-Seq"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 11
    ),
    plot.caption = element_text(
      hjust = 1,
      face = "italic",
      size = 9
    )
  )

print(pca34_plot)

ggsave(
  "PCA_PC3_PC4_AfterBatchCorrection.png",
  pca34_plot,
  width = 8,
  height = 6,
  dpi = 600
)

ggsave(
  "PCA_PC3_PC4_AfterBatchCorrection.pdf",
  pca34_plot,
  width = 8,
  height = 6
)

##############################################################
# Identify PCA Outlier
##############################################################

pca34_data[
  which.max(pca34_data$PC3),
]

##############################################################
# Library Size Check
##############################################################

library_sizes <- colSums(combat_counts)

summary(library_sizes)

library_sizes["GTEX-WFG8-1226-SM-5CHTQ"]

sort(
  library_sizes,
  decreasing = TRUE
)[1:10]

sort(
  library_sizes
)[1:10]

pca_corrected[
  rownames(pca_corrected) == "GTEX-WFG8-1226-SM-5CHTQ",
]

raw_library_sizes <- colSums(raw_counts)

raw_library_sizes[c(
  "GTEX-WFG8-1226-SM-5CHTQ",
  "GTEX-W5X1-1126-SM-5CHTH"
)]


cat("\n=============================================\n")
cat(" PREPROCESSING PIPELINE COMPLETED SUCCESSFULLY\n")
cat("=============================================\n")
cat("Genes retained :", nrow(combat_counts), "\n")
cat("Samples retained :", ncol(combat_counts), "\n")
cat("Tumor samples :", sum(metadata$Condition == "Tumor"), "\n")
cat("Normal samples :", sum(metadata$Condition == "Normal"), "\n")
cat("Batch correction :", "ComBat-Seq\n")
cat("PCA completed :", "Yes\n")
cat("=============================================\n")