##############################################################
# install_packages.R
#
# Installs the R packages used by the HNSC RNA-seq pipeline in src/.
#
# IMPORTANT: package versions used during the internship were NOT
# recorded (no sessionInfo() was captured; see report sections 5.7
# and 7.6). This script installs current versions, which may differ
# from those originally used. After installing, record your own
# environment with:  writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
#
# Usage:  source("requirements/install_packages.R")
##############################################################

cran_packages <- c(
  "data.table",   # 05: header-indexed column subsetting of Xena matrix
  "httr",         # 02: GDC files API (UUID -> barcode)
  "jsonlite",     # 02: parse GDC API response
  "ggplot2",      # 06, 07: PCA, volcano, Reactome plots
  "ggpubr",       # 06, 07: theme_pubr
  "ggrepel",      # 07: volcano labels
  "pheatmap",     # 07: top-100 DEG heatmap
  "scales"        # 07: log-scale legend breaks/labels
)

bioc_packages <- c(
  "TCGAbiolinks",    # 01: GDC manifest handling
  "DESeq2",          # 06, 07: VST, PCA, differential expression
  "sva",             # 06: ComBat_seq
  "apeglm",          # 07: lfcShrink(type = "apeglm")
  "clusterProfiler", # 07: bitr ID conversion
  "ReactomePA",      # 07: enrichPathway
  "org.Hs.eg.db",    # 07, exploratory: human gene annotation
  "AnnotationDbi"    # exploratory: mapIds
)

missing_cran <- setdiff(cran_packages, rownames(installed.packages()))
if (length(missing_cran) > 0) {
  install.packages(missing_cran)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

missing_bioc <- setdiff(bioc_packages, rownames(installed.packages()))
if (length(missing_bioc) > 0) {
  BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)
}

cat("All pipeline packages are installed.\n")
