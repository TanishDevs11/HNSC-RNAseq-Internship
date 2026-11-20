##############################################################
# 03_DESeq2_Differential_Expression.R
# DESeq2 -> lfcShrink -> Volcano -> Heatmap -> Reactome
#
# Starts directly from the ComBat-Seq corrected data produced by
# 01_Data_Harmonization.R and 02_Preprocessing_BatchCorrection.R.
# No harmonization, batch detection, or ComBat-Seq happens here -
# that work is already done; re-running it here would just waste
# time and risk producing a second, inconsistent copy of it.
##############################################################

setwd("D:/TANISH/COLLEGE/INTERNSHIP/HNSC/Tanish/checked files/Differential_Expression")

library(DESeq2)
library(apeglm)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(pheatmap)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)

##############################################################
# STEP 1: Load ComBat-Seq Corrected Counts and Metadata
##############################################################

counts <- read.csv(
  "Combined_Count_Matrix_ComBatSeq.csv",
  row.names = 1,
  check.names = FALSE
)

metadata <- read.csv(
  "Combined_Metadata.csv",
  stringsAsFactors = FALSE
)

rownames(metadata) <- metadata$Sample_ID

dim(counts)
dim(metadata)

stopifnot(
  ncol(counts) == nrow(metadata)
)

cat("Genes   :", nrow(counts), "\n")
cat("Samples :", ncol(counts), "\n")

head(counts[, 1:6])
head(metadata)

##############################################################
# STEP 2: Verify Sample Order
##############################################################
# Script 01/02 should have already produced counts and metadata in
# matching order - this is a lightweight verification, not a full
# re-harmonization.

stopifnot(all(rownames(metadata) %in% colnames(counts)))

counts <- counts[, rownames(metadata)]

stopifnot(all(colnames(counts) == rownames(metadata)))

cat("Sample order verified:", ncol(counts), "samples.\n")

metadata$Condition <- relevel(
  factor(metadata$Condition),
  ref = "Normal"
)

##############################################################
# STEP 3: Run DESeq2
##############################################################

dds_corrected <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts)),
  colData   = metadata,
  design    = ~ Condition
)

dds_corrected <- dds_corrected[rowSums(counts(dds_corrected)) >= 10, ]
dds_corrected <- DESeq(dds_corrected)

resultsNames(dds_corrected)

##############################################################
# STEP 4: Extract Shrunken Log2 Fold Changes
##############################################################
# lfcShrink gives more reliable log2FC estimates than raw results(),
# especially for low-count genes, by shrinking noisy fold changes
# toward zero using the apeglm method.

shrink_coef <- "Condition_Tumor_vs_Normal"

stopifnot(shrink_coef %in% resultsNames(dds_corrected))

res <- lfcShrink(
  dds_corrected,
  coef = shrink_coef,
  type = "apeglm"
)

res <- res[order(res$padj, na.last = TRUE), ]
deg_results <- as.data.frame(res)

head(deg_results)
summary(res)

##############################################################
# STEP 5: Normalized Counts + MA Plot
##############################################################

norm_counts <- counts(dds_corrected, normalized = TRUE)
write.csv(norm_counts, "DESeq2_Normalized_Counts.csv")

pdf("MA_Plot.pdf")
plotMA(res, main = "Tumor vs Normal (shrunken LFC)")
dev.off()

png(
  "MA_Plot.png",
  width  = 8,
  height = 6,
  units  = "in",
  res    = 600
)
plotMA(res, main = "Tumor vs Normal (shrunken LFC)")
dev.off()

##############################################################
# STEP 6: Save Full Results + Significant DEGs
##############################################################

deg_results$Gene <- rownames(deg_results)
deg_results <- deg_results[!is.na(deg_results$padj), ]

##############################################################
# STEP 6b: ENSEMBL -> Gene Symbol Mapping
##############################################################
# Done here, before the volcano plot and heatmap, because both of those
# steps need Gene_Symbol already present on deg_results. This is a
# separate lookup from the ENSEMBL -> ENTREZID mapping further down in
# the Reactome section, which only covers significant genes and is used
# for pathway enrichment, not for labeling.

# Ensembl IDs in this dataset have no version suffix (rownames look like
# "ENSG00000000003", not "ENSG00000000003.15"), but strip defensively in
# case that ever changes upstream.
all_gene_ids_stripped <- sub("\\..*$", "", deg_results$Gene)

symbol_mapping <- bitr(
  unique(all_gene_ids_stripped),
  fromType = "ENSEMBL",
  toType   = "SYMBOL",
  OrgDb    = org.Hs.eg.db
)
symbol_mapping <- symbol_mapping[!duplicated(symbol_mapping$ENSEMBL), ]

deg_results$Gene_Symbol <- symbol_mapping$SYMBOL[
  match(all_gene_ids_stripped, symbol_mapping$ENSEMBL)
]

# Fall back to the Ensembl ID itself for anything that didn't map - every
# row should still carry a usable label, never NA/blank.
unmapped <- is.na(deg_results$Gene_Symbol)
deg_results$Gene_Symbol[unmapped] <- deg_results$Gene[unmapped]

cat(
  "Genes mapped to a gene symbol:",
  sum(!unmapped), "of", nrow(deg_results),
  "\n"
)

write.csv(deg_results, "DESeq2_Results_All_Genes.csv")

sig_deg <- subset(
  deg_results,
  padj < 0.05 & abs(log2FoldChange) >= 2
)

write.csv(sig_deg, "DESeq2_Significant_DEGs.csv")

deg_results$Significance <- "Not Significant"

deg_results$Significance[
  deg_results$padj < 0.05 & deg_results$log2FoldChange >= 2
] <- "Upregulated"

deg_results$Significance[
  deg_results$padj < 0.05 & deg_results$log2FoldChange <= -2
] <- "Downregulated"

deg_results$Significance <- factor(
  deg_results$Significance,
  levels = c("Upregulated", "Downregulated", "Not Significant")
)

##############################################################
# STEP 7: Volcano Plot
##############################################################

volcano_plot <- ggplot(
  deg_results,
  aes(x = log2FoldChange, y = -log10(padj), color = Significance)
) +
  geom_point(size = 1.8, alpha = 0.8) +
  scale_color_manual(
    values = c(
      "Upregulated"     = "#D6604D",
      "Downregulated"   = "#4393C3",
      "Not Significant" = "grey70"
    )
  ) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_pubr(legend = "right") +
  labs(
    title = "Volcano Plot: Tumor vs Normal",
    x = "Log2 Fold Change (shrunken)",
    y = "-Log10 Adjusted P-value",
    color = "Significance"
  )

top_labeled_genes <- subset(
  deg_results,
  padj < 0.05 & abs(log2FoldChange) >= 2
)

top_labeled_genes <- top_labeled_genes[order(top_labeled_genes$padj), ][1:min(30, nrow(top_labeled_genes)), ]

volcano_plot_labeled <- volcano_plot +
  geom_text_repel(
    data = top_labeled_genes,
    aes(label = Gene_Symbol),
    size = 3,
    max.overlaps = 30,
    box.padding = 0.4,
    segment.color = "grey40",
    show.legend = FALSE
  )

ggsave(
  "Volcano_Plot.png",
  plot   = volcano_plot_labeled,
  width  = 9,
  height = 7,
  dpi    = 600
)

##############################################################
# STEP 8: Heatmap of Top 100 DEGs
##############################################################

vsd_corrected <- vst(dds_corrected, blind = FALSE)
vsd_mat <- assay(vsd_corrected)

sig_genes <- subset(
  deg_results,
  padj < 0.05 & abs(log2FoldChange) >= 2
)
sig_genes <- sig_genes[order(sig_genes$padj), ]

top100_genes <- head(rownames(sig_genes), 100)

if (length(top100_genes) == 0) {
  stop("No significant genes found for heatmap.")
}

heatmap_mat <- vsd_mat[top100_genes, ]

# Row names are Ensembl IDs at this point (same order as top100_genes) -
# swap in gene symbols for display, falling back to the Ensembl ID for
# anything that didn't map (same rule used everywhere else in this script).
rownames(heatmap_mat) <- sig_genes[top100_genes, "Gene_Symbol"]

# Quick sanity check on the matrix actually being plotted - catches an
# empty/malformed matrix before spending time on a device that might
# get stuck open, and is generally worth eyeballing before any heatmap.
cat("Heatmap matrix dimensions:", paste(dim(heatmap_mat), collapse = " x "), "\n")
cat("Heatmap value range      :", paste(round(range(heatmap_mat), 2), collapse = " to "), "\n")
cat("Significant DEGs (sig_deg):", nrow(sig_deg), "\n")

annotation_col <- data.frame(Condition = colData(dds_corrected)$Condition)
rownames(annotation_col) <- colnames(vsd_mat)

annotation_colors <- list(
  Condition = c(Normal = "#4393C3", Tumor = "#D6604D")
)

# Close any graphics devices left open from a previous run/session
# BEFORE opening a new one. This is what protects against the
# "+ prompt never finished, dev.off() closed the wrong device, PNG
# looks black" failure mode - a stray open device from earlier in the
# session is the most common cause of that, not a problem with the
# heatmap data itself.
graphics.off()

png(
  "Heatmap_Top100_DEGs.png",
  width  = 10,
  height = 14,
  units  = "in",
  res    = 600
)

# tryCatch(..., finally = dev.off()) guarantees the device closes even
# if pheatmap() errors out partway through - without this, an error
# here would leave the PNG device open and stuck, silently corrupting
# whatever plot runs next. (Note: on.exit() would NOT achieve this at
# the top level of a plain script - it only fires reliably inside a
# function call, so tryCatch/finally is the correct tool here.)
tryCatch({
  pheatmap(
    heatmap_mat,
    scale             = "row",
    cluster_rows      = TRUE,
    cluster_cols      = TRUE,
    show_rownames     = TRUE,    # row names are now gene symbols (see above),
    # readable at 100 rows unlike the raw Ensembl IDs
    show_colnames     = FALSE,   # ~792 sample columns - labels would overlap into
    # an unreadable solid block of text; Condition is already shown via annotation_col
    treeheight_col    = 0,       # with 792 leaves the column dendrogram itself renders
    # as a solid black band even with labels off - clustering still determines column
    # order, only the tree drawing is suppressed
    fontsize_row      = 6,
    annotation_col    = annotation_col,
    annotation_colors = annotation_colors,
    color             = colorRampPalette(c("#4393C3", "white", "#D6604D"))(100),
    main              = "Top 100 Differentially Expressed Genes"
  )
}, finally = {
  dev.off()
})

##############################################################
# STEP 9: Reactome Pathway Analysis
##############################################################

if (!exists("deg_results")) {
  stop("DESeq2 results not found.")
}

sig_genes <- subset(
  deg_results,
  padj < 0.05 & abs(log2FoldChange) >= 2
)

if (nrow(sig_genes) == 0) {
  stop("No significant genes found (padj < 0.05 & |log2FC| >= 2). Cannot run pathway enrichment.")
}

gene_ids <- rownames(sig_genes)
gene_ids <- sub("\\..*$", "", gene_ids)
gene_ids <- unique(gene_ids)

gene_mapping <- bitr(
  gene_ids,
  fromType = "ENSEMBL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

gene_mapping <- gene_mapping[!duplicated(gene_mapping$ENTREZID), ]

cat("Significant genes:", length(gene_ids), "\n")
cat("Successfully mapped to Entrez ID:", nrow(gene_mapping), "\n")

if (nrow(gene_mapping) == 0) {
  stop("No Ensembl IDs could be mapped to Entrez IDs. Check ID formatting/species.")
}

reactome_results <- enrichPathway(
  gene          = gene_mapping$ENTREZID,
  organism      = "human",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  readable      = TRUE
)

reactome_df <- if (is.null(reactome_results)) data.frame() else as.data.frame(reactome_results)

if (nrow(reactome_df) == 0) {
  
  message("No significant Reactome pathways found.")
  
} else {
  
  write.csv(reactome_df, "Reactome_Pathway_Results.csv")

  # clusterProfiler's built-in barplot() sorts by Count with no
  # straightforward way to sort by p.adjust instead, so this is a custom
  # ggplot2 bar chart built directly from reactome_df - most significant
  # pathway (lowest p.adjust) at the top.
  reactome_plot_data <- reactome_df[order(reactome_df$p.adjust), ]
  reactome_plot_data <- head(reactome_plot_data, 20)
  reactome_plot_data$Description <- factor(
    reactome_plot_data$Description,
    levels = rev(reactome_plot_data$Description)
  )

  # barplot() on an enrichResult returns a ggplot object (unlike pheatmap(),
  # which draws immediately via grid graphics) - it must be assigned and
  # passed to print() explicitly, or nothing renders to the PNG device even
  # though no error is raised. Same applies to a plain ggplot object here.
  # p.adjust spans an enormous range here (down to ~1e-125), so a linear
  # color scale would render almost every bar as visually identical and
  # would pack the legend's auto breaks into overlapping unreadable text.
  # log10-transforming the fill scale fixes both: it spaces out the color
  # gradient meaningfully across the whole range and gives the legend
  # sane, evenly-spaced breaks.
  reactome_barplot <- ggplot(
    reactome_plot_data,
    aes(x = Count, y = Description, fill = p.adjust)
  ) +
    geom_col() +
    scale_fill_gradient(
      low    = "#D6604D",   # lower p.adjust (more significant) - same warm
      high   = "#4393C3",   # color used for "Upregulated"/"Tumor" elsewhere
      trans  = "log10",
      breaks = scales::breaks_log(n = 5),
      labels = scales::label_scientific(digits = 0),
      name   = "p.adjust"
    ) +
    scale_y_discrete(
      labels = function(x) sapply(strwrap(x, width = 40, simplify = FALSE), paste, collapse = "\n")
    ) +
    theme_pubr(legend = "right") +
    theme(
      axis.text.y = element_text(size = 8),
      legend.text = element_text(size = 7),
      plot.margin = margin(t = 10, r = 15, b = 10, l = 10)
    ) +
    guides(fill = guide_colorbar(barwidth = unit(0.4, "cm"), barheight = unit(4, "cm"))) +
    labs(
      title = "Reactome Pathway Enrichment",
      x     = "Count",
      y     = NULL,
      fill  = "p.adjust"
    )

  # Close any stray open devices before opening this one - same rationale
  # as the heatmap step above.
  graphics.off()

  png(
    "Reactome_Barplot.png",
    width  = 11,
    height = 9,
    units  = "in",
    res    = 600
  )

  tryCatch({
    print(reactome_barplot)
  }, finally = {
    dev.off()
  })

  cat("Reactome pathway analysis completed.\n")
}

##############################################################
# STEP 10: Final Summary
##############################################################

cat("\n=====================================\n")
cat(" DESeq2 ANALYSIS COMPLETED\n")
cat("=====================================\n")
cat("Genes analysed    :", nrow(deg_results), "\n")
cat("Significant DEGs  :", nrow(sig_deg), "\n")
cat("Upregulated       :", sum(sig_deg$log2FoldChange > 0), "\n")
cat("Downregulated     :", sum(sig_deg$log2FoldChange < 0), "\n")
cat("FDR threshold     : padj < 0.05\n")
cat("Fold-change cutoff: |log2FC| >= 2\n")
cat("=====================================\n")

###################

#Toppgene preparation

####################

sum(grepl("^ENSG", sig_deg$Gene_Symbol))

sig_genes <- sig_deg

sig_genes <- sig_genes[
  order(sig_genes$padj),
]

top100_genes <- head(sig_genes, 100)

head(rownames(sig_genes), 100)

head(top100_genes[, c("X", "Gene_Symbol")])

sum(grepl("^ENSG", top100_genes$Gene_Symbol))

# Keep only valid gene symbols
toppgene_input <- unique(
  top100_genes$Gene_Symbol[
    !grepl("^ENSG", top100_genes$Gene_Symbol)
  ]
)

length(toppgene_input)

write.table(
  toppgene_input,
  "ToppGene_Input_Top100.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)