# Methodology

This summarises Section 5 (Materials and Methods) of the [internship report](../report/Tanish_N_Internship_Report_Corrected_01_11.pdf). Parameters are given exactly as applied in the [source code](../src/).

## 1. Data acquisition and cohort selection

### 1.1 TCGA-HNSC (tumour and adjacent normal)

- **Source:** NCI Genomic Data Commons (GDC). The manifest was generated on 18 June 2026 (`gdc_manifest.2026-06-18.143401.txt`).
- **Filter:** entries whose filename matches `rna_seq` and a `counts|quantification` pattern, written to `rna_manifest.txt`. The download used `gdc-client download -m rna_manifest.txt`.
- **Retrieved set:** 521 open-access files. All are TCGA-HNSC, RNA-Seq (Illumina), Transcriptome Profiling / Gene Expression Quantification, and produced by the GDC **STAR – Counts** workflow.
- **Provenance:** reference hg38, gene model **GENCODE v36**. Two workflow versions are present: `5d8c131b` (506 files; GDC releases 32.0–45.0) and `5a02c129` (15 files; releases 36.0–45.0). This is a possible within-TCGA source of technical variation.
- **Validation:** 521 present; 0 missing, 0 extra, 0 duplicated, 0 unreadable; 60,660 genes in every file.
- **Parsing:** the `unstranded` raw counts were used. TPM and FPKM are not appropriate inputs for count-based models. The 4 STAR summary rows were removed and Ensembl version suffixes stripped.
- **Identifier resolution:** file UUIDs were resolved to TCGA barcodes through a POST to `https://api.gdc.cancer.gov/files` (`httr`, `jsonlite`; fields `file_id`, `file_name`, `cases.submitter_id`, `cases.samples.submitter_id`; size 1000). 28-character aliquot barcodes came from the retained GDC metadata.
- **Cohort selection:** 482 Primary Tumor, 37 Solid Tissue Normal, 2 Metastatic. The metastatic samples were excluded, leaving **519**. Barcodes were truncated to 15 characters (sample level) and uniqueness was confirmed.

### 1.2 GTEx (surrogate normal)

- **Tissue choice:** GTEx samples no oral or oropharyngeal mucosa. Four candidates were considered: Esophagus – Mucosa, Esophagus – Muscularis, Esophagus – Gastroesophageal Junction and Minor Salivary Gland. **Esophagus – Mucosa** was selected because it is the only stratified squamous epithelium among them and has precedent in prior TCGA-HNSC work (Rauf et al., 2024).
- **Data route:** the GTEx v11 RNASeQC release (1,366 Esophagus – Mucosa samples) was **rejected** because its quantification pipeline differs from TCGA's. The **UCSC Xena TCGA-TARGET-GTEx** hub (TOIL RNA-seq recompute, gene-level expected counts, shared Ensembl IDs) was used instead. This traded sample size for quantification consistency.
- **Selection:** `_study == "GTEX"` and `primary disease or tissue == "Esophagus - Mucosa"`, giving **273 samples**.
- **Extraction:** the header was read alone, columns were indexed by set membership, and only those columns were read with `data.table::fread(select = …)`. Missing IDs were enumerated before extraction.

## 2. Metadata harmonisation and count-matrix assembly

- **Five-field metadata:** `Sample_ID`, `Condition` (Tumor/Normal), `Batch` (TCGA/GTEX), `Tissue` (HNSC / Adjacent_Normal / Esophagus_Mucosa), `Study`.
- **Gene identifiers:** versions were stripped from the first period onward. **44** duplicate TCGA Ensembl IDs were removed (first occurrence kept), and the absence of duplicates was asserted.
- **Intersection:** TCGA 60,660 ∩ GTEx 60,498 = **56,600 genes**.
- **Combination:** expression blocks were re-extracted immediately before `cbind`, and their dimensions were asserted against the metadata. Final alignment was asserted, giving **56,600 genes × 792 samples**.
- **Safeguards:** every join is guarded by `stopifnot()`. Two development failures motivated these checks: a variable-name typo that silently skipped de-duplication, and a stale 482-sample object that produced a 755-column matrix.

## 3. Quality control and exploratory analysis

- DESeq2 dataset with `design = ~1` (no testing at this stage). Expected counts were rounded to integers.
- `Condition` factor levels `Normal`, `Tumor`; `Batch` as a factor.
- **Low-count filter:** total count ≥ 10 across all 792 samples, retaining **52,568** genes. The filter is deliberately permissive so that genes expressed in only one condition are kept.
- **VST** with `blind = TRUE`; `plotPCA(intgroup = c("Batch", "Condition"), returnData = TRUE)`.
- The batch × condition cross-tabulation was computed at this stage, quantifying the confound.

## 4. Batch-effect correction

- `sva::ComBat_seq(counts = raw_counts, batch = Batch, group = Condition)` on the **raw combined counts**, not on transformed or filtered counts.
- Condition is supplied as `group` so that biological variation is protected. No other covariates were used. SVA was not used because the batch is known, not latent.
- **Assessment:** a new DESeq2 object and VST (`blind = TRUE`) → PC1/PC2 with the same encoding → full PCA with `prcomp(center = TRUE, scale. = FALSE)` → PC3/PC4.
- **Diagnostics:** the max-PC3 sample was named (`GTEX-WFG8-1226-SM-5CHTQ`), and raw and corrected library sizes were inspected.
- **Caveat:** the post-correction decomposition used all 56,600 genes, because the low-count filter was not re-applied. The before/after variance comparison is therefore approximate.

## 5. Differential expression

- Input: the ComBat-Seq matrix and combined metadata, re-read from disk with sample order re-verified.
- `design = ~ Condition`, reference level `Normal`; coefficient `Condition_Tumor_vs_Normal`. Batch was not modelled: it was already corrected, and it is nearly collinear with condition.
- Filter: row sum ≥ 10, leaving **52,568** genes.
- `DESeq()` defaults: median-of-ratios size factors, empirical-Bayes dispersion, negative-binomial **Wald** test.
- `lfcShrink(type = "apeglm")`, which places a heavy-tailed (Cauchy) prior on effect sizes.
- **Significance:** BH-adjusted **padj < 0.05** **and** **|shrunken log2FC| ≥ 2**. Rows with NA padj were dropped first.
- **Annotation:** Ensembl → SYMBOL via `clusterProfiler::bitr` / `org.Hs.eg.db`. Unmapped genes keep their Ensembl ID as the label.

## 6. Pathway enrichment (two tools, two scopes)

| | Reactome (full list) | ToppGene (high-confidence subset) |
|---|---|---|
| Input | All 20,618 significant DEGs → Ensembl → Entrez (1,934 Reactome-annotated) | Top 100 DEGs by padj → **98** HUGO symbols |
| Tool | `ReactomePA::enrichPathway` (clusterProfiler framework) | ToppGene Suite, ToppFun web interface, 14 July 2026 |
| Settings | organism human, `pvalueCutoff = 0.05`, `pAdjustMethod = "BH"`, `readable = TRUE`, default universe (11,091 Reactome-annotated genes) | GO MF/BP/CC, pathways, disease (DisGeNET), drug (CTD); BH FDR < 0.05 |
| Result | **127** pathways | Leading terms: ECM structural constituent, extracellular matrix, PID Integrin1, Reactome ECM organisation |

Running the same signal through two independent tools at two gene-list scopes tests whether an enriched term is **stable under a change of scope**. Real biology should strengthen when the list is restricted to the most confident genes, while artefacts near the significance boundary should disappear.

## 7. Visualisation

- PCA: `ggplot2` + `ggpubr::theme_pubr`, 600 dpi.
- MA plot: `DESeq2::plotMA` on the shrunken results.
- Volcano: `ggplot2`; dashed guides at ±2 and padj 0.05; top 30 labelled with `ggrepel`.
- Heatmap: `pheatmap` over the VST corrected matrix, top 100 DEGs, row-scaled, rows and columns clustered, column dendrogram hidden (792 samples).
- Reactome: custom `ggplot2` bar chart from the retained table, top 20 by `p.adjust`, log10 fill scale.
- Defensive practices: `graphics.off()` before opening devices, and `tryCatch(..., finally = dev.off())` around drawing.

## Thresholds at a glance (report Table 4)

| Stage | Parameter | Value |
|---|---|---|
| Cohort | TCGA sample types kept | Primary Tumor; Solid Tissue Normal |
| Cohort | Excluded | Metastatic (2) |
| Cohort | Barcode truncation | 28 → 15 characters |
| Cohort | GTEx tissue | Esophagus – Mucosa |
| Harmonisation | Ensembl version | Stripped from first period |
| Harmonisation | Duplicate IDs | First kept (44 removed, TCGA) |
| Harmonisation | Gene set | Intersection (56,600) |
| Preprocessing | Low-count filter | Row sum ≥ 10 across 792 samples |
| Preprocessing | Transformation | VST, `blind = TRUE` |
| Batch correction | Method / batch / group | ComBat_seq on raw counts / study / condition |
| DE | Design / reference | `~ Condition` / Normal |
| DE | Test / shrinkage | Wald / apeglm |
| DE | Significance | padj < 0.05 (BH) and \|log2FC\| ≥ 2 |
| Reactome | Cutoff / adjustment | 0.05 / BH |
| ToppGene | Input / FDR | Top 100 by padj (98 symbols) / 0.05 BH |
| Figures | Resolution | 600 dpi |
