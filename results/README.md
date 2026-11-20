# Results

These are the final result tables written by [`src/differential_expression/07_DESeq2_Differential_Expression_Analysis.R`](../src/differential_expression/07_DESeq2_Differential_Expression_Analysis.R), plus the combined sample sheet from script 05. They are copied unchanged from the analysis output directory. Every count below was re-checked against these files and matches the report.

| File | Rows | Description |
|---|---:|---|
| `DESeq2_Results_All_Genes.csv` | 52,568 | DESeq2 + apeglm results for every tested gene with a non-NA adjusted p-value |
| `DESeq2_Significant_DEGs.csv` | 20,618 | Subset with `padj < 0.05` and `abs(log2FoldChange) >= 2` (630 up, 19,988 down) |
| `Reactome_Pathway_Results.csv` | 127 | `ReactomePA::enrichPathway` output, BH-adjusted p < 0.05 |
| `ToppGene_Input_Top100.txt` | 98 | HUGO symbols of the 100 most significant DEGs (2 unresolved Ensembl IDs dropped), submitted to ToppGene on 14 July 2026 |
| `Combined_Metadata.csv` | 792 | Sample sheet for the combined cohort (482 TCGA tumour, 37 TCGA normal, 273 GTEx normal) |

## Column dictionary

**`DESeq2_Results_All_Genes.csv` / `DESeq2_Significant_DEGs.csv`**

| Column | Meaning |
|---|---|
| *(first, unnamed)* | Ensembl gene ID (row name) |
| `baseMean` | Mean of DESeq2-normalised counts across all 792 samples |
| `log2FoldChange` | **apeglm-shrunken** log2 fold change, Tumour vs Normal (positive = higher in tumour) |
| `lfcSE` | Posterior standard deviation of the shrunken log2FC |
| `pvalue` | Wald test p-value |
| `padj` | Benjamini–Hochberg adjusted p-value |
| `Gene` | Ensembl gene ID |
| `Gene_Symbol` | HUGO symbol from `org.Hs.eg.db`; falls back to the Ensembl ID when unmapped |

**`Reactome_Pathway_Results.csv`**: standard `clusterProfiler` `enrichResult` columns: `ID`, `Description`, `GeneRatio` (DEGs in pathway / mapped DEGs, e.g. `278/1934`), `BgRatio` (pathway size / universe, e.g. `387/11091`), `pvalue`, `p.adjust` (BH), `qvalue`, `geneID` (`/`-separated symbols), `Count`.

**`Combined_Metadata.csv`**: `Sample_ID` (15-character TCGA sample barcode or GTEx sample ID), `Condition` (Tumor/Normal), `Batch` (TCGA/GTEX), `Tissue` (HNSC / Adjacent_Normal / Esophagus_Mucosa), `Study`.

## Interpretation caveats

- Batch and condition are nearly collinear in this cohort. Treat `padj` as a **ranking**, not a calibrated error rate. See the [README limitations](../README.md#limitations).
- Very large negative fold changes (below about −10), especially at low `baseMean`, are more consistent with cross-cohort quantification differences than with biology.
- The olfactory-receptor pathways at the top of the Reactome table are interpreted in the report as a **cohort artefact**.
- Significant genes are **candidates**, not validated biomarkers.

## Not included (and why)

| Output | Reason |
|---|---|
| `Combined_Count_Matrix.csv`, `Combined_Count_Matrix_ComBatSeq.csv` | 118–143 MB each, above GitHub's file limit; regenerable from public data with scripts 05–06 |
| `DESeq2_Normalized_Counts.csv` | About 425 MB; regenerable with script 07 |
| Per-cohort harmonised count matrices | 63–86 MB each; intermediate |
| PNG/PDF plot files | The figures used in the report are in [`figures/report_figures/`](../figures/report_figures/) |

Data provenance: TCGA-HNSC (open access, NCI GDC) and GTEx (open access, via UCSC Xena). These tables contain gene-level summary statistics and public sample identifiers only. They contain no expression values per sample and no clinical data.
