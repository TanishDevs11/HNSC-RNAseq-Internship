# Software and packages

> [!NOTE]
> **Versions are not recorded.** As the report states (§5.7, §7.6), the scripts were run interactively and no `sessionInfo()` output was captured, so the R version and package versions used are unknown. The table below lists *what* was used and *why*. It does not list versions. If you re-run the pipeline, record `sessionInfo()` and commit it alongside your outputs.

## R and Bioconductor packages

| Package | Source | Script(s) | Role in this pipeline |
|---|---|:-:|---|
| **TCGAbiolinks** | Bioconductor | 01 | Reading and handling the GDC manifest |
| **httr** | CRAN | 02 | POST query to the GDC `files` endpoint for barcode resolution |
| **jsonlite** | CRAN | 02 | Parsing the GDC API JSON response |
| **data.table** | CRAN | 05 | Header-indexed column subsetting (`fread(select = …)`) of the Xena count matrix |
| **org.Hs.eg.db** | Bioconductor | 07, exploratory | Ensembl ↔ Entrez ↔ HUGO symbol mapping |
| **AnnotationDbi** | Bioconductor | exploratory | `mapIds` for Entrez → symbol |
| **DESeq2** | Bioconductor | 06, 07 | Dataset construction, low-count filtering, VST, `plotPCA`, size factors, dispersion, Wald test, `plotMA` |
| **sva** | Bioconductor | 06 | `ComBat_seq` negative-binomial batch adjustment with condition protected |
| **apeglm** | Bioconductor | 07 | Adaptive shrinkage of log2 fold changes via `lfcShrink` |
| **clusterProfiler** | Bioconductor | 07 | `bitr` identifier conversion; enrichment framework |
| **ReactomePA** | Bioconductor | 07 | `enrichPathway` Reactome over-representation analysis |
| **ggplot2** | CRAN | 06, 07 | PCA, volcano and Reactome plots |
| **ggpubr** | CRAN | 06, 07 | `theme_pubr` styling |
| **ggrepel** | CRAN | 07 | Non-overlapping volcano labels |
| **pheatmap** | CRAN | 07 | Clustered heatmap of the top 100 DEGs |
| **scales** | CRAN | 07 | Log-scale breaks and labels for the Reactome fill (called as `scales::`) |
| base R `stats` | R | 06 | `prcomp` for PC3/PC4 |

Install everything with [`requirements/install_packages.R`](../requirements/install_packages.R).

## External tools and resources

| Tool / resource | Use |
|---|---|
| **R** | Analysis environment (scripts run interactively in RStudio) |
| **GDC Data Transfer Tool** (`gdc-client`) | Manifest-driven download of the 521 STAR count files |
| **NCI Genomic Data Commons** API (`api.gdc.cancer.gov/files`) | File UUID → TCGA barcode resolution |
| **UCSC Xena** TCGA-TARGET-GTEx hub | GTEx (and TCGA) expected counts from the TOIL RNA-seq recompute; phenotype table |
| **GTEx Portal** (v11 annotations) | Tissue survey for surrogate selection |
| **Reactome** knowledgebase | Pathway annotation for over-representation analysis |
| **ToppGene Suite** (ToppFun) | Independent web-based enrichment of the top 100 DEGs (14 July 2026) |
| **Git** / **GitHub** | Version control and publication of this portfolio |

## Methods referenced (key citations from the report)

- Love, Huber & Anders (2014). DESeq2. *Genome Biology* 15:550.
- Zhu, Ibrahim & Love (2019). apeglm. *Bioinformatics* 35:2084–2092.
- Zhang, Parmigiani & Johnson (2020). ComBat-seq. *NAR Genomics and Bioinformatics* 2:lqaa078.
- Leek et al. (2012). The sva package. *Bioinformatics* 28:882–883.
- Yu & He (2016). ReactomePA. *Molecular BioSystems* 12:477–479.
- Wu et al. (2021). clusterProfiler 4.0. *The Innovation* 2:100141.
- Chen et al. (2009). ToppGene Suite. *Nucleic Acids Research* 37:W305–W311.
- Colaprico et al. (2016). TCGAbiolinks. *Nucleic Acids Research* 44:e71.
- Vivian et al. (2017). Toil. *Nature Biotechnology* 35:314–316.
- Goldman et al. (2020). UCSC Xena. *Nature Biotechnology* 38:675–678.
- Nygaard, Rødland & Hovig (2016). Batch correction and exaggerated confidence. *Biostatistics* 17:29–39.

The full reference list (69 entries) is in Section 9 of the [report](../report/Tanish_N_Internship_Report_Corrected_01_11.pdf).
