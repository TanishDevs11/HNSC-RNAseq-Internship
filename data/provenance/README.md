# Data provenance

These are small, open-access records that define exactly which TCGA files the analysis used. **No expression data and no clinical data are stored in this repository.**

| File | Rows | Produced by | Description |
|---|---:|---|---|
| `rna_manifest.txt` | 521 | script 01 | GDC manifest (`id`, `filename`, `md5`, `size`, `state`) of the TCGA-HNSC STAR – Counts files that were downloaded. Use it with `gdc-client download -m rna_manifest.txt` to re-fetch the same files. |
| `HNSC_metadata_521.csv` | 521 | GDC metadata, subset in script 02 | GDC **file-level** metadata for the same 521 files: file UUID, aliquot barcode (`cases`), sample type (482 Primary Tumor, 37 Solid Tissue Normal, 2 Metastatic), workflow type and version, GDC data release, md5 and size. Every row has `access = open`. Input to script 05. |

Source: NCI Genomic Data Commons, project TCGA-HNSC. Manifest generated 18 June 2026. These records fall under the GDC **open-access** tier.

The full portal manifest (`gdc_manifest.2026-06-18.143401.txt`, about 32,000 rows covering every TCGA-HNSC data type) is **not** included. Only its RNA-seq subset above was used.
