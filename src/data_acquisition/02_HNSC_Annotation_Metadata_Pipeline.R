###############################################################
# TCGA HNSC Gene Annotation, Barcode Conversion & Metadata Integration Pipeline
###############################################################
#
# Description:
# 1. Add HUGO gene symbols to the parsed count matrix
# 2. Convert UUIDs to TCGA sample barcodes (16-character)
# 3. Extract metadata from the TCGA GDC portal
# 4. Save metadata for downloaded samples
# 5. Convert UUIDs to full TCGA aliquot barcodes (28-character)
#
###############################################################

#############################
# Adding HUGO Gene Symbols
############################

# Read counts
counts <- read.csv(
  "hnsc_counts_ENSG.csv",
  check.names = FALSE
)

# Remove the first column (saved row names)
counts <- counts[, -1]

# Read annotation
gene_annotation <- read.csv(
  "gene_annotation.csv",
  stringsAsFactors = FALSE
)

# Add HUGO and Ensembl
counts <- cbind(
  HUGO = gene_annotation$gene_name,
  Ensembl = gene_annotation$gene_id,
  counts
)

dim(counts)

# Save
write.csv(
  counts,
  "hnsc_counts_HUGO_1.csv",
  row.names = FALSE
)

############################
#UUID TO TCGA BARCODES 
############################

rna_manifest <- read.delim(
  "rna_manifest.txt",
  stringsAsFactors = FALSE
)

dim(rna_manifest)

head(rna_manifest)

#extract UUID 
# Extract filename UUIDs from the RNA manifest
filename_uuid <- sub(
  "\\.rna_seq\\.augmented_star_gene_counts\\.tsv$",
  "",
  rna_manifest$filename
)

# Create the lookup table
lookup <- data.frame(
  filename_uuid = filename_uuid,
  file_id = rna_manifest$id,
  stringsAsFactors = FALSE
)

# Check it
head(lookup)

dim(lookup)

#compare with count matrix 

counts <- read.csv(
  "hnsc_counts_HUGO_1.csv",
  check.names = FALSE
)



uuid_cols <- colnames(counts)[3:ncol(counts)]

sum(rna_filename_uuid %in% uuid_cols)


################################################
#mapping uuid to 16 charac_TCGA BARCODES 
################################################

library(httr)
library(jsonlite)

body <- list(
  filters = list(
    op = "in",
    content = list(
      field = "file_id",
      value = lookup$file_id
    )
  ),
  fields = paste(
    c(
      "file_id",
      "file_name",
      "cases.submitter_id",
      "cases.samples.submitter_id"
    ),
    collapse = ","
  ),
  format = "JSON",
  size = 1000
)

#gdc API

response <- POST(
  "https://api.gdc.cancer.gov/files",
  body = body,
  encode = "json"
)

#response from api 

res <- fromJSON(
  content(response, "text", encoding = "UTF-8"),
  flatten = TRUE
)

#check if it works 

names(res)

str(res, max.level = 2)

cat(content(response, "text", encoding = "UTF-8"))

###extracting uuid 

tcga_barcode <- sapply(
  res$data$hits$cases,
  function(x) x$samples[[1]]$submitter_id
)

head(tcga_barcode)

length(tcga_barcode)


##extract api mapping 

mapping <- data.frame(
  file_id = res$data$hits$file_id,
  file_name = res$data$hits$file_name,
  TCGA_Barcode = tcga_barcode,
  stringsAsFactors = FALSE
)

head(mapping)

dim(mapping)

#recover the filename UUID 

mapping$filename_uuid <- sub(
  "\\.rna_seq\\.augmented_star_gene_counts\\.tsv$",
  "",
  mapping$file_name
)

head(mapping)

#merging with lookup table 

lookup2 <- merge(
  lookup,
  mapping,
  by = c("file_id", "filename_uuid")
)

dim(lookup2)

head(lookup2)

#renaming dictionary 
barcode_map <- setNames(
  lookup2$TCGA_Barcode,
  lookup2$filename_uuid
)

head(barcode_map)

#renaming count matrix 

old_names <- colnames(counts)[3:ncol(counts)]

new_names <- barcode_map[old_names]

sum(is.na(new_names))

#Rename

colnames(counts)[3:ncol(counts)] <- new_names

#to verify 

head(colnames(counts), 10)

#save the files 

write.csv(
  counts,
  "hnsc_counts_HUGO_2.csv",
  row.names = FALSE
)


##############################
#Metadata saver 
##############################

#save metadata 
write.csv(
  metadata,
  "HNSC_metadata.csv",
  row.names = FALSE
)

dim(metadata)

head(metadata)

#metadata only for 521 files 

metadata_521 <- metadata[
  metadata$file_name %in% rna_manifest$filename,
]

dim(metadata_521)

write.csv(
  metadata_521,
  "HNSC_metadata_521.csv",
  row.names = FALSE
)

#verify all downloaded files are present 

sum(metadata_521$file_name %in% rna_manifest$filename)

##############################################
#mapping uuid to 28 charac_TCGA BARCODES 
##############################################

full_barcode_map <- metadata[, c(
  "file_name",
  "cases"
)]

head(full_barcode_map)

#recover file name UUID 

full_barcode_map$filename_uuid <- sub(
  "\\.rna_seq\\.augmented_star_gene_counts\\.tsv$",
  "",
  full_barcode_map$file_name
)

head(full_barcode_map)

#keeping only 521 samples 
full_barcode_map <- full_barcode_map[
  full_barcode_map$filename_uuid %in% lookup$filename_uuid,
]

dim(full_barcode_map)

##lookup vector

barcode28 <- setNames(
  full_barcode_map$cases,
  full_barcode_map$filename_uuid
)

#read hugo count matrix 

counts28 <- read.csv(
  "hnsc_counts_HUGO_1.csv",
  check.names = FALSE
)

#renaming uuid 

old_names <- colnames(counts28)[3:ncol(counts28)]

new_names <- barcode28[old_names]

sum(is.na(new_names))

#replacing column names 
colnames(counts28)[3:ncol(counts28)] <- new_names

#verify 

head(colnames(counts28), 10)

#saving 
write.csv(
  counts28,
  "hnsc_counts_HUGO_3.csv",
  row.names = FALSE
)

