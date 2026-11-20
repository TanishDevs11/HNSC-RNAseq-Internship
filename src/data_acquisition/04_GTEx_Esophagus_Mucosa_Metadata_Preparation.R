#######################################
# Load GTEx Phenotype Data
#######################################

# Read phenotype file
phenotype <- read.delim(
  "TcgaTargetGTEX_phenotype.txt.gz",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Check dimensions
dim(phenotype)

# View first few rows
head(phenotype)

# View column names
colnames(phenotype)

# Check available studies
unique(phenotype$`_study`)

# Extract GTEx samples
gtex <- subset(
  phenotype,
  `_study` == "GTEX"
)

# Number of GTEx samples
dim(gtex)

# Explore available GTEx tissue types
table(gtex$`primary disease or tissue`)

#######################################
# Load Gene Expression Count Matrix
#######################################

counts_preview <- read.delim(
  "TcgaTargetGtex_gene_expected_count.gz",
  nrows = 5,
  check.names = FALSE
)

# Verify dimensions
dim(counts_preview)

# Inspect first few genes and samples
head(counts_preview[, 1:5])

# Select GTEx Esophagus Mucosa Samples

gtex_esophagus <- subset(
  gtex,
  `primary disease or tissue` == "Esophagus - Mucosa"
)

# Number of selected samples
dim(gtex_esophagus)

# Verify tissue type
table(gtex_esophagus$`primary disease or tissue`)

# View first few sample IDs
head(gtex_esophagus$sample)

# Match Sample IDs

sum(
  gtex_esophagus$sample %in%
    colnames(counts_preview)
)

#######################################
# Extract GTEx Esophagus Mucosa Samples
#######################################

gtex_esophagus <- subset(
  gtex,
  `primary disease or tissue` == "Esophagus - Mucosa"
)

gtex_ids <- gtex_esophagus$sample

length(gtex_ids)

#Read only header 

header <- colnames(
  read.delim(
    "TcgaTargetGtex_gene_expected_count.gz",
    nrows = 1,
    check.names = FALSE
  )
)

#Find the required columns 

gtex_cols <- which(header %in% gtex_ids)

length(gtex_cols)

#Read only 

gtex_counts <- read.delim(
  "TcgaTargetGtex_gene_expected_count.gz",
  check.names = FALSE,
  colClasses = c(
    "character",
    rep("NULL", min(gtex_cols) - 2),
    rep("numeric", length(gtex_cols))
  )
)

#######################################
# Save GTEx Esophagus Mucosa Metadata
#######################################

write.csv(
  gtex_esophagus,
  "GTEx_Esophagus_Mucosa_Metadata.csv",
  row.names = FALSE
)

cat("GTEx Esophagus Mucosa metadata saved successfully.\n")
