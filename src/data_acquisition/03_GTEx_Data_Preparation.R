#####################################################
# GTEx RNA-seq Data Preparation Pipeline
#> For TCGA-HNSC + GTEx Integration
#####################################################
#
#
######################
#Tissues Checker
######################

#Read file
gtex_meta <- read.delim(
  "GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt",
  stringsAsFactors = FALSE
)

dim(gtex_meta)  #check dimensions 

colnames(gtex_meta) #check columns 

head(gtex_meta) #view first rows 

sort(unique(gtex_meta$SMTS)) #listing broad range of tissues which are present 

sort(unique(gtex_meta$SMTSD)) #detailed tissue categories 

sort(table(gtex_meta$SMTSD), decreasing = TRUE) #sample per tissues 

########################
#note :
# Candidate healthy control tissues for TCGA-HNSC
#
# - Esophagus - Mucosa
# - Esophagus - Muscularis
# - Esophagus - Gastroesophageal Junction
# - Minor Salivary Gland

#########################

#selecting healthy control tissue 

selected_tissue <- "Esophagus - Mucosa" 

#extracting Metadata for Selected Tissue

gtex_selected_meta <- subset(
  gtex_meta,
  SMTSD == selected_tissue
)

# Check dimensions
dim(gtex_selected_meta)

# Preview
head(gtex_selected_meta)

# Number of selected samples
nrow(gtex_selected_meta)

# Extracting GTEx Sample IDs

gtex_sample_ids <- gtex_selected_meta$SAMPID

length(gtex_sample_ids)

head(gtex_sample_ids)

# Saving Selected Metadata

write.csv(
  gtex_selected_meta,
  "gtex_selected_metadata.csv",
  row.names = FALSE
)

# Inspect GTEx Gene Count Matrix

readLines(
  gzfile("GTEx_Analysis_2026-05-19_v11_RNASeQCv2.4.3_gene_reads.gct.gz"),
  n = 5
)


gtex_counts_preview <- read.delim(
  gzfile("GTEx_Analysis_2026-05-19_v11_RNASeQCv2.4.3_gene_reads.gct.gz"),
  skip = 2,
  nrows = 10,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

dim(gtex_counts_preview)

head(gtex_counts_preview[,1:6])

colnames(gtex_counts_preview)[1:10]

# Remove Ensembl Version Numbers

gtex_counts_preview$Name <- sub(
  "\\..*",
  "",
  gtex_counts_preview$Name
)

# Verify Ensembl IDs after removing version numbers
head(gtex_counts_preview$Name)