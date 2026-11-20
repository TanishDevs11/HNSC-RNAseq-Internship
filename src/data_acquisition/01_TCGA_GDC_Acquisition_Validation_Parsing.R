######################################################################
# TCGA HNSC RNA-Seq Data Acquisition & Validation & Parsing Pipeline
######################################################################

# STEP 1: Load Required Library


library(TCGAbiolinks)



# STEP 2: Read the GDC Manifest File


manifest <- read.delim(
  "gdc_manifest.2026-06-18.143401.txt",
  stringsAsFactors = FALSE
)



# STEP 3: Inspect the Manifest


head(manifest)


# STEP 4: Convert Manifest into a TCGAbiolinks Object


query <- getManifest(manifest)



# STEP 5: Filter Only RNA-Seq STAR Count Files


rna_manifest <- manifest[
  grepl("rna_seq", manifest$filename) &
    grepl("counts|quantification", manifest$filename, ignore.case = TRUE),
]


# STEP 6: Extract File UUIDs


file_ids <- manifest$id

length(file_ids)
head(file_ids)



# STEP 7: Save RNA-Seq Only Manifest


rna_manifest <- manifest[
  grepl("rna_seq", manifest$filename) &
    grepl("counts|quantification", manifest$filename, ignore.case = TRUE),
]

write.table(
  rna_manifest,
  "rna_manifest.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# STEP 8: Download RNA-Seq Files

# Run this command in Command Prompt
#
# gdc-client.exe download -m rna_manifest.txt
#
# NOTE:
# If download is interrupted:
# Press Ctrl + C
# Re-run the same command to resume.
#
# If TEMP is on another drive (WinError 17):
#
# set TEMP=D:\...\temp
# set TMP=D:\...\temp
#
# Then rerun:
#
# gdc-client.exe download -m rna_manifest.txt


# STEP 9: Locate Downloaded STAR Count Files


files <- list.files(
  recursive = TRUE,
  pattern = "augmented_star_gene_counts.tsv$",
  full.names = TRUE
)

length(files)



# STEP 10: Read One STAR Count File (Testing)


test <- read.delim(
  files[1],
  skip = 1
)

dim(test)

head(test)



# STEP 11: Compare Downloaded Files with Manifest


expected <- rna_manifest$filename

downloaded <- basename(files)

missing <- setdiff(expected, downloaded)

missing



# STEP 12: Identify Corrupted or Unreadable Files

for(i in seq_along(files)){
  
  cat("Reading:", basename(files[i]), "\n")
  
  x <- try(
    read.delim(
      files[i],
      skip = 1
    ),
    silent = TRUE
  )
  
  if(inherits(x, "try-error")){
    cat("ERROR FOUND:\n")
    print(files[i])
    break
  }
}



# STEP 13: Refresh File List After Re-downloading Missing Files


files <- list.files(
  recursive = TRUE,
  pattern = "augmented_star_gene_counts.tsv$",
  full.names = TRUE
)

length(files)



# STEP 14: Check for Extra Files


extra <- setdiff(downloaded, expected)

length(extra)

extra



# STEP 15: Check for Duplicate Files


sum(duplicated(downloaded))



# STEP 16: Verify Every STAR File Opens Correctly


bad_files <- c()

for(f in files){
  
  x <- try(
    read.delim(
      f,
      skip = 1,
      nrows = 5
    ),
    silent = TRUE
  )
  
  if(inherits(x, "try-error")){
    bad_files <- c(bad_files, basename(f))
  }
  
}

length(bad_files)

bad_files


# STEP 17: Verify Number of Genes in Every File

gene_counts <- sapply(files, function(f){
  
  nrow(
    read.delim(
      f,
      skip = 1
    )
  )
  
})

summary(gene_counts)

# END OF DATA ACQUISITION & VALIDATION
#
# Expected Results:
#
# • Total STAR files = 521 #would change according to the dataset 
# • Missing files = 0
# • Extra files = 0
# • Duplicate files = 0
# • Bad files = 0
# • Gene count ≈ 60,660 in every file #subject to change 
#
# Dataset Status:
# READY FOR PARSING



#######################
#Parser
#######################

#To Read the first file

first <- read.table(
  files[1],
  sep = "\t",
  header = TRUE,
  skip = 1,
  stringsAsFactors = FALSE
)

dim(first)

#Save gene annotation

gene_annotation <- first[, c(
  "gene_id",
  "gene_name",
  "gene_type"
)]

head(gene_annotation)

#Creating sample names

sample_names <- sub(
  "\\.rna_seq\\.augmented_star_gene_counts\\.tsv$",
  "",
  basename(files)
)

head(sample_names)
length(sample_names)

#Creating empty count matrix

count_matrix <- matrix(
  
  nrow = nrow(first),
  
  ncol = length(files)
  
)

#check matrix 
dim(count_matrix)

#matrix column names

colnames(count_matrix) <- sample_names

head(colnames(count_matrix))

#fill the first column 

count_matrix[,1] <- first$unstranded

head(count_matrix[,1]) #check it 

#complete the matrix

for(i in seq_along(files)){
  
  cat("Reading sample", i, "of", length(files), "\n")
  
  temp <- read.table(
    files[i],
    sep = "\t",
    header = TRUE,
    skip = 1,
    stringsAsFactors = FALSE
  )
  
  count_matrix[, i] <- temp$unstranded
  
}

#to get counts 

dim(count_matrix) 

#verify matrix

head(count_matrix[,1])
head(count_matrix[,2])
head(count_matrix[,3])

#Cleaning Matrix 

#Removing STAR summary rows

count_matrix <- count_matrix[-c(1:4), ]

gene_annotation <- gene_annotation[-c(1:4), ]

#counting 

dim(count_matrix)

dim(gene_annotation)

#verifying alingment 

all(
  gene_annotation$gene_id ==
    first$gene_id[-c(1:4)]
)


#Gene IDs as row names

rownames(count_matrix) <- gene_annotation$gene_id

#to check 

head(rownames(count_matrix))

#Removing Ensemble version number 

gene_annotation$gene_id <- sub(
  "\\..*$",
  "",
  gene_annotation$gene_id
)

rownames(count_matrix) <- gene_annotation$gene_id

#to check 

head(rownames(count_matrix))


#Saving outputs 

write.csv(
  count_matrix,
  "hnsc_counts_ENSG.csv"
)

write.csv(
  gene_annotation,
  "gene_annotation.csv",
  row.names = FALSE
)


#final check for counts 

dim(count_matrix)

head(count_matrix[,1:5])

dim(gene_annotation)

head(gene_annotation)