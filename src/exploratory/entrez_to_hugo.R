library(org.Hs.eg.db) #load library 

library(AnnotationDbi) #


#to check counts 

counts <- read.delim(
  "data_mrna_seq_read_counts.txt",
  check.names = FALSE
)


#convert entrez id to Hugo 

counts$HUGO <- mapIds(
  org.Hs.eg.db,
  keys = as.character(counts$Entrez_Gene_Id),
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

#get Hugo name 
head(counts[, c("Entrez_Gene_Id", "HUGO")])

#conversion of entrez to Hugo (Mapping)
#to check how many genes mapped 

sum(!is.na(counts$HUGO))


#how many genes didnt map 

sum(is.na(counts$HUGO))

#clean dataset 
counts_clean <- counts[!is.na(counts$HUGO), ]

#remove duplicate hugo symbols 

counts_clean <- counts_clean[!duplicated(counts_clean$HUGO), ]

#to move hugo id to the front column

counts_clean <- counts_clean[, c(
  "HUGO",
  "Entrez_Gene_Id",
  setdiff(colnames(counts_clean), c("HUGO", "Entrez_Gene_Id"))
)]

#preview 

head(counts_clean[, 1:5])

#save file
write.csv(
  counts_clean,
  "HNSC_counts_HUGO.csv",
  row.names = FALSE
)

