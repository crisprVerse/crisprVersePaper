library(SummarizedExperiment)
library(stringr)

# Annotation
ann <- read.csv("../data/Cas13d_GFP_library.final.fa", head=FALSE)
ann <- ann[,1]
ann <- matrix(ann, ncol=2, byrow=TRUE)
colnames(ann) <- c("ID", "spacer")
ann <- as.data.frame(ann)
ann$ID <- gsub(">","", ann$ID)
ann$class <- str_extract(ann$ID,
                         "randomDouble|consecTriple|consecDouble|FirstOrder")
ann$class[is.na(ann$class)] <- "perfectMatch"
ann$class[ann$class=="FirstOrder"] <- "singleMismatch"
ann$class[grepl("rc_",ann$ID)]<- "control"


# Let's get the crrna ID
xs <- strsplit(ann$ID, "_")
ids <- sapply(xs, function(x) x[[1]])
ann$spacer_id <- ids
PM <- ann[ann$class=="perfectMatch",]
wh <- match(ann$spacer_id, PM$spacer_id)
ann$protospacer <- PM[wh,"spacer"]

dist <- sapply(1:nrow(ann), function(i){
    adist(ann$spacer[i],
          ann$protospacer[i],
          cost=c(del=1000, ins=1000, sub=1)
          )
})
ann$n_mismatches <- dist
ann <- ann[,c("ID", "spacer", "protospacer", "n_mismatches")]
ann$ID_ontarget <- ann$ID[match(ann$protospacer, ann$spacer)]
rownames(ann) <- ann$ID


# Let's get the data:
counts <- read.csv("../data/GFP_count_matrix_raw.csv")
counts <- as.matrix(counts)


# Let's do the normalization now before removing
mode="median"
mode="mean"
mode="ntc"
# Mdian normalization:
if (mode=="median"){
  counts <- log2(counts+1)
  factors <- colMedians(counts)
  factors <- factors - median(factors)
  counts <- sweep(counts, 2,factors, "-")
  counts <- 2^counts-1
} else if (mode=="mean"){
  # Mean normalization:
  factors <- colSums(counts)
  factors <- factors/median(factors)
  counts <- sweep(counts, 2,factors, "/")
} else if (mode=="ntc"){
  counts <- log2(counts+1)

  ntc <- grepl("rc_", rownames(counts))
  factors <- colMedians(counts[ntc,])
  factors <- factors - median(factors)
  counts <- sweep(counts, 2,factors, "-")
  counts <- 2^counts-1
}



# Subsetting:
counts <- counts[rownames(counts) %in% ann$ID,]
counts <- counts[ann$ID,]



# Let's get pheno:
pheno <- data.frame(sample=colnames(counts))
pheno$Replicate <- str_extract(pheno$sample, "SC[1-4]")
pheno$Replicate <- gsub("SC", "Rep", pheno$Replicate)
pheno$Condition <- str_extract(pheno$sample, "BIN[1-4]|Input")
pheno$Condition <- gsub("BIN", "Bin", pheno$Condition)
rownames(pheno) <- colnames(counts)

se <- SummarizedExperiment(counts,
                           rowData=ann,
                           colData=pheno)
save(se, file="objects/se.rda")


