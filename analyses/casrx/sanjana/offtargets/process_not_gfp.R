library(SummarizedExperiment)
library(stringr)


addClass <- function(ann){
  ann <- ann[!grepl("intron", ann$ID),]
  ann <- ann[!grepl("LengthVariant", ann$ID),]
  ann <- ann[!grepl("LengthVariant", ann$ID),]
  ann <- ann[!grepl("RevComp", ann$ID),]
  
  ann$class <- str_extract(ann$ID,
                           "randomDouble|consecTriple|consecDouble|FirstOrder")
  ann$class[grepl("rc_",ann$ID)]<- "control"
  ann$class[is.na(ann$class)] <- "perfectMatch"
  ann$class[ann$class=="FirstOrder"] <- "singleMismatch"
  ann$class[ann$class=="randomDouble"] <- "doubleMismatch"
  ann
}



# Annotation
genes <- c("CD46", "CD55", "CD71")
files <- paste0(genes, "_library_final.fa")
anns <- lapply(files, function(file){
  out <- read.csv(file.path("../data", file),
                  head=FALSE)
  out <- out[,1]
  out <- matrix(out, ncol=2, byrow=TRUE)
  colnames(out) <- c("ID", "spacer")
  out <- as.data.frame(out)
  out$ID <- gsub(">","", out$ID)
  out
})
names(anns) <- genes
anns <- lapply(anns, addClass)


anns <- lapply(anns, function(ann){
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
  ann
})
annsTiling <- anns
save(annsTiling,
     file="objects/annsTiling.rda")



