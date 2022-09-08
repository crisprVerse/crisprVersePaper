file <- "raw/readcount-HCT116_1-lib1"
data <- read.table(file, header=TRUE)
rownames(data) <- data[,1]
data <- data[,-c(1:2)]
data1=data
colnames(data1) <- gsub("LIB1_", "", colnames(data1))

file <- "raw/readcount-HCT116_1-lib2"
data <- read.table(file, header=TRUE)
rownames(data) <- data[,1]
data <- data[,-c(1:2)]
data2=data
colnames(data2) <- gsub("LIB2_", "", colnames(data2))
data <- rbind(data1,data2)


guideMap <- data.frame(name=rownames(data))
grnas <- strsplit(guideMap$name, split="_")
grnas <- unlist(lapply(grnas, function(x) x[[2]]))
genes <- strsplit(guideMap$name, split="_")
genes <- unlist(lapply(genes, function(x) x[[1]]))
guideMap$spacer <- grnas
guideMap$gene <- genes

library(crisprDesignData)
library(S4Vectors)
key <- mcols(txdb_human$cds)[, c("gene_symbol", "gene_id")]
key <- as.data.frame(key)
key <- key[!duplicated(key),]
wh <- match(guideMap$gene, key$gene_symbol)
guideMap$gene_id <- key$gene_id[wh]
guideMap <- guideMap[complete.cases(guideMap),]



Y <- log2(data+1)
lfc <- rowMeans(Y[,c("T18_A", "T18_B")])-Y[,"T0"]
wh <- match(guideMap$name, names(lfc))
guideMap$lfc <- lfc[wh]
rownames(guideMap) <- NULL


# Normalization with egs/negs
library(matrixStats)
load("../data/egs.rda")
load("../data/negs.rda")
neg.guides <- which(guideMap$gene %in% negs)
eg.guides <- which(guideMap$gene %in% egs)
meds <- median(guideMap$lfc[neg.guides], na.rm=TRUE)
guideMap$lfc <- guideMap$lfc-meds
meds <- abs(median(guideMap$lfc[eg.guides], na.rm=TRUE))
guideMap$lfc  <- guideMap$lfc/meds
guideMap$name <- paste0(guideMap$gene_id, "_", guideMap$spacer)
save(guideMap, file="processed/guideMap.rda")








