library(readxl)
library(crisprDesignData)
library(S4Vectors)
library(readr)
ann <- read_excel("raw/mmc2.xlsx", sheet=2)
ann <- as.data.frame(ann)
ann <- ann[,c(1,2,3)]
colnames(ann) <- c("ID","gene", "spacer")
guideMap <- ann

key <- mcols(txdb_human$cds)[, c("gene_symbol", "gene_id")]
key <- as.data.frame(key)
key <- key[!duplicated(key),]
wh <- match(guideMap$gene, key$gene_symbol)
guideMap$gene_id <- key$gene_id[wh]
guideMap <- guideMap[complete.cases(guideMap),]
rownames(guideMap) <- NULL
guideMap$name <- paste0(guideMap$gene_id, "_", guideMap$spacer)


data <- read_table("raw/mmc6/grna_count_5-AML-lines.txt")
data <- as.data.frame(data)
rownames(data) <- data[,1]
data <- data[,-c(1,2)]
data <- data[rownames(data) %in% guideMap$ID,]
data <- data[match(guideMap$ID, rownames(data)),]
cols <- c("HumanV1-1", "HumanV1-2", "HL60_rep_1", "HL60_rep_2")
data <- data[,cols]


# Normalize:
factors <- colSums(data)
factors <- factors/median(factors)
data <- sweep(data, 2, factors, "/")
Y <- log2(data+1)
Y1 <- Y[,c("HL60_rep_1", "HL60_rep_2")]
Y2 <- Y[,c("HumanV1-1", "HumanV1-2")]
lfc <- rowMeans(Y1-Y2)
guideMap$lfc <- lfc
guideMap[guideMap$gene=="RAN",]
guideMap$ID <- NULL


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



