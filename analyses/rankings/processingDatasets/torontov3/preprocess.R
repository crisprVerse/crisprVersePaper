library(dplyr)
library(crisprDesignData)
library(S4Vectors)
library(pbapply)


data <- read.csv("raw/hap1.foldchange", sep="\t")
data <- dplyr::rename(data, gene=GENE)
seqs <- data[,1]
seqs <- strsplit(seqs, split="_")
data$spacer <- vapply(seqs, function(x) x[[2]], FUN.VALUE="a")
data[,1] <- NULL
samples <- c("HAP1_T18_A",
             "HAP1_T18_A2",
             "HAP1_T18_B",
             "HAP1_T18_C")
data <- data[, c("gene", "spacer", samples)]
guideMap <- data

# Getting essential genes:
library(crisprDesignData)
library(S4Vectors)
key <- mcols(txdb_human$cds)[, c("gene_symbol", "gene_id")]
key <- as.data.frame(key)
key <- key[!duplicated(key),]
wh <- match(guideMap$gene, key$gene_symbol)
guideMap$gene_id <- key$gene_id[wh]
guideMap <- guideMap[complete.cases(guideMap),]
guideMap$lfc <- rowMeans(guideMap[,samples])
guideMap[,samples] <- NULL


# Normalization with egs/negs
load("../data/egs.rda")
load("../data/negs.rda")
neg.guides <- which(guideMap$gene %in% negs)
eg.guides <- which(guideMap$gene %in% egs)
meds <- median(guideMap$lfc[neg.guides], na.rm=TRUE)
guideMap$lfc <- guideMap$lfc-meds
meds <- abs(median(guideMap$lfc[eg.guides], na.rm=TRUE))
guideMap$lfc  <- guideMap$lfc/meds
guideMap$name <- paste0(guideMap$gene_id, "_", guideMap$spacer)
save(guideMap,
     file="processed/guideMap.rda")



