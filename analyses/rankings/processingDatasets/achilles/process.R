library(crisprDesignData)
library(S4Vectors)
library(dplyr)
library(pbapply)
library(scales)
data(txdb_human)


load("processed/guideMap.rda")
rankings <- readRDS("../../processingRankings/objects/rankings.rds")
guideMap$name <- paste0(guideMap$ensembl_id, "_", guideMap$spacer)
common <- intersect(guideMap$name,rankings$name)
guideMap <- guideMap[match(common, guideMap$name),]
rankings <- rankings[match(common, rankings$name),]


# Only taking in common
good <- which(complete.cases(rankings))
rankings <- rankings[good,]
guideMap <- guideMap[good,]

# Ready for analysis:
df <- rankings
df$lfc <- guideMap$lfc

# Getting essential genes:
load("../../data/egs.rda")
load("../../data/negs.rda")
key <- mcols(txdb_human$cds)[, c("gene_symbol", "gene_id")]
key <- as.data.frame(key)
key <- key[!duplicated(key),]
egs <- unique(key$gene_id[key$gene_symbol %in% egs])
negs <- unique(key$gene_id[key$gene_symbol %in% negs])
rankings_achilles <- df[df$ensembl_id %in% egs,]
save(rankings_achilles,
     file="../objects/rankings_achilles.rda")
rankings_achilles_neg <- df[df$ensembl_id %in% negs,]
save(rankings_achilles_neg,
     file="../objects/rankings_achilles_neg.rda")

