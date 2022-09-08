library(crisprDesignData)
library(S4Vectors)
library(dplyr)
library(stringr)
data(txdb_human)

# Spacers for SANGER are 19mer
# so we need to convert our spacers:
rankings <- readRDS("../../processingRankings/objects/rankings.rds")
spacers <- str_extract(rankings$name, "_[A-Z]+")
spacers <- gsub("_", "", spacers)
spacers <- substr(spacers,2,20)
rankings$spacer <- spacers
rankings$name <- paste0(rankings$ensembl_id, "_", rankings$spacer)


load("processed/guideMap.rda")
guideMap$name <- paste0(guideMap$gene_id, "_", guideMap$spacer)
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
rankings_yusa <- df[df$ensembl_id %in% egs,]
save(rankings_yusa,
     file="../objects/rankings_yusa.rda")
rankings_yusa_neg <- df[df$ensembl_id %in% negs,]
save(rankings_yusa_neg,
     file="../objects/rankings_yusa_neg.rda")


