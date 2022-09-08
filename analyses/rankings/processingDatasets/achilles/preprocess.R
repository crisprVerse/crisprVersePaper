library(crisprDesignData)
library(S4Vectors)
library(readr)
key <- mcols(txdb_human$cds)[, c("gene_symbol", "gene_id")]
key <- as.data.frame(key)
key <- key[!duplicated(key),]


guideMap <- read.csv("raw/Achilles_guide_map.csv")
colnames(guideMap) <- c("spacer", "coord", "gene_symbol", "n")
guideMap <- guideMap[,c(1,3,4)]
gene_symbol <- guideMap$gene_symbol
gene_symbol <- gsub(" \\([0-9]+\\)", "", gene_symbol)
guideMap$gene_symbol <- gene_symbol
guideMap <- guideMap[guideMap$gene_symbol %in% key$gene_symbol,]
wh <- match(guideMap$gene_symbol, key$gene_symbol)
guideMap$ensembl_id <- key$gene_id[wh]
guideMap <- guideMap[guideMap$n==1,]
guideMap$n <- NULL
guideMap <- guideMap[order(guideMap$ensembl_id),]
rownames(guideMap) <- NULL

#Let's load efficacy:
df <- read.csv("raw/Achilles_guide_efficacy.csv")
colnames(df) <- c("spacer", "efficacy")
wh <- match(guideMap$spacer, df$spacer)
guideMap$efficacy <- df[wh, "efficacy"]



# Let's get the logFCs:
lfc <- read_csv("raw/Achilles_logfold_change.csv")
lfc <- as.data.frame(lfc)
rownames(lfc) <- lfc[,1]
lfc <- lfc[,-1]
wh <- match(guideMap$spacer, rownames(lfc))
lfc <- lfc[wh,]
lfc <- as.matrix(lfc)


samples <- colnames(lfc)
meta <- read.csv("raw/Achilles_replicate_map.csv")
meta <- meta[meta[,1] %in% samples,]
nsamples <- length(unique(meta[,2]))




# Normalization
library(matrixStats)
load("../data/egs.rda")
load("../data/negs.rda")
neg.guides <- guideMap$spacer[guideMap$gene_symbol %in% negs]
eg.guides <- guideMap$spacer[guideMap$gene_symbol %in% egs]
meds <- colMedians(lfc[neg.guides,], na.rm=TRUE)
lfc  <- sweep(lfc, 2, meds, "-")
meds <- colMedians(lfc[eg.guides,], na.rm=TRUE)
meds <- meds/median(meds)
lfc  <- sweep(lfc, 2, meds, "/")


# Final guideMap
score <- matrixStats::rowMedians(as.matrix(lfc), na.rm=TRUE)
guideMap$lfc <- score
save(guideMap, file="processed/guideMap.rda")




