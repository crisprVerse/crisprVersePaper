library(readxl)
ann <- read_excel("raw/aac7041_sm_table_s1.xlsx")
ann <- as.data.frame(ann)
ann <- ann[,c(1,2,6)]
colnames(ann) <- c("ID", "gene", "spacer")
rownames(ann) <- ann$ID
guideMap <- ann

data <- read_excel("raw/aac7041_sm_table_s2.xlsx")
data <- as.data.frame(data)
rownames(data) <- data[,1]
data <- data[,-1]
data <- data[rownames(ann),]
colnames(data) <- gsub(" ", "_", colnames(data))
data <- data[,c("K562_final", "K562_initial")]

# Library size normalization:
Y <- log2(data+1)
factors <- colSums(Y)
factors <- factors/mean(factors)
Y <- sweep(Y, 2, factors, "/")
lfc <- Y[,"K562_final"]-Y[,"K562_initial"]
names(lfc) <- rownames(data)


# Normalization with egs/negs
library(matrixStats)
load("../data/egs.rda")
load("../data/negs.rda")
neg.guides <- which(guideMap$gene %in% negs)
eg.guides <- which(guideMap$gene %in% egs)
meds <- median(lfc[neg.guides], na.rm=TRUE)
lfc <- lfc-meds
meds <- abs(median(lfc[eg.guides], na.rm=TRUE))
lfc  <- lfc/meds

guideMap$lfc <- lfc
guideMap$baseline <- data[,"K562_initial"]




library(crisprDesignData)
library(S4Vectors)
key <- mcols(txdb_human$cds)[, c("gene_symbol", "gene_id")]
key <- as.data.frame(key)
key <- key[!duplicated(key),]
wh <- match(guideMap$gene, key$gene_symbol)
guideMap$gene_id <- key$gene_id[wh]
guideMap <- guideMap[complete.cases(guideMap),]
rownames(guideMap) <- NULL
guideMap$name <- paste0(guideMap$gene_id, "_", guideMap$spacer)

guideMap <- guideMap[guideMap$baseline>=30,]
save(guideMap, file="processed/guideMap.rda")


