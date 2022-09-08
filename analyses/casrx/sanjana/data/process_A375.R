library(SummarizedExperiment)

# Feature annotation
data <- readLines("A375_library_final.fa")
data <- matrix(data, ncol=2, byrow=TRUE)
data <- as.data.frame(data)
colnames(data) <- c("ID", "spacer_27mer")
data$ID <- gsub(">","", data$ID)
ann <- data
rownames(ann) <- ann$ID
dfs <- strsplit(ann$ID, split="_")
genes <- sapply(dfs, function(x) x[[1]])
ann$gene_symbol <- genes

#Let's get genes:
genes <- read.table("TargetGenes.txt", head=TRUE)
genes <- genes[,c(1,2,9)]
colnames(genes) <- c("gene_symbol", "class", "txid")
genes$txid <- gsub("\\.[0-9]+","",genes$txid)
wh <- match(ann$gene_symbol, genes$gene_symbol)
ann$class <- genes$class[wh]
ann$txid <- genes$txid[wh]
ann$class[is.na(ann$class)] <- "ntc"
ann$txid[is.na(ann$txid)] <- NA



# Count matrix:
count <- read.csv("A375_count_matrix_raw.csv")
ann <- ann[rownames(count),]
se <- SummarizedExperiment(count,
                           rowData=ann)
se <- se[order(ann$class, ann$gene_symbol),]



# Normalization
Y <- assays(se)[[1]]
factors <- colSums(Y)
factors <- factors/median(factors)
Y <- sweep(Y, 2, factors, "/")
assays(se)[[1]] <- Y
save(se, file="../objects/se.rda")












#Let's calculate log-fold change
Y <- log2(assays(se)[[1]]+1)
lfc <- Y[,c("D14_R1","D14_R2")]-Y[,c("Input_R1","Input_R2")]
#lfc <- Y[,c("D7_R1","D7_R2")]-Y[,c("Input_R1","Input_R2")]
#lfc <- Y[,c("D14_R1","D14_R2")]-Y[,c("D7_R1","D7_R2")]
lfc <- rowMeans(lfc)
head(ann)
col <- as.numeric(as.factor(rowData(se)$class))
col <- as.numeric(as.factor(rowData(se)$gene_symbol))
plot(lfc, pch=20, cex=0.3,col=col)
abline(h=0, lty=3)




















