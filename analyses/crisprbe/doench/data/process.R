library(readxl)
library(SummarizedExperiment)
library(matrixStats)

# Processing counts:
data <- read_excel("mmc2.xlsx")
data <- as.data.frame(data)
counts <- data[-c(1:2), c(2,6,7,8)]
counts <- as.matrix(counts)
counts <- apply(counts,2,as.numeric)
guides <- data[-c(1:2),1]
rownames(counts) <- guides
colnames(counts)[2:4] <- paste0("Dropout",1:3)
se <- SummarizedExperiment(counts)


# Processing feature annotation:
ann <- read_excel("mmc2.xlsx", sheet=3)
ann <- as.data.frame(ann)
colnames(ann) <- gsub(" ", "_", colnames(ann))
ann <- dplyr::rename(ann, spacer_20mer=sgRNA_sequence)
ann <- dplyr::rename(ann, gene_symbol=Gene_symbol)
ann <- dplyr::rename(ann, mutation=Mutation_category)

cols <- c("spacer_20mer",
          "gene_symbol", 
          "mutation")
ann <- ann[,cols]
rownames(ann) <- ann$spacer_20mer
ann <- ann[rownames(se),]
rowData(se) <- ann
save(se,
     file="../objects/se.rda")


# Normalization
Y <- assays(se)[[1]]
Y <- log2(Y+1)
meds <- colMedians(Y)
meds <- meds-median(meds)
Y <- sweep(Y, 2, meds, "-")
Y <- 2^Y-1
assays(se)[[1]] <- Y
save(se,
     file="../objects/se.rda")

