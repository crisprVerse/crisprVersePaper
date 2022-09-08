library(SummarizedExperiment)
library(stringr)


extractClass <- function(names){
    class <- rep("PerfectMatch", length(names))
    class[grepl("rc", names)] <- "NTC"
    class[grepl("RevComp", names)] <- "RevComp"
    class[grepl("intron", names)] <- "Intron"
    class[grepl("LengthVariant", names)] <- "LengthVariant"
    class[grepl("randomDouble", names)] <- "RandomDouble"
    class[grepl("FirstOrder", names)] <- "FirstOrder"
    return(class)
}



extractPheno <- function(names){
    out <- data.frame(sample=names)
    out$Condition <- str_extract(names, "top|bot|input")
    out$Replicate <- str_extract(names, "R1|R2|R3")
    return(out)
}




genes <- c("CD46", "CD55", "CD71")
fafiles <- paste0(genes, "_library_final.fa")
countfiles <- paste0(genes, "_count_matrix_raw.csv")


ses <- list()
for (i in 1:3){
    file <- fafiles[i]
    data <- readLines(file)
    data <- matrix(data, ncol=2, byrow=TRUE)
    data <- as.data.frame(data)
    colnames(data) <- c("ID", "spacer")
    data$class <- extractClass(data$ID)
    ann <- data
    ann$ID <- gsub(">","", ann$ID)
    rownames(ann) <- ann$ID
    
    file <- countfiles[i]
    count <- read.csv(file)
    pheno <- extractPheno(colnames(count))
    rownames(pheno) <- colnames(count)


    ann <- ann[rownames(count),]
    ses[[i]] <- SummarizedExperiment(as.matrix(count),
                                     rowData=ann,
                                     colData=pheno)
}
names(ses) <- genes


getDensities <- function(se){
    Y <- assays(se)[[1]]
    Y <- log2(Y+1)
    pheno <- colData(se)
    plot(density(Y[,1]), col="white")
    col <- as.numeric(as.factor(pheno$Condition))
    for (j in 1:ncol(Y)){
        lines(density(Y[,j]), col=col[j])
    }
}

par(mfrow=c(1,3))
getDensities(ses[[1]])
getDensities(ses[[2]])
getDensities(ses[[3]])


# Normalization
ses_norm <- lapply(ses, function(se){
    Y <- assays(se)[[1]]
    factors <- colMedians(Y)
    factors <- factors/median(factors)
    Y <- sweep(Y, 2, factors, "/")
    assays(se)[[1]] <- Y
    se
})


par(mfrow=c(1,3))
getDensities(ses_norm[[1]])
getDensities(ses_norm[[2]])
getDensities(ses_norm[[3]])

ses_tiling <- ses_norm
save(ses_tiling,
     file="../objects/ses_tiling.rda")



