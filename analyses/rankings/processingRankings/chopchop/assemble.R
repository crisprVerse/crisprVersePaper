#ml R/dev
#R
library(pbapply)
library(dplyr)
genes <- list.files("outputs/")
files <- file.path("outputs", genes, "results.txt")
    cols <- c("rank", "protospacer", "orientation", "AggregateRankedScore_medianRank")                
dfs <- pblapply(files, function(file){
    data <- read.csv(file, header=TRUE, sep="\t")
    data <- data[,c(1,2),drop=FALSE]
    colnames(data) <- c("rank", "protospacer")
    data <- data[!duplicated(data$protospacer),,drop=FALSE]
    data <- data[order(data$rank),,drop=FALSE]
    data$rank <- seq_len(nrow(data))
    return(data)
})
for (i in 1:length(dfs)){
    gene <- genes[i]
    dfs[[i]]$ensembl_id <- gene
}
lib <- bind_rows(dfs)
saveRDS(lib, "chopchop.results.rds")

