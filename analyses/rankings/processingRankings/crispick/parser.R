files <- list.files("outputs",
                    pattern="sgrna-designs.txt",
                    full.names=TRUE)

parseOneChunk <- function(file){
    data <- read.csv(file, sep="\t")
    data <- data[,c("sgRNA.Sequence",
                    "Target.Gene.Symbol",
                    "Combined.Rank",
                    "Pick.Order")]
    colnames(data) <- c("spacer_20mer", 
                        "gene_symbol",
                        "rank",
                        "pickingOrder")
    return(data)
}
chunks <- lapply(files, parseOneChunk)
lib <- do.call(rbind, chunks)
dir.create("../../objects")
saveRDS(lib,
     file="../../objects/crispick.results.rds")

