#ml R/dev
#R
library(pbapply)
library(dplyr)

genes <- list.files("outputs", full.names=TRUE)
files <- file.path(genes, paste0(basename(genes), ".txt"))
good  <- file.exists(files)
files <- files[good]

readData <- function(file){
    a=readLines(file, n=2)
    if (length(a)>0){
        data <- read.table(file)
        data <- data[,c(2,3,6),drop=FALSE]
        data <- data[!duplicated(data),,drop=FALSE]
        colnames(data) <- c("protospacer", "rank", "CRISPRater")
    } else {
        data <- NA
    }
    return(data)
}
dfs <- pblapply(files, readData)
names(dfs) <- gsub(".txt","",basename(files))
bad <- vapply(dfs, function(x){
    x <- nrow(x)
    length(x)==0
}, FUN.VALUE=TRUE)
dfs <- dfs[!bad]
for (i in 1:length(dfs)){
    dfs[[i]]$ensembl_id <- names(dfs)[[i]]
}

lib <- dplyr::bind_rows(dfs)
saveRDS(lib, "cctop.results.rds")




