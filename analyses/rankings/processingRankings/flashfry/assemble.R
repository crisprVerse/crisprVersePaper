#ml R/dev
#R
library(pbapply)
library(dplyr)
files <- list.files("outputs/", pattern=".output.scored", full.names=TRUE)
cols <- c("contig", "target", "orientation", "AggregateRankedScore_medianRank")
dfs <- pblapply(files, function(file){
    data <- read.table(file, header=TRUE)[,cols,drop=FALSE]
    data <- data[!duplicated(data$target),,drop=FALSE]
    return(data)
})
lib <- bind_rows(dfs)
lib <- lib[,c(1,2,4)]
colnames(lib) <- c("ensembl_id", "protospacer", "rank")
saveRDS(lib, "flashfry.results.rds")
