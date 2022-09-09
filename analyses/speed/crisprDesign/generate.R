library(crisprDesign)
library(BSgenome.Hsapiens.UCSC.hg38)
bsgenome <- BSgenome.Hsapiens.UCSC.hg38 
bowtie <- "/Users/fortinj2/crisprIndices/bowtie/hg38/hg38"
bwa <- "/Users/fortinj2/crisprIndices/bwa/hg38/hg38"

##### crisprDesign bowtie search #####
load("../inputs/inputGrs.rda")

getDelta <- function(i){
    gr <- inputGrs[[i]]
    time1 <- proc.time()
    gs <- findSpacers(gr, canonical=TRUE, bsgenome=bsgenome)
    gs <- addSpacerAlignmentsIterative(gs,
                                       n_mismatches=2,
                                       aligner_index=bowtie,
                                       bsgenome=bsgenome)
    time2 <- proc.time()
    delta <- as.numeric(time2-time1)
    return(delta)
}


resultsBowtieIter <- lapply(1:6, function(i){
    print(i)
    out <- do.call(rbind,lapply(1:3, function(x) getDelta(i)))
    out
})
save(resultsBowtieIter, file="resultsBowtieIter.rda")


getDelta <- function(i){
    gr <- inputGrs[[i]]
    time1 <- proc.time()
    gs <- findSpacers(gr, canonical=TRUE, bsgenome=bsgenome)
    gs <- addSpacerAlignmentsIterative(gs,
                                       n_mismatches=2,
                                       aligner="bwa",
                                       aligner_index=bwa,
                                       bsgenome=bsgenome)
    time2 <- proc.time()
    delta <- as.numeric(time2-time1)
    return(delta)
}


resultsBwaIter <- lapply(1:6, function(i){
    print(i)
    out <- getDelta(i)
    out
})
save(resultsBwaIter, file="resultsBwaIter.rda")









