library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome)
library(multicrispr)
bsgenome <- BSgenome.Hsapiens.UCSC.hg38 


### Multicrispr:
load("../inputs/inputGrs.rda")
# This requires BSgenome.Hsapiens.UCSC.hg38 to be 
# a bowtie index in the bowtie folder below
bowtie <- "/Users/fortinj2/crisprIndices/bowtie"

getDelta <- function(i){
    gr <- inputGrs[[i]]
    time1 <- proc.time()
    gs <- find_spacers(gr,
                       plot=FALSE,
                       bsgenome=bsgenome,
                       ontargetmethod=NULL,
                       offtargetmethod="bowtie",
                       mismatches=2,
                       outdir=tempdir(),
                       indexedgenomesdir=bowtie)
    time2 <- proc.time()
    delta <- as.numeric(time2-time1)
    print(delta)
    return(delta)
}


results <- lapply(4:6, function(i){
    print(i)
    out <- getDelta(i)
    out[c(1,2,4,5)] <- NA # Removing unnecessary columns
    out
})
save(results, file="results.rda")

