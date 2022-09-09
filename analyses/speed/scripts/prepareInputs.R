library(crisprDesign)
library(crisprDesignData)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome)
library(Biostrings)
bsgenome <- BSgenome.Hsapiens.UCSC.hg38 
data(txdb_human)
gr <- txdb_human$exons
gr <- gr[seqnames(gr)=="chr1"]
gr <- unique(gr)
gr <- union(gr, gr)

# Creating sets for alignment:
set.seed(10)
set1 <- sample(gr, 100)
set2 <- sample(gr, 200)
set3 <- sample(gr, 400)
set4 <- sample(gr, 800)
set5 <- sample(gr, 1600)
set6 <- sample(gr, 3200)
sets <- list(set1,set2,set3,set4,set5,set6)
sets <- lapply(sets, function(gr){
    names(gr) <- paste0("exon_", 1:length(gr))
    gr
})
nbases <- sapply(sets, function(x) sum(width(x)))

#Now let's get the sequences
seqs <- lapply(sets, function(gr){
    x <- getSeq(bsgenome,gr)
    names(x) <- paste0("exon_", 1:length(x))
    x
})
inputDir <- "../inputs"
dir.create(inputDir, r=TRUE)
names <- 1:6
for (i in 1:length(seqs)){
    name <- paste0("../inputs/sequences_",names[i], ".fasta" )
    writeXStringSet(seqs[[i]], format="fasta", filepath=name) 
}
inputGrs <- sets
save(inputGrs, file="../inputs/inputGrs.rda")
