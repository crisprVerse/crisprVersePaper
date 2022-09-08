library(crisprDesign)
library(crisprBase)
library(SummarizedExperiment)
library(crisprDesignData)
library(BSgenome.Hsapiens.UCSC.hg38)
data(CasRx)
bsgenome <- BSgenome.Hsapiens.UCSC.hg38
txObject <- txdb_human

#genes <- c("CD46", "CD55", "CD71")
genes <- c("CD46", "CD55", "TFRC")

key <- crisprDesign:::.getTx2GeneTable(txObject)
key <- key[key$gene_symbol %in% genes,]
txids <- key$tx_id


mrnaSequences <- lapply(txids,
                        getMrnaSequences,
                        bsgenome=bsgenome,
                        txObject=txObject)
guides <- lapply(mrnaSequences,
                 findSpacers,
                 crisprNuclease=CasRx)


# Adding on-target scores:
for (i in 1:length(guides)){
    guides[[i]] <- addOnTargetScores(guides[[i]])
    print(i)
}
guidesTiling <- guides
save(guidesTiling,
     file="objects/guidesTiling.rda")




