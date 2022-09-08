library(crisprDesign)
library(crisprBase)
library(SummarizedExperiment)
library(crisprDesignData)
library(BSgenome.Hsapiens.UCSC.hg38)
data(CasRx)
bsgenome <- BSgenome.Hsapiens.UCSC.hg38
txObject <- txdb_human


load("../objects/se.rda")
txids <- unique(rowData(se)$txid)
txids <- txids[!is.na(txids)]
txids <- txids[txids %in% txObject$cds$tx_id]


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
save(guides, file="objects/guides.rda")



