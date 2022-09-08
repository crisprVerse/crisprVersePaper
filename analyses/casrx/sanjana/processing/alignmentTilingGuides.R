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
guidesTiling <- guides
guidesTiling <- do.call(c, guidesTiling)


### Bowtie alignment:
bowtie_index="/Users/fortinj2/crisprIndices/bowtie/ensembl_human_104/ensembl_human_104"
guidesTilingWithAlignments <- addSpacerAlignments(guidesTiling,
                            addSummary=TRUE,
                            txObject=txObject,
                            n_mismatches=3,
                            aligner="bowtie",
                            aligner_index=bowtie_index)
guidesTilingWithAlignments <- temp
save(guidesTilingWithAlignments, file="objects/guidesTilingWithAlignments.rda")



