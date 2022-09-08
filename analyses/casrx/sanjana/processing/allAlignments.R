library(crisprDesign)
library(crisprBase)
library(SummarizedExperiment)
library(crisprDesignData)
library(BSgenome.Hsapiens.UCSC.hg38)
data(CasRx)
bsgenome <- BSgenome.Hsapiens.UCSC.hg38
txObject <- txdb_human
load("objects/guides.rda")


# Adding on-target scores:
for (i in 1:length(guides)){
    guides[[i]] <- addOnTargetScores(guides[[i]])
    print(i)
}


### Bowtie alignment:
bowtie_index="/Users/fortinj2/crisprIndices/bowtie/ensembl_human_104/ensembl_human_104"
temp <- addSpacerAlignments(guides[[1]],
                            addSummary=FALSE,
                            txObject=txObject,
                            n_mismatches=2,
                            aligner="bowtie",
                            aligner_index=bowtie_index)




