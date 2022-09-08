library(crisprDesign)
library(crisprDesignData)
library(crisprDesignGne)
data(tss_human)
data(txdb_human)
library(BSgenome.Hsapiens.UCSC.hg38)
bsgenome <- BSgenome.Hsapiens.UCSC.hg38
bowtie_index <- "/Users/fortinj2/crisprIndices/bowtie/hg38/hg38"
gene <- "KRAS"
txObject <- txdb_human
tssObject <- tss_human
mart_dataset <- "hsapiens_gene_ensembl"
vcf <- "/Users/fortinj2/crisprIndices/snps/dbsnp151.grch38/00-common_all.vcf.gz"

gr <- queryTxObject(txObject=txObject,
                    queryValue=gene,
                    queryColumn="gene_symbol",
                    featureType="cds")
gs <- findSpacers(gr, bsgenome=bsgenome)
gs <- addSpacerAlignments(gs,
                          bsgenome=bsgenome,
                          bowtie_index=bowtie_index,
                          n_mismatches=3,
                          txObject=txObject)
gs <- addGeneAnnotation(gs,
                        addPfam=FALSE,
                        mart_dataset=mart_dataset,
                        txObject=txObject)
gs <- addOffTargetScores(gs)
gs <- addTssAnnotation(gs,
                       tssObject=tssObject)
gs <- addSNPAnnotation(gs,
                       vcf=vcf)

aln <- alignments(gs, unlist=FALSE)[[1]]
ann <- geneAnnotation(gs, unlist=FALSE)[[1]]
snps <- snps(gs, unlist=FALSE)[[1]]

gs <- unique(gs)



subbset <- aln[c(1,3,483)]

