#ml R/dev
library(crisprDesignData)
library(crisprDesign)
library(BiocGenerics)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
data(txdb_human, package="crisprDesignData")
bsgenome <- BSgenome.Hsapiens.UCSC.hg38
txObject <- txdb_human
exons <- txObject[["exons"]]
exons <- exons[order(exons$tx_id, exons$exon_rank)]
# Getting flanking regions:
start(exons) <- start(exons)-35
end(exons) <- end(exons)+35
S4Vectors::mcols(exons)$sequence <- Biostrings::getSeq(bsgenome, exons)
exons <- S4Vectors::split(exons, f = exons$gene_id)
sequences <- lapply(exons, function(x){
    paste0(x$sequence, collapse = "")
})
codingGenes <- unique(txObject[["cds"]]$gene_id)
sequences <- sequences[names(sequences) %in% codingGenes]

genes <- names(sequences)
write.table(genes,
            row.names=FALSE,
            col.names=FALSE,
            quote=FALSE,
            file="../genes.txt")


outdir <- "../inputs/"
for (k in seq_along(sequences)){
    geneid <- names(sequences)[k]
    filename <- paste0(geneid, ".fasta")
    filename <- file.path(outdir, filename)
    seq <- DNAStringSet(sequences[[k]])
    names(seq) <- geneid
    writeXStringSet(seq, filename)
    print(k)
}




