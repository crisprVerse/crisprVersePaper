library(crisprDesignData)
genes <- unique(txdb_human$cds$gene_symbol)
cats <- cut(seq_along(genes),
            breaks=c(seq(1, length(genes),500), length(genes)))
cats <- as.numeric(as.factor(cats))
genes <- split(genes, f=cats)
for (i in 1:length(genes)){
    write.table(genes[[i]],
                quote=FALSE,
                row.names=FALSE,
                col.names=FALSE,
                file=file.path("inputs",
                               paste0("genes_chunk_", i, ".txt")))
}
