library("Gviz")
library("crisprDesign")
library("crisprDesignData")

kras_cds <- queryTxObject(txdb_human, "cds", "gene_id",
                          "ENSG00000133703")
kras_5utr <- queryTxObject(txdb_human, "fiveUTRs", "gene_id",
                           "ENSG00000133703")
kras_3utr <- queryTxObject(txdb_human, "threeUTRs", "gene_id",
                           "ENSG00000133703")

idTrack <- IdeogramTrack(chromosome="chr12", genome="hg38", size=0.5)
gaTrack <- GenomeAxisTrack(range=kras_tx, add53=TRUE, add35=TRUE)

kras_cds <- data.frame(chromosome=seqnames(kras_cds),
                       start=start(kras_cds),
                       end=end(kras_cds),
                       strand=as.character(strand(kras_cds)),
                       feature="protein_coding",
                       gene=kras_cds$gene_id,
                       exon=kras_cds$exon_id,
                       transcript=kras_cds$tx_id,
                       symbol=kras_cds$tx_id)
kras_5utr <- data.frame(chromosome=seqnames(kras_5utr),
                        start=start(kras_5utr),
                        end=end(kras_5utr),
                        strand=as.character(strand(kras_5utr)),
                        feature="utr5",
                        gene=kras_5utr$gene_id,
                        exon=kras_5utr$exon_id,
                        transcript=kras_5utr$tx_id,
                        symbol=kras_5utr$tx_id)
kras_3utr <- data.frame(chromosome=seqnames(kras_3utr),
                        start=start(kras_3utr),
                        end=end(kras_3utr),
                        strand=as.character(strand(kras_3utr)),
                        feature="utr3",
                        gene=kras_3utr$gene_id,
                        exon=kras_3utr$exon_id,
                        transcript=kras_3utr$tx_id,
                        symbol=kras_3utr$tx_id)
kras <- rbind(kras_cds, kras_5utr, kras_3utr)
geneTrack <- GeneRegionTrack(kras, genome="hg38", chromosome="chr12",
                             name="KRAS", showId=TRUE,
                             background.title="darkgray")


pdf(file="~/crisprOven_paper/figures/krasplot.pdf", width=7, height=1.75)
plotTracks(list(idTrack, gaTrack, geneTrack))
dev.off()
