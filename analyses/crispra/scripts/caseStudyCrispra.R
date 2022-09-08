# Author: Luke Hoberecht
# Copyright 2022, Genentech, Inc.

library(crisprDesign)
library(crisprDesignData)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Gviz)
library(VariantAnnotation)
library(AnnotationHub)
library(RMariaDB)
library(rtracklayer)
data("SpCas9",
     "SpGCas9",
     "AsCas12a",
     "enAsCas12a",
     "gr.repeats.hg38")

# Plot color scheme:
colorScheme <- list(
    background="#FFFFFF",
    genomeAxis="darkgray",
    highlight="#FFEBED",
    geneBorder="#808080",
    geneFill="#ffd966",
    trackDivider="#000000",
    score0="#D2D2D2",    # light gray
    score1="#000080",    # navy blue
    snps="#000000",
    dnase="#CECECE",
    cage="#CECECE")


bsgenome <- BSgenome.Hsapiens.UCSC.hg38
grList <- TxDb2GRangesList(getTxDb(), genome="hg38")
seqlevels(grList)[seqlevels(grList) == "MT"] <- "M"
seqlevels(grList) <- paste0("chr", seqlevels(grList))

gene <- queryTxObject(grList, "exons", "gene_symbol", "MMP7")
promoter <- queryTss(grList[["tss"]],
                     "gene_symbol",
                     "MMP7",
                     tss_window = c(-150, -75))


## small functions ============================================================

guideSet2GRanges <- function(guideSet
){
    gr <- GRanges(seqnames=seqnames(guideSet),
                  ranges=ranges(guideSet),
                  strand=strand(guideSet))
    mcols(gr) <- mcols(guideSet)
    return(gr)
}


addScoreColor <- function(guideSet,
                          scoreCol,
                          score0,
                          score1
){
    colFunc <- colorRampPalette(c(colorScheme$score0, colorScheme$score1))
    colors <- colFunc(100)
    if (is.na(scoreCol)){
        intervals <- rep(1, length(guideSet))
    } else {
        intervals <- cut(mcols(guideSet)[[scoreCol]], breaks=seq(0, 1, by=0.01))
        intervals <- as.numeric(intervals)
    }
    guideSet$color <- colors[intervals]
    guideSet
}





# move to track definitions
makeGeneTrack <- function(gene
){
    df <- data.frame(chromosome=seqnames(gene),
                     start=start(gene),
                     end=end(gene),
                     width=width(gene),
                     strand=as.character(strand(gene)),
                     feature="protein_coding",
                     gene=gene$gene_id,
                     exon=gene$exon_id,
                     transcript=gene$tx_id,
                     symbol=gene$gene_symbol)
    GeneRegionTrack(df, genome="hg38",
                    chromosome=as.character(unique(seqnames(gene))),
                    background.title="transparent",
                    col=colorScheme$geneBorder,
                    fill=colorScheme$geneFill)
}
geneTrack <- makeGeneTrack(gene)



ah <- AnnotationHub()
chain <- import.chain("../objects/hg19ToHg38.over.chain")

makeDnaseTrack <- function(ah,
                           chain
){
    # DNase I hypersensitive sites
    # title: E116-DNase.macs2.narrowPeak.gz
    # description: Narrow DNasePeaks for consolidated
    # epigenomes from EpigenomeRoadMap Project
    dnase <- ah[['AH30743']]
    dnase <- rtracklayer::liftOver(dnase, chain)
    dnase <- unlist(dnase)
    ov <- suppressWarnings(findOverlaps(promoter, dnase))
    dnase <- dnase[subjectHits(ov)]
    AnnotationTrack(dnase,
                    shape="box",
                    fill=colorScheme$dnase,
                    color=colorScheme$dnase,
                    background.title="transparent",
                    cex.title=0,
                    size=0.8)
}
dnaseTrack <- makeDnaseTrack(ah, chain)


makeCageTrack <- function(ah,
                          chain
){
    # CAGE peaks
    # title: RIKEN CAGE Loc
    # description: GRanges object from UCSC track
    #'RIKEN CAGE Loc'
    cage <- ah[["AH5084"]]
    cage <- rtracklayer::liftOver(cage, chain)
    cage <- unlist(cage)
    ov <- suppressWarnings(findOverlaps(promoter, cage))
    cage <- cage[subjectHits(ov)]
    AnnotationTrack(cage,
                    shape="box",
                    fill=colorScheme$cage,
                    color=colorScheme$cage,
                    background.title="transparent",
                    cex.title=0,
                    size=0.8)
}
cageTrack <- makeCageTrack(ah, chain)



ranges(promoter) <- IRanges(start=end(cageTrack)+75, width=75)


getGuides <- function(promoter,
                      bsgenome
){
    vcf_file <- "https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-common_all.vcf.gz"
    
    spcas9 <- findSpacers(promoter,
                          crisprNuclease=SpCas9,
                          bsgenome=bsgenome)
    spcas9 <- addOnTargetScores(spcas9,
                                methods="deephf")
    spcas9 <- addRepeats(spcas9,
                         gr.repeats=gr.repeats.hg38)
    spcas9 <- addSNPAnnotation(spcas9,
                               vcf=vcf_file)
    
    spgcas9 <- findSpacers(promoter,
                           crisprNuclease=SpGCas9,
                           bsgenome=bsgenome)
    spgcas9 <- addRepeats(spgcas9,
                          gr.repeats=gr.repeats.hg38)
    spgcas9 <- addSNPAnnotation(spgcas9,
                                vcf=vcf_file)
    
    ascas12a <- findSpacers(promoter,
                            crisprNuclease=AsCas12a,
                            bsgenome=bsgenome, canonical=FALSE)
    ascas12a <- addOnTargetScores(ascas12a,
                                  methods="deepcpf1")
    ascas12a <- addRepeats(ascas12a,
                           gr.repeats=gr.repeats.hg38)
    ascas12a <- addSNPAnnotation(ascas12a,
                                 vcf=vcf_file)
    
    enascas12a <- findSpacers(promoter,
                              crisprNuclease=enAsCas12a,
                              bsgenome=bsgenome, canonical=FALSE)
    enascas12a <- addOnTargetScores(enascas12a,
                                    methods="enpamgb")
    enascas12a <- addRepeats(enascas12a, 
                             gr.repeats=gr.repeats.hg38)
    enascas12a <- addSNPAnnotation(enascas12a,
                                   vcf=vcf_file)
    
    guides <- list(spcas9, spgcas9, ascas12a, enascas12a)
    score <- list("score_deephf",
                  NA,
                  "score_deepcpf1",
                  "score_enpamgb")
    guides <- lapply(1:4, function(x){
        spacers <- getSpacerRanges(guides[[x]],
                                   nuclease=crisprNuclease(guides[[x]]))
        spacers <- addScoreColor(spacers, score[[x]], colorScheme$score0, colorScheme$score1)
        spacers <- guideSet2GRanges(spacers)
        ranges(spacers) <- IRanges(start=spacers$pam_site, width=1)
        spacers
    })
    return(guides)
}
# spacers
guides <- suppressWarnings(getGuides(promoter, bsgenome))






repeats <- lapply(guides, function(x){
    any(x$inRepeats)
})
any(unlist(repeats)) # no repeats

snps <- lapply(guides, function(x){
    any(x$hasSNP)
})
any(unlist(snps)) # has SNPs


makeSnpTrack <- function(guides
){
    chr <- unique(as.character(seqnames(guides[[1]])))
    snps <- lapply(guides, function(x){
        hits <- unlist(x$snps)
        data.frame(rs=unique(hits$rs),
                   rs_site=unique(hits$rs_site))
    })
    snps <- Reduce("rbind", snps)
    snps <- unique(snps)
    snps <- GRanges(seqnames=chr,
                    ranges=IRanges(start=snps$rs_site, width=1))
    AnnotationTrack(range=snps, genome="hg38", chromosome=chr,
                    shape="box", stackHeight=1, size=0.6,
                    min.height=1, cex.title=0,
                    background.title="transparent",
                    fill=colorScheme$snps, col=colorScheme$snps)
}
snpTrack <- makeSnpTrack(guides)







# gRNA tracks
spacerTrack <- lapply(seq_along(guides), function(x){
    guides[[x]] <- guides[[x]][!duplicated(ranges(guides[[x]]))]
    AnnotationTrack(range=guides[[x]], genome="hg38", chromosome="chr4",
                    shape="box",
                    stackHeight=0.8, size=1, min.width=1,
                    min.height=1, fill=guides[[x]]$color,
                    background.title="transparent",
                    col="black", cex.title=0)
})









makeGapTrack <- function(promoter
){
    range <- GRanges(seqnames=unique(as.character(seqnames(promoter))),
                     ranges=IRanges(start=1, width=1))
    AnnotationTrack(range=range,
                    genome="hg38",
                    name="", min.height=1, size=0.3,
                    fontcolor.title="transparent",
                    fontcolor="transparent",
                    col.frame="transparent",
                    col.border.title="transparent",
                    col.grid="transparent",
                    col.axis="transparent",
                    background.title="transparent",
                    fontcolor.item="transparent",
                    fontcolor.group="transparent",
                    fill="transparent",
                    col="transparent",
                    col.line="transparent")
}
gapTrack <- makeGapTrack(promoter)








# non-Gviz plotting elements
drawGridText <- function(text,
                         y
){
    grid.text(text, x=0.02, y=y, just=c("left", "center"),
              gp=gpar(fontsize=10))
}

renderLegend <- function(){
    grid.rect(x=unit(0.25, "npc"), y=unit(0.03, "npc"),
              width=unit(0.5, "npc"), height=unit(0.92, "npc"),
              just=c(0, 0), gp=gpar(col="black", lwd=1, fill="transparent"))
    # gRNA
    grid.rect(x=unit(0.32, "npc"), y=unit(0.46, "npc"), just="left",
              width=unit(0.007, "npc"), height=unit(0.4, "npc"),
              gp=gpar(fill="transparent", col="black"))
    grid.text(label="PAM site", x=unit(0.32, "npc"), y=unit(0.78, "npc"),
              just="center", gp=gpar(fontsize=10))
    # color scale
    colFunc <- colorRampPalette(c(colorScheme$score0, colorScheme$score1))
    colors <- colFunc(100)
    grid.rect(x=unit(seq(0.4 ,0.7, length=100), "npc"), y=unit(0.46, "npc"),
              width=unit(0.003, "npc"), height=unit(0.4, "npc"), 
              just="left", gp=gpar(col="transparent", fill=colors, (100), "y"))
    grid.text(label="On-target score", x=unit(0.55, "npc"), y=unit(0.78, "npc"),
              just="center", gp=gpar(fontsize=10))
    grid.text(label="0", x=unit(0.4, "npc"), y=unit(0.15, "npc"), just="center",
              gp=gpar(fontsize=10))
    grid.text(label="1", x=unit(0.703, "npc"), y=unit(0.15, "npc"), just="center",
              gp=gpar(fontsize=10))
}






makePlot <- function(){
    chr <- as.character(unique(seqnames(promoter)))
    tss <- end(cageTrack)
    grid.newpage()
    pushViewport(viewport(layout=grid.layout(nrow=4, ncol=2,
                                             heights=c(9, 1, 16, 4),
                                             widths=c(1, 10))))
    grid.rect(gp=gpar(fill=colorScheme$background, col=colorScheme$background))
    # top plot: gene, promoter region
    idTrack <- IdeogramTrack(chromosome=chr, genome="hg38", size=0.7)
    gaTrack <- GenomeAxisTrack(add53=TRUE, add35=TRUE)
    trackList <- list(gaTrack, gapTrack,
                      geneTrack, gapTrack,
                      dnaseTrack, gapTrack,
                      cageTrack, gapTrack)
    hTrack <- HighlightTrack(trackList=trackList, chromosome=chr,
                             range=promoter,
                             fill=colorScheme$highlight, col="transparent")
    pushViewport(viewport(layout.pos.col=c(1,2), layout.pos.row=1))
    plotTracks(list(idTrack, hTrack), chromosome=chr, add=TRUE,
               from=tss-250, to=tss+500)
    grid.text("MMP7", x=0.132, y=0.442, gp=gpar(fontface="italic", fontsize=10))
    grid.text("DNase I hypersensitive site", x=0.325, y=0.283,
              gp=gpar(fontsize=9, fontface="italic"))
    grid.text("CAGE peak", x=0.276, y=0.14,
              gp=gpar(fontsize=9, fontface="italic"))
    # TSS arrow
    grid.lines(x=c(0.225, 0.225), y=c(0.396, 0.55),
               gp=gpar(col=colorScheme$geneBorder, lwd=2))
    grid.lines(x=c(0.225, 0.18), y=c(0.55, 0.55),
               gp=gpar(col=colorScheme$geneBorder, lwd=2),
               arrow=arrow(angle=17, length=unit(0.05, "npc"), type="closed"))
    popViewport(1)
    
    # bottom plot: guides
    pushViewport(viewport(layout.pos.col=2, layout.pos.row=3))
    tracks <- c(gaTrack, spacerTrack[[1]], spacerTrack[[2]],
                spacerTrack[[3]], spacerTrack[[4]], gapTrack,
                snpTrack, gapTrack)
    hgTrack <- HighlightTrack(trackList=tracks,
                              range=promoter,
                              chromosome=chr,
                              fill=colorScheme$highlight, col="transparent")
    plotTracks(hgTrack, chromosome=chr, add=TRUE,
               from=start(promoter)-20, to=end(promoter)+20)
    popViewport(1)
    # legend
    pushViewport(viewport(layout.pos.col=c(1,2), layout.pos.row=4))
    renderLegend()
    popViewport(1)

    # additional drawings
    grid.polygon(x=c(0.455,0.545,0.845,0.282),
                 y=c(0.715,0.715,0.65,0.65),
                 gp=gpar(fill=colorScheme$highlight,
                         col=colorScheme$highlight))
    # labels
    ys <- c(0.195, seq(from=0.283, by=0.079, length.out=4))
    labels <- c("SNPs",
                "enAsCas12a",
                "AsCas12a",
                "SpGCas9",
                "SpCas9")
    for (i in seq_along(labels)){
        drawGridText(labels[i], ys[i])
    }
}

dev.new(width=10, height=6)
makePlot()
dev.off()
