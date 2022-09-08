library(crisprDesign)
library(crisprDesignAux)
library(crisprDesignGne)
library(crisprDesignData)
dir.create("figures")
dir.create("objects")
data("SpCas9", package="crisprBase")
species <- "human"
conservationFile <- getConservationFiles(species)
bsgenome <- getGenomePackage(species)
crisprNuclease <- SpCas9


# Getting protein-coding annotation:
data <- readRDS("../rankings/processingRankings/crisprverse/crisprverse.results.rds")
spacerCol="spacer"
gs <- crisprDesignAux::DataFrameToGuideSet(data,
                                           spacerCol=spacerCol,
                                           bsgenome=bsgenome,
                                           crisprNuclease=crisprNuclease)
gs$ensembl_id <- data$ensembl_id
gs$name <- paste0(gs$ensembl_id, "_", gs$protospacer)
gs$percentCDS <- data$percentCDS
gs$percentCodingIsoforms <- data$percentCodingIsoforms
gs$isCommonCodingExon <- data$isCommonCodingExon
saveRDS(gs,
        file="objects/gs.rds")




load("objects/egs.rda")
load("objects/negs.rda")

gs <- readRDS("objects/gs.rds")
old=gs
gs=old
load("../rankings/processingDatasets/achilles/processed/guideMap.rda")
guideMap$name <- paste0(guideMap$ensembl_id,"_", guideMap$spacer)
gs <- gs[gs$name %in% guideMap$name]
wh <- match(gs$name, guideMap$name)
gs$lfc <- guideMap$lfc[wh]



# Let's add gene symbol:
key <- crisprDesign:::.getTx2GeneTable(txdb_human)
wh <- match(gs$ensembl_id, key$gene_id)
gs$gene_symbol <- key$gene_symbol[wh]





# Sensitivity analysis:
ks <- seq(1,50,1)
cors <- c()
for (k in ks){
    temp <- gs[gs$gene_symbol %in% egs]
    temp <- addConservationScores(temp,
                                  nucExtension=k,
                                  fun="mean",
                                  conservationFile=conservationFile)
    cors[k] <- cor(temp$score_conservation, temp$lfc)
    print(k)
}




# Adding final conservation scores:
gs <- addConservationScores(gs,
                            nucExtension=9,
                            fun="mean",
                            conservationFile=conservationFile)
gs$binary <- ifelse(gs$score_conservation<=0,0,1)













pdf("figures/conservation.pdf", width=7, height=2.2)
par(mfrow=c(1,3), mar=c(2,2,2,2))

plot(ks*2,cors,
     xlab="",
     ylab="",
     yaxt="n",
     bty="L",
     pch=20,
     ylim=c(-0.2, -0.1))
abline(v=18, lty=3)
axis(side=2,
     at=c(-0.2, -0.15,-0.1))

final <- gs[gs$gene_symbol %in% egs]
xs <- split(final$lfc, f=final$binary)
lwd=1.3
plot(density(xs[[1]]),
     main="",
     yaxt="n",
     xaxt="n",
     lwd=lwd,
     ylim=c(0,1.3),
     xlim=c(-2.5,1.5),
     xlab="",
     ylab="",
     bty="L", col=1)
axis(side=2, at=c(0,0.5, 1))
axis(side=1, at=c(-2,-1, 0, 1))
lines(density(xs[[2]]),col="red",lwd=lwd)
abline(v=0, lty=3)
legend("topleft", col=c(1, "red"), bty="n",
       cex=0.75,
       lty=1, c("Low conservation", "High conservation"))

final <- gs[gs$gene_symbol %in% egs]
final <- gs[gs$gene_symbol %in% negs]
xs <- split(final$lfc, f=final$binary)
lwd=1.3
plot(density(xs[[1]]),
     main="",
     yaxt="n",
     xaxt="n",
     ylim=c(0,3),
     xlim=c(-2.5,1.5),
     lwd=lwd,
     xlab="",
     ylab="",
     bty="L", col=1)
axis(side=2, at=c(0,1, 2))
axis(side=1, at=c(-2,-1, 0, 1))
lines(density(xs[[2]]),col="red",lwd=lwd)
abline(v=0, lty=3)
legend("topleft", col=c(1, "red"), bty="n",
       cex=0.75,
       lty=1, c("Low conservation", "High conservation"))
dev.off()




# Now let's look at the %cds
gs <- readRDS("objects/gs.rds")
old=gs
gs=old
load("../rankings/processingDatasets/toronto/processed/guideMap.rda")
gs <- gs[gs$name %in% guideMap$name]
wh <- match(gs$name, guideMap$name)
gs$lfc <- guideMap$lfc[wh]


# Let's add gene symbol:
key <- crisprDesign:::.getTx2GeneTable(txdb_human)
wh <- match(gs$ensembl_id, key$gene_id)
gs$gene_symbol <- key$gene_symbol[wh]



lowessPlot <- function(final){
    x <- final$percentCDS
    y <- final$lfc
    good <- !is.na(x)&!is.na(y)
    x <- x[good]
    y <- y[good]
    plot(x,y,
         xlab="",
         ylab="",
         main="",
         bty="L",
         pch=20,
         cex=0.1,
         col="grey80",
         ylim=c(-4,2))
    abline(h=0)
    abline(h=-1, lty=3)
    a <- lowess(y~x, f=0.1)
    lines(a, col="red", lwd=2)
}


pdf("figures/cds.pdf", width=7, height=3)
par(mfrow=c(1,2),mar=c(2,2,2,2))
lowessPlot(gs[gs$gene_symbol %in% egs])
abline(v=85)
lowessPlot(gs[gs$gene_symbol %in% negs])
dev.off()





