library(crisprBase)
library(crisprDesign)
library(crisprDesignData)
library(BSgenome.Hsapiens.UCSC.hg38)
data(BE4max, package="crisprBase")
bsgenome <- BSgenome.Hsapiens.UCSC.hg38
txid="ENST00000357654"

editingWindow=c(-20,-8)

txTable <- getTxInfoDataFrame(tx_id=txid,
                              txObject=txdb_human,
                              bsgenome=bsgenome)
gr <- queryTxObject(txdb_human,
                    queryValue=txid,
                    queryColumn="tx_id",
                    featureType="cds")
gs <- findSpacers(gr,
                  bsgenome=bsgenome,
                  crisprNuclease=BE4max)
gs <- unique(gs)
#gs <- gs[100:110]
gs <- addEditedAlleles(gs,
                       baseEditor=BE4max,
                       txTable=txTable,
                       editingWindow=editingWindow)

gs <- addOnTargetScores(gs,
                        methods="deephf")

guideSet <- gs
save(guideSet,
     file="../objects/guideSet.rda")



plotScore <- function(gs){
     xlab <- ""
     ylab <- ""
     y <- gs$score_nonsense
     z <- gs$score_deephf
     pam <- pamSites(gs)
     x <- txTable$aa_number[match(pam, txTable$pos)]
     good <- !is.na(x) & gs$maxVariant!="not_targeting"

     col <- gs$maxVariant
     label <- col
     col_silent   <- "grey75"
     col_silent   <- "grey75"
     col_missense <- "deepskyblue2"
     col_missense <- "cadetblue4"
     col_nonsense <- "firebrick3"
     col_nonsense <- "2"
     col[col=="silent"] <- col_silent
     col[col=="nonsense"] <- col_nonsense
     col[col=="missense"] <- col_missense



     xx <- x[good]
     yy <- y[good]
     zz <- z[good]
     coll <- col[good]
     labell <- label[good]
     cexx=(zz-0.2)*3
     ylim=c(0,1)
     plot(xx,
          yy,
          ylim=ylim,
          yaxt="n",
          xlab=xlab,
          ylab=ylab,
          col="white",
          bty="L",
          cex=cexx, pch=20)
     axis(side=2, at=c(0, 0.2, 0.4, 0.6, 0.8, 1), cex.axis=0.75)
     abline(v=500, lty=3,col="grey75")
     abline(v=1000, lty=3,col="grey75")
     abline(v=1500, lty=3,col="grey75")
     wh_silent <- which(labell=="silent")
     wh_nonsense <- which(labell=="nonsense")
     wh_missense <- which(labell=="missense")
     points(xx[wh_missense],
            yy[wh_missense],
            cex=cexx[wh_missense],
            col=coll[wh_missense],
            pch=20)
     points(xx[wh_nonsense],
            yy[wh_nonsense],
            cex=cexx[wh_nonsense],
            col=coll[wh_nonsense],
            pch=20)
     points(xx[wh_silent],
            yy[wh_silent],
            cex=cexx[wh_silent],
            col=coll[wh_silent],
            pch=20)
}


pdf("../figures/scores.pdf", width=5, height=3.5)
plotScore(gs)
dev.off()


prepareData <- function(){
     grna <- as.character(spacers(gs))
     y <- gs$score_nonsense
     z <- gs$score_deephf
     pam <- pamSites(gs)
     x <- txTable$aa_number[match(pam, txTable$pos)]
     good <- !is.na(x) & gs$maxVariant!="not_targeting"
     label <- gs$maxVariant
     gg <- grna[good]
     xx <- x[good]
     yy <- y[good]
     zz <- z[good]
     labell <- label[good]
     out <- data.frame(grna=gg,
                       aminoAcid=xx,
                       scoreDeepHF=zz,
                       scoreNonsense=yy,
                       label=labell)
     return(out)   
}

aaData <- prepareData()
write.csv(aaData,
          file="nonsenseData.csv",
          row.names=FALSE)


### Looking at one example:
cond1 <- gs$score_nonsense<0.4
cond2 <- gs$score_nonsense>0.3
cond3 <- gs$maxVariant=="missense"
wh=which(cond1 & cond2 & cond3) #140
x=gs$editedAlleles[[wh]]
meta <- metadata(x)
seq <- getSeq(bsgenome, names=meta$chr,
              start=meta$start-10, end=meta$end+10)

