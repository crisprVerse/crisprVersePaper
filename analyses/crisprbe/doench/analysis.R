library(crisprDesign)
library(crisprDesignData)
library(SummarizedExperiment)
library(BSgenome.Hsapiens.UCSC.hg38)
bsgenome <- BSgenome.Hsapiens.UCSC.hg38
txObject <- txdb_human
bowtie_index="/Users/fortinj2/crisprIndices/bowtie/hg38/hg38"
load("../objects/guideSet.rda")
load("objects/se.rda")


# Removing guides with off-targets
guideSet <- addSpacerAlignments(guideSet,
                                txObject=txObject,
                                n_mismatches=2,
                                aligner_index=bowtie_index,
                                bsgenome=bsgenome)
good <- guideSet$n0==1 & guideSet$n1_c==0 & guideSet$n1_c==0
guideSet <- guideSet[good]



# Getting log-fold change:
Y <- assays(se)[[1]]


# Normalization
factors <- colSums(Y)
factors <- factors/median(factors)
Y <- sweep(Y, 2, factors, "/")
assays(se)[[1]] <- Y


# Filtering
Y <- log2(Y+1)
pdna <- Y[,1]
bad <- which(pdna<=(mean(pdna)-3*sd(pdna)))
se <- se[-bad,]


gs <- guideSet
spacers <- spacers(gs, as.character=TRUE)
common <- intersect(rownames(se), spacers)
wh <- match(common, spacers)
gs <- gs[wh]
se <- se[common,]
ann <- rowData(se)
muts <- strsplit(ann$mutation, split=";")
muts <- lapply(muts, unique)
muts <- vapply(muts, function(x){
    if ("Nonsense" %in% x){
        x <- "Nonsense"
    } else if("Missense" %in% x){
        x <- "Missense"
    } else if ("Intron" %in% x){
        x <- "Intron"
    }
    x
}, FUN.VALUE="a")
gs$mutation_doench <- muts



Y <- assays(se)[[1]]
Y <- log2(Y+1)
lfc <- rowMeans(Y[,2:4]-Y[,1])
mcols(gs)$lfc <- lfc

# Ordering
o <- order(-start(gs))
gs <- gs[o]


lfc <- split(gs$lfc, f=gs$maxVariant)
#boxplot(lfc)



pdf("../figures/doench.pdf", width=3.8, height=4)
library(pROC)
col <- 1:4
col <- c("orange","firebrick2", "firebrick3", "firebrick")
cutoffs <- c(0.2,0.4,0.6,0.8)
results <- lapply(cutoffs, function(cutoff){
    y <- ifelse(gs$score_nonsense>=cutoff,1,0)    
    roc(y~gs$lfc)
})
xlab="1 - specificity"
ylab="Sensitivity"
plot(1-results[[1]]$specificities,
     results[[1]]$sensitivities,
     xlab=xlab,
     ylab=ylab,
     col="white",
     type="l")
for (i in 1:4){
    lines(1-results[[i]]$specificities,
     results[[i]]$sensitivities,col=col[i])
}
abline(a=0,b=1, lty=3)
legend("bottomright",
       title="Nonsense score cutoff",
       c("0.2", "0.4", "0.6", "0.8"),
       col=col, lty=1, bty="n", cex=0.7)
dev.off()








gs <- addOnTargetScores(gs,
                        methods="azimuth")
gs <- addOnTargetScores(gs,
                        methods="ruleset1")



col_nonsense <- 2
col_other <- "grey95"
x <- gs$score_ruleset1
y <- gs$lfc



lfc_cutoff <- -0.5
myPlot <- function(x,
                   score_cutoff=0.5,
                   xlim=c(0,1)
){
    y <- gs$lfc
    plot(x,y,
         xlim=xlim,
         xlab="",
         ylab="",
         bty="L",
         col="white",
         xaxt="n",
         pch=20)
    wh <- which(gs$maxVariant!="nonsense")
    points(x[wh], y[wh],
           cex=gs$maxVariantScore[wh],
           col=col_other, pch=20)
    wh <- which(gs$maxVariant=="nonsense")
    cond1 <- x[wh]<score_cutoff & y[wh]>lfc_cutoff
    cond2 <- x[wh]>score_cutoff & y[wh]<lfc_cutoff
    cond3 <- x[wh]<score_cutoff & y[wh]<=lfc_cutoff
    cond4 <- x[wh]>score_cutoff & y[wh]>=lfc_cutoff
    goodCond <- cond1 | cond2
    badCond <- cond3 | cond4
    #wh1 <- which(x[wh]<score_cutoff)
    #points(x[wh], y[wh],
    #       cex=gs$score_nonsense[wh]*1.5,
    #       col=col_nonsense, pch=20)
    points(x[wh][goodCond], y[wh][goodCond],
           cex=gs$score_nonsense[wh][goodCond]*1.5,
           col="2", pch=20)
    points(x[wh][badCond], y[wh][badCond],
           cex=gs$score_nonsense[wh][badCond]*1.5,
           col="cadetblue4", pch=20)
    abline(v=score_cutoff, lty=3)
    abline(h=lfc_cutoff, lty=3)
    #abline(h=0, lty=1)
    lines(lowess(y[wh]~x[wh], f=0.7), col="grey55", lty=3)
}
pdf("../figures/doench_scores.pdf", width=5, height=2.8)
par(mfrow=c(1,3), mar=c(4,2,4,1))
myPlot(gs$score_ruleset1, xlim=c(0,0.8), score_cutoff=0.1)
axis(side=1, at=c(0,0.4, 0.8))
myPlot(gs$score_azimuth, xlim=c(0.3,0.7), score_cutoff=0.43)
axis(side=1, at=c(0.3,0.5, 0.7))
myPlot(gs$score_deephf, xlim=c(0.2, 0.8), score_cutoff=0.54)
axis(side=1, at=c(0.2,0.4, 0.6, 0.8))
dev.off()
method="pearson"
cor(gs$score_ruleset1, y, method=method)
cor(gs$score_azimuth, y, method=method)
cor(gs$score_deephf, y, method=method)



prepareData <- function(){
    grna <- as.character(spacers(gs))
    out <- data.frame(grna=grna)
    out$score_deephf <- gs$score_deephf
    out$score_ruleset1 <- gs$score_ruleset1
    out$score_azimuth <- gs$score_azimuth
    out$lfc <- gs$lfc
    out$variant <- gs$maxVariant
    rownames(out) <- NULL
    out
}

data <- prepareData()
write.csv(data, file="scoringData.csv", row.names=FALSE)

