library(crisprDesignData)
library(S4Vectors)
library(dplyr)
library(pbapply)
library(scales)
col_best  <- "firebrick3"
col_worst <- "grey30"

load("processingDatasets/objectsFinal/rankings_achilles.rda")
load("processingDatasets/objectsFinal/rankings_sabatini.rda")
load("processingDatasets/objectsFinal/rankings_toronto.rda")
load("processingDatasets/objectsFinal/rankings_toronto3.rda")
load("processingDatasets/objectsFinal/rankings_yusa.rda")
rankings <- list(achilles=rankings_achilles,
                 sabatini=rankings_sabatini,
                 toronto=rankings_toronto,
                 toronto3=rankings_toronto3,
                 yusa=rankings_yusa)
load("processingDatasets/objectsFinal/rankings_achilles_neg.rda")
load("processingDatasets/objectsFinal/rankings_sabatini_neg.rda")
load("processingDatasets/objectsFinal/rankings_toronto_neg.rda")
load("processingDatasets/objectsFinal/rankings_toronto3_neg.rda")
load("processingDatasets/objectsFinal/rankings_yusa_neg.rda")
rankings_neg <- list(achilles=rankings_achilles_neg,
                     sabatini=rankings_sabatini_neg,
                     toronto=rankings_toronto_neg,
                     toronto3=rankings_toronto3_neg,
                     yusa=rankings_yusa_neg)
rankings <- lapply(rankings, function(x){
    x$rank_crisprverse <- x$rank_crisprdesign_genic
    x
})
rankings_neg <- lapply(rankings_neg, function(x){
    x$rank_crisprverse <- x$rank_crisprdesign_genic
    x
})
n_datasets <- length(rankings)

# Some filtering on how-targets
rankings <- lapply(rankings, function(x){
    x <- x[x$n0==1,]
    x <- x[x$n1_c==0,]
    x <- x[x$n2_c==0,]
    x
})
rankings_neg <- lapply(rankings_neg, function(x){
    x <- x[x$n0==1,]
    x <- x[x$n1_c==0,]
    x <- x[x$n2_c==0,]
    x
})

allCols <- c("rank_cctop", 
             "rank_flash", 
             "rank_cpick", 
             "rank_chopchop", 
             "rank_crisprverse")

cutoff=15
results <- lapply(seq_len(n_datasets), function(k){
    df <- rankings[[k]]
    deltas <- vapply(allCols, function(col){
        x=df$lfc[which(df[[col]]<=cutoff)]
        z=df$lfc[which(df[[col]]>=cutoff)]
        return(mean(x)-mean(z))
    }, FUN.VALUE=1)
    deltas
})
results_neg <- lapply(seq_len(n_datasets), function(k){
    df <- rankings_neg[[k]]
    deltas <- vapply(allCols, function(col){
        x=df$lfc[which(df[[col]]<=cutoff)]
        z=df$lfc[which(df[[col]]>=cutoff)]
        return(mean(x)-mean(z))
    }, FUN.VALUE=1)
    deltas
})
pvals <- lapply(seq_len(n_datasets), function(k){
    df <- rankings[[k]]
    pvals <- vapply(allCols, function(col){
        x=df$lfc[which(df[[col]]<=cutoff)]
        z=df$lfc[which(df[[col]]>=cutoff)]
        t.test(x,z)$p.value
    }, FUN.VALUE=1)
    pvals
})
pvals_neg <- lapply(seq_len(n_datasets), function(k){
    df <- rankings_neg[[k]]
    pvals <- vapply(allCols, function(col){
        x=df$lfc[which(df[[col]]<=cutoff)]
        z=df$lfc[which(df[[col]]>=cutoff)]
        t.test(x,z)$p.value
    }, FUN.VALUE=1)
    pvals
})

results <- do.call(rbind, results)
results_neg <- do.call(rbind, results_neg)
pvals <- do.call(rbind, pvals)
pvals_neg <- do.call(rbind, pvals_neg)

# Ordering
meds <- colMeans(results)
results <- results[,order(meds)]
results_neg <- results_neg[,order(meds)]
rownames(results) <- names(rankings)
rownames(results_neg) <- names(rankings_neg)
rownames(pvals) <- names(rankings)
rownames(pvals_neg) <- names(rankings_neg)

o <- order(-colMeans(results))
results <- results[,o]
results_neg <- results_neg[,o]
o <- order(-rowMeans(results))
results <- results[o,]
results_neg <- results_neg[o,]
pvals <- pvals[rownames(results), colnames(results)]
pvals_neg <- pvals_neg[rownames(results_neg), colnames(results_neg)]


getPaperPalette <- function(){
    blue1   <- c(64, 116, 177)
    blue2 <- c(194,214, 236)
    green1  <- c(94, 128, 63)
    green2  <- c(197, 213, 175)
    yellow   <- c(249, 217, 120)
    red <- c(235, 50, 35)
    cols <- list(blue1, blue2, green1, green2, yellow)
    cols <- list(red, blue2, green1, green2, yellow)
    cols <- rev(cols)
    cols <- vapply(cols, function(col){
        rgb(col[1], col[2], col[3], maxColorValue=255)
    }, FUN.VALUE="a")
    cols[5] <- alpha(cols[5], 0.8)
    return(cols)
}
colors <- getPaperPalette()


labels <- c("sabatini",
            "yusa",
            "achilles",
            "toronto3",
            "toronto")
results <- results[labels,]
results_neg <- results_neg[labels,]
newLabels <- c(sabatini="Wang2015",
               achilles="Achilles",
               toronto="Hart2015",
               toronto3="Hart2017",
               yusa="Tzelepis2016")
newLabels <- newLabels[match(rownames(results), names(newLabels))]
rownames(results) <- newLabels
newLabels <- c(sabatini="Wang2015",
               achilles="Achilles",
               toronto="Hart2015",
               toronto3="Hart2017",
               yusa="Tzelepis2016")
newLabels <- newLabels[match(rownames(results_neg), names(newLabels))]
rownames(results_neg) <- newLabels



mets <- colnames(results)
mets <- gsub("rank_","", mets)
newMets <- c(cctop="CCTop",
             chopchop="CHOPCHOP",
             flash="FlashFry",
             cpick="CRISPick",
             crisprverse="crisprDesign")
mets <- newMets[match(mets, names(newMets))]             
colnames(results) <- mets
mets <- colnames(results_neg)
mets <- gsub("rank_","", mets)
newMets <- c(cctop="CCTop",
             chopchop="CHOPCHOP",
             flash="FlashFry",
             cpick="CRISPick",
             crisprverse="crisprDesign")
mets <- newMets[match(mets, names(newMets))]             
colnames(results_neg) <- mets

pdf("figures/rankings_barplots.pdf", width=5, height=3.5)
par(xaxt="n")
barplot(-(t(results)),
        xaxt="n",
        ylim=c(0,1),
        cex.axis=0.8,
        border="white",
        beside=TRUE, las=1, col=colors)
legend("topleft",
       cex=0.7,
       fill=colors,
       border=NA,
       legend=mets,
       bty="n")
dev.off()





pdf("figures/rankings_neg_barplots.pdf", width=5, height=3.5)
par(xaxt="n")
barplot(-(t(results_neg)),
        xaxt="n",
        ylim=c(0,1),
        cex.axis=0.8,
        border="white",
        beside=TRUE, las=1, col=colors)
legend("topleft",
       cex=0.7,
       fill=colors,
       border=NA,
       legend=mets,
       bty="n")
dev.off()


# Now we will illustrate with essential genes
cols <- c("rank_cctop",
          "rank_chopchop",
          "rank_flash",
          "rank_cpick",
          "rank_crisprverse")
names <- c("CCTop",
           "CHOPCHOP",
           "FlashFry",
           "CRISPick",
           "crisprDesign")


deltas <- c()
pvals  <- c()
densPlot <- function(df,
                     name,
                     col=NULL,
                     lwd=1.5,
                     ylim=c(0,1),
                     ...
){
    col_worst <- getPaperPalette()[3]
    col_best <- getPaperPalette()[5]
    xtop <- df$lfc[which(df[[col]]<=cutoff)]
    xbot <- df$lfc[which(df[[col]]>cutoff)]
    plot(density(xbot, na.rm=TRUE),
         bty="L",
         ylim=ylim,
         main=name,
         col=col_worst,
         lwd=lwd,
         yaxt="n",
         ...)
    axis(side=2, at=c(0,0.25,0.5))
    lines(density(xtop, na.rm=TRUE), col=col_best, lwd=lwd)
    abline(v=0, lty=3)
}
densPlot2 <- function(df,
                      df2,
                      name,
                      col=NULL,
                      lwd=1.5,
                      ylim=c(0,1),
                      ...
){
    col_worst <- getPaperPalette()[3]
    col_best <- getPaperPalette()[5]
    xtop <- df$lfc[which(df[[col]]<=cutoff)]
    xbot <- df$lfc[which(df[[col]]>cutoff)]
    xtop2 <- df2$lfc[which(df2[[col]]<=cutoff)]
    xbot2 <- df2$lfc[which(df2[[col]]>cutoff)]
    plot(density(xbot, na.rm=TRUE),
         bty="L",
         ylim=ylim,
         main=name,
         col=col_worst,
         lwd=lwd,
         yaxt="n",
         ...)
    axis(side=2, at=c(0,0.25,0.5))
    lines(density(xtop2, na.rm=TRUE), col=col_best, lwd=1.4, lty=3)
    lines(density(xbot2, na.rm=TRUE), col=col_worst, lwd=1.4, lty=3)
    lines(density(xtop, na.rm=TRUE), col=col_best, lwd=lwd)
    abline(v=0, lty=3)
}




pdf("figures/rankings_densities_both.pdf", width=8, height=1.7)
par(mfrow=c(1,5), mar=c(2,2,2,1))
final <- rankings[["toronto"]]
final_neg <- rankings_neg[["toronto"]]
for (i in 1:5){
    densPlot2(final,
              final_neg,
              col=cols[i],
             name="",
             lwd=2.5,
             ylim=c(0,0.6),
             xlim=c(-7,2))
}
dev.off()










