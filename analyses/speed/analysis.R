# Getting data
load("flashfry/resultsMm2.rda")
data <- list(flash=results)

load("crisprDesign/resultsBowtieIter.rda")
resultsBowtieIter <- lapply(resultsBowtieIter, matrixStats::colMedians)
data$bowtiei <- resultsBowtieIter

load("crisprDesign/resultsBwaIter.rda")
data$bwai <- resultsBwaIter

load("chopchop/results.rda")
data$chopchop <- results

load("multicrispr/results.rda")
results <- lapply(results, function(x){as.numeric(x[1,])})
data$multicrispr <- results

load("cctop/results.rda")
data$cctop <- results


data <- lapply(data, function(x){
    vapply(x, function(y) y[3], FUN.VALUE=1)
})
data <- do.call(rbind, data)

# Getting x-axis:
load("inputs/inputGrs.rda")
ns <- vapply(inputGrs, function(gr){
    sum(BiocGenerics::width(gr))
}, FUN.VALUE=1)


rownames(data) <- c("FlashFry",
                    "crisprDesign-bowtie",
                    "crisprDesign-bwa",
                    "CHOPCHOP",
                    "multicrispr",
                    "CCTop")
data <- data[6:1,]

cols <- 1:6
cols <- c("orange","olivedrab4",
          "grey50","firebrick2",
          "deepskyblue2", "black")



pdf("figures/speed.pdf", width=5, height=5)
plot(ns,data[1,],
     las=2,
     xlab="",
     ylab="",
     type="l",
     ylim=c(0,20000),
     bty="L", col="white")
for (i in 1:6){
    lines(ns,data[i,],
          lwd=1.2,
          type="l",
          col=cols[i])    
    points(ns,data[i,],
           col=cols[i],
           pch=3,
           cex=0.7)
}
legend("topright",
       bty="n",
       cex=0.75,
       lty=1,
       col=cols,
       legend=rownames(data))
dev.off()

prepareData <- function(){
    out <- data
    colnames(out) <- ns
    out
}

write.csv(prepareData(),file="timings.csv")


