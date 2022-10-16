load("../objects/results_kras.rda")
load("../objects/results_znf.rda")
load("../objects/results_egfr.rda")


labels <- c("bowtie", "bowtie-iter",
            "bwa", "bwa-iter")
myMatplot <- function(results,
                      ylim=c(0,320),
                      main="",
                      legend.cex=0.7
){

    colors <- c(rep("deepskyblue2",2),
            rep("firebrick2",2))
    lty <- c(3,1,3,1)
    ylab <- "Time (seconds)"
    xlab <- "Number of mismatches"
    matplot(results,
            main=main,
            ylab=ylab,
            xlab=xlab,
            ylim=ylim,
            xaxt="n",
            pch=20,
            cex=0.5,,
            bty="L",
            col=colors)
    for (k in 1:ncol(results)){
        x <- as.integer(rownames(results))+1
        y <- results[,k]
        good <- !is.na(y)
        lines(x[good],y[good],
              lwd=1.2,
              col=colors[k],
              lty=lty[k])
    }
    axis(side=1,
         at=c(1,2,3,4,5,6),
         labels= c(0,1,2,3,4,5))
    legend("topleft",
           col=colors,
           lty=lty,
           cex=legend.cex,
           legend=labels,
           bty="n")
}

pdf("../figures/offtarget_comparison.pdf", width=7, height=3)
par(mfrow=c(1,3))
ylim=c(0,130)
myMatplot(results_kras, main="KRAS", ylim=ylim)
ylim=c(0,700)
myMatplot(results_egfr, main="EGFR",ylim=ylim)
ylim=c(0,350)
myMatplot(results_znf, main="ZNF101",ylim=ylim)
dev.off()


write.csv(results_kras, file="results_kras.csv")
write.csv(results_egfr, file="results_egfr.csv")
write.csv(results_znf, file="results_znf.csv")

