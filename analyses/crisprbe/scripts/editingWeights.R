library(crisprBase)
data(BE4max)
baseEditor=BE4max
lwd=1.7

myPlotEditingWeights <- function(baseEditor,
                                 lwd=1.7
){
    choices <- crisprBase:::.getComboNames()  
    substitutions <- choices
    ws <- editingWeights(baseEditor,
                         substitutions=substitutions)
    ws <- crisprBase:::.getReducedEditingMatrix(ws)
    x <- as.numeric(colnames(ws))
    top <- max(ws, na.rm = TRUE)
    ylim <- c(0, top)
    plot(x, ws[1, ],
         xaxt="n",
         yaxt="n",
         col="white",
         xlab="",
         ylab="",
        ylim = ylim)
    ns <- nrow(ws)
    col <- c("gold3", 
             "2",
             "burlywood3",
             "darkolivegreen3",
             "cadetblue4")
    #col <- crisprBase:::.subColors()[seq_len(ns)]
    for (k in seq_len(ns)) {
        lines(x, ws[k, ], col = col[k], lwd=lwd)
    }
    legend("topleft",
           legend = rownames(ws),
           lwd=lwd,
           bty="n",
           col = col, lty = 1, 
           cex = 0.6)
    axis(side=1, at=c(0,-5,-10,-15,-20,-25,-30, -35), cex.axis=0.7)
    axis(side=2, at=c(0,0.2,0.4,0.6,0.8,1), cex.axis=0.7)
}

pdf("../figures/plotEditingWeights.pdf",
    height=3.5, width=3.5)
myPlotEditingWeights(BE4max)
dev.off()











