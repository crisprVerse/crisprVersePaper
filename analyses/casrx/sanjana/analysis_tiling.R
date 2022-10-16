library(SummarizedExperiment)
library(crisprDesign)
library(crisprDesignData)
library(RColorBrewer)
load("objects/ses_tiling.rda")
txObject <- txdb_human
key <- crisprDesign:::.getTx2GeneTable(txObject)
bowtie_index="/Users/fortinj2/crisprIndices/bowtie/ensembl_human_104/ensembl_human_104"
#CD46: ENST00000367042
#CD55: ENST00000367064
#CD71: ENST00000360110





colors <- brewer.pal(9,"YlOrRd")


# Calculate log-fold change
getLogFC <- function(se){
    se <- se[,order(colData(se)$Condition, colData(se)$Replicate)]
    pheno <- colData(se)
    Y <- log2(assays(se)[[1]]+1)
    Y1 <- Y[, pheno$Condition=="top"]
    Y2 <- Y[, pheno$Condition=="input"]
    rowMeans(Y1-Y2)
}



completeAnnotation <- function(k){
    gene <- names(ses_tiling)[k]
    se <- ses_tiling[[k]]
    if (gene=="CD71"){
        gene <- "TFRC"
    }

    # Only perfect match:
    se <- se[rowData(se)$class=="PerfectMatch",]
    ann <- rowData(se)
    keyGene <- key[key$gene_symbol==gene,]

    # Let's load guides annotation:
    load("processing/objects/guidesTiling.rda")
    names(guidesTiling) <- sapply(guidesTiling, function(x){
        as.character(seqnames(x)[1])
    })
    guidesTiling <- guidesTiling[keyGene$tx_id]
    spacers <- lapply(guidesTiling, spacers, as.character=TRUE)
    overlaps <- sapply(spacers, function(x){
        length(intersect(x, ann$spacer))
    })
    ns <- sapply(spacers, length)
    cbind(ns, overlaps)
    if (k==1){
        guides <- guidesTiling[["ENST00000367042"]]
    } else if (k==2){
        guides <- guidesTiling[["ENST00000367064"]]
    } else{
        guides <- guidesTiling[["ENST00000360110"]]
    }


    guides$spacer <- spacers(guides,as.character=TRUE)
    common <- intersect(guides$spacer, ann$spacer)
    guides <- guides[match(common, guides$spacer)]
    see <- se[match(common, rowData(se)$spacer),]
    guides$lfc <- getLogFC(see)


    aln <- addSpacerAlignments(guides,
                                addSummary=TRUE,
                                txObject=txObject,
                                n_mismatches=0,
                                aligner="bowtie",
                                aligner_index=bowtie_index)
    guides$ntx <- aln$n0_tx
    return(guides)
}

guidesCD46 <- completeAnnotation(1)
guidesCD55 <- completeAnnotation(2)
save(guidesCD46, file="objects/guidesCD48.rda")
save(guidesCD55, file="objects/guidesCD55.rda")


### Let's add DeepCas13 score:
deepcas13 <- read.csv("objects/deepcas13d_cd46.csv")
deepcas13 <- deepcas13[,c(2,4)]
wh <- match(guidesCD46$spacer, deepcas13[,1])
guidesCD46$deepcas13 <- deepcas13[,2][wh]


### Let's add DeepCas13 score:
deepcas13 <- read.csv("objects/deepcas13d_cd55.csv")
deepcas13 <- deepcas13[,c(2,4)]
wh <- match(guidesCD55$spacer, deepcas13[,1])
guidesCD55$deepcas13 <- deepcas13[,2][wh]
save(guidesCD46, file="objects/guidesCD48.rda")
save(guidesCD55, file="objects/guidesCD55.rda")



txs <- list(CD46="ENST00000367042",
            CD55="ENST00000367064",
            TFRC="ENST00000360110")

#Let's get UTRs and CDS coordinates
cds <- sapply(txs, function(txid){
    x <- txObject$cds
    x <- x[x$tx_id==txid]
    sum(width(x))
})
utr5 <- sapply(txs, function(txid){
    x <- txObject$fiveUTRs
    x <- x[x$tx_id==txid]
    sum(width(x))
})
utr3 <- sapply(txs, function(txid){
    x <- txObject$threeUTRs
    x <- x[x$tx_id==txid]
    sum(width(x))
})

metadata(guidesCD46)$cds_start <- utr5[["CD46"]]
metadata(guidesCD46)$cds_end <- utr5[["CD46"]]+cds[["CD46"]]
metadata(guidesCD55)$cds_start <- utr5[["CD55"]]
metadata(guidesCD55)$cds_end <- utr5[["CD55"]]+cds[["CD55"]]
save(guidesCD46, file="objects/guidesCD48.rda")
save(guidesCD55, file="objects/guidesCD55.rda")


######## Adding RNA-Seq  ########
#Let's get RA-Seq
library(readr)
rna <- read_tsv("data/transcript_rna_celline.tsv")
rna <- as.data.frame(rna)
wh <- which(grepl("HEK",colnames(rna)))
rna <- rna[,c(1:2, wh)]
wh <- grepl("TPM", colnames(rna))
rna$TPM <- rowMeans(rna[, wh])
colnames(rna)[1:2] <- c("geneid", 'txid')
rna <- rna[,c("txid","TPM")]
rownames(rna) <- NULL




#Let's get the list of spacers:
load("processing/objects/guidesTiling.rda")
names(guidesTiling) <- sapply(guidesTiling, function(x){
    as.character(seqnames(x)[1])
})


# CD46
key <- crisprDesign:::.getTx2GeneTable(txObject)
txids <- key[key$gene_symbol=="CD46","tx_id"]
guidesTiling <- guidesTiling[txids]
spacerList <- lapply(guidesTiling, spacers, as.character=TRUE)
spacers <- spacers(guidesCD46, as.character=TRUE)
member <- lapply(spacerList, function(x){
    spacers %in% x
})
member <- do.call(cbind, member)

temp <- rna[rna$txid %in% txids,]
temp <- temp[match(txids, temp$txid),]
guidesCD46$tpm <- apply(member, 1, function(x){
    sum(temp$TPM[x])
})




# CD55
load("processing/objects/guidesTiling.rda")
names(guidesTiling) <- sapply(guidesTiling, function(x){
    as.character(seqnames(x)[1])
})

key <- crisprDesign:::.getTx2GeneTable(txObject)
txids <- key[key$gene_symbol=="CD55","tx_id"]
guidesTiling <- guidesTiling[txids]
spacerList <- lapply(guidesTiling, spacers, as.character=TRUE)
spacers <- spacers(guidesCD55, as.character=TRUE)
member <- lapply(spacerList, function(x){
    spacers %in% x
})
member <- do.call(cbind, member)

temp <- rna[rna$txid %in% txids,]
temp <- temp[match(txids, temp$txid),]
guidesCD55$tpm <- apply(member, 1, function(x){
    sum(temp$TPM[x])
})

save(guidesCD46, file="objects/guidesCD48.rda")
save(guidesCD55, file="objects/guidesCD55.rda")



load("objects/guidesCD48.rda")
load("objects/guidesCD55.rda")









library(scales)
pdf("figures/cd46_tpm.pdf", width=10, height=1.7)
par(mar=c(4,2,2,2))
x <- pamSites(guidesCD46)
y <- guidesCD46$tpm
fit <- loess(y~x, span=0.01)
xx <- c(x, rev(x))
y <- fitted(fit)
yy <- c(rep(-10, length(y)), rev(y))
col <- alpha("firebrick2",0.1)
plot(x, y, type="l",
     col="white",
     lwd=2,
     bty="L",
     yaxt="n",
     xlab="",
     ylab="")
polygon(xx,yy,
        border=NA,
        col=col)
axis(side=2, at=c(35,60))
dev.off()

pdf("figures/cd55_tpm.pdf", width=10, height=1.7)
par(mar=c(4,2,2,2))
x <- pamSites(guidesCD55)
y <- guidesCD55$tpm
fit <- loess(y~x, span=0.01)
xx <- c(x, rev(x))
y <- fitted(fit)
yy <- c(rep(-10, length(y)), rev(y))
col <- alpha("firebrick2",0.1)
plot(x, y, type="l",
     ylim=c(5,max(y)),
     col="white",
     lwd=2,
     bty="L",
     yaxt="n",
     xlab="",
     ylab="")
polygon(xx,yy,
        border=NA,
        col=col)
axis(side=2, at=c(5,20))
dev.off()




load("objects/guidesCD48.rda")
load("objects/guidesCD55.rda")


### Getting cs files:
genes <- c("CD46", "CD55", "CD71")
cs_files <- paste0(genes, "_screen_crRNA_enrichments.csv")
cs <- lapply(cs_files, function(file){
    temp <- read.csv(file.path("data", file))
    temp <- as.data.frame(temp)
    return(temp)
})
names(cs) <- genes


ann <- rowData(ses_tiling[["CD46"]])
wh <- match(guidesCD46$spacer, ann$spacer)
guidesCD46$ID <- ann$ID[wh]
wh <- match(guidesCD46$ID, rownames(cs[["CD46"]]))
guidesCD46$cs <- cs[["CD46"]][wh,"meanCS.Top"]
guidesCD46 <- guidesCD46[!is.na(guidesCD46$cs)]

ann <- rowData(ses_tiling[["CD55"]])
wh <- match(guidesCD55$spacer, ann$spacer)
guidesCD55$ID <- ann$ID[wh]
wh <- match(guidesCD55$ID, rownames(cs[["CD55"]]))
guidesCD55$cs <- cs[["CD55"]][wh,"meanCS.Top"]
guidesCD55 <- guidesCD55[!is.na(guidesCD55$cs)]


save(guidesCD55, file="objects/guidesCD55.rda")
save(guidesCD46, file="objects/guidesCD48.rda")








load("objects/guidesCD48.rda")
load("objects/guidesCD55.rda")






pdf("figures/cd46_zigzag.pdf", width=10, height=4)
guides <- guidesCD46
gene <- "CD46"
ntx <- guides$ntx
col <- rep("grey55", length(guides))
col[ntx<=4] <- "deepskyblue"
col[ntx==5] <- "deepskyblue1"
col[ntx==6] <- "deepskyblue2"
col[ntx==8] <- "deepskyblue3"
col[ntx==9] <- "orange"
col[ntx==10] <- "firebrick2"
col[ntx>10] <- "firebrick"

x=pamSites(guides)
y=guides$cs
plot(x,
     y,
     col="white", pch=20,
     cex=guides$score_casrxrf,
     #cex=0.7,
     bty="L",
     xlab="",
     ylab="")

xleft=metadata(guides)[["cds_start"]]
xright=metadata(guides)[["cds_end"]]
xx <- seq(xleft, xright,1)
yy1 <- rep(-2.2, length(xx))
yy2 <- rep(20, length(xx))
polygon(c(xx, rev(xx)),
        c(yy1, rev(yy2)),
        border="white",
        col="grey90")
abline(h=0, lty=3)
points(x,
       y,
       col=col, pch=20,
       cex=guides$score_casrxrf)
lines(x,predict(loess(y~x, span=0.05)), lwd=1)
dev.off()

getData1 <- function(){
    guides <- guidesCD46
    grna <- as.character(spacers(guides))
    gene <- "CD46"
    ntx <- guides$ntx
    x=pamSites(guides)
    y=guides$cs
    out <- data.frame(grna=grna,
                      gene=gene,
                      pos=x,
                      lfc=y,
                      ntx=ntx,
                      score_casrxrf=guides$score_casrxrf)
    rownames(out) <- NULL
    out
}

getData2 <- function(){
    guides <- guidesCD55
    grna <- as.character(spacers(guides))
    gene <- "CD55"
    ntx <- guides$ntx
    x=pamSites(guides)
    y=guides$cs
    out <- data.frame(grna=grna,
                      gene=gene,
                      pos=x,
                      lfc=y,
                      ntx=ntx,
                      score_casrxrf=guides$score_casrxrf)
    rownames(out) <- NULL
    out
}

data1 <- getData1()
data2 <- getData2()
toSave <- rbind(data1,data2)
write.csv(toSave, file="tiling.csv", row.names=FALSE)

pdf("figures/cd55_zigzag.pdf", width=10, height=4)
guides <- guidesCD55
gene <- "CD55"
ntx <- guides$ntx
col <- ntx
palette(colors)

x=pamSites(guides)
y=guides$cs
ylim=c(-1, 3)
plot(x,
     y,
     ylim=ylim,
     col="white", pch=20,
     cex=guides$score_casrxrf,
     #cex=0.7,
     bty="L",
     xlab="",
     ylab="")
xleft=metadata(guides)[["cds_start"]]
xright=metadata(guides)[["cds_end"]]


xx <- seq(xleft, xright,1)
yy1 <- rep(-0.95, length(xx))
yy2 <- rep(20, length(xx))
polygon(c(xx, rev(xx)),
        c(yy1, rev(yy2)),
        border="white",
        col="grey95")
abline(h=0, lty=3)
points(x,
       y,
       col=col+1, pch=20,
       cex=guides$score_casrxrf)
lines(x,predict(loess(y~x, span=0.05)), lwd=1)
dev.off()





####################### Selection CD46 #######################
### Only retaining guides within CDS:
pdf("figures/densities_casrx.pdf", width=4, height=4)
guides <- guidesCD46
xcutoff <- metadata(guides)[["cds_end"]]
aa <- guides[pamSites(guides)<=xcutoff]

cutoff=0.5
txcutoff=11
xlim=c(-4,4)
lwd=2
plot(density(aa$cs),
     xlim=xlim,
     xaxt="n",
     yaxt="n",
     bty="L",
     main="",
     xlab="",
     ylab="",
     lwd=lwd,
     ylim=c(0,1))
#abline(v=0, lty=3)
segments(x0=0, x1=0, y0=-10, y1=0.7, lty=3)
wh=which(aa$score_casrxrf>=cutoff)
lines(density(aa$cs[wh]), col="orange", lwd=lwd)
wh=which(aa$ntx>=txcutoff)
lines(density(aa$cs[wh]), col="deepskyblue3",lwd=lwd)
wh=which(aa$score_casrxrf>=cutoff & aa$ntx>=txcutoff)
wh1=wh
lines(density(aa$cs[wh]), col="red", lwd=lwd)
axis(side=1, at=c(-4,-2,0,2,4))
axis(side=2, at=c(0,0.5,1))
legend("topleft", lty=1, bty="n",
       bg="white",
       cex=0.7,
       col=c("black","deepskyblue3","orange", "red"),
       c("No selection",
         "Common exon",
         "High score",
         "High score and common exon"))
dev.off()


prepareDensityData <- function(){
    guides <- guidesCD46
    xcutoff <- metadata(guides)[["cds_end"]]
    aa <- guides[pamSites(guides)<=xcutoff]
    grna <- as.character(spacers(aa))
    txcutoff=11
    cutoff=0.5
    common <- aa$ntx>=txcutoff
    highScore <- aa$score_casrxrf>=cutoff
    out <- data.frame(grna=grna,
                      score_casrxrf=aa$score_casrxrf,
                      lfc=aa$cs,
                      commonExon=common,
                      highScore=highScore)
    rownames(out) <- NULL
    out
}

write.csv(prepareDensityData(),
          file="densitydata.csv",
          row.names=FALSE)


pdf("figures/densities_deepcas13.pdf", width=4, height=4)
guides <- guidesCD46
xcutoff <- metadata(guides)[["cds_end"]]
aa <- guides[pamSites(guides)<=xcutoff]

cutoff=0.5
txcutoff=11
xlim=c(-4,4)
plot(density(aa$cs),
     xlim=xlim,
     xaxt="n",
     yaxt="n",
     bty="L",
     main="",
     xlab="",
     ylab="",
     lwd=lwd,
     ylim=c(0,1))
lwd=2
#abline(v=0, lty=3)
segments(x0=0, x1=0, y0=-10, y1=0.7, lty=3)
wh=which(aa$deepcas13>=cutoff)
lines(density(aa$cs[wh]), col="orange", lwd=lwd)
wh=which(aa$ntx>=txcutoff)
lines(density(aa$cs[wh]), col="deepskyblue3",lwd=lwd)
wh=which(aa$deepcas13>=cutoff & aa$ntx>=txcutoff)
wh1=wh
lines(density(aa$cs[wh]), col="red", lwd=lwd)
axis(side=1, at=c(-4,-2,0,2,4))
axis(side=2, at=c(0,0.5,1))
legend("topleft", lty=1, bty="n",
       bg="white",
       cex=0.7,
       col=c("black","deepskyblue3","orange", "red"),
       c("No selection",
         "Common exon",
         "High score",
         "High score and common exon"))
dev.off()
##############################################







#### On-target and isoform specific trends
load("objects/guidesCD48.rda")
load("objects/guidesCD55.rda")


guides <- list(guidesCD46,
               guidesCD55)
guides <- lapply(guides, function(x){
    pos <- pamSites(x)
    x <- x[pos>=metadata(x)[["cds_start"]]]
    x <- x[pos<=metadata(x)[["cds_end"]]]
    x
})
lfc <- lapply(guides, function(x) x$cs)
medians <- sapply(lfc, median)
medians <- medians/median(medians)
lfc[[1]] <- lfc[[1]]/medians[[1]]
lfc[[2]] <- lfc[[2]]/medians[[2]]
lfc <- do.call(c, lfc)
score <- lapply(guides, function(x) x$score_casrxrf)
score <- do.call(c, score)
grnas <- lapply(guides, function(x){
    as.character(spacers(x))
})
grnas <- do.call(c, grnas)


pdf("figures/on_target.pdf", width=4, height=4)
library(scales)
col1 <- "firebrick2"
col2 <- "firebrick3"
x=score
y=lfc
plot(x,y,
     bty="L",
     xaxt="n",
     yaxt="n",
     xlab="",
     ylab="",
     xlim=c(0,1),
     ylim=c(-1,2.5),
     pch=20,
     cex=0.3,
     col=alpha(col1,0.3))
axis(side=1, at=c(0, 0.5, 1))
axis(side=2, at=c(-1, 0, 1,2))
lines(lowess(y~x, f=0.5),lwd=2, col=col2)
abline(h=0, lty=3)
cor(x,y)
dev.off()



prepareOnTargetData <- function(){
    out <- data.frame(grna=grnas,
                      score_casrxrf=score,
                      lfc=lfc)
    out
}


write.csv(prepareOnTargetData(),
          file="ontarget.csv",
          row.names=FALSE)


pdf("figures/isoforms_relationship.pdf", 
    height=4, width=4.5)
par(mfrow=c(1,2), bty="L", mar=c(2,2,2,2))
col <- alpha("firebrick2", 0.5)
guides <- guidesCD46
ntx <- guides$ntx
tpm <- guides$tpm
x=pamSites(guides)
y=guides$cs
tpm <- guides$tpm
cat <- rep("CDS", length(x))
cat[x<=metadata(guides)$cds_start] <- "5utr"
cat[x>=metadata(guides)$cds_end] <- "3utr"
cat[cat=="CDS" & ntx>=10] <- "high"
cat[cat=="CDS" & ntx<10]  <- "low"
dfs <- split(y, f=cat)
names(dfs) <- NULL
boxplot(dfs[c(1,2,4,3)],
        col=col,
        ylim=c(-1.5, 3),
        xaxt="n",
        las=2,
        outline=FALSE)
abline(h=0, lty=3)

guides <- guidesCD55
ntx <- guides$ntx
tpm <- guides$tpm
x=pamSites(guides)
y=guides$cs
cat <- rep("CDS", length(x))
cat[x<=metadata(guides)$cds_start] <- "5utr"
cat[x>=metadata(guides)$cds_end] <- "3utr"
cat[cat=="CDS" & ntx>=4] <- "high"
cat[cat=="CDS" & ntx<4]  <- "low"
dfs <- split(y, f=cat)
boxplot(dfs[c(1,2,4,3)],
        col=col,
        xaxt="n",
        ylim=c(-1.5, 3),
        las=2, outline=FALSE)
abline(h=0, lty=3)
dev.off()


prepareIsoformData <- function(){
    guides <- guidesCD46
    ntx <- guides$ntx
    tpm <- guides$tpm
    x=pamSites(guides)
    y=guides$cs
    tpm <- guides$tpm
    cat <- rep("CDS", length(x))
    cat[x<=metadata(guides)$cds_start] <- "5utr"
    cat[x>=metadata(guides)$cds_end] <- "3utr"
    cat[cat=="CDS" & ntx>=10] <- "high"
    cat[cat=="CDS" & ntx<10]  <- "low"
    grna <- as.character(spacers(guides))
    out1 <- data.frame(gnra=grna,
                       nIsoforms=ntx,
                       TPM=tpm,
                       region=cat,
                       lfc=y,
                       gene="CD46")
   
    guides <- guidesCD55
    ntx <- guides$ntx
    tpm <- guides$tpm
    x=pamSites(guides)
    y=guides$cs
    cat <- rep("CDS", length(x))
    cat[x<=metadata(guides)$cds_start] <- "5utr"
    cat[x>=metadata(guides)$cds_end] <- "3utr"
    cat[cat=="CDS" & ntx>=4] <- "high"
    cat[cat=="CDS" & ntx<4]  <- "low"
    grna <- as.character(spacers(guides))
    out2 <- data.frame(gnra=grna,
                       nIsoforms=ntx,
                       TPM=tpm,
                       region=cat,
                       lfc=y,
                       gene="CD55")
    out <- rbind(out1,out2)
    rownames(out) <- NULL
    out  
}



write.csv(prepareIsoformData(),
          file="isoform.csv",
          row.names=FALSE)




