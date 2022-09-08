library(SummarizedExperiment)
library(crisprDesign)
load("objects/se.rda")
se <- se[!grepl("rc_", rownames(se)),]
ann <- rowData(se)
cs <- read.csv("data/GFP_screen_crRNA_enrichments.csv")
cs <- as.data.frame(cs)


# Annotating mismatches:
mm <- annotateMismatches(spacers=ann$spacer,
                         protospacers=ann$protospacer)
mm <- mm[,c("mm1", "mm2", "mm3")]
ann <- cbind(ann, mm)

# Adding logFC:
ann$lfc <- cs[,"meanCS.BIN1"][match(ann$ID, rownames(cs))]
ann$lfc_ontarget <- ann$lfc[match(ann$ID_ontarget, ann$ID)]
ann <- ann[which(ann$lfc_ontarget>=0),]
ann <- ann[!is.na(ann$lfc) & !is.na(ann$lfc_ontarget),]

#Let's calculate a delta log-fold changes:
ann$delta <-  ann$lfc-ann$lfc[match(ann$ID_ontarget, ann$ID)]
ann$dist <- ann$mm2-ann$mm1


pdf("figures/mm_casrx_weights.pdf", width=10, height=4)
####### Getting SM weights #######
library(scales)
par(mfrow=c(1,2))
df <- ann[ann$n_mismatches==1,]
dfs <- split(df$delta, f=df$mm1)
boxplot(rev(dfs),
        main="A",
        col=alpha("firebrick2",0.3),
        border="firebrick3",
        xlab="Mismatch position",
        ylab="Log2 fold change",
        outline=FALSE,
        cex.axis=0.7,
        las=2)
abline(h=0, lty=3)
max_delta <- median(ann$lfc[ann$n_mismatches==0])

xx <- 28-df$mm1
fit <- loess(df$delta~xx, span=0.3)
fitted.values <- predict(fit)
fitted.values <- fitted.values[match(27:1, xx)]
lines(27:1,fitted.values)
abline(h=-max_delta, lty=3)

# Obtaining weights:
weights <- 1-fitted.values/(-max_delta)
names(weights) <- 1:27
weights[weights<0] <- 0
weights[weights>1] <- 1
plot(1:27,
     main="B",
     weights[27:1],
     col="firebrick3",
     ylab="Mismatch tolerance score",
     xlab="Mismatch position",
     pch=20,
     xaxt="n",
     ylim=c(0,1))
axis(side=1, at=1:27, labels=27:1, las=2,cex.axis=0.7,)
dev.off()
######################################################



#To save
CasRxWeights <- data.frame(w=weights)
CasRxWeights$position <- 1:27
save(CasRxWeights, file="cfd_weights/CasRxWeights.rda")
write.csv(CasRxWeights,
          quote=FALSE,
          row.names=FALSE,
          file="cfd_weights/CasRxWeights.csv")








library(crisprDesign)
gene <- "CD55"
load("objects/annsTiling.rda")
ann <- annsTiling[[gene]]
ann <- ann[!is.na(ann$protospacer),]

mm <- annotateMismatches(spacers=ann$spacer,
                         protospacers=ann$protospacer)
mm <- mm[,c("mm1", "mm2", "mm3")]
ann <- cbind(ann, mm)


# Adding logFC:
file <- paste0("../data/",
               gene,
               "_screen_crRNA_enrichments.csv")
cs <- read.csv(file)
cs <- as.data.frame(cs)
ann$lfc <- cs[,"medianCS.Top"][match(ann$ID, rownames(cs))]
ann$lfc_ontarget <- ann$lfc[match(ann$ID_ontarget, ann$ID)]
ann <- ann[which(ann$lfc_ontarget>=0),]
ann <- ann[!is.na(ann$lfc) & !is.na(ann$lfc_ontarget),]
#Let's calculate a delta log-fold changes:
ann$delta <-  ann$lfc-ann$lfc[match(ann$ID_ontarget, ann$ID)]
ann$dist <- ann$mm2-ann$mm1
#df <- ann[ann$n_mismatches==2,]
df <- ann



load("cfd_weights/CasRxWeights.rda")
weights <- CasRxWeights
weights <- weights[,1]
names(weights) <-CasRxWeights[,2]
############ Adding scores ############
mm <- as.data.frame(df[, c("mm1", "mm2","mm3")])
scores <- c()
for (i in 1:nrow(mm)){
    pos <- mm[i,]
    pos <- pos[!is.na(pos)]
    if (length(pos)==0){
        scores[i] <- NA
    } else {
        scores[i] <- prod(weights[as.character(pos)])
    }
}
df$score <- scores
########################################


############ Testing ############
test <- df
test <- test[test$n_mismatches!=0,]
test$predicted_lfc <- test$lfc_ontarget*test$score

col <- test$n_mismatches
x <- test$predicted_lfc[col==1]
y <- test$lfc[col==1]
plot(x,y,
     xlim=c(0,1.2),
     ylim=c(-1,2),
     pch=20,cex=0.4,
     col=1)
lines(lowess(y~x, f=0.5),lwd=3)
cor(x,y)
abline(h=0, lty=3)




###### Off-target figures ######
pdf("figures/off_target_2mm.pdf", width=2.5, height=4)
library(scales)
x <- test$predicted_lfc[col==2]
y <- test$lfc[col==2]
cex <- test$score[col==2]
cor(x,y)
col1 <- "firebrick2"
col2 <- "firebrick3"
plot(x,y,
     bty="L",
     xaxt="n",
     yaxt="n",
     xlab="",
     ylab="",
     xlim=c(0,1.2),
     ylim=c(-1,2),
     pch=20,
     cex=0.3,
     col=alpha(col1,0.3))
axis(side=1, at=c(0, 0.5, 1))
axis(side=2, at=c(-1, 0, 1,2))
lines(lowess(y~x, f=0.9),lwd=2, col=col2)
abline(h=0, lty=3)
cor(x,y)
dev.off()


pdf("figures/off_target_1mm.pdf", width=2.5, height=4)
library(scales)
x <- test$predicted_lfc[col==1]
y <- test$lfc[col==1]
cex <- test$score[col==2]
cor(x,y)
col1 <- "firebrick2"
col2 <- "firebrick3"
plot(x,y,
     bty="L",
     xaxt="n",
     yaxt="n",
     xlab="",
     ylab="",
     xlim=c(0,1.2),
     ylim=c(-1,2),
     pch=20,
     cex=0.3,
     col=alpha(col1,0.3))
axis(side=1, at=c(0, 0.5, 1))
axis(side=2, at=c(-1, 0, 1,2))
lines(lowess(y~x, f=0.7),lwd=2, col=col2)
abline(h=0, lty=3)
cor(x,y)
dev.off()


###### Off-target figures ######









