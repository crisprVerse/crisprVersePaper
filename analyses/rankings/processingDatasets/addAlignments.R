library(crisprDesign)
library(crisprDesignGne)
gs <- readRDS("../processingRankings/crisprverse/crisprverse.results.rds")
cols <- c("ensembl_id", "spacer", "percentGC", "polyT", 
          "n0", "n0_c", "n1", "n1_c", "n2", "n2_c", "n3", "n3_c", 
          "hasSNP")
newCols <- c("percentGC", "polyT",  "hasSNP",
             "n0", "n0_c", "n1", "n1_c", 
             "n2", "n2_c", "n3", "n3_c")
gs <- gs[,cols]
gs$name20 <- paste0(gs$ensembl_id, "_", gs$spacer)
gs$name19 <- paste0(gs$ensembl_id, "_", substr(gs$spacer,2,20))


load("objects/rankings_achilles.rda")
wh <- match(rankings_achilles$name, gs$name20)
rankings_achilles <- cbind(rankings_achilles, gs[wh, newCols])
save(rankings_achilles, file="objectsFinal/rankings_achilles.rda")

load("objects/rankings_achilles_neg.rda")
wh <- match(rankings_achilles_neg$name, gs$name20)
rankings_achilles_neg <- cbind(rankings_achilles_neg, gs[wh, newCols])
save(rankings_achilles_neg, file="objectsFinal/rankings_achilles_neg.rda")


load("objects/rankings_sabatini.rda")
wh <- match(rankings_sabatini$name, gs$name20)
rankings_sabatini <- cbind(rankings_sabatini, gs[wh, newCols])
save(rankings_sabatini, file="objectsFinal/rankings_sabatini.rda")

load("objects/rankings_sabatini_neg.rda")
wh <- match(rankings_sabatini_neg$name, gs$name20)
rankings_sabatini_neg <- cbind(rankings_sabatini_neg, gs[wh, newCols])
save(rankings_sabatini_neg, file="objectsFinal/rankings_sabatini_neg.rda")

load("objects/rankings_toronto.rda")
wh <- match(rankings_toronto$name, gs$name20)
rankings_toronto <- cbind(rankings_toronto, gs[wh, newCols])
save(rankings_toronto, file="objectsFinal/rankings_toronto.rda")

load("objects/rankings_toronto_neg.rda")
wh <- match(rankings_toronto_neg$name, gs$name20)
rankings_toronto_neg <- cbind(rankings_toronto_neg, gs[wh, newCols])
save(rankings_toronto_neg, file="objectsFinal/rankings_toronto_neg.rda")

load("objects/rankings_toronto3.rda")
wh <- match(rankings_toronto3$name, gs$name20)
rankings_toronto3 <- cbind(rankings_toronto3, gs[wh, newCols])
save(rankings_toronto3, file="objectsFinal/rankings_toronto3.rda")


load("objects/rankings_toronto3_neg.rda")
wh <- match(rankings_toronto3_neg$name, gs$name20)
rankings_toronto3_neg <- cbind(rankings_toronto3_neg, gs[wh, newCols])
save(rankings_toronto3_neg, file="objectsFinal/rankings_toronto3_neg.rda")

load("objects/rankings_yusa.rda")
wh <- match(rankings_yusa$name, gs$name19)
rankings_yusa <- cbind(rankings_yusa, gs[wh, newCols])
save(rankings_yusa, file="objectsFinal/rankings_yusa.rda")

load("objects/rankings_yusa_neg.rda")
wh <- match(rankings_yusa_neg$name, gs$name19)
rankings_yusa_neg <- cbind(rankings_yusa_neg, gs[wh, newCols])
save(rankings_yusa_neg, file="objectsFinal/rankings_yusa_neg.rda")




