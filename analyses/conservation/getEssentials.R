egs <- read.csv("../rankings/achilles/raw/Achilles_common_essentials.csv")[,1]
egs <- gsub(" \\([0-9]+\\)","", egs)
save(egs, file="objects/egs.rda")
negs <- read.csv("../rankings/achilles/raw/nonessentials.csv")[,1]
negs <- gsub(" \\([0-9]+\\)","", negs)
save(negs, file="objects/negs.rda")
