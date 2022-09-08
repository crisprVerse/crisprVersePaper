egs <- read.table("essential_genes_hart2014.txt")[,1]
negs <- read.table("nonessential_genes_hart2014.txt")[,1]
save(egs, file="egs.rda")
save(negs, file="negs.rda")
