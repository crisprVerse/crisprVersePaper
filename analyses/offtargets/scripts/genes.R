library(crisprDesign)
library(crisprDesignData)
library(BSgenome.Hsapiens.UCSC.hg38)
bsgenome <- BSgenome.Hsapiens.UCSC.hg38
bowtie_index <- "/Users/fortinj2/crisprIndices/bowtie/hg38/hg38"
bwa_index    <- "/Users/fortinj2/crisprIndices/bwa/hg38/hg38"
txObject <- txdb_human



addBiostringAlignments <- function(gs,
                                   bsgenome,
                                   n_mismatches
){
    seqNames <- seqnames(bsgenome)
    k <- 1
    for (chr in seqNames){
        addSpacerAlignments(gs,
                            custom_seq=bsgenome[[chr]],
                            aligner="biostrings",
                            bsgenome=bsgenome,
                            n_mismatches=n_mismatches)
        k <- k+1
        print(k)
    }
}

getResults <- function(gs){
    ns_bowtie <- c(0,1,2,3)
    ns_bwa    <- c(0,1,2,3,4,5)
    ns_biostrings <- c(0,1,2)
    results <- matrix(NA, nrow=length(ns_bwa), ncol=4)
    rownames(results) <- ns_bwa
    colnames(results) <- c("bowtie", "bowtie_iter",
                           "bwa", "bwa_iter")
    
    for (mm in ns_bowtie){
        time <- system.time(addSpacerAlignments(gs,
                                                bsgenome=bsgenome,
                                                aligner_index=bowtie_index,
                                                n_mismatches=mm))[1]
        results[mm+1, 1] <- time
    }
    for (mm in ns_bowtie){
        time <- system.time(addSpacerAlignmentsIterative(gs,
                                                         bsgenome=bsgenome,
                                                         aligner_index=bowtie_index,
                                                         n_mismatches=mm))[1]
        results[mm+1, 2] <- time
    }
    for (mm in ns_bwa){
        time <- system.time(addSpacerAlignments(gs,
                                                bsgenome=bsgenome,
                                                aligner="bwa",
                                                aligner_index=bwa_index,
                                                n_mismatches=mm))[1]
        results[mm+1, 3] <- time
    }
    for (mm in ns_bwa){
        time <- system.time(addSpacerAlignmentsIterative(gs,
                                                         bsgenome=bsgenome,
                                                         aligner="bwa",
                                                         aligner_index=bwa_index,
                                                         n_mismatches=mm))[1]
        results[mm+1, 4] <- time
    }
    # for (mm in ns_biostrings){
    #     time <- system.time(addBiostringAlignments(gs=gs,
    #                                                bsgenome=bsgenome,
    #                                                n_mismatches=mm))[1]
    #     results[mm+1, 5] <- time
    # }
    return(results)
}


gene <- "KRAS"
gr <- queryTxObject(txObject=txObject,
                    queryValue=gene,
                    queryColumn="gene_symbol",
                    featureType="cds")
gs <- findSpacers(gr, bsgenome=bsgenome)
gs <- unique(gs)
nkras <- length(gs)
results_kras <- getResults(gs)
save(results_kras,
     file="../objects/results_kras.rda")


gene <- "ZNF101"
gr <- queryTxObject(txObject=txObject,
                    queryValue=gene,
                    queryColumn="gene_symbol",
                    featureType="cds")
gs <- findSpacers(gr, bsgenome=bsgenome)
gs <- unique(gs)
nznf <- length(gs)
results_znf <- getResults(gs)
save(results_znf,
     file="../objects/results_znf.rda")


gene <- "EGFR"
gr <- queryTxObject(txObject=txObject,
                    queryValue=gene,
                    queryColumn="gene_symbol",
                    featureType="cds")
gs <- findSpacers(gr, bsgenome=bsgenome)
gs <- unique(gs)
negfr <- length(gs)
results_egfr <- getResults(gs)
save(results_egfr,
     file="../objects/results_egfr.rda")






### Check identity of results between Bowtie and BWA:
checkAlignments <- function(gs){
    ns <- c(0,1,2,3)
    results_bowtie <- list()
    results_bwa    <- list()
    
    for (mm in ns){
        temp <- addSpacerAlignments(gs,
                                    bsgenome=bsgenome,
                                    aligner_index=bowtie_index,
                                    n_mismatches=mm)
        results_bowtie[[mm+1]] <- temp
    }
    for (mm in ns){
        temp<- addSpacerAlignments(gs,
                                   bsgenome=bsgenome,
                                   aligner="bwa",
                                   aligner_index=bwa_index,
                                   n_mismatches=mm)
        results_bwa[[mm+1]] <- temp
    }
    for (i in 1:4){
        print(all.equal(results_bowtie[[i]]$n0,results_bwa[[i]]$n0))
    }
    for (i in 1:4){
        print(all.equal(results_bowtie[[i]]$n1,results_bwa[[i]]$n1))
    }
    for (i in 1:4){
        print(all.equal(results_bowtie[[i]]$n2,results_bwa[[i]]$n2))
    }
    for (i in 1:4){
        print(all.equal(results_bowtie[[i]]$n3,results_bwa[[i]]$n3))
    }
    #return(list(bowtie=results_bowtie,
    #            bwa=results_bwa))
}


gene <- "KRAS"
gr <- queryTxObject(txObject=txObject,
                    queryValue=gene,
                    queryColumn="gene_symbol",
                    featureType="cds")
gs <- findSpacers(gr, bsgenome=bsgenome)
gs <- unique(gs)
checkAlignments(gs)




