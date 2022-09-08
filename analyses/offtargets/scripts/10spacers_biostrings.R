library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
targets <- c("ATTAAAGAAAATATCATCTT",
             "TCTGTATCTATATTCATCAT",
             "CATGGTGCATCTGACTCCTG",
             "GTAACGGCAGACTTCTCCTC",
             "TGTAGAAATCCTTCCAGTCA",
             "ATCCTTCCAGTCAGGGCCAT",
             "AGCAGCTGGGGCAGTGGTGG",
             "GCAGCTGGGGCAGTGGTGGG",
             "GGGGCAGTGGTGGGGGGCCT",
             "GCAGTGGTGGGGGGCCTTGG")

genome <- BSgenome::as.list(BSgenome.Hsapiens.UCSC.hg38)
genome <- DNAStringSet(genome)

off_target_counts <- lapply(seq_along(targets), function(i){
    target <- targets[i]
    counts <- vapply(0:3, function(ii){
        spacer <- target
        fwd <- sum(Biostrings::vcountPattern(spacer,
                                             genome,
                                             max.mismatch=ii))
        spacer_rev <- as.character(reverseComplement(DNAString(target)))
        rev <- sum(Biostrings::vcountPattern(spacer_rev,
                                             genome,
                                             max.mismatch=ii))
        fwd + rev
    }, FUN.VALUE=numeric(1))
    data.frame(n0=counts[1],
               n1=counts[2]-counts[1],
               n2=counts[3]-counts[2],
               n3=counts[4]-counts[3])
})
results <- Reduce(rbind, off_target_counts)
rownames(results) <- targets

