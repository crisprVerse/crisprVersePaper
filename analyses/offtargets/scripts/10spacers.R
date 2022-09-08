library(crisprDesign)
library(crisprBowtie)
library(crisprBwa)
library(BSgenome.Hsapiens.UCSC.hg38)
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
#bowtie_index <- "/gstore/data/omni/crispr/crisprIndices/bowtie/hg38/hg38"
#bwa_index <- "/gstore/data/omni/crispr/crisprIndices/bwa/hg38/hg38"
bowtie_index <- "/Users/fortinj2/crisprIndices/bowtie/hg38/hg38"
bwa_index <- "/Users/fortinj2/crisprIndices/bwa/hg38/hg38"
bsgenome <- BSgenome.Hsapiens.UCSC.hg38

aln <- getSpacerAlignments(targets,
                           bsgenome=bsgenome,
                           all_alignments=TRUE,
                           ignore_pam=TRUE,
                           standard_chr_only=FALSE,
                           aligner_index=bowtie_index,
                           n_mismatches=3)
aln <- as.data.frame(aln)

aln_bowtie <- aln
dfs <- split(aln, f=aln$spacer)[targets]
ns <- lapply(dfs, function(df){
    n0=sum(df$n_mismatches==0)
    n1=sum(df$n_mismatches==1)
    n2=sum(df$n_mismatches==2)
    n3=sum(df$n_mismatches==3)
    c(n0, n1, n2, n3)
})
ns_bowtie <- do.call(rbind, ns)




# BWA
aln <- runCrisprBwa(targets,
                    bsgenome=bsgenome,
                    ignore_pam=TRUE,
                    aligner_index=bwa_index,
                    n_mismatches=3)
aln_bwa <- aln
dfs <- split(aln, f=aln$spacer)[targets]
ns <- lapply(dfs, function(df){
    n0=sum(df$n_mismatches==0)
    n1=sum(df$n_mismatches==1)
    n2=sum(df$n_mismatches==2)
    n3=sum(df$n_mismatches==3)
    c(n0, n1, n2, n3)
})
ns_bwa <- do.call(rbind, ns)



aln1 <- aln_bowtie
aln1 <- aln1[aln1$n_mismatches==3,]
aln1$chr <- aln1$seqnames
aln1$pam_site <- aln1$start
aln1$id <- paste0(aln1$chr, ":",
                  aln1$pam_site,":",
                  aln1$spacer, ":",
                  aln1$protospacer)
#aln1$id <- paste0(aln1$seqnames, ":",
#                  aln1$start,":")

aln2 <- aln_bwa
aln2 <- aln2[aln2$n_mismatches==3,]
aln2$id <- paste0(aln2$chr, ":",
                  aln2$pam_site,":",
                  aln2$spacer, ":",
                  aln2$protospacer)
#aln2$id <- paste0(aln2$seqnames, ":",
#                  aln2$start,":")

setdiff(aln1$id, aln2$id)
setdiff(aln2$id, aln1$id)




# Without alternative loci:
aln <- getSpacerAlignments(targets,
                           bsgenome=bsgenome,
                           all_alignments=TRUE,
                           ignore_pam=TRUE,
                           standard_chr_only=TRUE,
                           aligner_index=bowtie_index,
                           n_mismatches=3)
aln <- as.data.frame(aln)
dfs <- split(aln, f=aln$query)[targets]
ns <- lapply(dfs, function(df){
    n0=sum(df$n_mismatches==0)
    n1=sum(df$n_mismatches==1)
    n2=sum(df$n_mismatches==2)
    n3=sum(df$n_mismatches==3)
    c(n0, n1, n2, n3)
})
ns <- do.call(rbind, ns)
print(ns)


