library(crisprDesignData)
library(S4Vectors)
library(dplyr)
library(pbapply)
key <- mcols(txdb_human$cds)[, c("gene_symbol", "gene_id")]
key <- as.data.frame(key)
key <- key[!duplicated(key),]
cols <- c("name",  "ensembl_id", "spacer", "rank")

getCCTOP <- function(){
    cctop <- readRDS("objects/cctop.results.rds")
    cctop$spacer <- substr(cctop$protospacer,1,20)
    cctop$name <- paste0(cctop$ensembl_id, "_", cctop$spacer)
    cctop <- cctop[,c("spacer", "name", "rank", "ensembl_id")]
    cctop$rank <- 1000-cctop$rank+1 
    cctop[,cols]
}

getFlash <- function(){
    flash <- readRDS("objects/flashfry.results.rds")
    flash$spacer <- substr(flash$protospacer,1,20)
    flash$name <- paste0(flash$ensembl_id, "_", flash$spacer)
    flash[,cols]
}

getCpick <- function(){
    cpick <- readRDS("objects/crispick.results.rds")
    cpick <- dplyr::rename(cpick, spacer=spacer_20mer)
    cpick <- cpick[cpick$gene_symbol %in% key$gene_symbol,] 
    cpick$ensembl_id <- key$gene_id[match(cpick$gene_symbol, key$gene_symbol)]
    cpick$name <- paste0(cpick$ensembl_id, "_", cpick$spacer)
    cpick[,cols]
}

getChopChop <- function(){
    chopchop <- readRDS("objects/chopchop.results.rds")
    chopchop$spacer <- substr(chopchop$protospacer,1,20)
    chopchop$name <- paste0(chopchop$ensembl_id, "_", chopchop$spacer)
    chopchop
}

# Getting the data:
cctop <- getCCTOP()
flash <- getFlash()
cpick <- getCpick()
chopchop <- getChopChop()


# Renaming:
cctop <- dplyr::rename(cctop, rank_cctop=rank)
flash <- dplyr::rename(flash, rank_flash=rank)
cpick <- dplyr::rename(cpick, rank_cpick=rank)
chopchop <- dplyr::rename(chopchop, rank_chopchop=rank)


# Merging everything together:
final <- cctop
final <- full_join(final, flash[,c("name", "rank_flash")], by="name")
final <- full_join(final, cpick[,c("name", "rank_cpick")], by="name")
final <- full_join(final, chopchop[,c("name", "rank_chopchop")], by="name")




# Adding crisprverse:
new   <- readRDS("objects/crisprverse.results.rds")
new <- new[, c("name", 
               "rank_crisprdesign",
               "rank_crisprdesign_genic")]
final <- full_join(final, new, by="name")
rankings <- final



# Cleaning up rankings
dfs <- split(rankings,
             f=rankings$ensembl_id)
dfs <- pblapply(dfs, function(df){
    df <- df[!is.na(df$rank_crisprdesign),,drop=FALSE]
    df$rank_flash <- rank(df$rank_flash, ties.method="first")
    df$rank_cctop <- rank(df$rank_cctop, ties.method="first")
    df$rank_chopchop <- rank(df$rank_chopchop, ties.method="first")
    df$rank_cpick <- rank(df$rank_cpick, ties.method="first")
    return(df)
})
rankings <- dplyr::bind_rows(dfs)
saveRDS(rankings,
        file="objects/rankings.rds")



