library(pbapply)
lib <- readRDS("crisprverse.results.rds")


# Removing useless columns:
cols <- c("spacer", "polyT", "ensembl_id",
          "n0", "n0_c", "n1", "n1_c",  "n2", "n2_c",  "n3", "n3_c", 
          "score_azimuth", "score_deephf", "score_deepspcas9", 
          "score_ruleset1", "score_crisprater",
          "score_crisprscan", "score_cfd", "score_mit", 
          "hasSNP", "percentGC", 
          "score_conservation", "score_composite",
          "percentCDS")
lib <- lib[,cols]
lib$name <- paste0(lib$ensembl_id, "_", lib$spacer)
lib$percentCDS[is.na(lib$percentCDS)] <- 100
lib$score_conservation_binary <- ifelse(lib$score_conservation>=0,1,0)
lib$score_cds <- ifelse(lib$percentCDS<=85,1,0)

# Creating bins:
lib$round <- NA
lib$round[lib$n0!=1 | (lib$n1_c + lib$n2_c)>5]   <- 3
lib$round[lib$n0==1 & (lib$n1_c + lib$n2_c)<=5]  <- 2
lib$round[lib$n0==1 & lib$n1_c==0 & lib$n2_c==0] <- 1
lib$round[lib$polyT | lib$hasSNP | lib$percentGC<20 | lib$percentGC>80 | is.na(lib$n1) | is.na(lib$n2)] <- 4


# Creating alternative rankings:
dfs <- split(lib, f=lib$ensembl_id)
dfs <- pblapply(dfs, function(df){
    df$rank_deephf  <- rank(-df$score_deephf, ties.method="random")
    df$rank_deepspcas9 <- rank(-df$score_deepspcas9, ties.method="random")
    df
})


# Creating final score:
dfs <- pblapply(dfs, function(df){
    df <- df[order(df$round,
                   -df$score_composite),,drop=FALSE]
    df$rank_crisprdesign <- seq_len(nrow(df))
    df <- df[order(df$round,
                   -df$score_cds,
                   -df$score_conservation_binary,
                   -df$score_composite),,drop=FALSE]
    df$rank_crisprdesign_genic <- seq_len(nrow(df))
    df
})


lib <- dplyr::bind_rows(dfs)
saveRDS(lib,
        file="../objects/crisprverse.results.rds")






