#module load Java/1.8.0_144
#database=/gstore/data/omni/crispr/flashfry/cas9_human
getDelta <- function(i){
    time1 <- proc.time()
    index <- "/Users/fortinj2/crisprIndices/bowtie/hg38/hg38"
    input <- paste0("inputs/sequences_",i,".fasta")
    output <- tempdir()
    cmd <- paste0("cctop --input ",
              input,
              " --index ",
              index,
              " --output ",
              output,
              " --totalMM 2 --coreMM 2 --maxOT 100000")
    cmd <- paste0("conda init bash;conda activate cctop; ", cmd)
    cmd <- paste0("source /Users/fortinj2/miniconda3/etc/profile.d/conda.sh; ", cmd)
    system(cmd)
    time2 <- proc.time()
    delta <- as.numeric(time2-time1)
    return(delta)
}
results <- lapply(1:6, function(i){
    print(i)
    out <- getDelta(i)
    print(out)
    out
})
save(results, file="results.rda")




1. 2468.863
2. 3193.461
3. 5452.507
4. 8394.946
5. 18551.761

library(readxl)
results <- read_excel("results.xlsx", col_names=FALSE)
results <- as.data.frame(results)
results <- split(results,f=1:6)
save(results, file="resultsUpTo5.rda")