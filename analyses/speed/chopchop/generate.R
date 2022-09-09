outdir <- tempdir()
getDelta <- function(i){
    cmd <- paste0("rm -rf ", outdir, "/*")
    system(cmd)
    time1 <- proc.time()
    input <- paste0("inputs/newsequences_",i,".fasta")
    outfile <- file.path(outdir, "results.txt")

    cmd <- paste0("./chopchop/chopchop.py -F -Target ",
              input,
              " -t WHOLE -G hg38 -v 2 --scoringMethod DOENCH_2014  -o ",
              outdir,
              " > ",
              outfile)
    cmd <- paste0("conda init bash;conda activate chopchop; ", cmd)
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
