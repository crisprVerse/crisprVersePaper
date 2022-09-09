#module load Java/1.8.0_144
#database=/gstore/data/omni/crispr/flashfry/cas9_human
getDelta <- function(i){
    time1 <- proc.time()
    database <- "cas9_human"
    input <- paste0("inputs/sequences_",i,".fasta")
    output <- paste0("outputs/output.txt")
    cmd <- paste0("java -Xmx4g -jar FlashFry-assembly-1.15.jar discover --database ",
     database,
     " --fasta ",input,
     " --maximumOffTargets ", 100000L,
     " --forceLinear ",
     " --maxMismatch ", 2,
     " --output ",output)
    system(cmd)
    time2 <- proc.time()
    delta <- as.numeric(time2-time1)
    return(delta)
}
results <- lapply(1:6, function(i){
    print(i)
    out <- getDelta(i)
    out
})
save(results, file="resultsMm2.rda")
