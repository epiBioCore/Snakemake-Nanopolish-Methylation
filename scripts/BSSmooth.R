log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(bsseq)
library(BiocParallel)

load(snakemake@input[["filt"]])



bss_smoothed <- BSmooth(
    BSseq = bss_filtered, 
    BPPARAM = MulticoreParam(workers = snakemake@resources["cpus"]), 
    verbose = TRUE)


save(bss_smoothed,file=snakemake@output[["s"]])