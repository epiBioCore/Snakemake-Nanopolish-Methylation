log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(DSS)
library(bsseq)
library(BiocParallel)






load(file = snakemake@input[["filt"]])

contrast <- snakemake@params[[1]]
col <- contrast[3]
treat <- contrast[1]
control <- contrast[2]

col
treat
control
p <- as.data.frame(pData(bss_filtered))
p
if (!col %in% colnames(p)) { 
    stop(paste("Variable",col,"is not in bsseq metadata."))
    }

if (!treat %in% p[,col]) { 
    stop(paste(treat,"is not a level of", col)) 
    }

if (!control %in% p[,col]) { 
    stop(paste(control,"is not a level of", col)) 
    }

g1 <- rownames(p[p[,col] == treat,])

g2 <- rownames(p[p[,col] == control,])

g1

g2
mParam = MulticoreParam(workers=snakemake@resources[["cpus"]], progressbar=TRUE)
mParam
dml <- DMLtest(BSobj = bss_filtered,
               group1 = g1,
               group2 = g2,
               smoothing = T,
               BPPARAM=mParam)

name <- contrast

save(dml,name,file=snakemake@output[[1]])


