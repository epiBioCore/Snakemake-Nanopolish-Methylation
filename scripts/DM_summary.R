log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


dml_files <- snakemake@input[["dml"]]

dmls <- lapply(dml_files,read.delim,header = T)
names(dmls) <- gsub("_SIG_Meth_CpGs.txt","",basename(dml_files))

dml_summary <- data.frame(
    Hypomethylated = sapply(dmls,function(x) nrow(x[x$diff<0,])),
    Hypermethylated = sapply(dmls,function(x) nrow(x[x$diff >0,])),
    Total = sapply(dmls,nrow)
)


## DMRs

dmr_files <- snakemake@input[["dmr"]]

dmrs <- lapply(dmr_files,read.delim,header = T)
names(dmrs) <- gsub("_SIG_Meth_regions.txt","",basename(dmr_files))

dmr_summary <- data.frame(
    Hypomethylated = sapply(dmrs,function(x) nrow(x[x$diff<0,])),
    Hypermethylated = sapply(dmrs,function(x) nrow(x[x$diff >0,])),
    Total = sapply(dmrs,nrow)
)


DM_summary <- merge(dml_summary,dmr_summary,by.x="row.names",by.y="row.names",suffixes = c(".dml",".dmr"))
colnames(DM_summary)[1] <- "Comparison"

write.table(DM_summary,file=snakemake@output[[1]],row.names=F,sep="\t",quote=F)