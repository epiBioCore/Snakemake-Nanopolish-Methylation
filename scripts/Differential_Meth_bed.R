 log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")
 
 sig <- read.table(snakemake@input[[1]],header = T)

 if ("pos" %in% colnames(sig)) {
     bed <- data.frame(chr=sig$chr,
                        start=sig$pos -1,
                        end=sig$pos)
 } else {
     bed <- data.frame(chr=sig$chr,
                        start=sig$start -1,
                        end=sig$end)
 }

 write.table(bed,file=snakemake@output[[1]],row.names=F,col.names=F,sep="\t",quote=F)