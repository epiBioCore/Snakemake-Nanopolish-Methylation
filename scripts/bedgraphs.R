log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


load(file=snakemake@input[["s"]])
sample <- gsub("\\.bdg","",basename(snakemake@output[[1]]))

bss_sub <- bss_smoothed[,sample]

meth <- getMeth(bss_sub,type = "smooth")

Ranges <- rowRanges(bss_sub)
Ranges$Meth <- meth

bdg <-  data.frame(seqnames=seqnames(Ranges),
starts=start(Ranges)-1,
ends=end(Ranges),
scores=elementMetadata(Ranges)[,1]
)

write.table(bdg,file=snakemake@output[[1]],row.names=F,col.names=F,sep="\t",quote=F)