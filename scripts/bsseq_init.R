library("tidyverse")
library("bsseq")
library("data.table")
#library(TxDb.Mmusculus.UCSC.mm10.knownGene)

##vars
# in
meta <- snakemake@input[["meta"]]
files <- snakemake@input[["meth"]]

#out
full <- snakemake@output[["full"]]
filt <- snakemake@output[["filt"]]


###inport methylation files
meth_data <- lapply(files,fread,header=T,stringsAsFactors=F)
names(meth_data) <- gsub("_frequency.tsv","",basename(files))

meth_data <- lapply(meth_data,function(x) {
  x[,start:=start-1][order(chromosome,start)]
})


## tidy
meth_data_tidied <- Map(function(x,name) {
  dat <- x[,.(chromosome,start,end,called_sites,called_sites_methylated)]
  setnames(dat,c("called_sites","called_sites_methylated"),c(paste0("Cov.",name),paste0("M.",name)))
  return(dat)
},x=meth_data,name=names(meth_data))

## merge dataframes together
merge_dfs <- function(x,y) merge(x,y,by=c("chromosome","start","end"),all=T)

Combined <- Reduce(function(x,y) merge(x,y,by=c("chromosome","start","end"),all=T),meth_data_tidied) 

Cov_mat <- dplyr::select(Combined,starts_with("Cov")) %>%
  rename_all(function(x) gsub("Cov.","",x)) %>%
  mutate_all(function(x) ifelse(is.na(x),0,x)) %>%
  as.matrix()


Meth_mat <- dplyr::select(Combined,starts_with("M")) %>%
  rename_all(function(x) gsub("M.","",x)) %>%
  mutate_all(function(x) ifelse(is.na(x),0,x)) %>%
  as.matrix()

## make BSS object

bss <- BSseq(chr=Combined$chromosome,pos = Combined$start,Cov=Cov_mat,M=Meth_mat,sampleNames = colnames(Cov_mat))
bss
## add metadata

metadata <- read.delim(meta)

rownames(metadata) <- metadata[,1]
## make sure metadata matches methylation data
metadata <- metadata[match(colnames(getCoverage(bss,type="M")),rownames(metadata)),]
colnames(getCoverage(bss,type="M"))
pData(bss) <- DataFrame(metadata)
## filter sites with no CpGs in any samples
loci.idx <- which(rowSums(getCoverage(bss, type="Cov")>0) >=1)

bss_filtered <- bss[loci.idx]

save(bss_filtered,file = filt)


save(bss,file=full)
