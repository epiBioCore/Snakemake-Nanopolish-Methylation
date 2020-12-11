log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(DSS)
library(bsseq)
library(BiocParallel)
library(tidyverse)

load(snakemake@input[[1]])

treat <- name[1]
control <- name[2]

## create Directories for output
dir.create(dirname(snakemake@output[["dml"]]),recursive=T,showWarnings=F)

## sig CpGs
Sig_DML <- callDML(dml,delta=snakemake@params[["thres"]])


Sig_DML <- rename(Sig_DML,!!paste(treat,"mean",sep="."):=mu1,
                        !!paste(control,"mean",sep="."):=mu2,
                        !!paste(treat,"dispersion",sep="."):=phi1,
                        !!paste(control,"dispersion",sep="."):=phi2)

write.table(Sig_DML,file=snakemake@output[["dml"]],sep="\t",row.names=F,quote=F)

## sig regions
Sig_DMR <- callDMR(dml,delta=snakemake@params[["thres"]])

Sig_DMR <- rename(Sig_DMR,!!paste(treat,"meanMethyl",sep="."):=meanMethy1,
                        !!paste(control,"meanMethyl",sep="."):=meanMethy2)

write.table(Sig_DMR,file=snakemake@output[["dmr"]],sep="\t",row.names=F,quote=F)

## all
all <- as.data.frame(dml)
class(all) <- "data.frame"
all <- rename(all,!!paste(treat,"mean",sep="."):=mu1,
                        !!paste(control,"mean",sep="."):=mu2,
                        !!paste(treat,"dispersion",sep="."):=phi1,
                        !!paste(control,"dispersion",sep="."):=phi2)

write.table(all,file=snakemake@output[["all"]],sep="\t",row.names=F,quote=F)