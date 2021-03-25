log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(bsseq)
library(tidyverse)
library(ggrepel)

load(snakemake@input[["s"]])

meth <- getMeth(bss_smoothed,type = "smooth")

head(meth)
pca <- prcomp(t(meth))

head(pca$x)
pca_data <- as.data.frame(pca$x)
pca_data <- bind_cols(as.data.frame(pData(bss_smoothed)),pca_data)
percentVar <- round((pca$sdev^2/sum(pca$sdev^2))*100)

dims <- list(c(1,2),c(2,3),c(1,3))

lapply(dims,function(d) {
    dim1<-paste0("PC",d[1])
    dim2<-paste0("PC",d[2])

    lapply(snakemake@params[["col"]],function(v) {
        p <- ggplot(data=pca_data,aes_string(dim1,dim2,color=v)) +
                geom_point(size=3) + 
                xlab(paste0(dim1,": ",percentVar[d[1]],"% variance")) +
                ylab(paste0(dim2,": ",percentVar[d[2]], "% variance")) +
                theme_bw()
       if (snakemake@params[["labs"]]) {
          p <-  p +
           geom_text_repel(aes(label = CoreNumber))
       }
       ggsave(p,file=paste0("results/Differential_Methylation/PCA_",dim1,"vs",dim2,"_colored_by_",v,".pdf"))         
    })
})   