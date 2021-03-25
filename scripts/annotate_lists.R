log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

## packages
library(ChIPseeker)
library(annotatr)

source(file.path(snakemake@scriptdir, 'Common.R') )

txdb_pkg <- snakemake@params[["txdb"]]
species_pkg <- snakemake@params[["species"]]

cat("txdb package")
txdb_pkg
snakemake@input[["txdb"]]


cat("species package")
species_pkg
snakemake@input[["species_anno"]]

load_bioconductor_package(snakemake@input[["txdb"]], txdb_pkg)

load_bioconductor_package(snakemake@input[["species_anno"]], species_pkg)


txdb <- get(txdb_pkg)

l <- read.delim(snakemake@input[["list"]],header=T)

if ("pos" %in% colnames(l))  {
    l$end <- l$pos
    colnames(l)[2] <- "start"
}

gr <- makeGRangesFromDataFrame(l,keep.extra.columns = T)
annotated_l <- annotatePeak(gr,TxDb=txdb,annoDb = species_pkg)

pdf(file=snakemake@output[["gene_annot_pie"]])
plotAnnoPie(annotated_l)
dev.off()

## CpG island annotations
annots <- paste0(snakemake@params[["build"]],'_cpgs')
annots
annot_ranges <- build_annotations(genome = snakemake@params[["build"]],annotations = annots)
annotated_2 <- as.data.frame(annotate_regions(annotations = annot_ranges,regions = annotated_l@anno))


write.table(annotated_2,file=snakemake@output[["list"]],sep="\t",row.names=F,quote=F)

### CpG island barplot

annotated_2 <- dplyr::count(annotated_2,annot.type)
annotated_2 <- mutate(annotated_2,annot.type=gsub("inter","open_sea",annot.type),
                annot.type=gsub("mm10_cpg_","",annot.type))
 
cpg_island_counts <- annotated_2$n
names(cpg_island_counts) <- annotated_2$annot.type
cpg_island_counts
cpg_island_counts <- cpg_island_counts[order(factor(names(cpg_island_counts),levels=c("sea","island","shores","shelvs")))]
cpg_island_percent <- round(cpg_island_counts*100/sum(cpg_island_counts),2)

labels <- paste0(names(cpg_island_percent)," (",cpg_island_percent,"%)")

pdf(file=snakemake@output[["cpg_annot_pie"]])
pie(cpg_island_counts,labels=labels,clockwise = T,main="Differentially Methylated CpGs Annotation",init.angle = 0)
dev.off()



