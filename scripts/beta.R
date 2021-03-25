log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(bsseq)
library(tidyverse)
library(dtplyr)

load(snakemake@input[["s"]])

p <- as.data.frame(pData(bss_smoothed))
meth <- getMeth(bss_smoothed,type = "smooth")

meth %>%
    t() %>%
  as.data.frame() %>%
  cbind(p,.) %>%
  rownames_to_column(var="Sample") %>%
  pivot_longer(cols = starts_with("V"),names_to = "site",values_to = "beta") %>% as.data.frame() -> plot_df




for (var in snakemake@params[["var"]]) {
 
if (! var %in% colnames(p)) {
  cat(paste(var,"is not a variable in metadata"))
  next
  } 

g <- ggplot(plot_df,aes_string(x="beta",color=var)) +
geom_density(aes(group=Sample)) +
theme_bw()

ggsave(g,file=paste0("results/Differential_Methylation/beta_distribution_colored_by_",var,".pdf"))
  
}