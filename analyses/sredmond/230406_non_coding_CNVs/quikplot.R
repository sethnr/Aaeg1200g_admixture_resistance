library(tidyverse)
library(patchwork)

gaps <- read.table("chrom_no_repeats_no_genes_10kb.bed.gz",sep="\t",col.names = c("chrom","start","end"))
genes <- read.table("gene_regions.bed",sep="\t",col.names = c("chrom","start","end","name","score","strand"))
genes$length <- genes$end-genes$start
gaps$length <- gaps$end-gaps$start

geneplot <- ggplot(genes,aes(x=length)) + 
  geom_histogram(bins=50) + 
  xlim(0,5e4) +  ylab("genes") +
  geom_vline(xintercept = quantile(genes$length,c(0.25,0.75)))

gapplot <- ggplot(subset(gaps,length>1e3),aes(x=length)) + 
  geom_histogram(bins=50) + 
  xlim(NA,5e4) + ylab("gaps") +
  geom_vline(xintercept = quantile(genes$length,c(0.25,0.75)))

geneplot / gapplot
ggsave("gene_gap_size_distribution.png",width=200,height=160,units="mm")
