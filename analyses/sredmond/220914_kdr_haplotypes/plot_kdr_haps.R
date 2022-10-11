library(ape)
library(pegas)
library(vcfR)

f <- "randthin_kdr_SNPs.coding.vcf"

snps <- VCFloci(f)
samples <- VCFlabels(f)

v <- pegas::read.vcf(f)
h <- pegas::haplotype(v,is.snp(snps))
dh <- dist.haplotype.loci(h)

x <- vcfR2DNAbin(read.vcfR(f, verbose = FALSE))
h <- haplotype(x)
d <- dist.dna(h,"n")

# #min spanning tree
# nt <- rmst(d)
#minimum spanning network
nt <- mst(d)


meta <- read.table("meta_Aaeg1kg_spp.txt",header=T,sep="\t")
rownames(meta) <- meta$sample

R <- haploFreq(x,fac = rep(meta[samples,"region"],each=2),h)

ntlabs <- attr(nt,"labels")

sz <- summary(h)


png(gsub(".vcf",".png",f),width=800,height=800)
plot(nt,fast=T,threshold=c(5,10),
     pie=R[ntlabs,],size=sz[ntlabs]/100,
     legend = c(0.15, 9.5))
dev.off()

