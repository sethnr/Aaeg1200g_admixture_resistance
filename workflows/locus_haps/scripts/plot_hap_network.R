
if(!require(ape)){install.packages("ape")
                          library(ape)}
if(!require(pegas)){install.packages("pegas")
                          library(pegas)}
if(!require(vcfR)){install.packages("vcfR")
                          library(vcfR)}

f <- "snpeff.kdr.nonsyn.vcf.gz"

library("getopt")

opttab <- matrix(c("vcf","i","1","character",
                   "out","o","1","character",
                   "meta","m","1","character",
),byrow=T,ncol=4)
opt <- getopt(opttab)

f <- opt$vcf
outprefix <- opt$out
metafile <- opt$meta

snps <- VCFloci(f)
samples <- VCFlabels(f)

v <- pegas::read.vcf(f)
h <- pegas::haplotype(v,is.snp(snps))
dh <- dist.haplotype.loci(h)

x <- vcfR2DNAbin(read.vcfR(f, verbose = FALSE))
h <- haplotype(x)
d <- dist.dna(h,"n")

# #min spanning tree
mstree <- rmst(d)
#minimum spanning network
msnet <- mst(d)
nt <- msnet

meta <- read.table(metafile,header=T,sep="\t")
rownames(meta) <- meta$sample

R <- haploFreq(x,fac = rep(meta[samples,"region"],each=2),h)

ntlabs <- attr(nt,"labels")
sz <- summary(h)

png(paste(prefix,".png",),width=800,height=800)
plot(nt,fast=T,
     pie=R[ntlabs,],size=sz[ntlabs],
     legend = T)
dev.off()
