library(ape)
library(pegas)
#library(vcfR)

library("getopt")

opttab <- matrix(c("vcf","i","1","character",
                   "out","o","1","character",
                   "meta","m","1","character"
                   ),byrow=T,ncol=4)
opt <- getopt(opttab)

f <- opt$vcf
outprefix <- opt$out
metafile <- opt$meta

write(f,file=stderr())

snps <- VCFloci(f)
samples <- VCFlabels(f)

v <- pegas::read.vcf(f)
write("is snp",file=stderr())
h <- pegas::haplotype(v,is.snp(snps))
write("makehaps",file=stderr())
dh <- dist.haplotype.loci(h)

write("make DNA bin",file=stderr())

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
     pie=R[ntlabs,],size=sz[ntlabs]/10,
     legend = T)
dev.off()
