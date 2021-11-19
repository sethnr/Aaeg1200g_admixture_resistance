# devtools::install_github("petrelharp/local_pca/lostruct")
library("lostruct")

args = commandArgs(trailingOnly=TRUE)

vcf = args[1]
chrom = as.numeric(args[2])
block = as.numeric(args[3])
sppfile = args[4]
spp = args[5]
outfile = args[6]

write(paste("vcf =",vcf),file=stderr())
write(paste("chrom =",chrom),file=stderr())
write(paste("block =",block),file=stderr())
write(paste("sppfile =",sppfile),file=stderr())
write(paste("spp =",spp),file=stderr())
write(paste("outfile =",outfile),file=stderr())

chromname <- c("NC_035107.1","NC_035108.1","NC_035109.1")

vcfblocksize <- 1e7
blockmb <- ((block-1)*vcfblocksize)
pcblocksize <- 1e5
regions <- data.frame("chrom"=rep(chromname[chrom],100),
                      "start"=(c(0:99)*pcblocksize)+1+blockmb,
                      "end"=c(1:100)*pcblocksize+blockmb)

#removed for custom windower (cuts of incomplete block at end)
# snps <- vcf_windower(vcf,size=pcblocksize,type='bp')
#pcs <- eigen_windows(snps,k=2)

spptab <- read.table(spptable,header=T)
samples <- spptab$sample[spptab$spp==spp]

get_set_size_blocks <- function(n) {vcf_query(vcf,regions=regions[n,],samples=samples)}
attr(get_set_size_blocks,"max.n") <- vcfblocksize/pcblocksize
#attr(get_set_size_blocks,"max.n") <- 10
attr(get_set_size_blocks,"samples") <- vcf_samples(vcf)

#get 1/2 principal components for all windows:
write(paste("getting 2PCs for ",nrow(regions),"blocks"),file=stderr())
pcs <- eigen_windows(get_set_size_blocks,k=2)

snps1 <- get_set_size_blocks(1)
write("=== SNPS ===",file=stderr())
write.table(snps1[1:10,],file=stderr(),col.names=T,quote=F,row.names=F)
write("=== REGIONS ===",file=stderr())
write.table(regions[1:10,],file=stderr(),col.names=T,quote=F,row.names=F)
write("=== PCS ===",file=stderr())
write.table(pcs[1:10,1:8],file=stderr(),col.names=T,quote=F,row.names=F)

outtable <- cbind(regions,pcs)

write.table(outtable[1:10,1:10],file=stderr(),col.names=T,quote=F,row.names=F)
write(paste("writing table",nrow(outtable),"blocks x",ncol(outtable),"inds to",outfile),file=stderr())

write.table(outtable,file=outfile,col.names=T,quote=F,row.names=F)
