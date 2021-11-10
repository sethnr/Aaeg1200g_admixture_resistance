# devtools::install_github("petrelharp/local_pca/lostruct")
library("lostruct")

args = commandArgs(trailingOnly=TRUE)

vcf = args[1]
chrom = as.numeric(args[2])
block = as.numeric(args[3])
outfile = args[4]
outposns = args[5]

chromname <- c("NC_035107.1","NC_035108.1","NC_035109.1")

vcfblocksize <- 1e7
blockmb <- ((block-1)*vcfblocksize)
pcblocksize <- 1e5
regions <- data.frame("chrom"=rep(chromname[chrom],100),
                      "start"=(c(0:99)*pcblocksize)+1,
                      "end"=c(1:100)*pcblocksize)

#vcf <- paste("data/snpable_no_outliers/aegy.wgs.aaa.aaf.mapQ20.GoodSites.chr",chrom,".10MB.blck.",block,".reg.comp.AC10.vcf.gz",sep="")

snps <- vcf_windower(vcf,size=pcblocksize,type='bp')

get_set_size_blocks <- function(n) {vcf_query(vcf,regions[n,])}
attr(get_set_size_blocks,"max.n") <- vcfblocksize/pcblocksize
attr(get_set_size_blocks,"samples") <- vcf_samples(vcf)

#get 1/2 principal components for all windows:
#pcs <- eigen_windows(snps,k=2)
pcs <- eigen_windows(get_set_size_blocks,k=2)

outtable <- cbind(regions,pcs)

write.table(outtable,file=outfile,col.names=T,quote=F,row.names=F)
