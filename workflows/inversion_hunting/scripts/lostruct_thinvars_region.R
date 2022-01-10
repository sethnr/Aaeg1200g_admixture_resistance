library("lostruct")
library("tidyverse")

args = commandArgs(trailingOnly=TRUE)

vcf <- args[1]
sppfile <- args[2]
chrom <- as.numeric(args[3])
cgroup <- args[4]
outtxt <- args[5]
outpng <- args[6]

write(paste("vcf: ",vcf),file=stderr())
write(paste("sppfile: ",sppfile),file=stderr())
write(paste("chrom: ",chrom),file=stderr())
write(paste("region: ",cgroup),file=stderr())
write(paste("pcs: ",outtxt),file=stderr())
write(paste("png: ",outpng),file=stderr())



chromname <- c("NC_035107.1","NC_035108.1","NC_035109.1")
chromlen <- c(310827022,474425716,409777670)
names(chromlen) <- chromname

pcblocksize <- 5e5

#get number of blocks in chrom
regionct <- ceiling(chromlen[chrom] / pcblocksize)
regions <- data.frame("chrom"=rep(chromname[chrom],regionct),
                      "start"=(0:(regionct-1)*pcblocksize) + 1,
                      "end"=1:regionct*pcblocksize)


spptab <- read.table(sppfile,header=T, sep="\t")
samples <- spptab$sample[tolower(spptab$contgroup)==cgroup]

get_set_size_blocks <- function(n) {vcf_query(vcf,regions=regions[n,],samples=samples)}
attr(get_set_size_blocks,"max.n") <- regionct
attr(get_set_size_blocks,"samples") <- samples #vcf_samples(vcf)


#get 1/2 principal components for all windows:
pcs <- eigen_windows(get_set_size_blocks,k=2,mc.cores=3)

outtable <- cbind(regions,pcs)
write.table(outtable,file=outtxt,col.names=T,quote=F,row.names=F)


pcdist <- pc_dist(pcs,npc=2,w=1)
pcdistdf <- as.data.frame(pcdist)



blocks <- paste(regions[,1],regions[,2],sep=":")
regions$block <- blocks
regions$pos <- regions$end / pcblocksize

colnames(pcdistdf) <- blocks
pcdistdf$x <- blocks
pcdistflat <- pivot_longer(pcdistdf,cols=all_of(blocks),names_to = "y")
pcdistflat <- merge(merge(pcdistflat,regions,by.x="x",by.y="block"),regions,by.x="y",by.y="block",suffixes = c(".x",".y"))


ggplot(pcdistflat,aes(x=pos.x,y=pos.y,fill=value)) + geom_raster() + coord_fixed() +
  ggtitle(paste(vcf,"\nchrom",chrom," ",cgroup," (",pcblocksize,"bp blocks)"))
ggsave(outpng)



