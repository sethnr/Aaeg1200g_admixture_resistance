
library("tidyverse")
library("patchwork")
library("lostruct")


args = commandArgs(trailingOnly=TRUE)

# dists <- "redmond-lab-aaeg1000g/results_lostruct/lostruct_all_chr/lostruct_chr1_Senegal.txt"
# snps <- "/Volumes/Mosquito_raw_data/Aedes/Aaeg1000g/thinrand/lostruct_chr1.vcf.gz"
# country <- "Senegal"
# chrom <- 1
# outtxt <- "./inv_candidates_chr1_Senegal.txt"
# outpng <- "./inv_candidates_chr1_Senegal.png"
# sppfile <- "resources/meta_Aaeg1kg_spp.txt"
# invfile <- "resources/redmond_2020_inversion_calls.txt"

dists <- args[1]
snps <- args[2]
country = args[3]
chrom <- as.numeric(args[4])
sppfile <- args[5]
#invfile <- args[6]
outtxt <- args[6]
outpng <- args[7]



chromname <- c("NC_035107.1","NC_035108.1","NC_035109.1")
chromlen <- c(310827022,474425716,409777670)
names(chromlen) <- chromname

pcblocksize <- 5e5

spptab <- read.table(sppfile,header=T, sep="\t")

samples <- spptab$sample[spptab$country==country]

#read table, parse out PC dists
pcdistdf <- read.table(dists,header=T)
pcdistdf <- pcdistdf[,c(1:dim(pcdistdf)[1])]
pos <- pcblocksize * 1:dim(pcdistdf)[1]
blocks <- paste(chrom,pos,sep=":")

regions <- data.frame(
                    "chrom"=rep(chrom,length(pos)),
                    "pos"=pos,
                    "block"=blocks)
colnames(pcdistdf) <- blocks

pcdistdf$x <- blocks
pcdistflat <- pivot_longer(pcdistdf,cols=all_of(blocks),names_to = "y")
pcdistflat <- merge(merge(pcdistflat,regions,by.x="x",by.y="block"),regions,by.x="y",by.y="block",suffixes = c(".x",".y"))



invcands <- data.frame("inversion"=character(),
                       "chrom"=numeric(),
                       "start"=numeric(),
                       "end"=numeric())


#remake indices from posn x/y
pcdistflat$y <- pcdistflat$pos.y/pcblocksize
pcdistflat$x <- pcdistflat$pos.x/pcblocksize

pcdistmat <- matrix(nrow=max(pcdistflat$x),ncol=max(pcdistflat$y))
for(i in c(1:nrow(pcdistflat))) {
  pcdistmat[pcdistflat$y[i],pcdistflat$x[i]] <- pcdistflat$value[i] 
}


mdsk <- 5
#set NAs to 1
pcdistmat[is.na(pcdistmat)] <- 1
mds <- cmdscale(pcdistmat,k=mdsk)

outliersT <- mds==NA

maxgap=10
minsize=4
#for(k in c(1:mdsk)) {

for(k in c(1:mdsk)) {
    koutliers <- which(mds[,k] %in% boxplot.stats(mds[,k])$out)
    if(length(koutliers) > minsize) {
      start <- koutliers[1]
      end <- koutliers[1]
      for(i in c(2:length(koutliers))) {
        if(koutliers[i] < end+maxgap) {
          end = koutliers[i-1]
        } else {
          if(end-start > minsize) {
            #write(paste(mdsk,k,start/2,end/2,end-start,end-start>minsize),file=stderr())
            ii <- nrow(invcands)+1
            invcands[ii,"inversion"] <- ii
            invcands[ii,"chrom"] <- chrom
            invcands[ii,"start"] <- start*pcblocksize
            invcands[ii,"end"] <- end*pcblocksize
            }
          start=koutliers[i]
          end=koutliers[i]
        }
      }
      if(end-start > minsize) {
        #write(paste(mdsk,k,start/2,end/2,end-start,end-start>minsize),file=stderr())
        ii <- nrow(invcands)+1
        invcands[ii,"inversion"] <- ii
        invcands[ii,"chrom"] <- chrom
        invcands[ii,"start"] <- start*pcblocksize
        invcands[ii,"end"] <- end*pcblocksize
      }
    }
}

axlab <- as_mapper(~ .x /1e6)

  
distplot <- ggplot(pcdistflat,aes(x=pos.x,y=pos.y,fill=value)) + geom_raster() + coord_fixed() +
    geom_rect(aes(xmin=start,xmax=end,ymin=start,ymax=end),
              data=invcands[invcands$chrom==chrom,],
              inherit.aes=F,fill=NA,color="orange",alpha=0.5) +
    scale_x_continuous(expand = c(0,0),labels=axlab) + scale_y_continuous(expand = c(0,0),labels=axlab) +
    theme(legend.position="none") + ggtitle(paste("chrom",chrom,"mdsk:",mdsk," gap:",maxgap," size",minsize))

#distplot


invcands$country=country
write.table(invcands,file=outtxt,sep="\t",quote=F,col.names=T,row.names=F)

invcands <- unique(invcands[order(invcands$chrom,invcands$start,invcands$end),c("chrom","start","end")])
invcands$chromname <- chromname[invcands$chrom]

rm("pcs")
for(ii in c(1:nrow(invcands))) {
    chr = invcands[ii,"chrom"]
    st = invcands[ii,"start"]/1e6
    en = invcands[ii,"end"]/1e6
    write(paste("PCA: ",ii,chr,st,en),file=stderr())
    invsnps <- vcf_query(paste("SNPs/lostruct_chr",chr,".vcf.gz",sep=""),
                         samples=samples,
                         regions=invcands[ii,c("chromname","start","end")])
    goodsnps <- invsnps[apply(invsnps,1,function(x) {!any(is.na(x))}),]
    pca <- prcomp(t(goodsnps))
    
    invpcs <- as.data.frame(pca$x[,c("PC1","PC2")])
    invpcs$sample <- samples
    invpcs <- merge(invpcs,spptab)
    invpcs$inv <- ii
    invpcs$chrom <- chr
    
    if(exists("pcs")) {
      pcs <- rbind(pcs,invpcs)
    } else {
      pcs <- invpcs
    }
    #spca <- summary(pca)
}

ncol=round(sqrt(length(unique(pcs$inv))))
invpca <- ggplot(pcs,aes(x=PC1,y=PC2)) + geom_point() + coord_fixed() +
  facet_wrap("inv ~ .",ncol=ncol)

distplot | invpca
ggsave(outpng,width=350,height=200,units="mm")
