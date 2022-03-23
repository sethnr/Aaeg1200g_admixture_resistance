
library("tidyverse")
library("patchwork")
library("lostruct")


args = commandArgs(trailingOnly=TRUE)

setwd("~/Gits/Aaeg1000g_analyses/analyses/sredmond/220118_inversion_region_calls/")
distsf <- "redmond-lab-aaeg1000g/results_lostruct/lostruct_all_chr/"
snpsf <- "/Volumes/Mosquito_raw_data/Aedes/Aaeg1000g/thinrand/"
sppfile <- "resources/meta_Aaeg1kg_spp.txt"
invfile <- "resources/redmond_2020_inversion_calls.txt"


chromname <- c("NC_035107.1","NC_035108.1","NC_035109.1")
chromlen <- c(310827022,474425716,409777670)
names(chromlen) <- chromname

pcblocksize <- 5e5

spptab <- read.table(sppfile,header=T, sep="\t")

rm("allpcdists")
allplots = list();

country <- "Senegal"
samples <- spptab$sample[spptab$country==country]

distfolderlist <- list.files(distsf,full.names = T)
alldistfiles <- vector()
for(chrom in c(1,2,3)) {
  #for(country in c(region,countries)) {
  cdist <- distfolderlist[grep(paste("chr",chrom,"_",country,".txt",sep=""),distfolderlist)]
  if(length(cdist)==1) {
    write(paste(country,"->",cdist[1]),file=stderr())
    alldistfiles[country] <- cdist[1]
  } else {
    write(paste("found ",length(cdist),"files for ",country),file=stderr())
  }

  
  write(paste("parsing",country,"::",alldistfiles[country]),file=stderr())
  regionfile <- alldistfiles[country]
  regiondf <- read.table(regionfile,header=T)
  pcdistdf <- regiondf[,c(1:dim(regiondf)[1])]

  pos <- pcblocksize * 1:dim(regiondf)[1]
  blocks <- paste(chrom,pos,sep=":")

  regions <- data.frame(
                      "chrom"=rep(chrom,length(pos)),
                      "pos"=pos,
                      "block"=blocks)
  colnames(pcdistdf) <- blocks

  pcdistdf$x <- blocks
  pcdistflat <- pivot_longer(pcdistdf,cols=all_of(blocks),names_to = "y")
  pcdistflat <- merge(merge(pcdistflat,regions,by.x="x",by.y="block"),regions,by.x="y",by.y="block",suffixes = c(".x",".y"))
  pcdistflat$set <- country
  
  if(exists("allpcdists")) {
    allpcdists <- rbind(allpcdists,pcdistflat)
  } else {
    allpcdists <- pcdistflat
  }

}
  

invcands <- data.frame("inversion"=character(),
                       "chrom"=numeric(),
                       "start"=numeric(),
                       "end"=numeric())


#find outliers in PC dists
pcdistplots <- list()
for(chrom in c(1:3)) {
  
  pcdistscnt <- allpcdists[allpcdists$set==country & allpcdists$chrom.x==chrom,]
  
  #remake indices from posn x/y
  pcdistscnt$y <- pcdistscnt$pos.y/pcblocksize
  pcdistscnt$x <- pcdistscnt$pos.x/pcblocksize
  
  pcdistmat <- matrix(nrow=max(pcdistscnt$x),ncol=max(pcdistscnt$y))
  for(i in c(1:nrow(pcdistscnt))) {
    pcdistmat[pcdistscnt$y[i],pcdistscnt$x[i]] <- pcdistscnt$value[i] 
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
  
  
  p <- ggplot(pcdistscnt,aes(x=pos.x,y=pos.y,fill=value)) + geom_raster() + coord_fixed() +
    geom_rect(aes(xmin=start,xmax=end,ymin=start,ymax=end),
              data=invcands[invcands$chrom==chrom,],
              inherit.aes=F,fill=NA,color="orange",alpha=0.5) +
    scale_x_continuous(expand = c(0,0),labels=axlab) + scale_y_continuous(expand = c(0,0),labels=axlab) +
    theme(legend.position="none") + ggtitle(paste("chrom",chrom,"mdsk:",mdsk," gap:",maxgap," size",minsize))
  pcdistplots[[chrom]] <- p
}




invcands <- unique(invcands[order(invcands$chrom,invcands$start,invcands$end),c(2:4)])
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

mychrom <- 3
ncol=round(sqrt(length(unique(pcs$inv[pcs$chrom==mychrom]))))
invpca <- ggplot(subset(pcs,chrom==mychrom),aes(x=PC1,y=PC2)) + geom_point() + coord_fixed() +
  facet_wrap("inv ~ .",ncol=ncol)
pcdistplots[[mychrom]] | invpca
