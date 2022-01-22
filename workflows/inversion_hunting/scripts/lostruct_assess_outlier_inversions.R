
library("tidyverse")
library("lostruct")

#library("patchwork")
library("gridExtra")
library("grid")

args = commandArgs(trailingOnly=TRUE)

# setwd("~/Gits/Aaeg1000g_analyses/analyses/sredmond/220118_inversion_region_calls/")
# dists <- "redmond-lab-aaeg1000g/results_lostruct/lostruct_all_chr/lostruct_chr1_Uganda.txt"
# vcffile <- "/Volumes/Mosquito_raw_data/Aedes/Aaeg1000g/thinrand/lostruct_chr1.vcf.gz"
# country <- "Uganda"
# chrom <- 2
# outtxt <- "./inv_candidates_chr2_Uganda.txt"
# outpng <- "./inv_candidates_chr2_Uganda.png"
# sppfile <- "resources/meta_Aaeg1kg_spp.txt"
# invfile <- "resources/redmond_2020_inversion_calls.txt"
# 

dists <- args[1]
vcffile <- args[2]
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

write(paste("found",length(samples),"samples for",country),file=stderr())

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
                       "end"=numeric(),
                       "country"=character())


#remake indices from posn x/y
pcdistflat$y <- pcdistflat$pos.y/pcblocksize
pcdistflat$x <- pcdistflat$pos.x/pcblocksize

pcdistmat <- matrix(nrow=max(pcdistflat$x),ncol=max(pcdistflat$y))
for(i in c(1:nrow(pcdistflat))) {
  pcdistmat[pcdistflat$y[i],pcdistflat$x[i]] <- pcdistflat$value[i]
}



#calculate MDS values
mdsk <- 5
#set NAs to 1
pcdistmat[is.na(pcdistmat)] <- 1
mds <- cmdscale(pcdistmat,k=mdsk)



#run through MDS to find outlier regions, assess for consistency (at least 4 blocks long and not separated by 10 blocks)
maxgap=10
minsize=4
outliersT <- mds==NA

for(k in c(1:mdsk)) {
    #write(paste("surveying MDS",k),file=stderr())
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
            invcands[ii,"country"] <- country
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
        invcands[ii,"country"] <- country
      }
    }
}


axlab <- as_mapper(~ .x /1e6)

write(paste("found",nrow(invcands),"inversion candidates on chrom",chrom),file=stderr())

if(nrow(invcands)>0) {
  distplot <- ggplot(pcdistflat,aes(x=pos.x,y=pos.y,fill=value)) + geom_raster() + coord_fixed() +
      geom_rect(aes(xmin=start,xmax=end,ymin=start,ymax=end),
                data=invcands[invcands$chrom==chrom,],
                inherit.aes=F,fill=NA,color="orange",alpha=0.5) +
      scale_x_continuous(expand = c(0,0),labels=axlab) + scale_y_continuous(expand = c(0,0),labels=axlab) +
      theme(legend.position="none") + ggtitle(paste("chrom",chrom,"mdsk:",mdsk," gap:",maxgap," size",minsize))
} else {
  distplot <- ggplot(pcdistflat,aes(x=pos.x,y=pos.y,fill=value)) + geom_raster() + coord_fixed() +
      scale_x_continuous(expand = c(0,0),labels=axlab) + scale_y_continuous(expand = c(0,0),labels=axlab) +
      theme(legend.position="none") + ggtitle(paste("chrom",chrom,"mdsk:",mdsk," gap:",maxgap," size",minsize))
}


invcands <- unique(invcands[order(invcands$chrom,invcands$start,invcands$end),c("chrom","start","end")])
invcands$chromname <- chromname[invcands$chrom]
invcands$valid <- logical()
invcands$name <- character()

rm("pcs")
#for(ii in c(1:nrow(invcands))) {
for(ii in rownames(invcands)) {
    chr = invcands[ii,"chrom"]
    st = invcands[ii,"start"]/1e6
    en = invcands[ii,"end"]/1e6
    
    invcands[ii,"name"] <- paste(chr,":",st,"-",en,sep="")
    
    write(paste("PCA: ",ii,chr,st,en),file=stderr())
    invsnps <- vcf_query(vcffile,
                         samples=samples,
                         regions=invcands[ii,c("chromname","start","end")])
    goodsnps <- invsnps[apply(invsnps,1,function(x) {!any(is.na(x))}),]
    pca <- prcomp(t(goodsnps))

    invpcs <- as.data.frame(pca$x[,c("PC1","PC2")])
    invpcs$sample <- samples
    invpcs <- merge(invpcs,spptab)
    invpcs$inv <- paste(chr,":",st,"-",en,sep="")
    invpcs$chrom <- chr

    if(exists("pcs")) {
      pcs <- rbind(pcs,invpcs)
    } else {
      pcs <- invpcs
    }
}

#### assess PCA clusters
if(exists("pcs")) {
  pcs$inv <- factor(pcs$inv,levels=invcands$name,ordered=T)
  pcs$valid<-factor(NA,levels=c("aa","ab","bb"))
  for(I in unique(pcs$inv)) {
    pcsinv <- subset(pcs,inv==I)
    kmpca <- kmeans(pcsinv[,c("PC1")],centers=3,nstart=20,iter.max=50)
    
    clusters <- factor(kmpca$cluster)
    clustorder <- order(kmpca$centers)
    
    #calculate deviation of middle cluster from center point
    centers <- kmpca$centers[clustorder]
    delta <- ((centers[2] - centers[1])-(centers[3] - centers[2])) / 
      ((centers[3]-centers[1])/2)
    
    #rename levels 0/1/2 based on pcorder
    levels(clusters)[clustorder] <- c('aa','ab','bb')
      
    n<-length(clusters)
    wss <- kmpca$withinss[clustorder]
    tot.wss <- kmpca$tot.withinss
    bss <- kmpca$betweenss
    tot.ss <- kmpca$totss
    
    
    
    MAXD <- 0.25
    MAXSS <- 20
    MINBSS <- 0.95
    
    invpass <- ((abs(delta)<=MAXD) && ((tot.wss/n) <= MAXSS) && (bss/tot.ss)>=MINBSS)
    write(paste("inv:",I,
                #"\tctr:",paste(round(centers,1),collapse="/"),
                #"\twss:",paste(round(wss,1),collapse="/"),
                "\tbss/tot:",round(bss/tot.ss,2),
                "\tss",round(tot.wss/n,2),
                "d",round(abs(delta),2),
                "\tpass:",invpass),file=stderr())
    
    invcands[invcands$name==I,"valid"] <- invpass
    if(invpass) {
      pcs$valid[which(pcs$inv==I)] <- clusters
    }
  }
}


write(paste("writing",nrow(invcands),"candidates on chrom",chrom),file=stderr())
write.table(invcands,file=outtxt,sep="\t",quote=F,col.names=T,row.names=F)





if(exists("pcs")) {
  invcols <- scale_color_manual(values=c("aa"="yellow","ab"="orange","bb"="red"),na.value = "dark grey")
  ncol=round(sqrt(length(unique(pcs$inv))))
  ncol=round(sqrt(length(unique(pcs$inv))))
  invpca <- ggplot(pcs,aes(x=PC1,y=PC2,color=valid)) + geom_point() + coord_fixed() + invcols +
    facet_wrap("inv ~ .",ncol=ncol) + theme(legend.position = "none")

  invpca
  
  png(filename = outpng,width=350,height=200,units="mm",res=400)
  grid.arrange(distplot, invpca, ncol=2)
  dev.off()
} else {
  png(filename = outpng,width=350,height=200,units="mm",res=400)
  grid.arrange(distplot,  grid.rect(gp=gpar(col="white")), ncol=2)
  dev.off()
}
