
library("tidyverse")
library("patchwork")

library("pvclust")


args = commandArgs(trailingOnly=TRUE)

setwd("~/Gits/Aaeg1000g_analyses/analyses/sredmond/220617_inversion_region_improvements/")
distsf <- "data/"
sppfile <- "resources/meta_Aaeg1kg_spp.txt"
invfile <- "resources/redmond_2020_inversion_calls.txt"


chromname <- c("NC_035107.1","NC_035108.1","NC_035109.1")
chromlen <- c(310827022,474425716,409777670)
names(chromlen) <- chromname

pcblocksize <- 5e5
chrom = 1

spptab <- read.table(sppfile,header=T, sep="\t")

rm("allpcdists")
allplots = list();

country <- "Senegal"
samples <- spptab$sample[spptab$country==country]

distfile <- "data/lostruct_chr1_Senegal.txt"
  

pcdists <- read.table(distfile,header=T)
#strip to just distance matrix
pcdists <- pcdists[,c(1:dim(pcdists)[1])]


pos <- pcblocksize * 1:dim(pcdists)[1]
blocks <- paste(chrom,pos,sep=":")
regions <- data.frame(
                      "chrom"=rep(chrom,length(pos)),
                      "pos"=pos,
                      "block"=blocks)
colnames(pcdists) <- blocks
rownames(pcdists) <- blocks


#make distance plot
pcdistflat <- pivot_longer(cbind(pcdistdf,blocks),cols=all_of(blocks),names_to = "y") %>% rename("block"="blocks")
pcdistflat <- merge(merge(pcdistflat,regions,by="block"),
                    regions,by.x="y",by.y="block",suffixes = c(".x",".y"))
distplot <- ggplot(pcdistflat,aes(x=pos.x,y=pos.y,fill=value)) + 
  geom_tile() + xlim(0,35e07)

cmethods <- c("average", "ward.D", "ward.D2", "single", "complete", "mcquitty", "median", "centroid")
dmethods <- c("correlation", "uncentered", "abscor")
for(cmeth in cmethods) {
  for (dmeth in dmethods) {
  distpv <- pvclust(pcdists,method.dist=dmeth,method.hclust=cmeth,nboot=50,r=1)
  distpvpk <- pvpick(distpv)
  
  clustdf<-data.frame("block"=numeric(),
             "cluster"=numeric())
  i<-0
  for(c in 1:length(distpvpk$clusters)){
    cblocks=distpvpk$clusters[[c]]
    if(length(cblocks)>5) {
      i=i+1
      cluster=rep(i,length(cblocks))
      clustdf <- rbind(clustdf,data.frame("block"=cblocks,
                                        "cluster"=cluster))
    }
  }
  clustdf <- merge(clustdf,regions,by="block")
  clustplot <- ggplot(clustdf,aes(x=pos,y=cluster)) + geom_tile() + 
    xlim(0,35e07) + ggtitle(paste(cmeth,"/",dmeth))
  
  clustplot / distplot + plot_layout(heights=c(2,5))
  ggsave(paste("senegal_pvclust_",cmeth,"_",dmeth,".png",sep=""))
  # mindistplot <- ggplot(pcdistflat,aes(x=pos.x,y=pos.y,fill=value<mindist)) + geom_tile() + xlim(0,35e07)
  # clustplot / mindistplot + plot_layout(heights=c(2,5))
  }
}
