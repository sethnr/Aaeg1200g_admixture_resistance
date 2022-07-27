
library("tidyverse")
library("patchwork")

library("pvclust")


args = commandArgs(trailingOnly=TRUE)

setwd("~/Gits/Aaeg1000g_analyses/analyses/sredmond/220617_inversion_region_improvements/")
distsf <- "data/"
sppfile <- "resources/meta_Aaeg1kg_spp.txt"
invfile <- "resources/redmond_2020_inversion_calls.txt"


#parse children from merges
getchildren <- function(x, mergetable=merges) {
  a = mergetable[x,1]
  b = mergetable[x,2]
  children = c()
  #print(paste(a,b),stderr())
  if(a < 0) {children = c(children,(a*-1))}
  else {children = c(children,getchildren(a, mergetable))}
  
  if(b < 0) {children = c(children,(b*-1))}
  else {children = c(children,getchildren(b, mergetable))}
  return(children)
}




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

cmeth <- "ward.D2"
dmeth <- "uncentered"
minP <- 0.99
minblocks <- 5
maxblocks <- 300 
nboots <- 1000

distpv <- pvclust(pcdists,method.dist=dmeth,r=c(1),method.hclust=cmeth,nboot=nboots,parallel=T)

#manually look at all nodes on tree with high Pvals
disttree <- distpv$hclust

edgetable <- distpv$edges
#edgetable <- subset(edgetable,au>minP)
goodedges <- as.numeric(rownames(edgetable[edgetable$au>minP,]))

#recurse through merges to get all children
merges <- disttree$merge
#goodmerges <- merges[goodedges,]

clust=list()
for(node in goodedges) {
  kids = sort(getchildren(node, merges))
  
  if(length(kids) >=minblocks & length(kids) < maxblocks){
    write(paste(node,length(kids),sep="::"),stderr())
    write(paste(kids,collapse="/"),stderr())
    clust[as.character(node)] = list(kids)
  }
}


clustdf<-data.frame("block"=numeric(),
                    "cluster"=character())
i=0
for(cid in names(clust)){
  i=i+1
  #get label IDs:
  cblocks=disttree$labels[clust[[cid]]]
  cluster=rep(cid,length(cblocks))
  clustdf <- rbind(clustdf,data.frame("block"=cblocks,
                                        "cluster"=cluster))
  }

clustdf <- merge(clustdf,regions,by="block",)
# clustplot <- ggplot(clustdf,aes(x=pos,y=as.factor(cluster))) + geom_tile() + 
#   xlim(0,35e07) + ggtitle(paste(cmeth,"/",dmeth,", ",
#                                 minblocks,"-",maxblocks," blocks ",
#                                 "P>",minP," ",nboots," boots",sep=""))
# 
# clustplot / distplot + plot_layout(heights=c(2,5))

edgetable$cluster = rownames(edgetable)

clustdf <- merge(clustdf,edgetable,by="cluster")
clustplotP <- ggplot(clustdf,aes(x=pos,fill=au,y=as.factor(cluster))) + geom_tile() + 
  xlim(0,35e07) + ylab("cluster") + ggtitle(paste(cmeth,"/",dmeth,", ",
                                minblocks,"-",maxblocks," blocks ",
                                "P>",minP," ",nboots," boots",sep=""))
clustplotP / distplot + plot_layout(heights=c(2,5))
ggsave(paste("senegal_pvclust_",cmeth,"_",dmeth,"_hier_P",minP,"_n",nboots,".png",sep=""))

