
library("tidyverse")
#library("patchwork")

library("pvclust")
library("getopt")
library("gridExtra")


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


chromlen <- c(310827022,474425716,409777670)

blocksize <- 5e5
chrom = 1

#distfile <- "data/lostruct_chr1_Trinidad.txt"
cmeth <- "single"
dmeth <- "uncentered"
minP <- 0.99
filterP <- 0.95
minblocks <- 5
maxblocks <- 300
nboots <- 1000
cores <- 0

opttab <- as.matrix(data.frame("long"=c("infile","outfile","chr",
                            "cluster_method","distance_method",
                            "cores",
                            "p_value","nboots",
                            "minblocks","maxblocks","blocksize"),
                   "short"=c("i","o","c",
                            "m","M",
                            "C",
                            "p","n",
                            "b","B","s"),
                   "argument"=as.character(c(1,1,1,
                             2,2,
                             2,
                             2,2,
                             2,2,2)),
                   "type"=c("character","character","integer",
                            "character","character",
                            "integer",
                            "double","integer",
                            "integer","integer","integer")))


opt <- getopt(opttab,)

infile <- opt$infile
outfile <- opt$outfile
chrom <- opt$chr

if (!is.null(opt$cluster_method) ) {cmeth <- opt$cluster_method}
if (!is.null(opt$distance_method) ) {dmeth <- opt$cluster_method}
if (!is.null(opt$p_value) ) {minP <- opt$p_value}
if (!is.null(opt$cores) ) {cores <- opt$cores}
if (!is.null(opt$minblocks) ) {minblocks <- opt$minblocks}
if (!is.null(opt$maxblocks) ) {maxblocks <- opt$maxblocks}
if (!is.null(opt$blocksize) ) {blocksize <- opt$blocksize}
if (!is.null(opt$blocksize) ) {blocksize <- opt$blocksize}



####
#read and parse PCA distance matrix
####

pcdists <- read.table(infile,header=T)
write(paste("found",dim(pcdists)[1],"blocks in ",infile),stderr())

#strip to just distance matrix
pcdists <- pcdists[,c(1:dim(pcdists)[1])]

pos <- blocksize * 1:dim(pcdists)[1]
blocks <- paste(chrom,pos,sep=":")
regions <- data.frame(
                      "chrom"=rep(chrom,length(pos)),
                      "pos"=pos,
                      "block"=blocks)
colnames(pcdists) <- blocks
rownames(pcdists) <- blocks

#######
# make distance plot
#######

pcdistflat <- pivot_longer(cbind(pcdists,blocks),cols=all_of(blocks),names_to = "y") %>% rename("block"="blocks")
pcdistflat <- merge(merge(pcdistflat,regions,by="block"),
                    regions,by.x="y",by.y="block",suffixes = c(".x",".y"))
distplot <- ggplot(pcdistflat,aes(x=pos.x,y=pos.y,fill=value)) +
  geom_tile() + xlim(0,chromlen[chrom])


######
# h-cluster, identify significant clusters
######

if(cores==0) {
  write(paste("running",
              dmeth,"/",cmeth," clustering with",
              nboots,"bootstraps on all cores"),stderr())
  distpv <- pvclust(pcdists,method.dist=dmeth,method.hclust=cmeth,nboot=nboots,parallel=T)
} else if(cores>1) {
  write(paste("running",
              dmeth,"/",cmeth," clustering with",
              nboots,"bootstraps on",cores,"cores"),stderr())
  distpv <- pvclust(pcdists,method.dist=dmeth,method.hclust=cmeth,nboot=nboots,parallel=cores)
} else {
  write(paste("running",
              dmeth,"/",cmeth," clustering with",
              nboots,"bootstraps on 1 core"),stderr())
  distpv <- pvclust(pcdists,method.dist=dmeth,method.hclust=cmeth,nboot=nboots,parallel=F)
}


#extract nodes on tree with Pvals > minP
edgetable <- distpv$edges
goodedges <- as.numeric(rownames(edgetable[edgetable$au>filterP,]))

#recurse through merges to get all children
disttree <- distpv$hclust
merges <- disttree$merge

#find all children of high-P nodes
clust=list()
for(node in goodedges) {
  kids = sort(getchildren(node, merges))

  if(length(kids) >= minblocks & length(kids) < maxblocks){
    clust[as.character(node)] = list(kids)
  }
}

#make data frame of blocks / clusters
clustdf<-data.frame("block"=numeric(),
                    "cluster"=character())
for(cid in names(clust)){
  #get label IDs from distance tree:
  cblocks=disttree$labels[clust[[cid]]]
  cluster=rep(cid,length(cblocks))
  clustdf <- rbind(clustdf,data.frame("block"=cblocks,
                                        "cluster"=cluster))
  }

#add back regions + p-values
clustdf <- merge(clustdf,regions,by="block",)
edgetable$cluster = rownames(edgetable)
clustdf <- merge(clustdf,edgetable,by="cluster")


clustplotP <- ggplot(subset(clustdf,au>=minP),aes(x=pos,fill=au,y=as.factor(cluster))) + geom_tile() +
  xlim(0,chromlen[chrom]) + ylab("cluster") + ggtitle(paste(cmeth,"/",dmeth,", ",
                                minblocks,"-",maxblocks," blocks ",
                                "P>",minP," ",nboots," boots",sep=""))

grid.arrange(clustplotP,distplot,ncol=1,heights=c(2,5))


#ggsave(paste("trinidad_chr1_pvclust_",cmeth,"_",dmeth,"_hier_P",minP,"_n",nboots,".png",sep=""))
ggsave(paste(outfile,".png",sep=""))

meandists <- data.frame("cluster"=integer(),
                        "meandist"=numeric(),
                        "lowdist"=logical())
#lower bound for outliers
distlim <-  boxplot.stats(pcdists[upper.tri(pcdists)])$stats[1]
for (C in unique(clustdf$cluster)) {
  i<- nrow(meandists)+1
  cblocks <- clustdf$block[clustdf$cluster==C]
  cblockdists <- pcdists[cblocks,cblocks]
  cblockdists <- cblockdists[upper.tri(cblockdists)]
  meandists[i,"cluster"] <- C
  meandists[i,"meandist"] <- mean(cblockdists)
  meandists[i,"lowdist"] <- (mean(cblockdists) < distlim)
}

clustdf <- merge(clustdf,meandists)
write.table(clustdf,file=paste(outfile,".txt",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
