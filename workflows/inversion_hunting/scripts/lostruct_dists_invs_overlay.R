library(tidyverse)
#library(patchwork)
library(gridExtra)

args = commandArgs(trailingOnly=TRUE)

dists = args[1]
pcs = args[2]
invs = args[3]
chrid = as.numeric(args[4])
spp = tolower(args[5])
outfile1 = args[6]
outfile2 = args[7]

write(paste("dists =",dists),file=stderr())
write(paste("pcs =",pcs),file=stderr())
write(paste("invs =",invs),file=stderr())
write(paste("chrom =",chrid),file=stderr())
write(paste("spp =",spp),file=stderr())
write(paste("outfile1 =",outfile1),file=stderr())
write(paste("outfile2 =",outfile2),file=stderr())

chromname <- c("NC_035107.1","NC_035108.1","NC_035109.1")

#print PCAs every <pcsep> bp along chromosome
pcsep = 1e7


#######
# get inversion calls
#######

invcalls <- read.table(invs,header=T,sep="\t")
colnames(invcalls) <- tolower(colnames(invcalls))
invcalls$chrom <- chromname[invcalls$chrom]

invcalls <- subset(invcalls,chrom==chromname[chrid])
#write.table(invcalls,stderr())
######
# get lostruct distances
######

pcdistdf <- read.table(dists)
blocks <- colnames(pcdistdf)
allposns <- t(data.frame(strsplit(blocks,"\\.")))
allposns <- data.frame("block"=blocks,
                       "chrom"=paste(allposns[,1],allposns[,2],sep="."),
                       "pos"=as.numeric(allposns[,3]))

pcdistdf$x <- blocks
pcdistflat <- pivot_longer(pcdistdf,cols=all_of(blocks),names_to = "y")
pcdistflat <- merge(merge(pcdistflat,allposns,by.x="x",by.y="block"),allposns,by.x="y",by.y="block",suffixes = c(".x",".y"))

#create dir if not exists
dir.create(dirname(outfile),recursive=T)
#plot distances with inversions overlaid
#write(paste("writing to",paste(outfile,"_invs_overlay.png",sep="")),stderr())
#png(paste(outfile,"_invs_overlay.png",sep=""),width=7,height=7,res=400)
write(paste("writing to",outfile1),stderr())
png(outfile1,sep=""),width=7,height=7,res=400)

pcdistsplot <- ggplot(pcdistflat,aes(x=pos.x,y=pos.y,fill=value)) + geom_raster() + 
  geom_rect(aes(xmin=start,xmax=end,ymin=start,ymax=end),
            data=invcalls,
            inherit.aes=F,fill=NA,color="orange") +
  ggtitle(paste(chrid,"lostruct invs",spp)) + 
  coord_fixed()
pcdistsplot
dev.off()
#ggsave(,plot=pcdistsplot)



# plot with regular PCAs
pcpos <- data.frame(x=seq(0,max(pcdistflat$pos.x),by=pcsep)+1,
                    y=seq(0,max(pcdistflat$pos.x),by=pcsep)+1,
                    ord=c(0:floor(max(pcdistflat$pos.x)/pcsep)))

pcs <- read.table(pcs,header=T)

#get every 10mth block
pcs <- pcs[pcs$start %% pcsep==1,]

pc1s <- colnames(pcs)[substr(colnames(pcs),0,4)=="PC_1"]
pc2s <- colnames(pcs)[substr(colnames(pcs),0,4)=="PC_2"]
headcols <- colnames(pcs)[1:6]
pc1tab <- pivot_longer(pcs[,c(headcols,pc1s)],cols=all_of(pc1s),names_to = "sample",values_to = "pc1")
pc1tab$sample <- gsub("PC_1_","",pc1tab$sample)
pc2tab <- pivot_longer(pcs[,c(headcols,pc2s)],cols=all_of(pc2s),names_to = "sample",values_to = "pc2")
pc2tab$sample <- gsub("PC_2_","",pc2tab$sample)
pcslong <- merge(pc1tab,pc2tab)
pcslong$label <- (pcslong$start-1)/pcsep

pcsplot <- ggplot(pcslong,aes(x=pc1,y=pc2)) + geom_point(size=0.2) + 
  facet_wrap("label",ncol=5) 

combplot <- arrangeGrob(
  pcdistsplot + 
    geom_text(data=pcpos,aes(x=x,y=y,label=ord),inherit.aes = F,color="red") + 
    theme(legend.position="bottom"),
  pcsplot,
  ncol=2,widths=c(5,5))

#write(paste("writing to",paste(outfile,"_pca_composite.png",sep="")),stderr())
#ggsave(paste(outfile,"_pca_composite.png",sep=""),plot=combplot, width=12, height=7,dpi = 400)
write(paste("writing to",outfile2),stderr())
ggsave(outfile2,plot=combplot, width=12, height=7,dpi = 400)
