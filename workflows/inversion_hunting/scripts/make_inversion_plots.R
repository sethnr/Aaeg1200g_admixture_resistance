
library("tidyverse")

#library("patchwork")
library("gridExtra")
library("grid")

library("getopt")

opttab <- matrix(c("inversions","i","1","character",
                   "aims","a","1","character",
                   "dists","d","1","character",
                   "meta","m","1","character",
                   "out","o","1","character"
),byrow=T,ncol=4)
opt <- getopt(opttab)

metafile <- opt$meta
distfile <- opt$dist
blockfile <- opt$inversions
aimsfile <- opt$aims
outprefix <- opt$out

#distfile <- "test_merged_chr1_Gabon_aims.txt"
# aimsfile <- "test_merged_chr1_Gabon_aims.txt"
# blockfile <- "test_merged_chr1_Gabon_blocks.txt"
# metafile <- "resources/meta_Aaeg1kg_spp.txt"
# outprefix <- "test_merged_aims_chr1_Gabon"

outaims <- paste(outprefix,"aims.png",sep="_")
outdists <- paste(outprefix,"dists.png",sep="_")

#blocksize<-5e05
# chromname <- c("NC_035107.1","NC_035108.1","NC_035109.1")
# chromlen <- c(310827022,474425716,409777670)
# names(chromlen) <- chromname

write(file.size(aimsfile),stderr())

if(file.size(aimsfile)==0L) {
  file.create(outaims)
  file.create(outdists)
  write(paste("no aims in file",aimsfile,"\n","writing empty files for",outaims,outdists),stderr())
  quit("no",0)
}


aims <- read.table(aimsfile,header=T)

metatab <- read.table(metafile,header=T, sep="\t")
samples <- metatab$sample
samples <- samples[samples %in% colnames(aims)]



chromlen <- c(310827022,474425716,409777670)
blocksize<-5e05

#######
# make inversion positions plots
#######


invblocks <- read.table(blockfile,header=T)
invblocks <- invblocks[order(invblocks$chrom,invblocks$pos),]

chrom <- unique(invblocks$chrom)[[1]]

invedges <- invblocks %>% group_by(inv) %>% summarize(st = min(pos), en = max(pos))

invplot <- ggplot(invblocks,aes(x=pos,y=as.factor(inv))) + geom_tile(height=0.9) +
  geom_text(data=subset(invedges,st<=(chromlen[chrom]/2)), aes(x=en,label=inv),hjust=-0.2) +
  geom_text(data=subset(invedges,st>(chromlen[chrom]/2)), aes(x=st,label=inv),hjust=1.2) +
  scale_x_continuous(limits=c(0,chromlen[chrom]),expand = c(0,0,0,0)) + 
  theme(axis.title=element_blank(),axis.text.y=element_blank())

aimplots=list()
for(invname in unique(invblocks$inv)) {
  invnamesafe <- paste("X",gsub("\\D",".",invname,perl=T),sep="")
  compinvids <- as.numeric(strsplit(invname,"/")[[1]])
  
  #get SNPs for inversion, remove duplicates, re-index
  invsnps <- subset(aims,inv == invname)
  meansnp <- apply(invsnps[,samples],2,FUN=function(x) {mean(na.omit(x))})
  
  cntinvorder <- metatab$sample[order(metatab$contgroup,metatab$region,metatab$country,meansnp[metatab$sample])]
  cntorder <- unique(metatab$country[order(metatab$contgroup,metatab$region)])
  
  aimsM <- pivot_longer(invsnps,all_of(samples),names_to = "sample") %>% rename("invcountry"="country")
  aimsM <- merge(aimsM,metatab,by="sample")
  aimsM$sample <- factor(aimsM$sample,levels = cntinvorder,ordered=T)
  aimsM$country <- factor(aimsM$country,levels = cntorder,ordered=T)
  aimsM$cncode <- aimsM$country
  levels(aimsM$cncode) <- substr(levels(aimsM$country),0,3)
  
  aimplot <- ggplot(aimsM,aes(x=i,y=as.numeric(sample),fill=as.factor(value))) + geom_raster() + 
    ylab("samples") + xlab("SNPs")+ theme(legend.position="none")+
    scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) +
    scale_fill_manual(values=c("0"="blue","1"="purple","2"="red")) +
    facet_grid("cncode ~ .",scale="free_y",space="free_y") +
    ggtitle(paste(invname,"aims","(",nrow(invsnps),")"))+
    theme(panel.spacing = unit(0.2, "mm"),
          axis.text.y=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks.y=element_blank())
  aimplots[[invnamesafe]] <- aimplot
  
  
  
}

# aimsM <- pivot_longer(aims,all_of(samples),names_to = "sample") %>% rename("invcountry"="country")
# aimsM <- merge(aimsM,metatab,by="sample")
# aimsM$sample <- factor(aimsM$sample,levels = cntinvorder,ordered=T)
# aimsM$country <- factor(aimsM$country,levels = cntorder,ordered=T)
# aimsM$cncode <- aimsM$country
# levels(aimsM$cncode) <- substr(levels(aimsM$country),0,3)
# 
# aimplotC <- ggplot(aimsM,aes(x=i,y=as.numeric(sample),fill=as.factor(value))) + geom_raster() +
#   ylab("samples") + xlab("SNPs")+ theme(legend.position="none")+
#   scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) +
#   scale_fill_manual(values=c("0"="blue","1"="purple","2"="red")) +
#   facet_grid("cncode ~ .",scale="free_y",space="free_y") +
#   theme(panel.spacing = unit(0.2, "mm"),
#         axis.text.y=element_blank(),
#         axis.title.y=element_blank(),
#         axis.ticks.y=element_blank()) +
#   facet_grid(cncode ~ inv,scale="free",space="free")
# aimplotC

aimplot <- do.call("arrangeGrob", c(aimplots, nrow=1))

png(outaims,res=400,width=400,height=200,units='mm')
grid.arrange(aimplot,invplot,heights=c(8,2))
dev.off()




#####
# make pc distance plot
#####

pcdists <- read.table(distfile,header=T)
write(paste("found",dim(pcdists)[1],"blocks in ",distfile),stderr())

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


invedges$mid <- invedges$st+(invedges$en-invedges$st)/2

pcdistflat <- pivot_longer(cbind(pcdists,blocks),cols=all_of(blocks),names_to = "y") %>% rename("block"="blocks")
pcdistflat <- merge(merge(pcdistflat,regions,by="block"),
                    regions,by.x="y",by.y="block",suffixes = c(".x",".y"))


invedges$label="L"
invedges$label[invedges$mid < chromlen[chrom]*0.1] <- "R"
invedges$label[invedges$mid > chromlen[chrom]*0.5 & invedges$mid < chromlen[chrom]*0.9] <- "R"
#invedges 

distplot <- ggplot(pcdistflat,aes(x=pos.x,y=pos.y,fill=value)) +
  geom_tile() + 
  scale_x_continuous(limits=c(0,chromlen[chrom]),expand = c(0,0,0,0)) +
  scale_y_continuous(limits=c(0,chromlen[chrom]),expand = c(0,0,0,0)) +
  geom_text(data=subset(invedges,label=="L"),aes(y=mid,x=mid-3e7,label=inv),inherit.aes=F,color="orange") +
  geom_text(data=subset(invedges,label=="R"),aes(y=mid,x=mid+3e7,label=inv),inherit.aes=F,color="orange") +
  geom_rect(data=invedges,aes(xmin=st,xmax=en,ymin=st,ymax=en),inherit.aes=F,fill=NA,color="orange") +
  coord_fixed()
distplot


png(outdists,res=400,width=220,height=200,units='mm')
distplot
dev.off()
