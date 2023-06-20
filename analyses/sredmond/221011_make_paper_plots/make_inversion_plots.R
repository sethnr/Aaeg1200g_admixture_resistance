library("tidyverse")

library("patchwork")
library("lostruct")
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


# chrom <- 1
# country <- "Gabon"
# invname <- "1_Gabon_1"
# clusters <- c("124","162")

chrom <- 1
country <- "Senegal"
invname <- "1_Senegal_1"
clusters <- c("291","194")



vcffile <- "resources/aegy.wgs.aaa.aaf.norep5-30x.ac10.thinrand01.chr1.vcf.gz"
distfile <- paste("data/lostruct_chr",chrom,"_",country,".txt.gz",sep="")
invcallsfile <- paste("data/invs_chr",chrom,"_inv_calls.txt",sep="")
aimsfile <- paste("data/invs_chr",chrom,"_aims.txt",sep="")
blockfile <- paste("data/regions_chr",chrom,"_",country,".txt",sep="")
metafile <- paste("resources/meta_Aaeg1kg_spp.txt",sep="")
outprefix <- paste("combined_plots_",country,sep="")

assesstab <- read.table(paste("data/assess_chr",chrom,"_",country,".txt",sep=""),header=T)
pccallstab <- read.table(paste("data/assess_chr",chrom,"_",country,"_calls.txt",sep=""),header=T)



metatab <- read.table(metafile,header=T, sep="\t")
samples <- metatab$sample[metatab$country==country]
#samples <- samples[samples %in% colnames(aims)]



chromlen <- c(310827022,474425716,409777670)
chromname <- c("NC_035107.1","NC_035108.1","NC_035109.1")
blocksize<-5e05

#######
# get hclust positions
#######


invblocks <- read.table(blockfile,header=T)
invblocks$cluster <- factor(invblocks$cluster)
invblocks$chromname <- chromname[invblocks$chrom]
invblocks$start <- invblocks$pos
invblocks$end <- invblocks$pos+blocksize

clusteredges <- invblocks %>% dplyr::group_by(cluster) %>% 
  dplyr::summarise(chrom=min(chrom),
                   start=min(pos),
                   end=max(pos)+blocksize,
                   size=round((blocksize+max(pos)-min(pos))/1e6,2)) %>% 
  mutate(region=paste(chrom,":",as.integer(start),"-",as.integer(end),sep="")) %>%
  dplyr::relocate(region,.before = size)


#####
# make pc distance plot
#####

pcdists <- read.table(distfile,header=T)
write(paste("found",dim(pcdists)[1],"blocks in ",distfile),stderr())

#strip to just distance matrix
pcdists <- pcdists[,c(1:dim(pcdists)[1])]

chrom <- 1
pos <- blocksize * 1:dim(pcdists)[1]
blocks <- paste(chrom,pos,sep=":")
regions <- data.frame(
  "chrom"=rep(chrom,length(pos)),
  "pos"=pos,
  "block"=blocks)
colnames(pcdists) <- blocks
rownames(pcdists) <- blocks


pcdistflat <- pivot_longer(cbind(pcdists,blocks),cols=all_of(blocks),names_to = "y") %>% dplyr::rename("block"="blocks")
pcdistflat <- merge(merge(pcdistflat,regions,by="block"),
                    regions,by.x="y",by.y="block",suffixes = c(".x",".y"))
pcdistflat$value[pcdistflat$value>1] <- 1

distplot <- ggplot(pcdistflat,aes(x=pos.x,y=pos.y,fill=value)) +
  geom_tile() +
  scale_x_continuous(limits=c(0,chromlen[chrom]),expand = c(0,0,0,0),breaks = seq(0,3e8,by=0.5e8),labels=seq(0,3e2,by=0.5e2)) +
  scale_y_continuous(limits=c(0,chromlen[chrom]),expand = c(0,0,0,0)) +
  scale_fill_continuous(limits=c(0,1),name="pc distance")+
  geom_rect(data=subset(clusteredges,size>=5),aes(xmin=start,xmax=end,ymin=start,ymax=end),inherit.aes=F,fill=NA,color="orange",size=0.5,alpha=0.6) +
  coord_fixed() + theme(axis.title = element_blank(),axis.text.y = element_blank(),legend.position = "bottom")
#distplot




################
#do PCA plots for selected inversion candidates
###############

# Run PCAs for all candidate inversion regions

pcaplots <- list()

if(exists("allpcs")){rm("allpcs")}
for(C in clusters) {
  invsnps <- vcf_query(vcffile,
                       samples=samples,
                       regions=invblocks[invblocks$cluster==C,c("chromname","start","end")])
  goodsnps <- invsnps[apply(invsnps,1,function(x) {!any(is.na(x))}),]
  write(paste("found",nrow(goodsnps),"snps for cluster",C),stderr())
  
  pca <- prcomp(t(goodsnps))
  
  pcs <- as.data.frame(pca$x[,c("PC1","PC2")])
  pcs$sample <- samples
  pcs <- merge(pcs,metatab,all.x=T)
  pcs$calls <- pccallstab[,paste("X",C,sep="")]
  pcs$angle <- assesstab$angle[assesstab$cluster==C]
  pcs$cluster<-paste("cluster",C)
  
  #plot according to assess table:
  if(exists("allpcs")) {
    allpcs <- rbind(allpcs,pcs)
  } else {
    allpcs <- pcs
  }
}

invcols <- scale_color_manual(values=c("aa"="blue2","ab"="purple1","bb"="red2"),na.value = "dark grey")


pcaplots <- ggplot(allpcs,aes(x=PC1,y=PC2,color=calls)) +
  geom_abline(aes(slope=(angle/45)*-1,intercept=0),linetype=2)+
  geom_point() + coord_fixed() + invcols +
  theme(legend.position = "none",axis.title=element_blank(),axis.text=element_blank()) +
  facet_wrap("cluster ~ .",ncol=1) 
#pcaplots




######
# plot aims
######
countries <- metatab$country
names(countries) <- metatab$sample
cntorder <- rev(unique(metatab$country[order(metatab$contgroup,metatab$region)]))

aims <- read.table(aimsfile,header=T)
aims <- aims[aims$cluster==invname,]
aimsamples <- colnames(aims)[8:1140]
aimsM <- pivot_longer(aims,all_of(aimsamples),
                      names_to = "sample")   #%>% rename("invcountry"="country")

write("sorting countries",stderr())
aimsM$country <- countries[aimsM$sample]
aimsM <- aimsM[!is.na(aimsM$country),]
aimsM$country <- factor(aimsM$country,levels = cntorder,ordered=T)

#order samples by country, then mean SNP call
meansnp <- apply(aims[,aimsamples],2,FUN=function(x) {mean(na.omit(x))})
snporder <- metatab$sample[order(metatab$country,meansnp[metatab$sample])]

aimsM$sample <- factor(aimsM$sample,levels = snporder,ordered=T)
aimsM$cncode <- aimsM$country
levels(aimsM$cncode) <- substr(levels(aimsM$country),0,3)
aimsM$y <- as.numeric(aimsM$sample)
aimsM$aim <- NA
aimsM$aim[aimsM$country==country] <- aimsM$value[aimsM$country==country]



  write("plotting AIMs",stderr())
aimplot <- ggplot(aimsM[aimsM$country==country,],aes(x=i,y=y,fill=as.factor(aim))) + geom_tile() +
  ylab("samples") + xlab("SNPs")+ theme(legend.position="none")+
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) +
  scale_fill_manual(values=c("0"="blue2","1"="purple2","2"="red2")) +
  facet_grid(". ~ 'assoc SNPs'",scale="free",space="free") +
  theme(panel.spacing = unit(0.2, "mm"),
        axis.text.y=element_blank(),
        axis.title=element_blank(),
        axis.ticks.y=element_blank())
#aimplot

snpplot <- ggplot(aimsM,aes(x=i,y=y,fill=as.factor(value))) + geom_tile() +
  ylab("samples") + xlab("SNPs")+ theme(legend.position="none")+
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) +
  scale_fill_manual(values=c("0"="blue2","1"="purple2","2"="red2")) +
  facet_grid("country ~ 'marker SNPs'",scale="free",space="free") +
  theme(panel.spacing = unit(0.2, "mm"),
        axis.text.y=element_blank(),
        axis.title=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text.y = element_blank())
#snpplot


callstab <- read.table(invcallsfile,header=T,sep="\t")
callstab <- callstab[,c("sample",paste("X",invname,sep=""),"country")] %>% dplyr::rename("call"=paste("X",invname,sep=""))
callstab$country <- factor(callstab$country,levels = cntorder,ordered=T)
callsorder <- callstab$sample[order(callstab$country,callstab$call)]
callstab$sample <- factor(callstab$sample,levels=callsorder,ordered=T)
callstab$cncode <- callstab$country
levels(callstab$cncode) <- substr(levels(callstab$country),0,3)

callplot <- ggplot(callstab,aes(x=1,y=as.numeric(sample),fill=as.factor(call))) + geom_tile() +
  ylab("samples") + xlab("SNPs")+ theme(legend.position="none")+
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) +
  scale_fill_manual(values=c("0"="blue2","1"="purple2","2"="red2")) +
  facet_grid("cncode ~ 'calls'",scale="free",space="free") +
  theme(panel.spacing = unit(0.2, "mm"),
        axis.text=element_blank(),
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        strip.text.y.right = element_text(angle = 0))
#callplot


design <- c("  1111112233445
               1111112277445
               1111116666445  ")

combplot <- distplot + pcaplots + aimplot + snpplot + callplot + plot_spacer() + plot_spacer() + plot_layout(design=design) + plot_annotation(tag_level=c("a"))
combplot

ggsave(paste(outprefix,".png",sep=""),combplot,width=12,height=7,dpi=400)
ggsave(paste(outprefix,".pdf",sep=""),combplot,width=12,height=7,dpi=400)
