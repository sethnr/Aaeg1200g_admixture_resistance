
library("tidyverse")
library("lostruct")

#library("patchwork")
library("gridExtra")
library("grid")

library("getopt")

opttab <- matrix(c("infile","i","1","character",
                   "vcf","v","1","character",
                   "country","c","1","character",
                   "samples","s","1","character",
                   "outfile","o","1","character",
                   "blocksize","b","2","integer",
                   "maxd","D","2","double",
                   "maxss","S","2","double",
                   "minbss","B","2","double"
                      ),byrow=T,ncol=4)
opt <- getopt(opttab)

invcandfile <- opt$infile
outfile <- opt$outfile
#chrom <- opt$chr
vcffile <- opt$vcf
country <- opt$country
sppfile <- opt$samples

#default inversion validation criteria
MAXD <- 0.25
MAXSS <- 20
MINBSS <- 0.95
blocksize <- 5e5
if (!is.null(opt$blocksize)  ) {blocksize <- opt$blocksize}
if (!is.null(opt$maxd)  ) {MAXD <- opt$maxd}
if (!is.null(opt$maxss) ) {MAXSS <- opt$maxss}
if (!is.null(opt$minbss)) {MINBSS <- opt$minbss}
#invcandfile <- "uganda_chr1_pvclust_single_uncentered_hier_P0.99_n1000.txt"

chromname <- c("NC_035107.1","NC_035108.1","NC_035109.1")
chromlen <- c(310827022,474425716,409777670)
names(chromlen) <- chromname

spptab <- read.table(sppfile,header=T, sep="\t")
samples <- spptab$sample[spptab$country==country]
write(paste("found",length(samples),"samples for",country),stderr())


invcands <- read.table(invcandfile,header=T)
invcands$end <- as.numeric((as.data.frame(strsplit(invcands$block,":"))[2,]))
invcands$start <- invcands$end-blocksize
invcands$chromname <- chromname[invcands$chrom]

for(C in unique(invcands$cluster)) {
  #write.table(invcands[invcands$cluster==C,c("chromname","start","end")],stderr())
  #write(samples,stderr())
  invsnps <- vcf_query(vcffile,
                        samples=samples,
                        regions=invcands[invcands$cluster==C,c("chromname","start","end")])
  goodsnps <- invsnps[apply(invsnps,1,function(x) {!any(is.na(x))}),]
  write(paste("found",nrow(goodsnps),"snps for",country,"cluster",C),stderr())

  pca <- prcomp(t(goodsnps))

  invpcs <- as.data.frame(pca$x[,c("PC1","PC2")])
  invpcs$sample <- samples
  invpcs <- merge(invpcs,spptab)
  invpcs$inv <- C
  #invpcs$chrom <- chr

  if(exists("pcs")) {
    pcs <- rbind(pcs,invpcs)
  } else {
    pcs <- invpcs
  }
}

write("assessing PCAs as inversions",file=stderr())

invsummary <- invcands %>% 
                    group_by(cluster) %>% 
                    mutate(size = length(block)*blocksize)  %>% 
                    select(c("cluster","au","bp","meandist","lowdist","size")) %>%
                    unique()

#### assess PCA clusters
if(exists("pcs")) {
  pcs$inv <- factor(pcs$inv,levels=invsummary$cluster,ordered=T)
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

    invpass <- ((abs(delta)<=MAXD) && ((tot.wss/n) <= MAXSS) && (bss/tot.ss)>=MINBSS)
    write(paste("inv:",I,
                #"\tctr:",paste(round(centers,1),collapse="/"),
                #"\twss:",paste(round(wss,1),collapse="/"),
                "\tbss/tot:",round(bss/tot.ss,2),
                "\tss",round(tot.wss/n,2),
                "d",round(abs(delta),2),
                "\tpass:",invpass),file=stderr())

    invsummary[invsummary$cluster==I,"valid"] <- invpass
    invsummary[invsummary$cluster==I,"d"] <- round(abs(delta),2)
    invsummary[invsummary$cluster==I,"mean_wss"] <- round(tot.wss/n,2)
    invsummary[invsummary$cluster==I,"bss_tot"] <- round(bss/tot.ss/n,2)


    if(invpass) {
      pcs$valid[which(pcs$inv==I)] <- clusters
    }
  }
}


write(paste("writing",nrow(invsummary),"candidates"),file=stderr())
write.table(invsummary,file=paste(outfile,"txt",sep="."),sep="\t",quote=F,col.names=T,row.names=F)

write(paste("plotting",nrow(invsummary),"PCs"),file=stderr())

clustplot <- ggplot(invcands,aes(x=pos,fill=meandist,y=as.factor(cluster))) + geom_tile() + ylab("cluster")

if(exists("pcs")) {
  saveRDS(pcs,file=paste(outfile,"pcs.Rds",sep="_"))

  invcols <- scale_color_manual(values=c("aa"="yellow","ab"="orange","bb"="red"),na.value = "dark grey")
  ncol=round(sqrt(length(unique(pcs$inv))))

  invpca <- ggplot(pcs,aes(x=PC1,y=PC2,color=valid)) + geom_point() + coord_fixed() + invcols +
    facet_wrap("inv ~ .",ncol=ncol) + theme(legend.position = "none")

  png(filename = paste(outfile,"png",sep="."),width=350,height=200,units="mm",res=400)
  grid.arrange(clustplot, invpca, ncol=2)
  invpca
  dev.off()
} else {
  file.create(paste(outfile,"pcs.Rds",sep="_"))
  #file.create(paste(outfile,"txt",sep="."))
  png(filename = paste(outfile,"png",sep="."),width=350,height=200,units="mm",res=400)
  grid.arrange(clustplot, grid.rect(gp=gpar(col="white")), ncol=2)
  dev.off()
}
