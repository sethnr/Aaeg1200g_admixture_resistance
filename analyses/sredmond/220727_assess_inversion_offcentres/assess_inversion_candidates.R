
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
MAXWSS <- 2
MINBSS <- 10
MINBSP <- 0.95
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
#good: 194, 181, 129
#bad: 30, 50 
invsummary <- invcands %>%
                    group_by(cluster) %>%
                    mutate(size = length(block)*blocksize)  %>%
                    select(c("cluster","au","bp","meandist","lowdist","size")) %>%
                    unique()

assessInvK <- function(pcs,maxd=MAXD,maxwss=MAXWSS,minbss=MINBSS,minbsp=MINBSP) {
  
  kmpca <- kmeans(pcs,centers=3,nstart=50,iter.max=100)
  clusters <- factor(kmpca$cluster)
  clustorder <- order(kmpca$centers)
  
  #calculate deviation of middle cluster from center point
  centers <- kmpca$centers[clustorder]
  delta <- ((centers[2] - centers[1])-(centers[3] - centers[2])) /
    ((centers[3]-centers[1])/2)
  
  #rename levels 0/1/2 based on pcorder
  levels(clusters)[clustorder] <- c('aa','ab','bb')
  
  #get sum squares
  n<-length(clusters)
  npairs <- n*(n-1)/2 
  wss <- kmpca$withinss[clustorder]
  tot.wss <- kmpca$tot.withinss
  nwss <- sum(apply(table(clusters),1,FUN=function(x) {x*(x-1)/2}))
  
  tot.bss <- kmpca$betweenss
  ncl <- table(clusters)
  nbss <- sum(ncl[1]*ncl[2] + ncl[2]*ncl[3] + ncl[1]*ncl[3])
  
  tot.ss <- kmpca$totss
  c(nwss,nbss,npairs)
  
  #assess pass
  invpass <- (
    (tot.bss/nbss   >= minbss && 
       tot.wss/nwss   <= maxwss )
    && 
      (abs(delta)     <= maxd && 
         tot.bss/tot.ss >= minbsp ) )
  
  list("valid"=invpass,
       "d"=abs(delta),
       "mean_ss"=tot.ss/npairs,
       "mean_wss"=tot.wss/nwss,
       "mean_bss"=tot.bss/nbss,
       "prop_bss"=tot.bss/tot.ss,
       "clusters"=clusters)
}


rotateXY <- function(coords,a,origin=c(0,0)) {
  arad <- a*(pi/180)
  rotm <- matrix(c(cos(arad),sin(arad),-sin(arad),cos(arad)),ncol=2)
  t(rotm %*% t(coords))
}


angles <- c(0,rep(seq(5,45,5),each=2)*c(1,-1))
#angles <- c(0,rep(seq(5,45,5)))

#### assess PCA clusters via kmeans with angle rotate
if(exists("pcs")) {
  pcs$inv <- factor(pcs$inv,levels=invsummary$cluster,ordered=T)
  pcs$valid<-factor(NA,levels=c("aa","ab","bb"))
  
  for(I in unique(pcs$inv)) {
    pcsinv <- subset(pcs,inv==I)
    
    anygood = FALSE
    for(a in angles) {
      rotPCs <- rotateXY(as.matrix(pcsinv[,c("PC1","PC2")]),a)
      colnames(rotPCs) <- c("PC1","PC2")
      assk <- assessInvK(rotPCs[,"PC1"])
      write(paste(I,
                  a,
                  round(assk$d,2),
                  round(assk$mean_wss,2),
                  round(assk$mean_bss,2),
                  assk$valid),
            stderr())
      if(assk$valid) {
        anygood=T
        break}
    }
    if(anygood) {
      pcs$valid[which(pcs$inv==I)] <- assk$clusters
    } else {
      assk <- assessInvK(pcsinv[,c("PC1","PC2")])
      a <- 0
    }
    
    invsummary[invsummary$cluster==I,"valid"]    <- assk$valid
    invsummary[invsummary$cluster==I,"d"]        <- round(assk$d,2)
    invsummary[invsummary$cluster==I,"mean_wss"] <- round(assk$mean_wss,2)
    invsummary[invsummary$cluster==I,"mean_bss"]  <- round(assk$mean_bss,3)
    invsummary[invsummary$cluster==I,"angle"]    <- a
    
    
  }
}


write(paste("writing",nrow(invsummary),"candidates"),file=stderr())
if(nrow(invsummary) > 0) {
    write.table(invsummary,file=paste(outfile,"txt",sep="."),sep="\t",quote=F,col.names=T,row.names=F)
} else {
    file.create(paste(outfile,"txt",sep="."))
}


write(paste("plotting",nrow(invsummary),"PCs"),file=stderr())

invsummary$inv <- factor(invsummary$cluster)
invcandsV <- merge(invcands,invsummary[,c("cluster","valid")],all.x=T)
clustplot <- ggplot(invcandsV,aes(x=pos,fill=valid,y=as.factor(cluster))) + geom_tile() + ylab("cluster")

if(exists("pcs")) {
  saveRDS(pcs,file=paste(outfile,"pcs.Rds",sep="_"))

  invcols <- scale_color_manual(values=c("aa"="yellow","ab"="orange","bb"="red"),na.value = "dark grey")
  ncol=round(sqrt(length(unique(pcs$inv))))
  invpca <- ggplot(pcs,aes(x=PC1,y=PC2,color=valid)) + 
    geom_abline(aes(slope=(angle/45)*-1,intercept=0),linetype=2,invsummary)+
    geom_point() + coord_fixed() + invcols +
    facet_wrap("inv ~ .",ncol=ncol) + 
    theme(legend.position = "none")
  
  png(filename = paste(outfile,"png",sep="."),width=350,height=200,units="mm",res=400)
  grid.arrange(clustplot, invpca, ncol=2)
  invpca
  dev.off()
} else {
  file.create(paste(outfile,"pcs.Rds",sep="_"))
  png(filename = paste(outfile,"png",sep="."),width=350,height=200,units="mm",res=400)
  grid.arrange(clustplot, grid.rect(gp=gpar(col="white")), ncol=2)
  dev.off()
}
