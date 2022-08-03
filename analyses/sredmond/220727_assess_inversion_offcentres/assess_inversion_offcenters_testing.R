library("tidyverse")
library("patchwork")
library("gridExtra")
library("grid")

library("getopt")

invcandfile <- "regions_chr1_Senegal.txt"
vcffile <- "data/aegy.wgs.aaa.aaf.norep5-30x.ac10.thinrand01.chr1.vcf.gz"
country <- "Senegal"
sppfile <- "resources/meta_Aaeg1kg_spp.txt"

# invcandfile <- "regions_chr1_Kenya.txt"
# vcffile <- "data/aegy.wgs.aaa.aaf.norep5-30x.ac10.thinrand01.chr1.vcf.gz"
# country <- "Kenya"
# sppfile <- "resources/meta_Aaeg1kg_spp.txt"


# invcandfile <- "regions_chr1_Trinidad.txt"
# vcffile <- "data/aegy.wgs.aaa.aaf.norep5-30x.ac10.thinrand01.chr1.vcf.gz"
# country <- "Trinidad"
# sppfile <- "resources/meta_Aaeg1kg_spp.txt"



#default inversion validation criteria
MAXD <- 0.25
MAXWSS <- 2
MINBSS <- 10
MINBSP <- 0.95
blocksize <- 5e5

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


#invcands <- invcands[invcands$cluster %in% c(194,181,50,217,271),]
invcands <- invcands[invcands$cluster %in% c(17),]

rm("pcs")
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


invsummary$inv <- factor(invsummary$cluster)

invcandsV <- merge(invcands,invsummary[,c("cluster","valid")],all.x=T)
clustplot <- ggplot(invcandsV,aes(x=pos,fill=valid,y=as.factor(cluster))) + geom_tile() + ylab("cluster")

invcols <- scale_color_manual(values=c("aa"="yellow","ab"="orange","bb"="red"),na.value = "dark grey")

ncol=round(sqrt(length(unique(pcs$inv))))

invpca <- ggplot(pcs,aes(x=PC1,y=PC2,color=valid)) + 
            geom_abline(aes(slope=(angle/45)*-1,intercept=0),linetype=2,invsummary)+
            geom_point() + coord_fixed() + invcols +
            facet_wrap("inv ~ .",ncol=ncol) + 
            theme(legend.position = "none")
#invpca
#grid.arrange(clustplot, invpca, ncol=2)
clustplot | invpca

ggsave(paste("inversions",country,MINBSS,MAXWSS,MAXD,MINBSP,"tests.png",sep="_"),
       width=10,height=6,dpi=400)




#for(I in unique(pcs$inv)) {
for(I in c(50,194,172,29,217,120,161)) {
  pcsinv <- subset(pcs,inv==I)


  rm("allpcs")
  rotsummary <- data.frame("I"=numeric(),
                           "angle"=numeric(),
                           "mean_wss"=numeric(),
                           "mean_bss"=numeric(),
                           "prop_bss"=numeric(),
                           "d"=numeric(),
                           "valid"=logical())
  
  for(a in angles) {
    rotPCs <- as.data.frame(rotateXY(as.matrix(pcsinv[,c("PC1","PC2")]),a))
    colnames(rotPCs) <- c("PC1","PC2")
    assk <- assessInvK(rotPCs[,"PC1"])
    write(paste(I,
                a,
                round(assk$mean_wss,2),
                round(assk$mean_bss,2),
                round(assk$prop_bss,2),
                round(assk$d,2),
                assk$valid),
          stderr())
    i<-nrow(rotsummary)+1
    rotsummary[i,c("I","angle","mean_wss","mean_bss","prop_bss","d")] <- c(I,a,assk$mean_wss,assk$mean_bss,assk$prop_bss,assk$d)
    rotsummary[i,"valid"] <- assk$valid
    rotPCs$cluster <- assk$clusters
    rotPCs$angle <- a
    if(exists("allpcs")) {
      allpcs <- rbind(allpcs,rotPCs)
    } else {
      allpcs <- rotPCs
    }
  }
  rotsummary <- pivot_longer(rotsummary,cols =c("mean_wss","mean_bss","prop_bss","d"))
  rotsummary$value <- as.numeric(rotsummary$value)
  rotsummary$angle <- as.numeric(rotsummary$angle)
  ncol=ceiling(sqrt(length(angles)))
  
  
  allrotplots <- ggplot(allpcs,aes(x=PC1,y=PC2,color=cluster)) + 
    geom_abline(aes(slope=angle/45,intercept=0),linetype=2,allpcs)+
    geom_point() + coord_fixed() + invcols +
    facet_wrap("angle ~ .",ncol=ncol) + 
    theme(legend.position = "none") +
    ggtitle(paste("all rotations, inv",I))
  
  varlimits <- data.frame(intercept=c(MINBSS,MAXWSS,MAXD,MINBSP),name=c("mean_bss","mean_wss","d","prop_bss"))
  
  dpropplots <- ggplot(subset(rotsummary,name %in% c("d","prop_bss")),aes(x=angle,y=value,color=name,group=name)) + 
    geom_line() + 
    geom_point(data=subset(rotsummary,name %in% c("d","prop_bss") & valid)) +
    geom_hline(data=subset(varlimits,name %in% c("d","prop_bss")),
               aes(yintercept=intercept,color=name),linetype="dashed")+
    ggtitle(paste("D / propbss",I))
  
  valplots <- ggplot(subset(rotsummary,name %in% c("mean_wss","mean_bss")),aes(x=angle,y=value,color=name,group=name)) + 
    geom_line() + 
    geom_hline(data=subset(varlimits,name %in% c("mean_wss","mean_bss")),
               aes(yintercept=intercept,color=name),linetype="dashed")+
    ggtitle(paste("mean wss , bss",I))
  
  allrotcombplots <- allrotplots | (dpropplots / valplots)
  ggsave(paste("all_rotations_",country,"_",I,".png",sep=""),
         allrotcombplots, width=10,height=6,dpi=400)
}

