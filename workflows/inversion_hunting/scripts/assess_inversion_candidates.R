
library("tidyverse")
library("lostruct")

#library("patchwork")
library("gridExtra")
library("grid")

library("getopt")

#######
# default inversion validation criteria
#######
#MAXD <- 0.2
#MAXWSS <- 5
#MINBSS <- 20
#MINBSP <- 0.95
blocksize <- 5e5



#######
# functions
#######

vcf_positions <- function (file, regions) {
  bcf.sites <- data.table::fread(cmd = paste("bcftools query -f '%CHROM\\t%POS\\n'",
                                             file,"-r",region_string(regions)),
                                 header = FALSE, sep = "\t", data.table = FALSE)
  colnames(bcf.sites) <- c("chrom", "pos")
  bcf.sites
}

vcf_genotypes <- function (file, regions, samples) {
  txtgenos <- data.table::fread(cmd = paste("bcftools query -f '[ %GT]\\n'",
                                            "-r",region_string(regions),
                                            "-s",shQuote(paste(samples, collapse = ",")),
                                            file),
                                header = FALSE, sep = "\t", data.table = FALSE)
  txtgenos
}

meanF3 <- function(invsnps,p3,p1,p2) {
  #calculate f3 stat for 3 'pops'
  getmaf = function(x) {sum(na.omit(x)) / sum(!is.na(x))*2}
  getmafs <- function(inds,snps) {
    if(sum(inds)==1) {  #if only one ind
      if (is.null(nrow(invsnps))) {   #if only one SNP
        maf <- snps[inds]/2
      } else {
        maf <- snps[,inds]/2
      }
    } else {
      if (is.null(nrow(invsnps))) {   #if only one SNP
        maf <- getmaf(snps[inds])
      } else {
      maf <- apply(snps[,inds],1,getmaf)
      }
    }
    maf
  }

  p1maf <- getmafs(p1,invsnps)
  p2maf <- getmafs(p2,invsnps)
  p3maf <- getmafs(p3,invsnps)

  f3=(p3maf-p1maf)*(p3maf-p2maf)
  mean(na.omit(f3))
}

jackknifeF3se <- function(snps,posns,p3,p1,p2,blocksize=1e5) {
  starts <- seq(floor(min(posns$pos)/blocksize)*blocksize, floor(max(posns$pos)/blocksize)*blocksize, by=blocksize)
  snpcounts <- c()
  for(s in starts) {
    snpcounts = c(snpcounts,sum(posns$pos>=s & posns$pos < s+blocksize))}
  starts <- starts[snpcounts>0]
  snpcounts <- snpcounts[snpcounts>0]


  f3blocks <- c()
  for(s in starts) {
    #write(paste("  getting F3 for block",s,snpcounts[starts==s],"SNPs"),stderr())
    blocksnps <- invsnps[posns$pos>=s & posns$pos<s+blocksize,]
    blockf3 <- meanF3(blocksnps,p3,p2,p1)
    f3blocks = c(f3blocks,blockf3)}

  f3jacks <- c()
  for(s in starts) {
    f3jacks = c(f3jacks,mean(f3blocks[starts != s]))}

  # compute mean of jackknife values
  #m <- weighted.mean(f3jacks,snpcounts)
  m <- mean(f3jacks)
  n <- length(starts)
  # compute standard error
  sv = ((n - 1) / n) * sum((f3jacks - m)^2)
  se = sqrt(sv)
  se
}



assessInvK <- function(pcs,invsnps,maxd,maxwss,minbss,minbsp) {

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
       "prop_wss"=tot.wss/tot.ss,
       "prop_bss"=tot.bss/tot.ss,
       "clusters"=clusters)
}

rotateXY <- function(coords,a,origin=c(0,0)) {
  arad <- a*(pi/180)
  rotm <- matrix(c(cos(arad),sin(arad),-sin(arad),cos(arad)),ncol=2)
  t(rotm %*% t(coords))
}








opttab <- matrix(c("infile","i","1","character",
                   "vcf","v","1","character",
                   "country","c","1","character",
                   "samples","s","1","character",
                   "outfile","o","1","character",
                   "blocksize","b","2","integer",
                   "maxd","D","2","double",
                   "maxwss","W","2","double",
                   "minbss","B","2","double",
                   "minbsp","P","2","double"
                      ),byrow=T,ncol=4)
opt <- getopt(opttab)

invcandfile <- opt$infile
outfile <- opt$outfile
#chrom <- opt$chr
vcffile <- opt$vcf
country <- opt$country
sppfile <- opt$samples


if (!is.null(opt$blocksize)  ) {blocksize <- opt$blocksize}
if (!is.null(opt$maxd)  ) {MAXD <- opt$maxd}
if (!is.null(opt$maxwss) ) {MAXWSS <- opt$maxwss}
if (!is.null(opt$minbss)) {MINBSS <- opt$minbss}
if (!is.null(opt$minbsp)) {MINBSP <- opt$minbsp}

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

#invcands <- invcands[invcands$cluster %in% c(485),]


# Run PCAs for all candidate inversion regions
for(C in unique(invcands$cluster)) {
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

  if(exists("pcs")) {
    pcs <- rbind(pcs,invpcs)
  } else {
    pcs <- invpcs
  }
}



#### assess PCA clusters via kmeans with angle rotation

write("assessing PCAs as inversions",file=stderr())
write(paste("D<=",MAXD,"WSS<=",MAXWSS,"BSS>=",MINBSS,"BSP>=",MINBSP),file=stderr())
invsummary <- invcands %>%
                    group_by(cluster) %>%
                    mutate(size = length(block)*blocksize)  %>%
                    select(c("cluster","au","bp","meandist","lowdist","size")) %>%
                    unique()

angles <- c(0,rep(seq(5,45,5),each=2)*c(1,-1))
if(exists("pcs")) {
  pcs$inv <- factor(pcs$inv,levels=invsummary$cluster,ordered=T)
  pcs$valid<-factor(NA,levels=c("aa","ab","bb"))

  for(I in unique(pcs$inv)) {
    pcsinv <- subset(pcs,inv==I)

    anygood = FALSE
    for(a in angles) {
      rotPCs <- rotateXY(as.matrix(pcsinv[,c("PC1","PC2")]),a)
      colnames(rotPCs) <- c("PC1","PC2")
      assk <- assessInvK(rotPCs[,"PC1"],
                        maxd=MAXD,
                        maxwss=MAXWSS,
                        minbss=MINBSS,
                        minbsp=MINBSP)

      if(assk$valid) {
        anygood=T
        break
      }
    } #all angles tested


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
    invsummary[invsummary$cluster==I,"prop_bss"]  <- round(assk$prop_bss,3)
    # invsummary[invsummary$cluster==I,"f3"]  <- round(assk$f3,3)

    invsummary[invsummary$cluster==I,"angle"]    <- a


  }
}

########
# do f3 test on candidates
########
invsummary$admixed=NA
invsummary$f3=NA
invsummary$f3se=NA


for(I in unique(invsummary$cluster[invsummary$valid])) {
  invsnps <- vcf_query(vcffile,
                       samples=samples,
                       regions=invcands[invcands$cluster==I,c("chromname","start","end")])

  invposns <- vcf_positions(vcffile,
                            invcands[invcands$cluster==I,c("chromname","start","end")])

  goodloci <- apply(invsnps,1,function(x) {!any(is.na(x))})

  invsnps <- invsnps[goodloci,]
  invposns <- invposns[goodloci,]
  clusters <- pcs$valid[which(pcs$inv==I)]
  p1 <- clusters=='aa'
  p3 <- clusters=='ab'
  p2 <- clusters=='bb'

  f3 <- meanF3(invsnps,p3,p1,p2)
  f3se <- jackknifeF3se(invsnps,invposns,p3,p1,p2)
  invsummary$f3[invsummary$cluster==I] <- f3
  invsummary$f3se[invsummary$cluster==I] <- f3se
  if(f3 < 0-(2*f3se)) {
    invsummary$admixed[invsummary$cluster==I] <- T
  } else {
    invsummary$admixed[invsummary$cluster==I] <- F
  }

}







write(paste("writing",nrow(invsummary),"candidates"),file=stderr())
if(nrow(invsummary) > 0) {
    write.table(invsummary,file=paste(outfile,"txt",sep="."),sep="\t",quote=F,col.names=T,row.names=F)
} else {
    file.create(paste(outfile,"txt",sep="."))
}


write(paste("plotting",nrow(invsummary),"PCs"),file=stderr())

invsummary$inv <- factor(invsummary$cluster,levels=sort(unique(invsummary$cluster)),ordered=T)
invcandsV <- merge(invcands,invsummary[,c("cluster","valid","admixed")],all.x=T)

invcandsV$validation <- "fail"
invcandsV$validation[invcandsV$valid] <- "pass-d"
invcandsV$validation[invcandsV$admixed] <- "pass-F3"

clustplot <- ggplot(invcandsV,aes(x=pos,fill=validation,y=as.factor(cluster))) + geom_tile() + ylab("cluster")

outcalls <- paste(outfile,"calls.txt",sep="_")

if(exists("pcs") & !all(is.na(pcs$valid))) {
  saveRDS(pcs,file=paste(outfile,"pcs.Rds",sep="_"))

  callstable <- pivot_wider(pcs[,c("sample","valid","inv")],names_from = inv,values_from = c("valid"))
  write.table(callstable,file=outcalls,sep="\t",quote=F,row.names = F)

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
  file.create(outcalls)
  file.create(paste(outfile,"pcs.Rds",sep="_"))
  png(filename = paste(outfile,"png",sep="."),width=350,height=200,units="mm",res=400)
  grid.arrange(clustplot, grid.rect(gp=gpar(col="white")), ncol=2)
  dev.off()
}
