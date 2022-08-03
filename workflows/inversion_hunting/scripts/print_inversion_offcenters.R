library("tidyverse")
library("patchwork")
library("lostruct")
library("gridExtra")
library("grid")

library("getopt")

######
# functions
######

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
  maf = function(x) {sum(na.omit(x)) / sum(!is.na(x))*2}
  getmafs <- function(inds,snps) {
    if(sum(inds)>1) {
      maf <- apply(snps[,inds],1,maf)
    } else {
      maf <- snps[,inds]/2
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
    # jacksnps <- invsnps[posns$pos<s | posns$pos>=s+blocksize,]
    # jackf3 <- meanF3(jacksnps,p3,p2,p1)
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



assessInvK <- function(pcs,invsnps,maxd=MAXD,maxwss=MAXWSS,minbss=MINBSS,minbsp=MINBSP) {
  
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
    #(tot.bss/nbss   >= minbss &&
    # tot.wss/nwss   <= maxwss )
    # &&
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







#default inversion validation criteria
MAXD <- 0.25
MAXWSS <- 2
MINBSS <- 10
MINBSP <- 0.95
blocksize <- 5e5


opttab <- matrix(c("infile","i","1","character",
                   "vcf","v","1","character",
                   "country","c","1","character",
                   "samples","s","1","character",
                   "outfile","o","1","character",
                   "inversion","n","1","integer",
                   "maxd","D","2","double",
                   "maxss","S","2","double",
                   "minbss","B","2","double"
),byrow=T,ncol=4)
opt <- getopt(opttab)

invcandfile <- opt$infile
outfile <- opt$outfile
vcffile <- opt$vcf
I <- as.integer(opt$inversion)
country <- opt$country
sppfile <- opt$samples

if (!is.null(opt$blocksize)  ) {blocksize <- opt$blocksize}
if (!is.null(opt$maxd)  ) {MAXD <- opt$maxd}
if (!is.null(opt$maxss) ) {MAXSS <- opt$maxss}
if (!is.null(opt$minbss)) {MINBSS <- opt$minbss}



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
invcands <- invcands[invcands$cluster == I,]

write(paste("found",nrow(invcands),"blocks for",I),stderr())


invsnps <- vcf_query(vcffile,
                      samples=samples,
                      regions=invcands[invcands$cluster==I,c("chromname","start","end")])

invposns <- vcf_positions(vcffile,
                          invcands[invcands$cluster==I,c("chromname","start","end")])

goodloci <- apply(invsnps,1,function(x) {!any(is.na(x))})

invsnps <- invsnps[goodloci,]
invposns <- invposns[goodloci,]

write(paste("found",nrow(invsnps),"snps for",country,"cluster",I),stderr())

pca <- prcomp(t(invsnps))

invpcs <- as.data.frame(pca$x[,c("PC1","PC2")])
invpcs$sample <- samples
invpcs <- merge(invpcs,spptab)
invpcs$inv <- I
#invpcs$chrom <- chr

if(exists("pcs")) {
  pcs <- rbind(pcs,invpcs)
} else {
  pcs <- invpcs
}


write(paste("assessing inversion",I,"in all rotations"),stderr())
pcsinv <- subset(pcs,inv==I)


rm("allpcs")
rotsummary <- data.frame("I"=numeric(),
                         "angle"=numeric(),
                         "mean_wss"=numeric(),
                         "mean_bss"=numeric(),
                         "prop_wss"=numeric(),
                         "prop_bss"=numeric(),
                         "d"=numeric(),
                         "valid"=logical())

angles <- c(0,rep(seq(5,45,5),each=2)*c(1,-1))
for(a in angles) {
  rotPCs <- as.data.frame(rotateXY(as.matrix(pcsinv[,c("PC1","PC2")]),a))
  colnames(rotPCs) <- c("PC1","PC2")
  assk <- assessInvK(rotPCs[,"PC1"],invsnps)
  write(paste(I,
              a,
              round(assk$mean_wss,2),
              round(assk$mean_bss,2),
              round(assk$prop_bss,2),
              round(assk$d,2),
              assk$valid),
        stderr())
  
  p1 <- assk$clusters=='aa'
  p3 <- assk$clusters=='ab'
  p2 <- assk$clusters=='bb'
  f3 <- meanF3(invsnps,p3,p1,p2)
  assk$f3 <- f3
  if(assk$valid) {
    f3se <- jackknifeF3se(invsnps,invposns,p3,p1,p2)
    assk$f3se <- f3se
    #if fails F3 test, set valid to false
    if(f3 > 0-(2*f3se)) {
      assk$valid <- F
    }
  }
  
  i<-nrow(rotsummary)+1
  rotsummary[i,c("I","angle","mean_wss","mean_bss","prop_wss","prop_bss","d","f3")] <- c(I,a,
                                                                              assk$mean_wss,assk$mean_bss,
                                                                              assk$prop_wss,assk$prop_bss,
                                                                              assk$d,assk$f3)
  rotsummary[i,"valid"] <- assk$valid
  
  rotPCs$cluster <- assk$clusters
  rotPCs$angle <- a
  if(exists("allpcs")) {
    allpcs <- rbind(allpcs,rotPCs)
  } else {
    allpcs <- rotPCs
  }
}
rotsummary <- pivot_longer(rotsummary,cols =c("mean_wss","mean_bss","prop_wss","prop_bss","d","f3"))
rotsummary$value <- as.numeric(rotsummary$value)
rotsummary$angle <- as.numeric(rotsummary$angle)




################
# do all plots
################

ncol=ceiling(sqrt(length(angles)))

invcols <- scale_color_manual(values=c("aa"="yellow","ab"="orange","bb"="red"),na.value = "dark grey")

allrotplots <- ggplot(allpcs,aes(x=PC1,y=PC2,color=cluster)) + 
  geom_abline(aes(slope=angle/45,intercept=0),linetype=2,allpcs)+
  geom_point() + coord_fixed() + invcols +
  facet_wrap("angle ~ .",ncol=ncol) + 
  theme(legend.position = "none") +
  ggtitle(paste("all rotations, inv",I))

varlimits <- data.frame(intercept=c(MINBSS,MAXWSS,MAXD,MINBSP),name=c("mean_bss","mean_wss","d","prop_bss"))

dpropplots <- ggplot(subset(rotsummary,name %in% c("d","prop_bss","prop_wss")),aes(x=angle,y=value,color=name,group=name)) + 
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

f3plots <- ggplot(subset(rotsummary,name %in% c("f3")),aes(x=angle,y=value,color=name,group=name)) + 
  geom_line() + 
  geom_hline(yintercept=0,linetype="dashed")+
  ggtitle(paste("f3",I))

allrotcombplots <- allrotplots | (dpropplots / valplots / f3plots)

allrotcombplots

ggsave(outfile,
       allrotcombplots, width=10,height=6,dpi=400)

