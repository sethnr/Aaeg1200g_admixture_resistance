library("tidyverse")

library("patchwork")
library("gridExtra")

library("zoo")
library("plyr")



indir <- "inv_regions"

#make function to read/parse all vars and create table
parse_inv <- function(x) {
  fnspl <- strsplit(basename(x),"_")[[1]]
  chr = as.integer(gsub("chr","",fnspl[2]))
  pop <- gsub(".txt","",fnspl[3])
  read.table(x,header=T) %>% 
    add_column("pop"=pop,"chr"=chr,.before = 1)
}

#read in block files:
write(paste("getting blocks from ",indir),stderr())
invfiles <- list.files(indir,pattern="inversions.*txt",full.names = T)
invtab <- ldply(invfiles, parse_inv)


invtab <- invtab %>% rename(c("cluster"="inversion","pop"="country"))

mantab <- read.table("200810_visual_inspect_good_inversions.tsv",header=T)

invtab <- merge(invtab,mantab,all=T)
invtab$quality[is.na(invtab$quality)] <- 0

#invtab$quality <- factor(invtab$quality)

ggplot(invtab,aes(x=d,y=mean_wss,color=as.factor(quality))) + geom_point(alpha=0.5) + facet_wrap("quality")
ggplot(invtab,aes(x=d,y=mean_wss,color=as.factor(quality))) + geom_point(alpha=0.5) + facet_wrap("quality") + xlim(0,1) + ylim(0,200)

ggplot(invtab,aes(x=mean_bss,y=d,color=as.factor(quality))) + geom_point(alpha=0.5) + facet_wrap("quality")
ggplot(invtab,aes(x=mean_bss,y=d,color=as.factor(quality))) + geom_point(alpha=0.5) + facet_wrap("quality") + xlim(0,200)

ggplot(invtab,aes(x=mean_wss,y=mean_bss,color=as.factor(quality))) + geom_point(alpha=0.5) + facet_wrap("quality")
ggplot(invtab,aes(x=mean_wss,y=mean_bss,color=as.factor(quality))) + geom_point(alpha=0.5) + facet_wrap("quality") + xlim(0,100)

ggplot(invtab,aes(x=d,y=prop_bss,color=as.factor(quality))) + geom_point(alpha=0.5) + facet_wrap("quality")

#bss proportion lower limit
bsplim <- 0.95
#mean wss upper limit
wsslim <- 10 # (10?)
#d upper limit
dlim <- 0.2

invtab$hardlim <- invtab$d < dlim & invtab$mean_wss < wsslim & invtab$prop_bss > bsplim


ggplot(invtab,aes(x=mean_wss,y=mean_bss,color=hardlim)) + geom_point(alpha=0.5) + facet_wrap("quality")

FP = sum(invtab$hardlim & invtab$quality==0)
TP = sum(invtab$hardlim & invtab$quality!=0)
FN = sum(!invtab$hardlim & invtab$quality!=0)
TN = sum(!invtab$hardlim & invtab$quality==0)

sens <- TP/(TP+FN)
spec <- TN/(TN+FP)
c("sens"=sens,"spec"=spec)


speclims <- function(dlim,bsplim,wsslim,invtab) {
  hardlim <- invtab$d < dlim & invtab$mean_wss < wsslim & invtab$prop_bss > bsplim
  
  
  FP = sum(hardlim & invtab$quality==0)
  TP = sum(hardlim & invtab$quality!=0)
  FN = sum(!hardlim & invtab$quality!=0)
  TN = sum(!hardlim & invtab$quality==0)
  sens <- TP/(TP+FN)
  spec <- TN/(TN+FP)
  
  FP2 = sum(hardlim & invtab$quality==0)
  TP2 = sum(hardlim & invtab$quality!=0)
  FN2 = sum(!hardlim & invtab$quality>=2)
  TN2 = sum(!hardlim & invtab$quality==0)
  
  sens2 <- TP/(TP+FN2)
  spec2 <- TN/(TN+FP)
  
  
  list("d"=dlim,"bsp"=bsplim,"wss"=wsslim,
    "sens1"=sens,"spec1"=spec,
    "sens2"=sens2
    )
}


ds <- seq(0.05,0.5,by=0.05)
wsss <- seq(0,50,by=5)
bsps <- seq(0.8,0.99,by=0.05)

roctab <- data.frame(d=numeric(),
                     bsp=numeric(),
                     wss=numeric(),
          sens1=numeric(),
           sens2=numeric(),
           spec=numeric())

for(d in ds) {
  for(wss in wsss) {
    for (bsp in bsps) {
      res = speclims(dlim=d,
               bsplim=bsp,
               wsslim=wss,
               invtab)
      i=nrow(roctab)+1
      roctab[i,] <- c(res$d,
                      res$bsp,
                      res$wss,
                      res$sens1,
                      res$sens2,
                      res$spec
      )
    }
  }
}


# ggplot(roctab,aes(x=spec,y=sens2,color=as.factor(bsp),group=as.factor(bsp))) + geom_line()
# ggplot(roctab,aes(x=spec,y=sens2,color=as.factor(wss),group=as.factor(wss))) + geom_line()


ggplot(subset(roctab,d==0.2),aes(x=wss,y=sens2,color=as.factor(bsp),group=bsp)) + 
  geom_line() + geom_line(aes(y=spec),linetype="dashed") + ylim(0,1)

ggplot(subset(roctab,bsp=="0.9"),aes(x=d,y=sens2,color=as.factor(wss),group=wss)) + 
  geom_line() + geom_line(aes(y=spec),linetype="dashed") + ylim(0,1)

ggplot(subset(roctab,wss==5),aes(x=bsp,y=sens2,color=as.factor(d),group=d)) + 
  geom_line() + geom_line(aes(y=spec),linetype="dashed") + ylim(0,1)



#invtab$hardlim <- invtab$d <= 0.2 & invtab$mean_wss < 5 & invtab$prop_bss > 0.9 &  invtab$mean_bss > 20
#best guess limits based on curves above
invtab$hardlim <- invtab$d <= 0.2 & invtab$mean_wss < 5 & invtab$mean_bss > 25
ggplot(invtab,aes(x=mean_wss,y=mean_bss,color=hardlim)) + geom_point(alpha=0.5) + facet_wrap("quality")

invtab[invtab$quality==3,c(1:3,10:19)]
invtab[invtab$quality==2,c(1:3,10:19)]
invtab[invtab$country=="Senegal" & invtab$chr==1,c(1:3,10:19)]

