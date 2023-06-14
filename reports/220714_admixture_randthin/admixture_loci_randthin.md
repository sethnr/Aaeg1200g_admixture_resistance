---
output: html_document
editor_options: 
  chunk_output_type: console
---

```r
library("tidyverse")

library("patchwork")
library("gridExtra")

library("zoo")
library("plyr")
library("getopt")

knitr::opts_chunk$set(fig.width=10,fig.height=8,dpi=300)
```


```r
indir <- "admix_tests"
outprefix <- "3pop_admix_loci"

#make function to read/parse all vars and create table
parse_3pop <- function(x) {
  fnspl <- strsplit(basename(x),"_")[[1]]
  chr = fnspl[2]
  mixpop <- fnspl[5]
  pop1 <- fnspl[3]
  pop2 <- fnspl[4]
  read.table(x,col.names = c("start","mid","end","f3")) %>% 
    add_column("p1"=pop1,"p2"=pop2,"p3"=mixpop,"chr"=chr)
}

#read in block files:
write(paste("getting blocks from ",indir),stderr())
pop3files <- list.files(indir,pattern="3pop.*_blocks.txt.gz",full.names = T)
pop3tab <- ldply(pop3files, parse_3pop)
pop3tab[pop3tab$p2>pop3tab$p1,c("p1","p2")] <- pop3tab[pop3tab$p2>pop3tab$p1,c("p2","p1")]
pop3tab$set <- paste(pop3tab$p3,pop3tab$p1,pop3tab$p2,sep="/")

write(paste("found",nrow(pop3tab),"block lines from",length(pop3files),"files"),stderr())


#read in summary (significance) files
pop3files <- list.files(indir,pattern="3pop.*_summaries.txt",full.names = T)
pop3sum <- ldply(pop3files, function(x) {read.table(x,col.names=c("p3","p1","p2","chr","f3","f3sd","f3z","f3sig"))})
pop3sum[pop3sum$p2>pop3sum$p1,c("p1","p2")] <- pop3sum[pop3sum$p2>pop3sum$p1,c("p2","p1")]
pop3sum$set <- paste(pop3sum$p3,pop3sum$p1,pop3sum$p2,sep="/")
pop3sum$parent <- paste(pop3sum$p1,pop3sum$p2,sep="/")
pop3sum$f3sig <- as.logical(pop3sum$f3sig)
write(paste("found",nrow(pop3sum),"summary lines from",length(pop3files),"files"),stderr())


#get sign admixed pops based on whole genome f3
sigpops <- unique(pop3sum$p3[pop3sum$f3sig])
#get those that show sig f3 but not from admixed parents
sigsets <- unique(pop3sum$set[pop3sum$f3sig & 
                                ((!pop3sum$p1 %in% sigpops) & (!pop3sum$p2 %in% sigpops) & (pop3sum$p3 %in% sigpops))])
#get those that show non-sig f3
nonsigsets <- subset(pop3sum,!p1 %in% sigpops & !p2 %in% sigpops &
                                !p3 %in% sigpops & !set %in% sigsets)$set



#calculate z-scores
pop3results <- merge(pop3tab,pop3sum[,c("set","parent","chr","f3sd","f3sig")],by=c("set","chr"),all=T)
pop3results$f3z <- pop3results$f3/pop3results$f3sd
```


```r
sigmixes <- subset(pop3sum,p3 %in% sigpops) %>% group_by(set) %>% 
                        dplyr::summarise(p1=first(p1),p2=first(p2),p3=first(p3),
                                         f3mean=mean(f3),f3zmean=mean(f3z),
                                         f3zmin=min(f3z),sig=max(f3sig))


minmixes <- subset(sigmixes,!p1 %in% sigpops & !p2 %in% sigpops) %>% 
                  group_by(p3) %>% 
                  slice_min(order_by = f3zmin)
```



```r
#######
# plot z scores
####### 

ggplot(subset(pop3results,set %in% sigsets),aes(x=mid,y=f3z,group=set,color=parent)) + 
  geom_line() + facet_grid(p3 ~ chr,scale="free_x",space="free_x") + 
  scale_y_continuous(limits=c(-100,0)) +
  theme(legend.position = "bottom")
```

```
## Warning: Removed 2044 row(s) containing missing values (geom_path).
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png)

```r
ggplot(subset(pop3results,!p1 %in% sigpops & !p2 %in% sigpops & !p3 %in% sigpops & !set %in% sigsets),aes(x=mid,y=f3z,group=set,color=parent)) + 
  geom_line() + 
  facet_grid(p3 ~ chr,scale="free_x",space="free_x") + 
  scale_y_continuous(limits=c(-100,0)) +
  theme(legend.position = "bottom")
```

```
## Warning: Removed 528624 row(s) containing missing values (geom_path).
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-2.png)

```r
ggplot(subset(pop3results,set %in% minmixes$set),aes(x=mid,y=f3z,group=set,color=parent)) + 
  geom_line() + facet_grid(p3 ~ chr,scale="free_x",space="free_x") + 
  scale_y_continuous(limits=c(-100,0)) +
  theme(legend.position = "bottom")
```

```
## Warning: Removed 55 row(s) containing missing values (geom_path).
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-3.png)


#maps of admixture populations(no arrows)

```r
library(maptools)
library(raster)
library(ggmap)

world <- map_data("world")
poptotals <- read.table("resources/aegy.wgs.pops.list.csv",sep=",",header=T,stringsAsFactors = F)
colnames(poptotals) <- tolower(colnames(poptotals))
poptotals$country <- gsub(" ","",poptotals$country)


countryfile <- "country_geocoords.txt"
if(file.exists(countryfile)){
  locations <- read.table(countryfile,header=T) 
} else {
  register_google(key="AIzaSyBFbT5X6MjrepPAlzIj3zz1HS58hnIPtEM", write=T)
  countries <- unique(poptotals$country)
  locations <- geocode(countries)
  locations$country <- countries
  write.table(locations,countryfile,sep="\t",row.names=F,col.names=T)
}



countrytotals <- unique(merge(aggregate(num_bams ~ country,poptotals,FUN=sum),
                              aggregate(f3sig ~ p3,pop3sum,FUN=max),by.x="country",by.y="p3",all.x=T))

countrytotals <- merge(countrytotals,locations,all.x=T) %>% dplyr::rename(n=num_bams)
countrytotals$f3sig[is.na(countrytotals$f3sig)] <- -1
countrytotals$admixed <- c("not tested","not admixed","admixed")[countrytotals$f3sig+2]

ggplot() +
    geom_polygon(data = world, aes(x=long, y = lat, group=group), fill='grey',color="black",size=0.2) +
    geom_point(data=countrytotals, aes(x=lon, y=lat, fill=admixed,size=n), shape=21,inherit.aes=F) +
    #geom_text(data=sample_table, aes(x=lon, y=lat, fill=set,label=Population), shape=21,inherit.aes=F) +
    scale_fill_manual(values=c("red","green","white")) +
    ggtitle(paste("sampling sites",sep="")) +
    coord_fixed(ylim=c(-48,48),xlim=c(-155,140))+
    theme(axis.text=element_blank(),
          panel.border=element_rect(fill=NA, color="black"),
          panel.background=element_rect(fill='light blue', color=NA),
          panel.grid = element_blank(),
          axis.title=element_blank(),
          legend.position="bottom")
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png)


#heatmaps of Z scores

```r
corder <- c("Kenya","Uganda","Gabon","Ghana","BurkinaFaso","Senegal","Brazil" ,"PuertoRico","Trinidad","USA","SaudiArabia","Philippines","Vietnam")

pop3sum$p2 <-factor(pop3sum$p2,levels=corder,ordered=T)
pop3sum$p1 <-factor(pop3sum$p1,levels=corder,ordered=T)

uppertri <- sigmixes$p1>sigmixes$p2
sigmixes[uppertri,c("p1","p2")] <- sigmixes[uppertri,c("p2","p1")]

sigmixes$f3zmean[sigmixes$f3zmean>0] <- 0
sigmixes$f3zmin[sigmixes$f3zmin>0] <- 0

ggplot(sigmixes,aes(x=p1,y=p2,fill=f3zmin)) + 
  geom_raster() + facet_wrap("p3") + coord_fixed() +
  theme(axis.text.x=element_text(angle=45,hjust=1),axis.title=element_blank())
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-1.png)

```r
ggplot(subset(sigmixes,!p1 %in% sigpops & !p2 %in% sigpops),aes(x=p1,y=p2,fill=f3zmin)) + 
  geom_raster() + facet_wrap("p3") + coord_fixed() +
  theme(axis.text.x=element_text(angle=45,hjust=1),axis.title=element_blank())
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-2.png)


#map admixed populations, w arrows

```r
minmixes <- merge(minmixes,locations,by.x="p3",by.y="country")
minmixes <- merge(minmixes,locations,by.x="p1",by.y="country",suffixes=c("","1"))
minmixes <- merge(minmixes,locations,by.x="p2",by.y="country",suffixes=c("","2"))


ggplot() +
    geom_polygon(data = world, aes(x=long, y = lat, group=group), fill='grey',color="black",size=0.2) +
    geom_point(data=countrytotals, aes(x=lon, y=lat, fill=admixed,size=n), shape=21,inherit.aes=F) +
    geom_curve(data=minmixes, 
                 aes(x=lon1, y=lat1, xend=lon,yend=lat),
                 arrow=arrow(length = unit(0.3,"cm")),inherit.aes=F) +
    geom_curve(data=minmixes, aes(x=lon2, y=lat2, xend=lon,yend=lat),
                 arrow=arrow(length = unit(0.3,"cm")),inherit.aes=F) +
    geom_point(data=subset(countrytotals,admixed=="not admixed"), 
               aes(x=lon, y=lat, fill=admixed,size=n), shape=21,inherit.aes=F) +
    ggtitle(paste("3-pop results (unadmixed parents, lowest z-score)",sep="")) +
    scale_fill_manual(values=c("red","green","white")) +
    coord_fixed(ylim=c(-48,48),xlim=c(-155,140))+
    theme(axis.text=element_blank(),
          panel.border=element_rect(fill=NA, color="black"),
          panel.background=element_rect(fill='light blue', color=NA),
          panel.grid = element_blank(),
          axis.title=element_blank(),
          legend.position="bottom")
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png)





```r
aldersum <- read.table("admix_tests/alder_all_summaries.txt",sep="\t", header=T)
aldersum$test.result <- gsub("\\s+.*","",aldersum$test.status,perl=T)
aldersum$inconsistent <- grepl("inconsistent",aldersum$test.status)

#get tested pops
aldertested <- unique(aldersum$test.pop)

aldersum <- aldersum[grep("success*",aldersum$test.status),] %>%      
                     rename(c("test.pop"="p3","ref.B"="p2","ref.A"="p1")) %>% 
                     mutate("set"=paste(p3,p1,p2,sep="/")) %>% 
                     mutate("parent"=paste(p1,p2,sep="/")) %>% 
                     select(c("test.result","inconsistent","p.value","p1","p2","p3",
                        "X2.ref.z.score", #"max.decay.diff..",
                        "X2.ref.decay","X2.ref.amp_exp",
                        "set","parent")) %>% 
                     mutate("parent"=paste(p1,p2,sep="/"))
```

```
## Error in (function (classes, fdef, mtable) : unable to find an inherited method for function 'select' for signature '"data.frame"'
```

```r
sigpopsAlder <- unique(aldersum$p3)
aldersum$p2 <-factor(aldersum$p2,levels=corder,ordered=T)
```

```
## Error in `$<-.data.frame`(`*tmp*`, p2, value = structure(integer(0), .Label = c("Kenya", : replacement has 0 rows, data has 858
```

```r
aldersum$p1 <-factor(aldersum$p1,levels=corder,ordered=T)
```

```
## Error in `$<-.data.frame`(`*tmp*`, p1, value = structure(integer(0), .Label = c("Kenya", : replacement has 0 rows, data has 858
```

```r
#select (discard 1-ref results)

uppertri <- aldersum$p1>aldersum$p2
aldersum[uppertri,c("p1","p2")] <- aldersum[uppertri,c("p2","p1")]
```

```
## Error in `[.data.frame`(aldersum, uppertri, c("p2", "p1")): undefined columns selected
```

```r
minalder <- subset(aldersum,!p2 %in% sigpopsAlder & !p1 %in% sigpopsAlder & !inconsistent) %>% 
                  group_by(p3) %>% 
                  slice_min(order_by = p.value)
```

```
## Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function '%in%': object 'p2' not found
```



```r
minalder <- merge(
                merge(
                  merge(minalder,locations,by.x="p3",by.y="country"),
                                locations,by.x="p1",by.y="country",suffixes=c("","1")),
                                locations,by.x="p2",by.y="country",suffixes=c("","2"))
```

```
## Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'merge': error in evaluating the argument 'x' in selecting a method for function 'merge': error in evaluating the argument 'x' in selecting a method for function 'merge': object 'minalder' not found
```

```r
countrytotals$alder <- "not tested"
countrytotals$alder[countrytotals$country %in% aldertested] <- "not admixed"
countrytotals$alder[countrytotals$country %in% sigpopsAlder] <- "admixed"


ggplot() +
    geom_polygon(data = world, aes(x=long, y = lat, group=group), fill='grey',color="black",size=0.2) +
    geom_point(data=countrytotals, aes(x=lon, y=lat, fill=alder,size=n), shape=21,inherit.aes=F) +
    geom_curve(data=minalder, 
                 aes(x=lon1, y=lat1, xend=lon,yend=lat),
                 arrow=arrow(length = unit(0.3,"cm")),inherit.aes=F) +
    geom_curve(data=minalder, aes(x=lon2, y=lat2, xend=lon,yend=lat),
                 arrow=arrow(length = unit(0.3,"cm")),inherit.aes=F) +
    geom_point(data=subset(countrytotals,alder=="not admixed"), 
               aes(x=lon, y=lat, fill=alder,size=n), shape=21,inherit.aes=F) +
    ggtitle(paste("Alder results (unadmixed parents, lowest p-val)",sep="")) +
    scale_fill_manual(values=c("red","green","white")) +
    coord_fixed(ylim=c(-48,48),xlim=c(-155,140))+
    theme(axis.text=element_blank(),
          panel.border=element_rect(fill=NA, color="black"),
          panel.background=element_rect(fill='light blue', color=NA),
          panel.grid = element_blank(),
          axis.title=element_blank(),
          legend.position="bottom")
```

```
## Error in fortify(data): object 'minalder' not found
```
