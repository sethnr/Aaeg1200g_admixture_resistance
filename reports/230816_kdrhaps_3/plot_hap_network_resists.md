---
output: html_document
editor_options: 
  chunk_output_type: console
---

```r
library("geneHapR")
library("vcfR")
library("pegas")
library("getopt")
library("tidyverse")
```


```r
vcffile <- "maf.kdr3recall.vcf.gz"
prefix <- "kdr3_resist"
metafile <- "resources/meta_Aaeg1kg_spp.txt"

includes=c("315999297")

resists=c("315939224","315983762","315983763","315999297","316014588","316080722")
names(resists) = c("F1534C","V1016G","V1016I","I915K","S723T","V410L")
```


```r
meta <- read.table(metafile,sep="\t",header=T)

metahap <- rbind(meta %>% mutate(sample=paste(sample,"0",sep="_")),
                 meta %>% mutate(sample=paste(sample,"1",sep="_")))
rownames(metahap) <- metahap$sample

poplocs <- read.table("resources/aegy.wgs.pops.list.csv",sep=",",header=T) %>% 
                    select("Pop","Lat","Long") %>% 
                    rename("pop"="Pop","poplat"="Lat","poplong"="Long")
cntlocs <- read.table("resources/country_geocoords.txt",sep="\t",header=T) %>% 
                    rename("cntlat"="lat","cntlong"="lon")


metahaploc <- merge(merge(metahap,cntlocs,by="country"),
                    poplocs,by="pop") %>% rename("Hap"="sample")
rownames(metahaploc) <- metahaploc$Hap


#vcf <- pegas::read.vcf(vcffile)
loci <- VCFloci(vcffile)
```

```
## Scanning file maf.kdr3recall.vcf.gz 
##  0.574225 Mb
## Done.
```

```r
vcfr <- read.vcfR(vcffile)
```

```
## Scanning file to determine attributes.
## File attributes:
##   meta lines: 18
##   header_line: 19
##   variant count: 78
##   column count: 1212
## Meta line 18 read in.
## All meta lines processed.
## gt matrix initialized.
## Character matrix gt created.
##   Character matrix gt rows: 78
##   Character matrix gt cols: 1212
##   skip: 0
##   nrows: 78
##   row_num: 0
## Processed variant: 78
## All variants processed
```

```r
gts <- extract.gt(vcfr)
haps <- extract.haps(vcfr)
```

```
## Variant 0 processedVariant 78 processed
```

```r
#only include samples in haps file (removed related, etc)
haps <- haps[,colnames(haps) %in% metahap$sample]

#paste(rep(VCFlabels(vcffile),each=2),c("0","1"),sep="_")
haptab <- cbind(loci[,c(1,2,4,5,8)],haps)

write.table(haptab,paste(prefix,"gttab.txt",sep="_"),sep="\t",quote=F,row.names=F,col.names=T)

if(is.null(includes)) {
  hapResult <- table2hap(haptab[sort(c(grep("MODERATE",haptab$INFO),
                                grep("HIGH",haptab$INFO))),])
} else {
  hapResult <- table2hap(haptab[sort(c(grep("MODERATE",haptab$INFO),
                                       grep("HIGH",haptab$INFO),
                                       grep(paste(includes,collapse="|"),haptab$POS) )),])
}

write.table(haptab,paste(prefix,"gttab.txt",sep="_"),sep="\t",quote=F,row.names=F,col.names=T)


hapSummary <- hap_summary(hapResult)
```


```r
#get hap GTS as data frame
hapgtsdf <- hapSummary[grepl(paste0("H", "[0-9]{1,}"), hapSummary[,1]), ] %>% as.data.frame()


# #get tree order
#   nhaps = length(hapgtsdf$Hap)
#   hdist = matrix(nrow=nhaps,ncol=nhaps)
#   for(h1i in seq(1,nhaps)){
#     h1 = hapgtsdf$Hap[h1i]
#     h1gt = hapgtsdf[hapgtsdf$Hap==h1,alleletab$POS]
#     for(h2i in seq(1,nhaps)) {
#       h2 = hapgtsdf$Hap[h2i]
#       h2gt = hapgtsdf[hapgtsdf$Hap==h2,alleletab$POS]  
#       hdist[h1i,h2i] = sum(h1gt != h2gt)
#     }
#   }
#   
#   rownames(hdist) = hapgtsdf$Hap
#   colnames(hdist) = hapgtsdf$Hap
#   haptree = ladderize(root(nj(hdist),"H001"))
#   treeorder <- haptree$tip.label[haptree$edge[haptree$edge[,2] < length(haptree$tip.label)+1,2]]
# 
# hapgtsdf$Hap <- factor(hapgtsdf$Hap,levels=treeorder,ordered=T)


alleletab = data.frame(t(hapSummary[c(1,2,4), grepl("[0-9]",hapSummary[c(2), ])])) %>% 
  rename("CHROM"="X1",
        "POS"="X2",
        "ALLELE"="X4") %>% 
  separate_wider_delim("ALLELE","/",names=c("REF","ALT")) %>% 
  mutate("RESIST"=POS %in% resists)

hapgtsdf$resisthap=F
for(pos in alleletab$POS[alleletab$RESIST]){
  pi = grep(pos,alleletab$POS)
  ra = alleletab$ALT[pi]
  hapgtsdf$resisthap[hapgtsdf[as.character(pos)]==ra] = T
}

hapgtsdfm <-  pivot_longer(hapgtsdf,
                           cols = alleletab$POS,
                           names_to = "posn",
                           values_to = "call")



hapgtsdfm$calltype=NA
for(pos in alleletab$POS){
  pi = grep(pos,alleletab$POS)
  alt = alleletab$ALT[pi]
  ref = alleletab$REF[pi]
  if(alleletab$RESIST[pi]) {
    hapgtsdfm$calltype[hapgtsdfm$posn==pos & hapgtsdfm$call==ref]="ref"
    hapgtsdfm$calltype[hapgtsdfm$posn==pos & hapgtsdfm$call==alt]="res"
  } else {
    hapgtsdfm$calltype[hapgtsdfm$posn==pos & hapgtsdfm$call==ref]="ref"
    hapgtsdfm$calltype[hapgtsdfm$posn==pos & hapgtsdfm$call==alt]="alt"
  }
}

gtplot <- ggplot(hapgtsdfm,aes(x=as.factor(posn),y=Hap,fill=calltype,label=call)) + 
  geom_raster() + 
  geom_text() + 
  scale_y_discrete(limits=rev)+
  scale_fill_manual(values=c("res"="red",
                            "ref"="lightgrey",
                            "alt"="darkgrey"),
                    guide=F) +
  theme(axis.text.x=element_text(angle=45,hjust=1),axis.title=element_blank())

hapfreq <- ggplot(hapgtsdf,aes(x=freq,y=Hap,fill=resisthap,label=freq)) + 
  geom_bar(stat="identity") + 
  geom_text(hjust=-0.2) +
  theme(axis.text.y=element_blank(),axis.title=element_blank()) +
  scale_y_discrete(limits=rev) +
  scale_fill_manual(values=c("darkgrey","red"),guide=F)

gtplot + hapfreq + plot_layout(ncol=2,widths=c(9,1))
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png)

```r
W=300
H=200
R=400
png(paste(prefix,"haptable_resist.png",sep="_"),width=W,height=H,res=R,units="mm")
#plot haplotype table
gtplot + hapfreq + plot_layout(ncol=2,widths=c(9,1))
dev.off()
```

```
## RStudioGD 
##         2
```

```r
plotHapTable(hapSummary)
```

```
## Warning in plotHapTable(hapSummary): using first 'angle': 0
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-2.png)

```
## Warning: Removed 1 rows containing missing values (`geom_text()`).
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-3.png)

```r
png(paste(prefix,"haptable.png",sep="_"),width=W,height=H,res=R,units="mm")
#plot haplotype table
plotHapTable(hapSummary)
```

```
## Warning in plotHapTable(hapSummary): using first 'angle': 0
```

```
## Warning: Removed 1 rows containing missing values (`geom_text()`).
```

```r
dev.off()
```

```
## RStudioGD 
##         2
```



```r
#plot haplotype network
hapNet <- get_hapNet(hapSummary,
                     AccINFO = metahap,
                     groupName = "region")

# hapcols = c("West Africa"="dark blue","East Africa"="light blue",
#             "Carribean"="Purple","South America"="Dark Red","North America"="Light Red",
#             "Middle East"="Light Green","Asia"="Green","Pacific"="Dark Green")
plotHapNet(hapNet,
           show.mutation = 2,
           legend=c(-35,-10),
           show_size_legend=F,
           scale='log2',
           threshold=0) #+ scale_fill_manual(values=hapcols)
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png)

```r
png(paste(prefix,"hapnet.png",sep="_"),width=W,height=H,res=R,units="mm")
plotHapNet(hapNet,
           show.mutation = 2,
           legend=c(-35,-10),
           show_size_legend=F,
           scale='log2',
           threshold=0) #+ scale_fill_manual(values=hapcols)
dev.off()
```

```
## RStudioGD 
##         2
```

```r
hapDistribution(hapResult,metahaploc,"poplong","poplat",
                c("H001","H002","H003","H004","H005","H006"),
                legend=TRUE)
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-2.png)![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-3.png)

```r
png(paste(prefix,"geodist_pop.png",sep="_"),width=W,height=H,res=R,units="mm")
hapDistribution(hapResult,metahaploc,"poplong","poplat",
                c("H001","H002","H003","H004","H005","H006"),
                legend=TRUE)
dev.off()
```

```
## RStudioGD 
##         2
```

```r
png(paste(prefix,"geodist_cnt.png",sep="_"),width=W,height=H,res=R,units="mm")
hapDistribution(hapResult,metahaploc,"cntlong","cntlat",
                c("H001","H002","H003","H004","H005",
                  "H006","H007","H008","H009","H010"),
                legend=TRUE)
dev.off()
```

```
## RStudioGD 
##         2
```




```r
hapResultRes <- table2hap(haptab[grep(paste(resists,collapse="|"),haptab$POS),])

plotHapTable(hapResultRes)
```

```
## Warning in plotHapTable(hapResultRes): using first 'angle': 0
```

```
## Warning: Removed 1 rows containing missing values (`geom_text()`).
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-1.png)

```r
hapResultResDF <- as.data.frame(hapResultRes[grepl(paste0("H", "[0-9]{1,}"),hapResultRes$Hap)
                           ,c("Hap","Accession")])
hapResultResDF$Hap <- factor(hapResultResDF$Hap,ordered=T)

levels(hapResultResDF$Hap) <- c("S",
                                sprintf('R%0.2d',
                                        c(1:(length(levels(hapResultResDF$Hap))-1))))


hapResultResDF$Sample <- gsub("_[0:1]$","",hapResultResDF$Accession)

sampleResultRes <- merge(hapResultResDF[grep("*_0",
                                             hapResultResDF$Accession),
                                        c("Hap","Sample")],
                        hapResultResDF[grep("*_1",
                                            hapResultResDF$Accession),
                                        c("Hap","Sample")],by="Sample")

sampleResultRes <- merge(sampleResultRes,meta,by.x="Sample",by.y="sample")

swapindex <- sampleResultRes$Hap.y<sampleResultRes$Hap.x
newx = sampleResultRes$Hap.y[swapindex]
newy = sampleResultRes$Hap.x[swapindex]
sampleResultRes$Hap.y[swapindex]=newy
sampleResultRes$Hap.x[swapindex]=newx

sampleResultRes$geno = paste(sampleResultRes$Hap.x,sampleResultRes$Hap.y,sep="/")

haptotals <- sampleResultRes[c("Hap.y","Hap.x")] %>% 
                              group_by(Hap.y,Hap.x,.drop = FALSE) %>% 
                              summarise(count=n()) %>%
                              mutate(count = replace(count, count == 0, NA)) %>%
                              filter(Hap.x <= Hap.y) %>%
                              filter(Hap.x < "R06",Hap.y < "R06")
```

```
## `summarise()` has grouped output by 'Hap.y'. You can override using the `.groups` argument.
```

```r
ggplot(haptotals,aes(x=Hap.x,y=Hap.y,fill=count)) + geom_raster()
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-2.png)

```r
countrytotals <- sampleResultRes[c("country","Hap.y","Hap.x")] %>%   
                              group_by(country,Hap.y,Hap.x,.drop = FALSE) %>% 
                              summarise(count=n()) %>%
                              mutate(count = replace(count, count == 0, NA)) %>%
                              filter(Hap.x <= Hap.y) %>%
                              filter(Hap.x < "R06",Hap.y < "R06")
```

```
## `summarise()` has grouped output by 'country', 'Hap.y'. You can override using the `.groups`
## argument.
```

```r
ggplot(countrytotals,aes(x=Hap.x,y=Hap.y,fill=count)) + geom_raster() + facet_wrap("country")
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-3.png)

```r
regiontotals <- sampleResultRes[c("region","Hap.y","Hap.x")] %>% 
                              group_by(region,Hap.y,Hap.x,.drop = FALSE) %>% 
                              summarise(count=n()) %>%
                              mutate(count = replace(count, count == 0, NA)) %>%
                              filter(Hap.x <= Hap.y) %>%
                              filter(Hap.x < "R06",Hap.y < "R06")
```

```
## `summarise()` has grouped output by 'region', 'Hap.y'. You can override using the `.groups`
## argument.
```

```r
ggplot(regiontotals,aes(x=Hap.x,y=Hap.y,fill=count)) + geom_raster() + facet_wrap("region")
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-4.png)



```r
sampleResultRes$resistF=1
sx = sampleResultRes$Hap.x=="S"
sampleResultRes$resistF[sx] = sampleResultRes$resistF[sx] - 0.5
sy = sampleResultRes$Hap.y=="S"
sampleResultRes$resistF[sy] = sampleResultRes$resistF[sy] - 0.5

world <- map_data("world")

popresists <- sampleResultRes[c("pop","resistF")] %>%   
                              group_by(pop,.drop = FALSE) %>% 
                              summarise(resistance_freq=mean(resistF)) %>% 
                              mutate("pop" = gsub("_","",pop)) %>%
                              mutate(pop=gsub("Virhembe","Virembe",pop)) %>% 
                              mutate(pop=gsub("Cuanda","Luanda",pop))

metapops <- read.table("aegy.wgs.pops.list_display.csv",sep=",",header=T,stringsAsFactors = F) %>% 
                rename_with(tolower) %>% 
                mutate("pop" = gsub("_","",pop)) %>%
                mutate(pop=gsub("Cuanda","Luanda",pop)) %>% 
                mutate(displaced = !is.na(latdisp))
metapops$latdisp[!metapops$displaced] <- metapops$lat[!metapops$displaced]
metapops$longdisp[!metapops$displaced] <- metapops$long[!metapops$displaced]

popresists = merge(popresists,metapops,by="pop",all=T)


ggplot(popresists,aes(x=longdisp, y=latdisp, label=pop,xend=long, yend=lat)) +
    geom_polygon(data = world, aes(x=long, y=lat, group=group), fill='#C7ECD3',color="black",size=0.2,inherit.aes=F) +
    #geom_point(data=subset(poptotals,is.na(latdisp)), aes(x=long, y=lat, fill=f3zdisp,size=num_bams), shape=21,inherit.aes=F) +
    geom_segment(data=subset(popresists,displaced)) +
    geom_point(data=popresists, aes(fill=resistance_freq,size=num_bams), shape=21) +
    #geom_text(data=subset(popresists,pop %in% llab), hjust=1.3) +
    #geom_text(data=subset(popresists,pop %in% rlab), hjust=-0.3) +
    coord_fixed(ylim=c(-48,48),xlim=c(-155,140))+
    ggtitle(paste("VGSC resistance haplotype prevalence",sep="")) +
    scale_fill_gradient(low="white",high="red",na.value="grey") +
    theme(axis.text=element_blank(),
          panel.border=element_rect(fill=NA, color="black"),
          panel.background=element_rect(fill='#C7E6F1', color=NA),
          panel.grid = element_blank(),
          axis.title=element_blank(),
          axis.ticks=element_blank(),
          legend.position="bottom")
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png)



```r
F1534C <- haps["NC_035109.1_315939224",]
F1534C <- data.frame("hap"=names(F1534C),"call"=F1534C)
rownames(F1534C) <- c(1:nrow(F1534C))

F1534C <- F1534C %>% 
  mutate("freq" = as.numeric(call=="C")/2) %>% 
  mutate("sample" = gsub("_[0|1]$","",hap)) %>% 
  group_by(sample) %>% 
  summarize("freq"=sum(freq))

F1534C <- merge(F1534C,meta,by="sample") %>% 
  mutate("pop" = gsub("_","",pop)) %>%
  mutate(pop=gsub("Virhembe","Virembe",pop)) %>% 
  mutate(pop=gsub("Cuanda","Luanda",pop)) %>%
  group_by(pop) %>%
  summarise(resistance_freq=mean(freq)) %>% 
  merge(metapops,by="pop")


ggplot(popresists,aes(x=longdisp, y=latdisp, label=pop,xend=long, yend=lat)) +
    geom_polygon(data = world, aes(x=long, y=lat, group=group), fill='#C7ECD3',color="black",size=0.2,inherit.aes=F) +
    geom_segment(data=subset(popresists,displaced)) +
    geom_point(data=F1534C, aes(fill=resistance_freq,size=num_bams), shape=21) +
    coord_fixed(ylim=c(-48,48),xlim=c(-155,140))+
    ggtitle(paste("VGSC F1534C prevalence",sep="")) +
    scale_fill_gradient(low="white",high="red",na.value="grey") +
    theme(axis.text=element_blank(),
          panel.border=element_rect(fill=NA, color="black"),
          panel.background=element_rect(fill='#C7E6F1', color=NA),
          panel.grid = element_blank(),
          axis.title=element_blank(),
          axis.ticks=element_blank(),
          legend.position="bottom")
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-1.png)




```r
V1016I <- haps["NC_035109.1_315983763",]
V1016I <- data.frame("hap"=names(V1016I),"call"=V1016I)
rownames(V1016I) <- c(1:nrow(V1016I))

V1016I <- V1016I %>% 
  mutate("freq" = as.numeric(call=="T")/2) %>% 
  mutate("sample" = gsub("_[0|1]$","",hap)) %>% 
  group_by(sample) %>% 
  summarize("freq"=sum(freq))

V1016I <- merge(V1016I,meta,by="sample") %>% 
  mutate("pop" = gsub("_","",pop)) %>%
  mutate(pop=gsub("Virhembe","Virembe",pop)) %>% 
  mutate(pop=gsub("Cuanda","Luanda",pop)) %>%
  group_by(pop) %>%
  summarise(resistance_freq=mean(freq)) %>% 
  merge(metapops,by="pop")


ggplot(popresists,aes(x=longdisp, y=latdisp, label=pop,xend=long, yend=lat)) +
    geom_polygon(data = world, aes(x=long, y=lat, group=group), fill='#C7ECD3',color="black",size=0.2,inherit.aes=F) +
    geom_segment(data=subset(popresists,displaced)) +
    geom_point(data=V1016I, aes(fill=resistance_freq,size=num_bams), shape=21) +
    coord_fixed(ylim=c(-48,48),xlim=c(-155,140))+
    ggtitle(paste("VGSC V1016I prevalence",sep="")) +
    scale_fill_gradient(low="white",high="red",na.value="grey") +
    theme(axis.text=element_blank(),
          panel.border=element_rect(fill=NA, color="black"),
          panel.background=element_rect(fill='#C7E6F1', color=NA),
          panel.grid = element_blank(),
          axis.title=element_blank(),
          axis.ticks=element_blank(),
          legend.position="bottom")
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-1.png)




```r
V1016G <- haps["NC_035109.1_315983762",]
V1016G <- data.frame("hap"=names(V1016G),"call"=V1016G)
rownames(V1016G) <- c(1:nrow(V1016G))

V1016G <- V1016G %>% 
  mutate("freq" = as.numeric(call=="C")/2) %>% 
  mutate("sample" = gsub("_[0|1]$","",hap)) %>% 
  group_by(sample) %>% 
  summarize("freq"=sum(freq))

V1016G <- merge(V1016G,meta,by="sample") %>% 
  mutate("pop" = gsub("_","",pop)) %>%
  mutate(pop=gsub("Virhembe","Virembe",pop)) %>% 
  mutate(pop=gsub("Cuanda","Luanda",pop)) %>%
  group_by(pop) %>%
  summarise(resistance_freq=mean(freq)) %>% 
  merge(metapops,by="pop")


ggplot(popresists,aes(x=longdisp, y=latdisp, label=pop,xend=long, yend=lat)) +
    geom_polygon(data = world, aes(x=long, y=lat, group=group), fill='#C7ECD3',color="black",size=0.2,inherit.aes=F) +
    geom_segment(data=subset(popresists,displaced)) +
    geom_point(data=V1016G, aes(fill=resistance_freq,size=num_bams), shape=21) +
    coord_fixed(ylim=c(-48,48),xlim=c(-155,140))+
    ggtitle(paste("VGSC V1016G prevalence",sep="")) +
    scale_fill_gradient(low="white",high="red",na.value="grey") +
    theme(axis.text=element_blank(),
          panel.border=element_rect(fill=NA, color="black"),
          panel.background=element_rect(fill='#C7E6F1', color=NA),
          panel.grid = element_blank(),
          axis.title=element_blank(),
          axis.ticks=element_blank(),
          legend.position="bottom")
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10-1.png)
