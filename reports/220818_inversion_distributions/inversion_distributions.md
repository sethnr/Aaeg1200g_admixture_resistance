---
output: html_document
editor_options: 
  chunk_output_type: console
---



```r
callsfiles <- list.files("lostruct_merge",pattern="invs_chr.*_inv_calls.txt",full.names = T)
callsfiles <- callsfiles[laply(callsfiles, file.size)>0]
rm("calltab")

for(x in callsfiles) {
  icalls <- read.table(x,header=T,sep="\t")
  invnames <- colnames(icalls)[!colnames(icalls) %in% c("sample","spp","pop","region","contgroup","country")]
  sample <- icalls$sample
  icallmat <- icalls[,invnames]
  if(!exists("calltab")) {
    calltab <- cbind(sample,icallmat)
  } else {
    calltab <- cbind(calltab,icallmat)
  }
  
}

meta <- icalls[,c("sample","pop","region","country")]
calltab <- merge(calltab,meta)


countries <- unique(calltab$country)
invnames <- colnames(calltab)[!colnames(calltab) %in% c("sample","spp","pop","region","contgroup","country")]
maf <- function(x) {round(sum(x) / (length(x)*2),2)}

countryfreqs <- aggregate(. ~ country,data=calltab[,c("country",as.character(invnames))],FUN=maf)
contfreqs <- aggregate(. ~ region,data=calltab[,c("region",as.character(invnames))],FUN=maf)
popfreqs <- aggregate(. ~ pop,data=calltab[,c("pop",as.character(invnames))],FUN=maf)

countryN <- as.data.frame(table(calltab[,c("country")]))
colnames(countryN) <- c("country","n")
contN <- as.data.frame(table(calltab[,c("region")]))
colnames(contN) <- c("cont","n")
popN <- as.data.frame(table(calltab[,c("pop")]))
colnames(popN) <- c("pop","n")



poplocs <- read.table("../../resources/aegy.wgs.pops.list.csv",sep=",",header=T)
```



```r
countryfile <- "country_geocoords.txt"
if(file.exists(countryfile)){
  countrylocs <- read.table(countryfile,header=T) 
} else {
  register_google(key="AIzaSyBFbT5X6MjrepPAlzIj3zz1HS58hnIPtEM", write=T)
  countries <- unique(invfreqs$country)
  countrylocs <- geocode(countries)
  countrylocs$country <- countries
  write.table(countrylocs,countryfile,sep="\t",row.names=F,col.names=T)
}

world <- map_data("world")

countryfqM <- pivot_longer(countryfreqs,cols = invnames,names_to="inv",values_to="maf")
countryfqM <- merge(merge(countryfqM,countrylocs),countryN)
ncols <- floor(sqrt(length(invnames)))
ggplot() +
  geom_polygon(data = world, aes(x=long, y = lat, group=group), fill='grey',color="black",size=0.2) +
  geom_point(data=countryfqM, aes(x=lon, y=lat, fill=maf,size=n), shape=21,inherit.aes=F) +
  ggtitle(paste("inversion distributions country",sep="")) +
  scale_fill_distiller(limits=c(0,1),palette="Reds",direction = 1)+
  coord_fixed(ylim=c(-48,48),xlim=c(-155,140))+
  theme(axis.text=element_blank(),
        panel.border=element_rect(fill=NA, color="black"),
        panel.background=element_rect(fill='light blue', color=NA),
        panel.grid = element_blank(),
        axis.title=element_blank(),
        legend.position="bottom") +
  facet_wrap(inv ~ .,ncol=ncols)
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2-1.png)


#calculate HWE

```r
hwetab <- data.frame("inversion"=factor(levels=invnames),
           "country"=factor(levels=unique(calltab$country)),
           "aa"=numeric(),
           "ab"=numeric(),
           "bb"=numeric(),
           "HWE"=numeric()
           )

for (C in unique(calltab$country)) {
      countrycalls <- calltab[calltab$country==C,]
      for(I in invnames) {
        geno <- genotype(c("a/a","a/b","b/b")[calltab[calltab$country==C,I]+1])
        if(nallele(geno)==2) {
          test <- HWE.test(geno)
          Pval <- test$test$p.value
        } else{
          Pval <- 1
        }
        i <- nrow(hwetab)+1
        hwetab[i,c("aa","ab","bb")] <- c(sum(geno=="a/a"),sum(geno=="a/b"),sum(geno=="b/b"))
        hwetab[i,c("country","inversion")] <- c(C,I)
        hwetab[i,c("HWE")] <- Pval
      }
}
```

```
## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call
```

```r
hwetabP <- subset(hwetab,HWE<0.05)
```

#get and order countries / inversions 

```r
#order countries in metatable by region (w african, eafrican, american, asian)
regorder <- c("East Africa","West Africa","Caribean","South America","North America","Middle East","Asia","Pacific")
meta$region <- factor(meta$region,levels=regorder,ordered=T)
cntorder <- unique(meta$country[order(meta$region)])

blocks <- ldply(as.list(list.files("./lostruct_merge","*blocks.txt",full.names = T)),read_tsv)
```

```
## Rows: 233 Columns: 7
## ── Column specification ────────────────────────────────────────────────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (5): inv, block, region, id, cluster
## dbl (2): chrom, pos
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
## Rows: 129 Columns: 7
## ── Column specification ────────────────────────────────────────────────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (5): inv, block, region, id, cluster
## dbl (2): chrom, pos
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
## Rows: 0 Columns: 0
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
blocks$y=1

blocks$inversion <- paste("X",blocks$cluster,sep="")
blocks <- blocks[blocks$inversion %in% unique(hwetab$inversion),]
blockextents <- blocks %>% dplyr::group_by(inversion) %>% dplyr::summarise(min=min(pos),mid=mean(pos),max=max(pos),chrom=min(chrom),y=1)
blockextents$size <- blockextents$max-blockextents$min

write("calculating display levels",stderr())
    clustersL2S <- blockextents$inversion[order(blockextents$size,decreasing=T)]
    clustersL2R <- blockextents$inversion[order(blockextents$min,decreasing=F)]

    for(C1 in clustersL2S) {
      for(C2 in clustersL2S) {
        chrC1 = blockextents$chrom[blockextents$inversion==C1]
        lenC1 = blockextents$size[blockextents$inversion==C1]
        minC1 = blockextents$min[blockextents$inversion==C1]
        maxC1 = blockextents$max[blockextents$inversion==C1]
        yC1   = min(blocks$y[blocks$inversion==C1])

        chrC2 = blockextents$chrom[blockextents$inversion==C2]
        lenC2 = blockextents$size[blockextents$inversion==C2]
        minC2 = blockextents$min[blockextents$inversion==C2]
        maxC2 = blockextents$max[blockextents$inversion==C2]
        yC2   = min(blocks$y[blocks$inversion==C2])

        if((minC2<maxC1 & maxC2>minC1 & chrC1==chrC2)) {
          if(lenC2<lenC1 & yC1==yC2 & chrC1==chrC2) {
            newYC2 <- yC2+1
            blocks$y[blocks$inversion==C2] <- newYC2
            blockextents$y[blockextents$inversion==C2] <- newYC2
            write(paste(C1,":",yC1," ",C2,":",yC2,"->",yC2+1,sep=""),stderr())
          }
        }
      }
    }





clustcols <- sample(rainbow(nrow(blockextents),s=0.6,v=0.9))
names(clustcols) <- blockextents$inversion
```


```r
hwetab$country <- factor(hwetab$country,levels=cntorder,ordered=T)
hwetab$inversion <- factor(hwetab$inversion,levels=clustersL2R,ordered=T)
hwetab <- merge(hwetab,unique(blocks[,c("inversion","chrom")]))
hwetabM <- pivot_longer(hwetab,cols=c("aa","ab","bb"),names_to = "call",values_to="count")

callplot <- ggplot(hwetabM,aes(x=inversion,y=count,fill=call)) + 
  geom_bar(stat="identity",width=0.99) + 
  geom_bar(data=subset(hwetab,HWE<0.01),aes(y=aa+ab+bb),stat="identity",fill=NA,color="cyan",size=1.5,width=0.99) + 
  facet_grid(country ~ chrom,scale="free",space="free") +
  scale_fill_manual(values=c("aa"="blue","ab"="purple","bb"="red")) +
  scale_x_discrete(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  theme(panel.spacing.x = unit(2, "mm"),
        panel.spacing.y = unit(0.2, "mm"),
              axis.text.y=element_blank(),
              axis.title=element_blank(),
              axis.ticks.y=element_blank(),
        plot.margin = unit(c(1,1,1,1), "mm"))

chromlen <- data.frame("min"=0,
                       "chrom"=c(1,2),
                       "max"=c(310827022,474425716))
blockclustplot <- ggplot(blocks,aes(x=pos,y=y,fill=as.factor(inversion))) +
      geom_raster() +
      geom_segment(data=chromlen,aes(x=min,xend=max,y=0,yend=0,color=NA),inherit.aes=F) +
      geom_segment(data=blockextents,aes(x=min,xend=max,y=y,yend=y,color=as.factor(inversion)),inherit.aes=F) +
      geom_text(data=blockextents,aes(x=mid,y=y,label=inversion),inherit.aes=F) +
      scale_fill_manual(values = clustcols)+
      scale_color_manual(values = clustcols)+
      facet_grid(. ~ chrom,scale="free_x",space="free_x") +
      theme(axis.title=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
            legend.position="none",legend.title.align = 1) +
      scale_x_continuous(expand = c(0,0,0,0))
#blockclustplot

callplot / blockclustplot + plot_layout(heights=c(9,1))
```

```
## Warning: Raster pixels are placed at uneven horizontal intervals and will be shifted. Consider using geom_tile() instead.
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png)



```r
popfqM <- pivot_longer(popfreqs,cols = invnames,names_to="inv",values_to="maf")
popfqM <- merge(merge(popfqM,poplocs,by.x="pop",by.y="Pop"),popN)

popfqM$pop[popfqM$pop == "Lope_village"] <- "LopeVillage"
popfqM$pop[popfqM$pop == "Virhembe"] <- "Virembe"

popfqM$Country <- gsub(" ","",popfqM$Country)


cmaps <- list()
bigcountries <- c("Brazil","SaudiArabia","Argentina","USA")
for(C in unique(hwetabP$country)) {
  lon <- countrylocs$lon[countrylocs$country==C]
  lat <- countrylocs$lat[countrylocs$country==C]
  if(C %in% bigcountries) {
    cmap <- get_map(c(lon,lat),zoom=4)
  } else {
    cmap <- get_map(c(lon,lat),zoom=6)}
  cmaps[[C]] <- cmap
}
```

```
## Warning in gzfile(file, "rb"): cannot open compressed file '/var/folders/33/h3nj351j6fgd7_jh5_cjk60r0000gn/T//RtmpopYILS/
## ggmap/index.rds', probable reason 'No such file or directory'
```

```
## Error in gzfile(file, "rb"): cannot open the connection
```

```r
hwIs <- unique(hwetabP$inversion)
for(I in hwIs) {
  ggmaps <- list()
  for(C in hwetabP$country[hwetabP$inversion==I]) {
  
  ccfreqs <- popfqM[popfqM$Country==C & popfqM$inv==I,]
  ggmaps[[C]] <- ggmap(cmaps[[C]]) +
      geom_point(data=ccfreqs, aes(x=Long, y=Lat, fill=maf,size=n), shape=21,inherit.aes=F) +
    #   #geom_text(data=sample_table, aes(x=lon, y=lat, fill=set,label=Population), shape=21,inherit.aes=F) +
      scale_fill_distiller(limits=c(0,1),palette="Reds",direction = 1)+
      scale_size_continuous(limits=c(0,50))+
      ggtitle(I)+
      theme(axis.text=element_blank(),
            panel.border=element_rect(fill=NA, color="black"),
            panel.background=element_rect(fill='light blue', color=NA),
            axis.title=element_blank(),
            legend.position="none")
  }
  do.call(grid.arrange,ggmaps)
}
```

```
## Error: ggmap plots objects of class ggmap, see ?get_map
```
