

``` r
library("geneHapR")
library("vcfR")
library("pegas")
library("getopt")
library("tidyverse")
library("patchwork")
library("plyr")
```


``` r
vcffile <-"./snpeff.kdr3recall.maf0.005.vcf.gz"
#vcffile <-"./locus.kdr3recall.fixref.vcf.gz"
prefix <- "kdr3nj"
metafile <- "resources/meta_Aaeg1kg_spp.txt"


resists=c(315939224,315983762,315983763,315999297,316014588,316080722)

includes=c(315939224,315983762,315983763,315999297,316014588,316080722)

read.table("resist_loci.txt")
```

```
##       V1        V2
## 1 F1534C 315939224
## 2 V1016G 315983762
## 3 V1016I 315983763
## 4  I915K 315999297
## 5  S723T 316014588
## 6  V410L 316080722
```

``` r
regorder <- c("East Africa","West Africa","Caribean","North America","South America","Middle East","Asia","Pacific")
meta <- read.table(metafile,sep="\t",header=T) %>%
            mutate(region=factor(region,levels=regorder,ordered=T))


metahap <- rbind(meta %>% mutate(sample=paste(sample,"0",sep="_")),
                 meta %>% mutate(sample=paste(sample,"1",sep="_")))
rownames(metahap) <- metahap$sample

poplocs <- read.table("resources/aegy.wgs.pops.list.csv",sep=",",header=T) %>% 
                    dplyr::select("Pop","Lat","Long") %>% 
                    dplyr::rename("pop"="Pop","poplat"="Lat","poplong"="Long")
cntlocs <- read.table("resources/country_geocoords.txt",sep="\t",header=T) %>% 
                    dplyr::rename("cntlat"="lat","cntlong"="lon")


metahaploc <- merge(merge(metahap,cntlocs,by="country"),
                    poplocs,by="pop") %>% dplyr::rename("Hap"="sample")
rownames(metahaploc) <- metahaploc$Hap
```


``` r
loci <- VCFloci(vcffile)
```

```
## Scanning file ./snpeff.kdr3recall.maf0.005.vcf.gz 
##  0.78469 Mb
## Done.
```

``` r
vcfr <- read.vcfR(vcffile)
```

```
## Scanning file to determine attributes.
## File attributes:
##   meta lines: 19
##   header_line: 20
##   variant count: 109
##   column count: 1212
## Meta line 19 read in.
## All meta lines processed.
## gt matrix initialized.
## Character matrix gt created.
##   Character matrix gt rows: 109
##   Character matrix gt cols: 1212
##   skip: 0
##   nrows: 109
##   row_num: 0
## Processed variant: 109
## All variants processed
```

``` r
gts <- extract.gt(vcfr)
haps <- extract.haps(vcfr)
```

```
## Variant 0 processedVariant 109 processed
```

``` r
#only include samples in haps file (removed related, etc)
hapnames <-  colnames(haps)[colnames(haps) %in% metahap$sample]
haps <- haps[,hapnames]
metahap <- metahap[hapnames,]

#paste(rep(VCFlabels(vcffile),each=2),c("0","1"),sep="_")
haptab <- cbind(loci[,c(1,2,4,5,8)],haps)

#check if 297 still present
haptab[haptab$POS==315999297,1:10]
```

```
##                             CHROM       POS REF ALT
## NC_035109.1_315999297 NC_035109.1 315999297   G   T
##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  INFO
## NC_035109.1_315999297 AC=1271;AF=0.528263;CM=0.068417;ANN=T|synonymous_variant|LOW|AAEL023266|AAEL023266.2|transcript|AAEL023266-RC|protein_coding|24/37|c.2712C>A|p.Ile904Ile|4429/12626|2712/6390|904/2129||,T|synonymous_variant|LOW|AAEL023266|AAEL023266|transcript|AAEL023266-RE|protein_coding|25/38|c.2733C>A|p.Ile911Ile|4449/12577|2733/6342|911/2113||,T|synonymous_variant|LOW|AAEL023266|AAEL023266|transcript|AAEL023266-RG|protein_coding|26/39|c.2742C>A|p.Ile914Ile|4458/12586|2742/6351|914/2116||,T|synonymous_variant|LOW|AAEL023266|AAEL023266|transcript|AAEL023266-RH|protein_coding|25/38|c.2709C>A|p.Ile903Ile|4425/12553|2709/6318|903/2105||,T|synonymous_variant|LOW|AAEL023266|AAEL023266|transcript|AAEL023266-RK|protein_coding|24/37|c.2712C>A|p.Ile904Ile|4428/12625|2712/6390|904/2129||,T|synonymous_variant|LOW|AAEL023266|AAEL023266|transcript|AAEL023266-RM|protein_coding|24/37|c.2712C>A|p.Ile904Ile|4428/12523|2712/6288|904/2095||;AN=2406
##                       Amacuzac_Mexico-_12_0 Amacuzac_Mexico-_12_1 Amacuzac_Mexico-_15_0
## NC_035109.1_315999297                     T                     G                     T
##                       Amacuzac_Mexico-_15_1 Amacuzac_Mexico-_2_0
## NC_035109.1_315999297                     T                    T
```



``` r
#filter out non-synonymous & manually re-add known resistants
if(is.null(includes)) {
  hapResult <- table2hap(haptab[sort(c(grep("MODERATE",haptab$INFO),
                                grep("HIGH",haptab$INFO))),])
} else {
  hapResult <- table2hap(haptab[unique(sort(c(grep("MODERATE",haptab$INFO),
                                       grep("HIGH",haptab$INFO),
                                       grep(paste(includes,collapse="|"),haptab$POS) ))),])
}

haptab <- haptab[unique(sort(c(grep("MODERATE",haptab$INFO),
                                       grep("HIGH",haptab$INFO),
                                       grep(paste(includes,collapse="|"),haptab$POS) ))),]

haptab[haptab$POS==315999297,1:10]
```

```
##                             CHROM       POS REF ALT
## NC_035109.1_315999297 NC_035109.1 315999297   G   T
##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  INFO
## NC_035109.1_315999297 AC=1271;AF=0.528263;CM=0.068417;ANN=T|synonymous_variant|LOW|AAEL023266|AAEL023266.2|transcript|AAEL023266-RC|protein_coding|24/37|c.2712C>A|p.Ile904Ile|4429/12626|2712/6390|904/2129||,T|synonymous_variant|LOW|AAEL023266|AAEL023266|transcript|AAEL023266-RE|protein_coding|25/38|c.2733C>A|p.Ile911Ile|4449/12577|2733/6342|911/2113||,T|synonymous_variant|LOW|AAEL023266|AAEL023266|transcript|AAEL023266-RG|protein_coding|26/39|c.2742C>A|p.Ile914Ile|4458/12586|2742/6351|914/2116||,T|synonymous_variant|LOW|AAEL023266|AAEL023266|transcript|AAEL023266-RH|protein_coding|25/38|c.2709C>A|p.Ile903Ile|4425/12553|2709/6318|903/2105||,T|synonymous_variant|LOW|AAEL023266|AAEL023266|transcript|AAEL023266-RK|protein_coding|24/37|c.2712C>A|p.Ile904Ile|4428/12625|2712/6390|904/2129||,T|synonymous_variant|LOW|AAEL023266|AAEL023266|transcript|AAEL023266-RM|protein_coding|24/37|c.2712C>A|p.Ile904Ile|4428/12523|2712/6288|904/2095||;AN=2406
##                       Amacuzac_Mexico-_12_0 Amacuzac_Mexico-_12_1 Amacuzac_Mexico-_15_0
## NC_035109.1_315999297                     T                     G                     T
##                       Amacuzac_Mexico-_15_1 Amacuzac_Mexico-_2_0
## NC_035109.1_315999297                     T                    T
```

``` r
hapSummary <- hap_summary(hapResult)
#get hap GTS as data frame
hapgtsdf <- hapSummary[grepl(paste0("H", "[0-9]{1,}"), hapSummary[,1]), ] %>% as.data.frame()
```

                                


``` r
# #get allele table
# 
# 
alleletab = data.frame(t(hapSummary[c(1,2,4), grepl("[0-9]",hapSummary[c(2), ])])) %>%
  dplyr::rename("CHROM"="X1",
        "POS"="X2",
        "ALLELE"="X4") %>%
  separate_wider_delim("ALLELE","/",names=c("REF","ALT")) %>%
  mutate("RESIST"=POS %in% resists)
```

    


``` r
#filter for hap alleles in final haptab (nonsyn or known resist)
goodpos <- gsub("NC_035109.1_","",rownames(haps)) %in% alleletab$POS 
haps <- haps[goodpos,]
haptab <- cbind(loci[goodpos,c(1,2,4,5,8)],haps)
haptab$resist = haptab$POS %in% resists

#get nj tree
hamdist <- dist.hamming(t(haps))

# njtree <- nj(hdist)
# 
# njtree <- midpoint_root(njtree)
# ggtree(njtree)


library(ggdendro)

hapclust <- hclust(hamdist)


#get order of clusters, then...
#make a new-to-old mapping so clusters run 1-47 left-to-right on graph
#clusters <- cutree(hapclust,h=0)
#treeorder <- hapclust$labels[hapclust$order]

#get treeorder and clusters when ladderized
treeorder <- rev(get_taxa_name(ggtree(hapclust)))
clusters <- cutree(ape::ladderize(as.phylo(hapclust)),h=0)

clorder <- clusters[treeorder][!duplicated(clusters[treeorder])]
clsort <- c(1:length(clorder))[order(clorder)]
newclusters <- clsort[clusters]
names(newclusters) <- hapclust$labels

#add to metahap table
metahap$clusters <- newclusters[metahap$sample]



newclusters[metahap$sample][is.na(newclusters[metahap$sample])]
```

```
## named integer(0)
```

``` r
#get colors for tree tips
treelabels <- rep(meta$region,each=2)
names(treelabels) <- hapnames
treelabels <- treelabels[treeorder]
tipcols <- rainbow(length(unique(treelabels)))[as.factor(treelabels)]


treeplot <- ggtree(hapclust, branch.length='none',ladderize = FALSE) + 
  layout_dendrogram() + 
  theme(axis.text=element_blank(),
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        panel.background=element_blank())




treeplot2 <- ggtree(hapclust,right=F) +
  layout_dendrogram() + 
  theme(axis.text=element_blank(),
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        panel.background=element_blank())

treeplot / treeplot2
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-1.png)


``` r
#make individual hap table


haptabm <- haptab %>% 
  pivot_longer(all_of(hapnames),names_to="sample",values_to="allele") %>% 
  mutate(genotype=as.numeric(allele==ALT))

haptabm <- merge(haptabm,metahap,by="sample")

  

haptabm$calltype = "ref"
haptabm$calltype[haptabm$genotype > 0] = "alt"
haptabm$calltype[haptabm$resist & haptabm$genotype > 0] = "res"

#reorder haps and samples
haptabm$posf <- factor(haptabm$POS,ordered=T)
haptabm$sample <- factor(haptabm$sample,levels=treeorder,ordered=T)
haptabm$region <- treelabels[haptabm$sample]




igtplot <- ggplot(haptabm,aes(x=sample,y=posf,fill=calltype)) + 
  geom_raster() + 
  scale_y_discrete(limits=rev)+
  scale_fill_manual(values=c("res"="red",
                            "ref"="lightgrey",
                            "alt"="darkgrey"),
                    labels=c("reference",
                             "alternate",
                             "resistant"),
                    name="") +
  theme(axis.text.x=element_blank(),axis.title=element_blank(),axis.ticks=element_blank())
igtplot
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png)


``` r
mhsorder <- metahap$sample[order(metahap$clusters,metahap$region)]
metahap$sample <- factor(metahap$sample,levels=mhsorder,ordered=T)

regcols <- c("East Africa"="#fff7bc",
             "West Africa"="#fec44f",
             "Caribean"="#DE2D26",
             "North America"="#FB6A4A",
             "South America"="#A50F15",
             "Middle East"="#3182BD",
             "Asia"="#6BAED6",
             "Pacific"="#08519C")


metaregplot <- ggplot(metahap,aes(x=sample,y=1,fill=region)) + 
  geom_raster() + scale_fill_manual(values=regcols) +
  theme(axis.text=element_blank(),axis.title=element_blank(),axis.ticks=element_blank())

metaclustplot <- ggplot(metahap,aes(x=sample,y=1,fill=factor(clusters%%2))) + 
  geom_raster() + scale_fill_manual(values=c("#444444","#BBBBBB"),guide=FALSE) +
  theme(axis.text=element_blank(),axis.title=element_blank(),axis.ticks=element_blank())



ly<-'AAAE
BBBE
CCCE
DDDE'

# treeplot + 
#   metaclustplot + 
#   metaregplot + 
#   igtplot + guide_area() + 
#   plot_layout(design=ly,guides="collect",
#                      heights=c(2,0.5,0.5,5),
#                      widths=c(9,2))


treeplot2 + 
  metaclustplot +
  metaregplot + 
  igtplot + guide_area() + 
  plot_layout(design=ly,guides="collect",
                     heights=c(2,0.5,0.5,5),
                     widths=c(9,2))
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-1.png)

``` r
ggsave("hclust_haplo_plot.png",width=300,height=200,dpi=400,units="mm")
ggsave("hclust_haplo_plot.svg",width=300,height=200,dpi=400,units="mm")
```



