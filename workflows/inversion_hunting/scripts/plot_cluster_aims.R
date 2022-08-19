
library("tidyverse")
library("lostruct")

library("patchwork")
library("gridExtra")
library("grid")

library("getopt")

opttab <- matrix(c("infile","i","1","character",
                   "blocks","b","1","character",
                   "samples","s","1","character",
                   "outfile","o","1","character"
),byrow=T,ncol=4)
opt <- getopt(opttab)

invaimsfile <- opt$infile
blocksfile <- opt$blocks
metafile <- opt$samples
outsnpspng <- opt$outfile


write("gathering meta",file=stderr())

metatab <- read.table(metafile,header=T, sep="\t")
metatab$contgroup <- factor(metatab$contgroup,levels=c("Wafrica","Safrica","Eafrica","Americas","Asia"),ordered=T)
metatab$region <- factor(metatab$region,levels=c("East Africa","West Africa","South America","Carribean","North America","Middle East","Asia","Pacific"),ordered=T)

countrysorttab <- unique(metatab[,c("country","contgroup","region")])
metatab$country <- factor(metatab$country,levels=unique(metatab$country[order(metatab$region)]),ordered=T)

chromlen <- c(310827022,474425716,409777670)


write("reading inv snps",file=stderr())

if (file.size(invaimsfile)>0) {
  allinvsnps <- read.table(invaimsfile,header=T,check.names = F)
  write(paste(length(unique(allinvsnps$inv)),"invs found in",invaimsfile),file=stderr())
  allinvnames <- unique(allinvsnps$inv)
} else {
  allinvnames <- c()
  write(paste("no AIMs found to plot in ",invaimsfile,"\n",
              "writing empty plot for",outsnpspng),stderr())
  file.create(outsnpspng)
  q("no",0,F)
}

nsamples <- nrow(metatab)
metatab <- metatab[metatab$sample %in% colnames(allinvsnps),]
write(paste("found",nrow(metatab),"samples of",nsamples),stderr())
samples <- metatab$sample

#order countries in metatable by region (w african, eafrican, american, asian)
cntorder <- unique(metatab$country[order(metatab$contgroup,metatab$region)])

meansnp <- apply(allinvsnps[,samples],2,FUN=function(x) {mean(na.omit(x))})
write("got means",stderr())
samples <- samples[samples %in% colnames(allinvsnps)]


allinvsnps <- allinvsnps[,c("chrom","pos","i","id","inv","cluster",samples)]
allinvsnps <- unique(allinvsnps)
write(colnames(allinvsnps)[!colnames(allinvsnps) %in% samples],stderr())

#removing invs with too few AIMs
allinvsnps <- allinvsnps[allinvsnps$include,]
goodclusters <- allinvsnps$cluster

write(paste(goodclusters),stderr())
write(length(goodclusters),stderr())
if(length(goodclusters) < 1) {
    file.create(outsnpspng)
    q("no",0)
}

aimsM <- pivot_longer(allinvsnps,all_of(samples),
			names_to = "sample")   #%>% rename("invcountry"="country")
write("pivoted",stderr())

countries <- metatab$country
names(countries) <- metatab$sample


write("sorting SNPs",stderr())
#order SNPs by country, then mean call in this inversion
snporder <- metatab$sample[order(metatab$country,meansnp[metatab$sample])]


write("sorting countries",stderr())
    aimsM$country <- countries[aimsM$sample]
    aimsM$country <- factor(aimsM$country,levels = cntorder,ordered=T)
    aimsM$sample <- factor(aimsM$sample,levels = snporder,ordered=T)
    aimsM$cncode <- aimsM$country
    levels(aimsM$cncode) <- substr(levels(aimsM$country),0,3)


write("reading blocks",stderr())
    blocks <- read.table(blocksfile,header=T)
    blocks$y=1
    blocks <- blocks[blocks$cluster %in% goodclusters,]

write("calculating block extents",stderr())
    blockextents <- blocks %>% group_by(cluster) %>% summarise(min=min(pos),mid=mean(pos),max=max(pos),y=min(y))
    clustcols <- sample(rainbow(nrow(blockextents),s=0.6,v=0.9))
    names(clustcols) <- blockextents$cluster
    blockextents$size <- blockextents$max-blockextents$min

#calculate level to put each cluster on graph
write("calculating display levels",stderr())
    clustersL2S <- blockextents$cluster[order(blockextents$size,decreasing=T)]
    clustersL2R <- blockextents$cluster[order(blockextents$min,decreasing=F)]

    for(C1 in clustersL2S) {
      for(C2 in clustersL2S) {
        lenC1 = blockextents$size[blockextents$cluster==C1]
        minC1 = blockextents$min[blockextents$cluster==C1]
        maxC1 = blockextents$max[blockextents$cluster==C1]
        yC1   = min(blocks$y[blocks$cluster==C1])

        lenC2 = blockextents$size[blockextents$cluster==C2]
        minC2 = blockextents$min[blockextents$cluster==C2]
        maxC2 = blockextents$max[blockextents$cluster==C2]
        yC2   = min(blocks$y[blocks$cluster==C2])

        if((minC2<maxC1 & maxC2>minC1)) {
          if(lenC2<lenC1 & yC1==yC2) {
            newYC2 <- yC2+1
            blocks$y[blocks$cluster==C2] <- newYC2
            blockextents$y[blockextents$cluster==C2] <- newYC2
            write(paste(C1,":",yC1," ",C2,":",yC2,"->",yC2+1,sep=""),stderr())
          }
        }
      }
    }


write("plotting AIMs",stderr())
    aimsM$cluster <- factor(aimsM$cluster,levels=clustersL2R,ordered=T)
    aimplot <- ggplot(aimsM,aes(x=i,y=as.numeric(sample),fill=as.factor(value))) + geom_raster() +
       ylab("samples") + xlab("SNPs")+ theme(legend.position="none")+
         scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) +
        scale_fill_manual(values=c("0"="blue","1"="purple","2"="red")) +
        facet_grid("cncode ~ cluster",scale="free",space="free") +
        #ggtitle(paste(invname,"aim calls","( top",nrow(invsnps),")"))+
        theme(panel.spacing = unit(0.2, "mm"),
              axis.text.y=element_blank(),
              axis.title=element_blank(),
              axis.ticks.y=element_blank())




write("plotting blocks",stderr())
    chrom <- blocks$chrom[1]
    blockclustplot <- ggplot(blocks,aes(x=pos,y=y,fill=as.factor(cluster))) +
      geom_raster() +
      geom_segment(data=blockextents,aes(x=min,xend=max,y=y,yend=y,color=as.factor(cluster)),inherit.aes=F) +
      geom_text(data=blockextents,aes(x=mid,y=y,label=cluster),inherit.aes=F) +
      scale_fill_manual(values = clustcols)+
      scale_color_manual(values = clustcols)+
      theme(axis.title=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
            legend.position="none",legend.title.align = 1) +
      scale_x_continuous(limits=c(0,chromlen[chrom]),expand = c(0,0,0,0))

write("plotting aims / blocks",stderr())
    combplot <- arrangeGrob(
      aimplot, blockclustplot,
      ncol=1,heights=c(8,2))

write("saving aims / blocks",stderr())
    ggsave(outsnpspng,plot=combplot,dpi = 300,width=250,height=175,units="mm")
