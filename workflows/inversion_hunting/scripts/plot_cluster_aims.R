
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

# invcands <- read.table(invfile,header=T)
# invcands <- subset(invcands,as.logical(valid))

write("reading inv snps",file=stderr())

if (file.size(invaimsfile)>0) {
  allinvsnps <- read.table(invaimsfile,header=T,check.names = F)
  write.table(table(allinvsnps$inv),file=stderr(),row.names = F,quote=F,col.names = F)
  write(paste(length(unique(allinvsnps$inv)),"invs found in",invaimsfile),file=stderr())
  allinvnames <- unique(allinvsnps$inv)
} else {
  allinvnames <- c()
}

nsamples <- nrow(metatab)
metatab <- metatab[metatab$sample %in% colnames(allinvsnps),]
write(paste("found",nrow(metatab),"samples of",nsamples),stderr())
samples <- metatab$sample
#write(samples,stderr())

#order countries in metatable by region (w african, eafrican, american, asian)
cntorder <- unique(metatab$country[order(metatab$contgroup,metatab$region)])


# mlpreds <- read.table(mlcallsfile,header=F,col.names = unique(allinvsnps$inv))
# meanpreds <- read.table(meancallsfile,header=T)

meansnp <- apply(allinvsnps[,samples],2,FUN=function(x) {mean(na.omit(x))})

write("got means",stderr())
write(dim(allinvsnps),stderr())
write(length(samples),stderr())
samples <- samples[samples %in% colnames(allinvsnps)]

allinvsnps <- allinvsnps[,c("chrom","pos","i","id","inv",samples)]
allinvsnps <- unique(allinvsnps)
write(colnames(allinvsnps)[!colnames(allinvsnps) %in% samples],stderr())

aimsM <- pivot_longer(allinvsnps,all_of(samples),
			names_to = "sample")   #%>% rename("invcountry"="country")
write("pivoted",stderr())

write.table(head(aimsM),stderr())

#aimsM$sample <- factor(aimsM$sample,levels=samples)
#metatab$sample <- factor(metatab$sample,levels=samples)


countries <- metatab$country
names(countries) <- metatab$sample
#aimsM <- merge(aimsM,metatab,by="sample")
#write("merged",stderr())


write("sorting countries",stderr())
aimsM$country <- countries[aimsM$sample]
aimsM$country <- factor(aimsM$country,levels = cntorder,ordered=T)

write("sorting SNPs",stderr())
#order SNPs by country, then mean call in this inversion
snporder <- metatab$sample[order(metatab$country,meansnp[metatab$sample])]

    #write(snporder,stderr())
aimsM$sample <- factor(aimsM$sample,levels = snporder,ordered=T)
aimsM$cncode <- aimsM$country
levels(aimsM$cncode) <- substr(levels(aimsM$country),0,3)

#write.table(head(aimsM),stderr())
#write(summary(aimsM),stderr())

write("making AIM plot",stderr())
aimplot <- ggplot(aimsM,aes(x=i,y=as.numeric(sample),fill=as.factor(value))) + geom_raster() +
   ylab("samples") + xlab("SNPs")+ theme(legend.position="none")+
     scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) +
    scale_fill_manual(values=c("0"="blue","1"="purple","2"="red")) +
    facet_grid("cncode ~ inv",scale="free",space="free") +
    #ggtitle(paste(invname,"aim calls","( top",nrow(invsnps),")"))+
    theme(panel.spacing = unit(0.2, "mm"),
          axis.text.y=element_blank(),
          axis.title=element_blank(),
          axis.ticks.y=element_blank())
aimplot

write("reading blocks",stderr())
blocks <- read.table(blocksfile,header=T)

aimcounts <- aimsM %>% group_by(inv) %>% summarise(aims=n_distinct(pos))
blocks <- blocks[blocks$cluster %in% aimcounts$inv,]

# csizeorder <- blocks %>% group_by(cluster) %>% summarise("size"=n_distinct(block)) %>% arrange(desc(size))
# csizeorder <- csizeorder$cluster

write("calculating display levels",stderr())
#calculate level to put each cluster on graph
blocks$y=1
for(C1 in unique(sort(blocks$cluster))) {
  for(C2 in unique(sort(blocks$cluster))) {
    lenC1 = sum(blocks$cluster==C1)
    minC1 = min(blocks$pos[blocks$cluster==C1])
    maxC1 = max(blocks$pos[blocks$cluster==C1])
    yC1   = min(blocks$y[blocks$cluster==C1])

    lenC2 = sum(blocks$cluster==C2)
    minC2 = min(blocks$pos[blocks$cluster==C2])
    maxC2 = max(blocks$pos[blocks$cluster==C2])
    yC2   = min(blocks$y[blocks$cluster==C2])


    if((minC2<maxC1 & maxC2>minC1)) {
      if(C2>C1 & yC1==yC2) {
        blocks$y[blocks$cluster==C2] <- blocks$y[blocks$cluster==C2]+1
        #write(paste(C1,":",yC1," ",C2,":",yC2,"->",yC2+1,sep=""),stderr())
      }
    }
  }
}


write("calculating block extents",stderr())
blockextents <- blocks %>% group_by(cluster) %>% summarise(min=min(pos),mid=mean(pos),max=max(pos),y=min(y))


clustcols <- sample(rainbow(nrow(blockextents),s=0.6,v=0.9))
names(clustcols) <- blockextents$cluster

chrom <- blocks$chrom[1]

write("plotting blocks",stderr())

blockclustplot <- ggplot(blocks,aes(x=pos,y=y,fill=as.factor(cluster))) +
  geom_raster() +
  geom_segment(data=blockextents,aes(x=min,xend=max,y=y,yend=y,color=as.factor(cluster)),inherit.aes=F) +
  geom_text(data=blockextents,aes(x=mid,y=y,label=cluster),inherit.aes=F) +
  scale_fill_manual(values = clustcols)+
  scale_color_manual(values = clustcols)+
  theme(axis.title=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        legend.position="none",legend.title.align = 1) +
  scale_x_continuous(limits=c(0,chromlen[chrom]),expand = c(0,0,0,0))

blockheight <- 0.5+(max(blocks$y)*0.25)

write("plotting aims / blocks",stderr())

#aimplot / blockclustplot + plot_layout(heights=c(10-blockheight,blockheight))

combplot <- arrangeGrob(
  aimplot, blockclustplot,
  ncol=1,heights=c(8,2))

write("saving aims / blocks",stderr())

ggsave(outsnpspng,plot=combplot,dpi = 300,width=250,height=175,units="mm")
