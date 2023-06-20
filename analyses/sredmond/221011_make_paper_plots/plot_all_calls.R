
library("tidyverse")
library("lostruct")

library("patchwork")
library("gridExtra")
library("grid")

library("plyr")
library("genetics")
# opttab <- matrix(c("infile","i","1","character",
#                    "blocks","b","1","character",
#                    "samples","s","1","character",
#                    "minaims","N","1","numeric",
#                    "outfile","o","1","character"
# ),byrow=T,ncol=4)
# opt <- getopt(opttab)
# 
# invaimsfile <- opt$infile
# blocksfile <- opt$blocks
# metafile <- opt$samples
# outsnpspng <- opt$outfile
# MINAIMS <- opt$minaims

invcalls <- "invs_chr.*_inv_calls.txt"
blocksfile <- "invs_chr.*_blocks.txt"
metafile <- "resources/meta_Aaeg1kg_spp.txt"
outprefix <- "all_inversion_calls"


write("gathering meta",file=stderr())

metatab <- read.table(metafile,header=T, sep="\t")
metatab$contgroup <- factor(metatab$contgroup,levels=c("Wafrica","Safrica","Eafrica","Americas","Asia"),ordered=T)
metatab$region <- factor(metatab$region,levels=c("East Africa","West Africa","South America","Carribean","North America","Middle East","Asia","Pacific"),ordered=T)

countrysorttab <- unique(metatab[,c("country","contgroup","region")])
metatab$country <- factor(metatab$country,levels=unique(metatab$country[order(metatab$region)]),ordered=T)

regorder <- c("East Africa","West Africa","Caribean","South America","North America","Middle East","Asia","Pacific")
metatab$region <- factor(metatab$region,levels=regorder,ordered=T)
cntorder <- unique(metatab$country[order(metatab$region)])

chromlen <- c(310827022,474425716,409777670)


#####
# read in all calls 
#####
callsfiles <- list.files("data",pattern=invcalls,full.names = T)
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

invnames <- colnames(calltab)[2:11]

#####
# calculate HWE
#####

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
hwetabP <- subset(hwetab,HWE<0.05)
hwetabP
#######
# get inv postions and arrange 
#######
#order countries in metatable by region (w african, eafrican, american, asian)

blocks <- ldply(as.list(list.files("data",blocksfile,full.names = T)),read_tsv)
blocks$y=0.5

blocks$inversion <- paste("X",blocks$cluster,sep="")
blocks <- blocks[blocks$inversion %in% unique(hwetab$inversion),]
blockextents <- blocks %>% dplyr::group_by(inversion) %>% dplyr::summarise(min=min(pos),mid=mean(pos),max=max(pos),chrom=min(chrom),y=0.5)
blockextents$size <- blockextents$max-blockextents$min

write("calculating display levels",stderr())
clustersL2S <- blockextents$inversion[order(blockextents$size,decreasing=T)]
clustersL2R <- blockextents$inversion[order(blockextents$chrom,blockextents$min,decreasing=F)]

clustersL2R

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

#swap names across
hwetab$country <- factor(hwetab$country,levels=cntorder,ordered=T)

hwetab$cncode <- hwetab$country
levels(hwetab$cncode) <- substr(levels(hwetab$country),0,3)

hwetab$inversion <- factor(hwetab$inversion,levels=clustersL2R,ordered=T)
hwetab <- merge(hwetab,unique(blocks[,c("inversion","chrom")]))
hwetabM <- pivot_longer(hwetab,cols=c("aa","ab","bb"),names_to = "call",values_to="count")


newnames <- c("1a","1b","1c","1d","1e","1f","1g","2a","2b","2c")
names(newnames) <- clustersL2R

levels(hwetabM$inversion) <- newnames
levels(hwetab$inversion) <- newnames


callplot <- ggplot(hwetabM,aes(x=inversion,y=count,fill=call)) + 
  geom_bar(stat="identity",width=0.99) + 
  geom_bar(data=subset(hwetab,HWE<0.01),aes(y=aa+ab+bb),stat="identity",fill=NA,color="cyan",size=1.5,width=0.99) + 
  facet_grid(cncode ~ chrom,scale="free",space="free") +
  scale_fill_manual(values=c("aa"="blue","ab"="purple","bb"="red")) +
  scale_x_discrete(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  theme(panel.spacing.x = unit(10, "mm"),
        panel.spacing.y = unit(0.2, "mm"),
        axis.text.y=element_blank(),
        axis.title=element_blank(),
        axis.ticks.y=element_blank(),
        plot.margin = unit(c(1,1,1,1), "mm"),
        strip.text.y.right = element_text(angle = 0))
callplot








chromlen <- data.frame("min"=0,
                       "chrom"=c(1,2),
                       "max"=c(310827022,474425716))

blocks$inversion <- newnames[blocks$inversion]
blockextents$inversion <- newnames[blockextents$inversion]
# clustcols <- sample(rainbow(nrow(blockextents),s=0.6,v=0.9))
# names(clustcols) <- blockextents$inversion
clustcols <- c("#CCCCCC","#555555")[0:length(blockextents$inversion) %% 2 + 1]
names(clustcols) <- newnames

blockclustplot <- ggplot(blocks,aes(x=pos,y=y,fill=as.factor(inversion))) +
  geom_raster() +
  geom_segment(data=chromlen,aes(x=min,xend=max,y=0,yend=0,color=NA),inherit.aes=F) +
  geom_segment(data=blockextents,aes(x=min,xend=max,y=y,yend=y,color=as.factor(inversion)),inherit.aes=F) +
  geom_text(data=blockextents,aes(x=mid,y=y+0.5,label=inversion),inherit.aes=F,vjust=-0.2) +
  scale_fill_manual(values = clustcols)+
  scale_color_manual(values = clustcols)+
  facet_grid(. ~ chrom,scale="free_x",space="free_x") +
  theme(axis.title=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        legend.position="none",legend.title.align = 1,panel.spacing.x = unit(10, "mm"),
        panel.grid.minor = element_blank(),panel.background = element_blank()) +
  scale_x_continuous(expand = c(0,0,0,0))+
  ylim(0,max(blockextents$y)+1)
blockclustplot

combplot <- callplot / blockclustplot + plot_layout(heights=c(9,1))
combplot

ggsave(paste(outprefix,".png",sep=""),combplot,width=12,height=7,dpi=400)
ggsave(paste(outprefix,".pdf",sep=""),combplot,width=12,height=7,dpi=400)


