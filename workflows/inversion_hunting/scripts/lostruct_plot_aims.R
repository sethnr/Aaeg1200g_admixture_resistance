
library("tidyverse")
library("lostruct")

#library("patchwork")
library("gridExtra")
# library("grid")

args = commandArgs(trailingOnly=TRUE)

# setwd("~/Gits/Aaeg1000g_analyses/analyses/sredmond/220118_inversion_region_calls/")
# metafile <- "resources/meta_Aaeg1kg_spp.txt"
# invaimsfile <- "inv_Senegal_n100_aim_snps.txt"
# outsnpspng <- "inv_Senegal_n100_calls_snps.png"


metafile <- args[1]
invfile <- args[2]
invaimsfile <- args[3]
outsnpspng <- args[4]


write("gathering meta",file=stderr())

metatab <- read.table(metafile,header=T, sep="\t")
metatab$contgroup <- factor(metatab$contgroup,levels=c("Wafrica","Safrica","Eafrica","Americas","Asia"),ordered=T)
metatab$region <- factor(metatab$region,levels=c("East Africa","West Africa","South America","Carribean","North America","Middle East","Asia","Pacific"),ordered=T)

countrysorttab <- unique(metatab[,c("country","contgroup","region")])
metatab$country <- factor(metatab$country,levels=unique(metatab$country[order(metatab$region)]),ordered=T)
samples <- metatab$sample


# invcands <- read.table(invfile,header=T)
# invcands <- subset(invcands,as.logical(valid))

allinvsnps <- read.table(invaimsfile,header=T,check.names = F)           

# mlpreds <- read.table(mlcallsfile,header=F,col.names = unique(allinvsnps$inv))           
# meanpreds <- read.table(meancallsfile,header=T)           

aimplots=list()
for(invname in unique(allinvsnps$inv)) {
    invnamesafe <- paste("X",gsub("\\D",".",invname,perl=T),sep="")
    invsnps <- subset(allinvsnps,inv==invname)
    
    meansnp <- apply(invsnps[,samples],2,FUN=function(x) {mean(na.omit(x))})
    
    cntinvorder <- metatab$sample[order(metatab$contgroup,metatab$region,metatab$country,meansnp[metatab$sample])]
    cntorder <- unique(metatab$country[order(metatab$contgroup,metatab$region)])
    
    aimsM <- pivot_longer(invsnps,all_of(samples),names_to = "sample") %>% rename("invcountry"="country")
    aimsM <- merge(aimsM,metatab,by="sample")
    aimsM$sample <- factor(aimsM$sample,levels = cntinvorder,ordered=T)
    aimsM$country <- factor(aimsM$country,levels = cntorder,ordered=T)
    aimsM$cncode <- aimsM$country
    levels(aimsM$cncode) <- substr(levels(aimsM$country),0,3)

    aimplot <- ggplot(aimsM,aes(x=i,y=as.numeric(sample),fill=as.factor(value))) + geom_raster() + 
       ylab("samples") + xlab("SNPs")+ theme(legend.position="none")+
       scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) +
      scale_fill_manual(values=c("0"="blue","1"="purple","2"="red")) +
      facet_grid("cncode ~ .",scale="free_y",space="free_y") +
      ggtitle(paste(invname,"aim calls","( top",nrow(invsnps),")"))+
      theme(panel.spacing = unit(0.2, "mm"),
            axis.text.y=element_blank(),
            axis.title.y=element_blank(),
            axis.ticks.y=element_blank())
    aimplots[[invname]] <- aimplot

    
    
}

if(length(aimplots>0)) {
png(outcallspng,res=400,width=200,height=200,units='mm')
  do.call("grid.arrange", c(aimplots, nrow=1))
dev.off()
} else {file.create(outcallspng)}