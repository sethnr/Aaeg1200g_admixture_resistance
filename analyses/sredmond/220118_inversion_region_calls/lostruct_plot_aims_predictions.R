
library("tidyverse")
library("lostruct")

#library("patchwork")
library("gridExtra")
# library("grid")

args = commandArgs(trailingOnly=TRUE)

# setwd("~/Gits/Aaeg1000g_analyses/analyses/sredmond/220118_inversion_region_calls/")
# country <- "Senegal"
# chrom <- 1
# metafile <- "resources/meta_Aaeg1kg_spp.txt"
# invfile <- "inv_candidates_chr1_Senegal.txt"
# invaimsfile <- "inv_Senegal_n100_aim_snps.txt"
# meancallsfile <- "inv_Senegal_n100_mean_calls.txt"
# mlcallsfile <- "inv_Senegal_n100_nbayes_predictions.txt"
# outcallspng <- "inv_Senegal_n100_calls_dist.png"
# outsnpspng <- "inv_Senegal_n100_calls_snps.png"


country = args[1]
chrom = args[2]
metafile <- args[3]
invfile <- args[4]
invaimsfile <- args[5]
meancallsfile <- args[6]
mlcallsfile <- args[7]

outcallspng <- args[8]
outsnpspng <- args[9]


write("gathering meta",file=stderr())

metatab <- read.table(metafile,header=T, sep="\t")
metatab$contgroup <- factor(metatab$contgroup,levels=c("Wafrica","Safrica","Eafrica","Americas","Asia"),ordered=T)
metatab$region <- factor(metatab$region,levels=c("East Africa","West Africa","South America","Carribean","North America","Middle East","Asia","Pacific"),ordered=T)

countrysorttab <- unique(metatab[,c("country","contgroup","region")])
metatab$country <- factor(metatab$country,levels=unique(metatab$country[order(metatab$region)]),ordered=T)
samples <- metatab$sample


invcands <- read.table(invfile,header=T)
invcands <- subset(invcands,as.logical(valid))

allinvsnps <- read.table(invaimsfile,header=T,check.names = F)           

mlpreds <- read.table(mlcallsfile,header=F,col.names = unique(allinvsnps$inv))           
meanpreds <- read.table(meancallsfile,header=T)           

aimplots=list()
for(invname in unique(allinvsnps$inv)) {
    invnamesafe <- paste("X",gsub("\\D",".",invname,perl=T),sep="")
    invsnps <- subset(allinvsnps,inv==invname)
    
    predictions = cbind(meanpreds[,c("sample",invnamesafe)] %>% rename("mean" = invnamesafe),
          "ml"=mlpreds[,invnamesafe])
    predictions <- merge(predictions,metatab,by="sample")
    
    meansnp <- apply(invsnps[,samples],2,FUN=function(x) {mean(na.omit(x))})
    
    cntinvorder <- predictions$sample[order(predictions$contgroup,predictions$region,predictions$country,meansnp[predictions$sample])]
    cntorder <- unique(predictions$country[order(predictions$contgroup,predictions$region)])
    
    aimsM <- pivot_longer(invsnps,all_of(samples),names_to = "sample") %>% rename("invcountry"="country")
    aimsM <- merge(aimsM,metatab,by="sample")
    aimsM$sample <- factor(aimsM$sample,levels = cntinvorder,ordered=T)
    aimsM$country <- factor(aimsM$country,levels = cntorder,ordered=T)
    
    aimplot <- ggplot(aimsM,aes(x=i,y=as.numeric(sample),fill=as.factor(value))) + geom_raster() + 
       ylab("samples") + xlab("SNPs")+ theme(legend.position="none")+
       scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) +
      scale_fill_manual(values=c("0"="blue","1"="purple","2"="red")) +
      facet_grid("country ~ .",scale="free_y",space="free_y") +
      ggtitle(paste(invname,country,"aim calls","( top",nrow(invsnps),")"))+
      theme(panel.spacing = unit(0.2, "mm"),
            strip.background = element_blank(),
            strip.text = element_blank(),
            axis.text.y=element_blank(),
            axis.title.y=element_blank(),
            axis.ticks.y=element_blank())
    aimplot
    aimplots[[invname]] <- aimplot

    predsM <- pivot_longer(predictions,all_of(c("mean","ml")),names_to = "pred")
    predsM$sample <- factor(predsM$sample,levels = cntinvorder,ordered=T)
    predsM$country <- factor(predsM$country,levels = cntorder,ordered=T)
    
    predplot <- ggplot(predsM,aes(x=pred,y=as.numeric(sample),fill=as.factor(value))) + geom_raster() + 
      ylab("samples") + xlab("calls")+ theme(legend.position="none")+
      scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) +
      scale_fill_manual(values=c("0"="blue","1"="purple","2"="red")) +
      facet_grid("country ~ .",scale="free_y",space="free_y") +
      ggtitle(paste(invname,country,"aim calls","( top",nrow(invsnps),")"))+
      theme(panel.spacing = unit(0.2, "mm"),
            axis.text.y=element_blank(),
            axis.title.y=element_blank(),
            axis.ticks.y=element_blank())
    predplot
    
    aimplots[[paste(invname,"preds")]] <- predplot
    
    
}

do.call("grid.arrange", c(aimplots, nrow=1))

if(exists("allinvsnps")) {
  write("found invsnps",file=stderr())
  write.table(allinvsnps,outsnps,col.names=T,quote=F,row.names=F,sep="\t")
  plotwidths <- rep(c(8,2),length(aimplots)/2)
  
  png(outsnpspng,res=400,width=200,height=200,units='mm')
    do.call("grid.arrange", c(aimplots, nrow=1))
  dev.off()
  
  png(outcallspng,res=400,width=200,height=200,units='mm')
    do.call("grid.arrange", c(callplots, ncol=1))
  dev.off()
} else {
  write("writing null tables",file=stderr())
  
  file.create(outsnps)
  #file.create(outsnpspng)
  #file.create(outcallspng)
}
