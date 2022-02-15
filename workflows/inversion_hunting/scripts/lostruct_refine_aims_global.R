
library("tidyverse")
library("lostruct")

#library("patchwork")
# library("gridExtra")
# library("grid")

args = commandArgs(trailingOnly=TRUE)

# setwd("~/Gits/Aaeg1000g_analyses/analyses/sredmond/220118_inversion_region_calls/")
# vcffile <- "/Volumes/Mosquito_raw_data/Aedes/Aaeg1000g/thinrand/lostruct_chr1.vcf.gz"
# country <- "Senegal"
# chrom <- 1
# metafile <- "resources/meta_Aaeg1kg_spp.txt"
# invfile <- "inv_candidates_chr1_Senegal.txt"
# aimsfile <- "./AIM_candidates_chr1_Senegal.txt"


vcffile <- args[1]
country = args[2]
chrom = args[3]
metafile <- args[4]
invfile <- args[5]
aimsfile <- args[6]
outcalls <- args[7]
outcallspng <- args[8]
outsnps <- args[9]
outsnpspng <- args[10]

plotaims <- F


#AIM criteria
MAXCHISQ <- 1e-9

vcf_positions <- function (file, regions) {
  bcf.sites <- data.table::fread(cmd = paste("bcftools query -f '%CHROM\\t%POS\\n'", 
                                             file,"-r",regions), 
                                  header = FALSE, sep = "\t", data.table = FALSE)
  colnames(bcf.sites) <- c("chrom", "pos")
  bcf.sites
}


vcf_genotypes <- function (file, regions, samples) {
  txtgenos <- data.table::fread(cmd = paste("bcftools query -f '[ %GT]\\n'", 
                                             "-r",region_string(regions),
                                             "-s",shQuote(paste(samples, collapse = ",")),
                                             file), 
                               header = FALSE, sep = "\t", data.table = FALSE)
  txtgenos
}

chromname <- c("NC_035107.1","NC_035108.1","NC_035109.1")
chromlen <- c(310827022,474425716,409777670)
names(chromlen) <- chromname

metatab <- read.table(metafile,header=T, sep="\t")
metatab$contgroup <- factor(metatab$contgroup,levels=c("Wafrica","Eafrica","Americas","Asia"),ordered=T)
metatab$region <- factor(metatab$region,levels=c("East Africa","West Africa","South America","Carribean","North America","Middle East","Asia","Pacific"),ordered=T)


countrysorttab <- unique(metatab[,c("country","contgroup","region")])

metatab$country <- factor(metatab$country,levels=unique(metatab$country[order(metatab$region)]),ordered=T)


samples <- metatab$sample

invcands <- read.table(invfile,header=T)
invcands <- subset(invcands,as.logical(valid))


allaims <- read.table(aimsfile)           

invcalls <- data.frame("sample"=samples,"country"=country)
if(exists("allinvsnps")) {rm("allinvsnps")}
aimplots <- list()
callplots <- list()

for(i in as.character(row.names(invcands))) {
    chr = invcands[i,"chrom"]
    chrname = invcands[i,"chromname"]
    
    st = invcands[i,"start"]
    en = invcands[i,"end"]
    invname = invcands[i,"name"]
    
    invaims <- allaims[allaims$inv==invname,]
    # if more than <maxaims> aims, only take top <maxaims> sorted by P-value
    maxaims <- 100
    if(nrow(invaims) > maxaims) {
      invaims <- invaims[order(invaims$assoc)[1:maxaims],]
      invaims <- invaims[order(invaims$i),]
    }
    invaims$i <- c(1:nrow(invaims))
    
    invsnps <- vcf_query(vcffile,
                         regions=data.frame("chrom"=invaims$chrom,
                                            "start"=invaims$pos,
                                            "end"=invaims$pos),
                         samples=samples)
    
    
    colnames(invsnps) <- samples
    
    
    nocall_rate <- round(apply(invsnps,1,function(x) {sum(is.na(x))/length(x)}),2)
    
    #if inversely correlated with modal value, flip call
    modecall <- apply(invsnps,2,function(x) {as.numeric(names(sort(table(na.omit(x)),decreasing = T))[1])})
    modecorr <- apply(invsnps,1,function(x) {cor(modecall[!is.na(x)],x[!is.na(x)])})
    invsnps[modecorr<0,] <- abs(invsnps[modecorr<0,]-2)
    
    #remove SNPs not in LD across whole dataset (mean -1x sd)
    #calculate mean r2 for each SNP
    r2s <- matrix(rep(-1,nrow(invsnps)^2),nrow=nrow(invsnps))
    for(si in c(1:nrow(invsnps))){
      for(sj in c(1:nrow(invsnps))){
        r2s[si,sj] = cor(invsnps[si,],invsnps[sj,],use = "pairwise.complete.obs")^2
      }
    }
    meanr2s <- apply(r2s,1,mean)
    ldinclude <- meanr2s >= (mean(r2s)-sd(r2s))
    invsnps <- invsnps[ldinclude,]
    invaims <- invaims[ldinclude,]
    invaims$i <- c(1:nrow(invaims))
    
    #order all SNPs by country, then mean inv call of high LD SNPs
    meancall <- apply(invsnps,2,function(x) {mean(na.omit(x))})
    cntinvorder <- metatab$sample[order(metatab$contgroup,metatab$country,meancall)]
    invaims$qual <- mean(abs(modecorr))

    
    invcall = rep(NA,length(meancall))
    invcall[meancall < 0.25] <- 0
    invcall[meancall < 1.25 & meancall > 0.75 ] <- 1
    invcall[meancall > 1.75 ] <- 2
    invcalls[[invname]] <- invcall
    
    invsnps <- cbind(invaims,invsnps)
    if(!exists("allinvsnps")) {allinvsnps <- invsnps} else {allinvsnps <- rbind(allinvsnps,invsnps)}
    
    
    aimsM <- pivot_longer(invsnps,all_of(samples),names_to = "sample") %>% rename("invcountry"="country")
    aimsM <- merge(aimsM,metatab,by="sample")
    aimsM$sample <- factor(aimsM$sample,levels = cntinvorder,ordered=T)
    
    aimplot <- ggplot(aimsM,aes(x=i,y=as.numeric(sample),fill=as.factor(value))) + geom_raster() + 
       ylab("samples") + xlab("SNPs")+ theme(legend.position="none")+
       scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) +
      scale_fill_manual(values=c("0"="blue","1"="purple","2"="red")) +
      facet_grid("country ~ .",scale="free_y",space="free_y") +
      ggtitle(paste(invname,country,"aim calls","( top",nrow(invaims),")"))+
      theme(panel.spacing = unit(0.2, "mm"),
            axis.text.y=element_blank(),
            axis.title.y=element_blank(),
            axis.ticks.y=element_blank())
    aimplots[[invname]] <- aimplot
    
      # ggsave(paste("inv_",gsub(":","_",invname),"_",country,"_n",maxaims,"_aims_LD_pruned.png",sep=""),
      #        aimplot,
      #        height=250,width=200,units="mm")
      # 
      
      
    callplot <- ggplot(data.frame("sample"=names(meancall),
                      "mean"=meancall,
                      "call"=invcall),aes(x=mean,group=call,fill=as.factor(call))) + 
      scale_fill_manual(values=c("0"="blue","1"="purple","2"="red"),na.value = "grey",guide="none") +
      geom_histogram(bins=100)+
      ggtitle(invname)
    callplots[[invname]] <- callplot
       # ggsave(paste("inv_",gsub(":","_",invname),"_",country,"_n",maxaims,"_calls_LD_pruned.png",sep=""),
       #       callplot,
       #       height=250,width=200,units="mm")
    
}

# write.table(invcalls,paste("inv_",chrom,"_",country,"_n",maxaims,"_mean_calls.txt",sep=""),col.names=T,quote=F,row.names=F,sep="\t")
# write.table(allinvsnps,paste("inv_",chrom,"_",country,"_n",maxaims,"_aim_snps.txt",sep=""),col.names=T,quote=F,row.names=F,sep="\t")
write.table(invcalls,outcalls,col.names=T,quote=F,row.names=F,sep="\t")

if(!exists("allinvsnps")) {
  write.table(allinvsnps,outsnps,col.names=T,quote=F,row.names=F,sep="\t")

  png(outsnpspng,res=400,width=200,height=200,units='mm')
    do.call("grid.arrange", c(aimplots, nrow=1))
  dev.off()
  
  png(outcallspng,res=400,width=200,height=200,units='mm')
    do.call("grid.arrange", c(callplots, ncol=1))
  dev.off()
} else {
  file.create(allinvsnps)
  file.create(outsnpspng)
  file.create(outcallspng)
}
