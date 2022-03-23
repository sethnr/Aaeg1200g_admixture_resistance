library("tidyverse")
library("lostruct")

setwd("~/Gits/Aaeg1000g_analyses/analyses/sredmond/220118_inversion_region_calls/")
vcffile <- "/Volumes/Mosquito_raw_data/Aedes/Aaeg1000g/thinrand/lostruct_chr1.vcf.gz"
country <- "Senegal"
metafile <- "resources/meta_Aaeg1kg_spp.txt"
invfile <- "inv_candidates_chr1_Senegal.txt"
pcsfile <- "lostruct_chr1_Senegal_pcs.Rds"
outtxt <- "./AIM_candidates_chr1_Senegal.txt"


#AIM criteria
MAXCHISQ <- 1e-9

vcf_positions <- function (file, regions) 
{
  bcf.sites <- data.table::fread(cmd = paste("bcftools query -f '%CHROM\\t%POS\\n'", 
                                             file,"-r",regions), 
                                  header = FALSE, sep = "\t", data.table = FALSE)
  colnames(bcf.sites) <- c("chrom", "pos")
  bcf.sites
}


chromname <- c("NC_035107.1","NC_035108.1","NC_035109.1")
chromlen <- c(310827022,474425716,409777670)
names(chromlen) <- chromname

metatab <- read.table(metafile,header=T, sep="\t")
samples <- metatab$sample[metatab$country==country]
write(paste("found",length(samples),"samples for",country),file=stderr())

pcs <- readRDS(pcsfile)

invcands <- read.table(invfile,header=T)
invcands <- subset(invcands,as.logical(valid))


allaims <- data.frame(chrom=character(),
                      pos=numeric(),
                      i=numeric(),
                      inv=character(),
                      country=character(),
                      assoc=numeric())
                      


for(i in as.character(row.names(invcands))) {
    chr = invcands[i,"chrom"]
    chrname = invcands[i,"chromname"]
    
    st = invcands[i,"start"]
    en = invcands[i,"end"]
    invname = invcands[i,"name"]
    invsnps <- vcf_query(vcffile,
                         samples=samples,
                         regions=invcands[i,c("chromname","start","end")])
    posns <- vcf_positions(vcffile,paste(chrname,":",st,"-",en,sep=""))
    
    
    
    #parse out inversion calls from PC file, 
    invcall <- as.numeric(pcs$valid[pcs$inv==invname])-1
    names(invcall) <- pcs$sample[pcs$inv==invname]
    invcall <- invcall[samples]
    invorder <- order(invcall)
    
    
    #get SNPs that can be used for chisq
    goodsnpsi <- apply(invsnps,1,function(x) {!any(is.na(x)) & length(unique(x))>1})
    goodsnps <- invsnps[goodsnpsi,]
    goodposns <- posns[goodsnpsi,]
    
    #do chisq with inversion calls for each
    assoc <- apply(goodsnps,1,FUN=function(x)(chisq.test(x,invcall)$p.value))
    goodposns$assoc <- assoc
    
    #get associated SNPS & posns
    write(paste(sum(goodposns$assoc < MAXCHISQ)," potential AIMs over",MAXCHISQ,"for",invname),file=stderr())
    if(sum(goodposns$assoc < MAXCHISQ)>0) {
      realgoodposns <- goodposns[goodposns$assoc < MAXCHISQ,] %>% add_column(i=c(1:sum(goodposns$assoc < MAXCHISQ)),
                                                                              "inv"=invname,
                                                                              "country"=country,
                                                                              .after="pos")
    
      allaims <- rbind(allaims,realgoodposns)
    }
}

write.table(allaims,outtxt)


for(i in as.character(row.names(invcands))) {
  chr = invcands[i,"chrom"]
  chrname = invcands[i,"chromname"]
  
  st = invcands[i,"start"]
  en = invcands[i,"end"]
  invname = invcands[i,"name"]
  
  write(paste(chr,st,en,invname),file=stderr())
  
  invsnps <- vcf_query(vcffile,
                       samples=samples,
                       regions=invcands[i,c("chromname","start","end")])
  posns <- vcf_positions(vcffile,paste(chrname,":",st,"-",en,sep=""))
  
  invcall <- as.numeric(pcs$valid[pcs$inv==invname])-1
  names(invcall) <- pcs$sample[pcs$inv==invname]
  invcall <- invcall[samples]
  invorder <- order(invcall)
  
  
  
  aimsinv <- subset(allaims,inv==invname)
  
  aimi <- posns$pos %in% aimsinv$pos
  aimsnps <- invsnps[aimi,]
  colnames(aimsnps) <- samples
  
  aims <- cbind(aimsinv,aimsnps)
  aimsM <- pivot_longer(aims,all_of(samples),names_to = "sample")
  
  
  
  
  aimsM$sample <- factor(aimsM$sample,levels=samples[invorder],ordered=T)
  aimsplot <- ggplot(aimsM,aes(x=i,y=as.numeric(sample),fill=as.factor(value))) + geom_raster() + 
    scale_fill_manual(values=c("0"="yellow","1"="orange","2"="red")) +
    ylab("sample") + xlab("SNPs")+ ggtitle(paste("AIMS",MAXCHISQ,country,invname)) + theme(legend.position="none")+
    scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))
  aimsplot
  png(gsub(":","-",paste("aims_",country,"_",invname,".png",sep="")),res=400,width=200,height=200,units='mm')
    print(aimsplot)
  dev.off()
}

