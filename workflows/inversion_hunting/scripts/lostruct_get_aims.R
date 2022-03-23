library("tidyverse")
library("lostruct")

args = commandArgs(trailingOnly=TRUE)

# setwd("~/Gits/Aaeg1000g_analyses/analyses/sredmond/220118_inversion_region_calls/")
# vcffile <- "/Volumes/Mosquito_raw_data/Aedes/Aaeg1000g/thinrand/lostruct_chr1.vcf.gz"
# country <- "Kenya"
# metafile <- "resources/meta_Aaeg1kg_spp.txt"
# invfile <- "inv_candidates_chr1_Kenya.txt"
# pcsfile <- "inv_candidates_chr1_Kenya.Rds"
# outtxt <- "./AIM_candidates_chr2_Uganda.txt"


vcffile <- args[1]
country = args[2]
metafile <- args[3]
invfile <- args[4]
pcsfile <- args[5]
outtxt <- args[6]


#AIM criteria
MAXCHISQ <- 1e-9


######
# helper functions
######
vcf_positions <- function (file, regions)
{
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


write("gathering meta",file=stderr())
  metatab <- read.table(metafile,header=T, sep="\t")
  metatab$contgroup <- factor(metatab$contgroup,levels=c("Wafrica","Eafrica",
                      "Americas","Asia"),ordered=T)
  metatab$region <- factor(metatab$region,levels=c("East Africa","West Africa",
                      "South America","Carribean","North America",
                      "Middle East","Asia","Pacific"),ordered=T)

  countrysorttab <- unique(metatab[,c("country","contgroup","region")])
  metatab$country <- factor(metatab$country,levels=unique(metatab$country[order(metatab$region)]),ordered=T)
  samples <- metatab$sample

  csamples <- metatab$sample[metatab$country==country]
  write(paste("found",length(csamples),"samples for",country),file=stderr())


write("reading PCs from lostruct analysis",file=stderr())
  pcs <- readRDS(pcsfile)

write("loading candidate regions",file=stderr())
  invcands <- read.table(invfile,header=T)
  write(paste(" ",nrow(invcands),"inversion candidates"),file=stderr())
  invcands <- subset(invcands,valid)
  write(paste(" ",nrow(invcands),"valid inversions"),file=stderr())


allaims <- data.frame(chrom=character(),
                      pos=numeric(),
                      i=numeric(),
                      inv=character(),
                      country=character(),
                      assoc=numeric())
######
# get all potential aims - any SNP associated with PCA in region
######

for(i in as.character(row.names(invcands))) {
    chr = invcands[i,"chrom"]
    chrname = invcands[i,"chromname"]
    st = invcands[i,"start"]
    en = invcands[i,"end"]
    invname = invcands[i,"name"]
    write(paste("finding aims for",invname),file=stderr())

   #get SNPs in inverted region
    invsnps <- vcf_query(vcffile,
                         samples=csamples,
                         regions=invcands[i,c("chromname","start","end")])
    posns <- vcf_positions(vcffile,paste(chrname,":",st,"-",en,sep=""))

   #parse out inversion calls from PC file,
    invcall <- as.numeric(pcs$valid[pcs$inv==invname])-1
    names(invcall) <- pcs$sample[pcs$inv==invname]
    invcall <- invcall[csamples]
    invorder <- order(invcall)


   #get SNPs that can be used for chisq
    goodsnpsi <- apply(invsnps,1,function(x) {!any(is.na(x)) & length(unique(x))>1})
    goodsnps <- invsnps[goodsnpsi,]
    goodposns <- posns[goodsnpsi,]

   #do chisq with inversion calls for each
    write(paste(length(invcall),dim(goodsnps)),file=stderr())
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



########
# refine aims based on LD Across all samples
########

for(i in as.character(row.names(invcands))) {
    chr = invcands[i,"chrom"]
    chrname = invcands[i,"chromname"]
    invname = invcands[i,"name"]


   #get calculated aim posns for this inversion
    write(paste("refining aim",invname),file=stderr())
    invaims <- allaims[allaims$inv==invname,]
    write(paste(" ",nrow(invaims),"aims found"),file=stderr())

   # if more than <maxaims> aims, only take top <maxaims> sorted by P-value
    maxaims <- 100
    if(nrow(invaims) > maxaims) {
      invaims <- invaims[order(invaims$assoc)[1:maxaims],]
      invaims <- invaims[order(invaims$i),]
    }
    if(nrow(invaims)==0) {next}


   #pull out only those SNPs from file for ALL samples
    invsnps <- vcf_query(vcffile,
                         regions=data.frame("chrom"=invaims$chrom,
                                            "start"=invaims$pos,
                                            "end"=invaims$pos),
                         samples=samples)
    colnames(invsnps) <- samples

   #if inversely correlated with modal value, flip call
    modecall <- apply(invsnps,2,function(x) {as.numeric(names(sort(table(na.omit(x)),decreasing = T))[1])})
    modecorr <- apply(invsnps,1,function(x) {cor(modecall[!is.na(x)],x[!is.na(x)])})
    invsnps[modecorr<0,] <- abs(invsnps[modecorr<0,]-2)


   #remove SNPs not in LD across whole dataset (mean -1x sd)
    write(paste("  LD filtering",nrow(invsnps),"SNPs"),file=stderr())
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
    write(paste("  -->",nrow(invsnps)),file=stderr())



   #order all SNPs by country, then mean inv call of high LD SNPs
    meancall <- apply(invsnps,2,function(x) {mean(na.omit(x))})
    cntinvorder <- metatab$sample[order(metatab$contgroup,metatab$country,meancall)]
    invaims$qual <- mean(abs(modecorr))

    invsnps <- cbind(invaims,invsnps)
    if(!exists("allinvsnps")) {allinvsnps <- invsnps} else {allinvsnps <- rbind(allinvsnps,invsnps)}


}


if(exists("allinvsnps")) {
  write.table(allinvsnps,outtxt)
} else {
  create.file(outtxt)
}
