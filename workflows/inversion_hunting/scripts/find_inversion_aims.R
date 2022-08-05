library("tidyverse")
library("lostruct")
library("getopt")
opttab <- matrix(c("inversions","i","1","character",
                   "assessment","a","1","character",
                   "vcf","v","1","character",
                   "pcs","p","1","character",
                   "country","c","1","character",
                   "samples","s","1","character",
                   "outfile","o","1","character"
),byrow=T,ncol=4)
opt <- getopt(opttab)

vcffile <- opt$vcf
country <- opt$country
metafile <- opt$samples
invfile <- opt$inversions
assessfile <- opt$assessment
pcsfile <- opt$pcs
outtxt <- opt$outfile


#AIM criteria
MAXCHISQ <- 1e-9
blocksize<-5e05
chromname <- c("NC_035107.1","NC_035108.1","NC_035109.1")
chromlen <- c(310827022,474425716,409777670)
names(chromlen) <- chromname


######
# helper functions
######
vcf_positions <- function (file, regions)
{
  bcf.sites <- data.table::fread(cmd = paste("bcftools query -f '%CHROM\\t%POS\\n'",
                                             file,"-r",region_string(regions)),
                                  header = FALSE, sep = "\t", data.table = FALSE)
  colnames(bcf.sites) <- c("chrom", "pos")
  bcf.sites
}

vcf_genotypes <- function (file, regions, samples) {

  write(paste("bcftools query -f '[ %GT]\\n'",
              "-r",regions,
              "-s",shQuote(paste(samples, collapse = ",")),
              file),stderr())
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

  invass <- read.table(assessfile,header=T)
  goodinvs <- invass$cluster[invass$valid & invass$lowdist]
  write(paste(" ",length(goodinvs),"valid candidates"),file=stderr())

  invcands <- subset(invcands,invcands$cluster %in% goodinvs)
  write(paste(" ",nrow(invcands),"valid blocks"),file=stderr())

  invcands$chromname <- chromname[invcands$chrom]
  invcands$end <- invcands$pos
  invcands$start <- invcands$end-blocksize


allaims <- data.frame(chrom=character(),
                      pos=numeric(),
                      i=numeric(),
                      inv=character(),
                      country=character(),
                      assoc=numeric())
######
# get all potential aims - any SNP associated with PCA in region
######

for(C in unique(invcands$cluster)) {
    write(paste("finding aims for cluster",C),file=stderr())
    #write.table(invcands[invcands$cluster==C,c("chromname","start","end")],stderr())
   #  regions <- paste(invcands$chromname[invcands$cluster==C],":",
   #                   invcands$start[invcands$cluster==C],"-",
   #                   invcands$end[invcands$cluster==C],
   #                   sep="",collapse=",")
    #get SNPs in inverted region

    invcands[invcands$cluster==C,c("chromname","start","end")]

    write("  get SNPs",stderr())
    invsnps <- vcf_query(vcffile,
                         samples=csamples,
                         regions=invcands[invcands$cluster==C,c("chromname","start","end")])

    write("  get posns",stderr())
    posns <- vcf_positions(vcffile,invcands[invcands$cluster==C,c("chromname","start","end")])

    write(paste(" ",nrow(posns),"SNPs in region"),stderr())

    #parse out inversion calls from PC file,
    invcall <- as.numeric(pcs$valid[pcs$inv==C])-1
    names(invcall) <- pcs$sample[pcs$inv==C]
    invcall <- invcall[csamples]
    invorder <- order(invcall)


   #get SNPs that can be used for chisq
    goodsnpsi <- apply(invsnps,1,function(x) {!any(is.na(x)) & length(unique(x))>1})
    goodsnps <- invsnps[goodsnpsi,]
    goodposns <- posns[goodsnpsi,]

   #do chisq with inversion calls for each
    #write(paste(length(invcall),dim(goodsnps)),file=stderr())
    assoc <- apply(goodsnps,1,FUN=function(x)(chisq.test(x,invcall)$p.value))
    goodposns$assoc <- assoc

   #get associated SNPS & posns
    write(paste(" ",sum(goodposns$assoc < MAXCHISQ)," potential AIMs over",MAXCHISQ,"for",C),file=stderr())
    if(sum(goodposns$assoc < MAXCHISQ)>0) {
      realgoodposns <- goodposns[goodposns$assoc < MAXCHISQ,] %>% add_column(i=c(1:sum(goodposns$assoc < MAXCHISQ)),
                                                                              "inv"=C,
                                                                              "country"=country,
                                                                              .after="pos")

      allaims <- rbind(allaims,realgoodposns)
      }
}



########
# refine aims based on LD Across all samples
########

for(C in unique(invcands$cluster)) {
    chr = invcands[invcands$cluster==C,"chrom"][1]
    chrname = chromname[chr]



   #get calculated aim posns for this inversion
    write(paste("refining aim",C),file=stderr())
    invaims <- allaims[allaims$inv==C,]
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
    write(length(modecall),stderr())
    write(paste(modecall,sep="",collapse="."),stderr())
    modecall <- as.numeric(modecall)
    write(paste(modecall,sep="",collapse="."),stderr())
    modecorr <- apply(invsnps,1,function(x) {if(sum(!is.na(x))>0) {cor(modecall[!is.na(x)],x[!is.na(x)])} else {0}})
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
  write.table(allinvsnps,outtxt,quote=F,row.names=F)
} else {
  file.create(outtxt)
}
