library("tidyverse")
library("lostruct")
library("getopt")
opttab <- matrix(c("blocks","b","1","character",
                   "calls", "c","1","character",
                   "vcf",   "v","1","character",
                   "chisq","P","1","numeric",
                   "maxaims","N","1","numeric",
                   "outfile","o","1","character"
),byrow=T,ncol=4)
opt <- getopt(opttab)

blockfile <- opt$blocks
callsfile <- opt$calls
vcffile <- opt$vcf
outtxt <- opt$outfile

#AIM criteria
maxaims <- opt$maxaims
MAXCHISQ <- as.numeric(opt$chisq)

blocksize<-5e05


if(file.size(callsfile)==0L) {
  file.create(outtxt)
  write(paste("no calls in file",callsfile,"\n","writing empty files for",outtxt),stderr())
  quit("no",0)
}

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
  txtgenos <- data.table::fread(cmd = paste("bcftools query -f '[ %GT]\\n'",
                                             "-r",region_string(regions),
                                             "-s",shQuote(paste(samples, collapse = ",")),
                                             file),
                               header = FALSE, sep = "\t", data.table = FALSE)
  txtgenos
}



write("loading inversion blocks",file=stderr())
  invblocks <- read.table(blockfile,header=T)
  write(paste(" ",length(unique(invblocks$inv)),"inversion candidates"),file=stderr())

  invblocks$chromname <- chromname[invblocks$chrom]
  invblocks$end <- invblocks$pos
  invblocks$start <- invblocks$end-blocksize


write("reading calls from lostruct merge",file=stderr())
calls <- read.table(callsfile,header=T,stringsAsFactors=T)
colnames(calls) <- gsub("X","",colnames(calls))

samples <- calls$sample

allaims <- data.frame(chrom=character(),
                      pos=numeric(),
                      i=numeric(),
                      inv=character(),
#                      country=character(),
                      assoc=numeric())
######
# get all potential aims - any SNP associated with PCA in region
######

for(C in as.character(unique(invblocks$inv))) {
    write(paste("finding aims for cluster",C),file=stderr())

    #get SNPs in inverted region
    write("  get SNPs",stderr())
    invsnps <- vcf_query(vcffile,
                         samples=samples,
                         regions=invblocks[invblocks$inv==C,c("chromname","start","end")])
    write("  get posns",stderr())
    posns <- vcf_positions(vcffile,invblocks[invblocks$inv==C,c("chromname","start","end")])
    write(paste(" ",nrow(posns),"SNPs in region"),stderr())

   #get SNPs that can be used for chisq (none missing)
    compsnpsi <- apply(invsnps,1,function(x) {!any(is.na(x)) & length(unique(x))>1})
    compsnps <- invsnps[compsnpsi,]
    compposns <- posns[compsnpsi,]

    #parse out inversion calls from PC file,
     invcall <- as.numeric(calls[,C])-1
     names(invcall) <- calls$sample

   #do chisq with inversion calls for each
    #write(paste(length(invcall),dim(compsnps)),file=stderr())
    assoc <- apply(compsnps,1,FUN=function(x)(chisq.test(x,invcall)$p.value))
    compposns$assoc <- assoc

   #get associated SNPS & posns
    write(paste(" ",sum(compposns$assoc < MAXCHISQ)," potential AIMs under",MAXCHISQ,"for",C,"(min",min(compposns$assoc),")"),file=stderr())
    if(sum(compposns$assoc < MAXCHISQ)>0) {
      assocposns <- compposns[compposns$assoc < MAXCHISQ,] %>% add_column(i=c(1:sum(compposns$assoc < MAXCHISQ)),
                                                                              "inv"=C,
                                                                              .after="pos")
   #if more than [maxaims] posns, take top 100 by chisq p-value
    if(nrow(assocposns)>maxaims) {
        assocposns <- assocposns[order(assocposns$assoc)[1:maxaims],]
        assocposns <- assocposns[order(assocposns$i),]
    }
    allaims <- rbind(allaims,assocposns)
    }
}




########
# refine aims based on LD Across all samples
########

for(C in unique(invblocks$inv)) {
    chr = invblocks[invblocks$inv==C,"chrom"][1]
    chrname = chromname[chr]

   #get calculated aim posns for this inversion
    write(paste("refining aim",C),file=stderr())
    invaims <- allaims[allaims$inv==C,]
    write(paste(" ",nrow(invaims),"aims found"),file=stderr())

    if(nrow(invaims)==0) {
        next
    } else if (nrow(invaims)==1) {
        #if only one aim, no LD filtering, just pull out AIM and add to set
        invsnps <- vcf_query(vcffile,
                         regions=data.frame("chrom"=invaims$chrom,
                                            "start"=invaims$pos,
                                            "end"=invaims$pos),
                         samples=samples)

	colnames(invsnps) <- samples
        invsnps <- cbind(invaims,as.data.frame(invsnps))

    } else {
       #pull out only those SNPs from file for ALL samples
        invsnps <- vcf_query(vcffile,
                             regions=data.frame("chrom"=invaims$chrom,
                                                "start"=invaims$pos,
                                                "end"=invaims$pos),
                             samples=samples)
        colnames(invsnps) <- samples

       #if inversely correlated with modal value, flip call
        modecall <- apply(invsnps,2,function(x) {as.numeric(names(sort(table(na.omit(x)),decreasing = T))[1])})
        modecall <- as.numeric(modecall)
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

        meancall <- apply(invsnps,2,function(x) {mean(na.omit(x))})
        #order all SNPs by country, then mean inv call of high LD SNPs
        #cntinvorder <- metatab$sample[order(metatab$contgroup,metatab$country,meancall)]

        invsnps <- cbind(invaims,invsnps)
        }


    if(!exists("allinvsnps")) {
        allinvsnps <- invsnps
    } else {
        allinvsnps <- rbind(allinvsnps,invsnps)}

}


if(exists("allinvsnps")) {
  write.table(allinvsnps,outtxt,quote=F,row.names=F)
} else {
  file.create(outtxt)
}
