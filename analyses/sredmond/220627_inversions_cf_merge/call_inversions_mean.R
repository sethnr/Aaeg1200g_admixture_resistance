
library("tidyverse")

library("patchwork")
library("gridExtra")
library("grid")

library("genetics")
library("getopt")

opttab <- matrix(c("infile","i","1","character",
                   "meta","m","1","character",
                   "outfile","o","1","character"
),byrow=T,ncol=4)
opt <- getopt(opttab)

invaimsfile <- opt$infile
metafile <- opt$meta
outprefix <- opt$outfile

# invaimsfile <- "invs_ldcf_chr1_r0.25_aims.txt"
# metafile <- "resources/meta_Aaeg1kg_spp.txt"


write("gathering meta",file=stderr())

metatab <- read.table(metafile,header=T, sep="\t")
metatab$contgroup <- factor(metatab$contgroup,levels=c("Wafrica","Safrica","Eafrica","Americas","Asia"),ordered=T)
metatab$region <- factor(metatab$region,levels=c("East Africa","West Africa","South America","Carribean","North America","Middle East","Asia","Pacific"),ordered=T)

countrysorttab <- unique(metatab[,c("country","contgroup","region")])
metatab$country <- factor(metatab$country,levels=unique(metatab$country[order(metatab$region)]),ordered=T)


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

invnames <- sort(unique(allinvsnps$inv))

callmatrix <- matrix(rep(NA,length(samples)*length(invnames)),
                     ncol = length(invnames),
                     dimnames = list(samples,invnames) )

for(i in c(1:length(invnames))) {
  I <- invnames[i]
  meansnpcalls <- apply(allinvsnps[allinvsnps$inv==I,samples],2,FUN=function(x) {as.integer(median(na.omit(x)))})
  callmatrix[samples,i] <- meansnpcalls
}


calldf <- as.data.frame(callmatrix)
calldf$sample <- row.names(calldf)



calldf <- merge(calldf,metatab,by="sample")


getCallFreqs <- function(x) {
  ab = sum(x==1)
  aa = sum(x==0)
  bb = sum(x==2)
  n = (ab+aa+bb)
  maf = (ab+(2*bb))/(2*n)
  list("aa"=aa,
       "ab"=ab,
       "bb"=bb,
       "n"=n,
       "maf"=maf)
}



freqtable <- data.frame(inv=character(),
           country=character(),
           pop=character(),
           aa=numeric(),
           ab=numeric(),
           bb=numeric(),
           n=numeric(),
           maf=numeric(),
           P=numeric(),
           HWE=logical())

for(I in as.character(invnames)) {
  #write(I,stderr())
  for(C in unique(calldf$country)) {
  #for(C in c("Kenya")) {
    #write(paste(" ",C),stderr())
    calls <- calldf[calldf$country==C,I]
    cfreqs <- getCallFreqs(calls)
    geno <- genotype(c("a/a","a/b","b/b")[calls+1])
    if(nallele(geno)==2) {
      test <- HWE.test(geno)
      Pval <- test$test$p.value
    } else{
      Pval <- 1
    }
    
    i = nrow(freqtable)+1
    freqtable[i,c("inv","country","pop")] <- c(I,C,"all")
    freqtable[i,c("aa","ab","bb","n","maf","P","HWE")] <- c(cfreqs[["aa"]],
                                               cfreqs[["ab"]],
                                               cfreqs[["bb"]],
                                               cfreqs[["n"]],
                                               cfreqs[["maf"]],
                                               Pval,
                                               Pval > 0.01)
    
    for (P in unique(calldf$pop[calldf$country==C])) {
      #write(paste(C,P),stderr())
      calls <- calldf[calldf$country==C & calldf$pop==P,I]
      cfreqs <- getCallFreqs(calls)
      geno <- genotype(c("a/a","a/b","b/b")[calls+1])
      if(nallele(geno)==2) {
        test <- HWE.test(geno)
        Pval <- test$test$p.value
      } else{
        Pval <- 1
      }
      
      i = nrow(freqtable)+1
      freqtable[i,c("inv","country","pop")] <- c(I,C,P)
      freqtable[i,c("aa","ab","bb","n","maf","P","HWE")] <- c(cfreqs[["aa"]],
                                                              cfreqs[["ab"]],
                                                              cfreqs[["bb"]],
                                                              cfreqs[["n"]],
                                                              cfreqs[["maf"]],
                                                              Pval,
                                                              Pval > 0.01)
      
    }
  }
}

write.table(subset(freqtable,HWE==0),stderr())
write.table(freqtable,paste(outprefix,"inv_freqs.txt",sep="_"),sep="\t",quote=F,row.names=F,col.names=T)

write.table(calldf,paste(outprefix,"inv_calls.txt",sep="_"),
            sep="\t",quote=F,row.names=F,col.names=T)





