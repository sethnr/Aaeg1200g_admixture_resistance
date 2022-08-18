
library("tidyverse")

library("patchwork")
library("gridExtra")
library("grid")

library("getopt")

opttab <- matrix(c("infile","i","1","character",
                   "meta","m","1","character",
                   "outfile","o","1","character"
),byrow=T,ncol=4)
opt <- getopt(opttab)

invaimsfile <- opt$infile
metafile <- opt$meta
outprefix <- opt$outfile

callsfile <- paste(outprefix,"inv_calls.txt",sep="_")

write("gathering meta",file=stderr())
metatab <- read.table(metafile,header=T, sep="\t")


write("reading inv snps",file=stderr())
if (file.size(invaimsfile)>0) {
  allinvsnps <- read.table(invaimsfile,header=T,check.names = F)
  write.table(table(allinvsnps$inv),file=stderr(),row.names = F,quote=F,col.names = F)
  write(paste(length(unique(allinvsnps$inv)),"invs found in",invaimsfile),file=stderr())
  allinvnames <- unique(allinvsnps$inv)
} else {
  allinvnames <- c()
  write(paste("no AIMs found in ",invaimsfile,"\n",
              "writing empty files for",callsfile),stderr())
  file.create(callsfile)

  q("no",0,F)

}

nsamples <- nrow(metatab)
metatab <- metatab[metatab$sample %in% colnames(allinvsnps),]

write(paste("found",nrow(metatab),"samples of",nsamples),stderr())
samples <- metatab$sample



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

write.table(calldf,callsfile,sep="\t",quote=F,row.names=F,col.names=T)
