
library("tidyverse")

library("patchwork")
library("gridExtra")
library("grid")

library("getopt")

opttab <- matrix(c("infile","i","1","character",
                   "meta","m","1","character",
                   "minaims","N","1","numeric",
                   "outfile","o","1","character"
),byrow=T,ncol=4)
opt <- getopt(opttab)

invaimsfile <- opt$infile
metafile <- opt$meta
outprefix <- opt$outfile
MINAIMS <- opt$minaims

callsfile <- paste(outprefix,"inv_calls.txt",sep="_")

write("gathering meta",file=stderr())
metatab <- read.table(metafile,header=T, sep="\t")


write("reading inv snps",file=stderr())
if (file.size(invaimsfile)>0) {
  allinvsnps <- read.table(invaimsfile,header=T,check.names = F)
  write.table(table(allinvsnps$cluster),file=stderr(),row.names = F,quote=F,col.names = F)
  write(paste(length(unique(allinvsnps$cluster)),"invs found in",invaimsfile),file=stderr())
  allinvnames <- unique(allinvsnps$cluster)
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


#removing invs with too few AIMs, stop if no clusters left
    write(paste("removing invs with <",MINAIMS,"aims"),stderr())
    aimcounts <- as.data.frame(table(allinvsnps$cluster))
    colnames(aimcounts) <- c("cluster","n")
    write.table(aimcounts,sep="\t",quote=F,stderr())
    goodclusters <- aimcounts$cluster[aimcounts$n>MINAIMS]
    if(length(goodclusters < 1)) {
        file.create(callsfile)
        q("no",0)
    }
    allinvsnps <- allinvsnps[allinvsnps$cluster %in% goodclusters,]

#invnames <- sort(unique(allinvsnps$cluster))
invnames <- goodclusters

callmatrix <- matrix(rep(NA,length(samples)*length(invnames)),
                     ncol = length(invnames),
                     dimnames = list(samples,invnames) )

for(i in c(1:length(invnames))) {
  I <- invnames[i]
  meansnpcalls <- apply(allinvsnps[allinvsnps$cluster==I,samples],2,FUN=function(x) {as.integer(median(na.omit(x)))})
  callmatrix[samples,i] <- meansnpcalls
}


calldf <- as.data.frame(callmatrix)
calldf$sample <- row.names(calldf)

calldf <- merge(calldf,metatab,by="sample")

write(paste("writing calls to",callsfile),stderr())
write.table(calldf,callsfile,sep="\t",quote=F,row.names=F,col.names=T)
