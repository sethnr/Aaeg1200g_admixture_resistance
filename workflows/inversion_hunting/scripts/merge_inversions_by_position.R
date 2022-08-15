library("tidyverse")
library("gridExtra")
library("grid")

library("getopt")

opttab <- matrix(c("blocks","b","1","character",
                   "calls","c","1","character",
                   "out","o","1","character"
),byrow=T,ncol=4)
opt <- getopt(opttab)

blockfile <- opt$blocks
callsfile <- opt$calls
outprefix <- opt$out

# aimsfile <- "data/lostruct_aims/aims_chr1_Kenya.txt"
# invfile <- "data/lostruct_inversion_candidates/regions_chr1_Kenya.txt"
# assessfile <- "data/lostruct_inversions/inversions_chr1_Kenya.txt"
# metafile <- "resources/meta_Aaeg1kg_spp.txt"
# outprefix <- "merged_aims_chr1_Kenya"

outblocks <- paste(outprefix,"blocks.txt",sep="_")
outcalls <- paste(outprefix,"calls.txt",sep="_")
outmerges <- paste(outprefix,"merges.txt",sep="_")


if(file.size(callsfile)==0L) {
  file.create(outblocks)
  file.create(outcalls)
  file.create(outmerges)
  write(paste("no calls in file",callsfile,"\n","writing empty files for",outcalls,outblocks,outmerges),stderr())
  quit("no",0)
}


invblocks <- read.table(blockfile,header=T)
calls <- read.table(callsfile,header=T)
colnames(calls) <- gsub("X","",colnames(calls))

#get only inversions with calls (i.e. no NA columns)
calls <- calls[,apply(calls,2,function(x) {!all(is.na(x))})]

#inverse order inversions by block size
csizes = c()

invids <- colnames(calls)[2:ncol(calls)]
#invids <- unique(invblocks$cluster)
for(C1 in invids) {
  B1 <- invblocks$block[invblocks$cluster==C1]
  csizes <- c(csizes,length(B1))}
names(csizes) <- invids
invids <- invids[order(csizes,decreasing = T)]


#merge all inversions that are >= 10% overlap and have identical calls by PCA/kmeans
blocksim <- 0.1 #minimum 2-way block overlap for merge

rawnames <- invids
newblocks <- list()

mergetab <- data.frame("cluster"=character(),
                       "inversion"=character())
ci=0
while(length(rawnames)>0) {
  ci<-ci+1
  C1 <- rawnames[1]
  B1 <- invblocks$block[invblocks$cluster==C1]
  ols=c(C1)
  olbs <- B1
  for(C2 in rawnames) {
    if (C1==C2) {next}
    B2 <- invblocks$block[invblocks$cluster==C2]
    if (length(B1)<length(B2)) {next}

    B1ol <- sum(B1 %in% B2)/length(B1)
    B2ol <- sum(B2 %in% B1)/length(B2)

    C1calls <- as.character(calls[,as.character(C1)])
    C2calls <- as.character(calls[,as.character(C2)])

    flipcall <- c('aa'='bb','ab'='ab','bb'='aa')
    callsmatch <- all(C1calls==C2calls) | all(flipcall[C1calls]==C2calls)

    if(B1ol>=blocksim & B2ol>=blocksim & callsmatch) {
      #write(paste(C1,"<-",C1,C2,length(B1),length(B2)),stderr())
      #write(paste(C1,"<-",C1,C2),stderr())
      ols <- c(ols,C2)
      olbs <- c(olbs,B2)
    }
  }
  rawnames <- rawnames[!rawnames %in% ols]
  write(paste(" ",C1,"<-",paste(ols,collapse="/")),stderr())

  #clustername <- paste(ols,collapse="/")
  clustername <- paste(ci,length(ols),sep="_")
  newblocks[[clustername]] = unique(olbs)
  mergetab <- rbind(mergetab,data.frame("cluster"=rep(clustername,length(ols)),
                         "inversion"=ols))

  calls[,clustername] <- calls[,ols[1]]
}

write.table(mergetab,outmerges,col.names=T,quote=F,row.names=F,sep="\t")

#for(invname in names(newblocks)) {
#  write(paste(" ",invname),file=stderr())
#  invnamesafe <- paste("X",gsub("\\D",".",invname,perl=T),sep="")
#  compinvids <- as.numeric(strsplit(invname,"/")[[1]])
#  #get SNPs for inversion, remove duplicates, re-index
#  #aims$inv[aims$inv %in% compinvids] <- invname
#}

newcalls <- calls[,c("sample",names(newblocks))]
write.table(newcalls,outcalls,col.names=T,quote=F,row.names=F,sep="\t")


#write blocks file
block = c()
inv = c()
for(invname in names(newblocks)) {
  block <- c(block,newblocks[[invname]])
  inv <- c(inv,rep(invname,length(newblocks[[invname]])))

}
chrom <- as.numeric(as.data.frame(strsplit(block,":"))[1,])
posn <- as.numeric(as.data.frame(strsplit(block,":"))[2,])
write.table(data.frame("inv"=inv,"block"=block,"chrom"=chrom,pos=posn),
            outblocks,col.names=T,quote=F,row.names=F,sep="\t")
