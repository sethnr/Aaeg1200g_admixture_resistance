
library("tidyverse")

#library("patchwork")
library("gridExtra")
library("grid")

library("getopt")

opttab <- matrix(c("inversions","i","1","character",
                   "assessment","a","1","character",
                   "vcf","v","1","character",
                   "aims","A","1","character",
                   "country","c","1","character",
                   "meta","m","1","character",
                   "out","o","1","character"
),byrow=T,ncol=4)
opt <- getopt(opttab)

vcffile <- opt$vcf
country <- opt$country
metafile <- opt$meta
invfile <- opt$inversions
assessfile <- opt$assessment
aimsfile <- opt$aims
outprefix <- opt$out

# aimsfile <- "data/lostruct_aims/aims_chr1_Kenya.txt"
# invfile <- "data/lostruct_inversion_candidates/regions_chr1_Kenya.txt"
# assessfile <- "data/lostruct_inversions/inversions_chr1_Kenya.txt"
# metafile <- "resources/meta_Aaeg1kg_spp.txt"
# outprefix <- "merged_aims_chr1_Kenya"

outaims <- paste(outprefix,"aims.txt",sep="_")
outblocks <- paste(outprefix,"blocks.txt",sep="_")
#outpng <- paste(outprefix,"png",sep=".")

#blocksize<-5e05
# chromname <- c("NC_035107.1","NC_035108.1","NC_035109.1")
# chromlen <- c(310827022,474425716,409777670)
# names(chromlen) <- chromname

write(file.size(aimsfile),stderr())

if(file.size(aimsfile)==0L) {
  file.create(outaims)
  file.create(outblocks)
  write(paste("no aims in file",aimsfile,"\n","writing empty files for",outaims,outblocks),stderr())
  quit("no",0)
}


aims <- read.table(aimsfile,header=T)

metatab <- read.table(metafile,header=T, sep="\t")
samples <- metatab$sample
samples <- samples[samples %in% colnames(aims)]


#samples <- metatab$sample[metatab$country==country]
#write(paste("found",length(samples),"samples for",country),stderr())


invcands <- read.table(invfile,header=T)
# invcands$end <- as.numeric((as.data.frame(strsplit(invcands$block,":"))[2,]))
# invcands$start <- invcands$end-blocksize
# invcands$chromname <- chromname[invcands$chrom]

#get inversions that pass PCA assessment:
invass <- read.table(assessfile,header=T)
invids <- invass$cluster[invass$valid]



#inverse order inversions by block size
csizes = c()
for(C1 in invids) {
  B1 <- invcands$block[invcands$cluster==C1]
  csizes <- c(csizes,length(B1))}
names(csizes) <- invids
invids <- invids[order(csizes,decreasing = T)]


blocksim <- 0.75 #minimum 2-way block overlap for merge

rawnames <- invids
newblocks <- list()
while(length(rawnames)>0) {
  C1 <- rawnames[1]
  B1 <- invcands$block[invcands$cluster==C1]
  ols=c(C1)
  olbs <- B1
  for(C2 in rawnames) {
    if (C1==C2) {next}
    B2 <- invcands$block[invcands$cluster==C2]
    if (length(B1)<length(B2)) {next}

    B1ol <- sum(B1 %in% B2)/length(B1)
    B2ol <- sum(B2 %in% B1)/length(B2)

    write(paste(B1ol,sum(B1 %in% B2),length(B1),
                B2ol,sum(B2 %in% B1),length(B2),
                blocksim),stderr())
                
    if(B1ol>=blocksim & B2ol>=blocksim) {
      #write(paste(C1,"<-",C1,C2,length(B1),length(B2)),stderr())
      #write(paste(C1,"<-",C1,C2),stderr())
      ols <- c(ols,C2)
      olbs <- c(olbs,B2)
    }
  }
  rawnames <- rawnames[!rawnames %in% ols]
  write(paste(" ",C1,"<-",paste(ols,collapse="/")),stderr())
  #newnames <- c(newnames,paste(ols,collapse="/"))

  newblocks[[paste(ols,collapse="/")]] = unique(olbs)
  #write(paste(allBlocks),stderr())

}
#length(newblocks)




for(invname in names(newblocks)) {
  write(paste(" ",invname),file=stderr())
  invnamesafe <- paste("X",gsub("\\D",".",invname,perl=T),sep="")
  compinvids <- as.numeric(strsplit(invname,"/")[[1]])
  #get SNPs for inversion, remove duplicates, re-index
  aims$inv[aims$inv %in% compinvids] <- invname
}

aims <- aims[!duplicated(aims[c("chrom","pos","country","inv")]),]
aims <- aims[order(aims$pos),]
for(invname in names(newblocks)) {
  aims$i[aims$inv %in% compinvids] <- c(1:sum(aims$inv %in% compinvids))
}
write.table(aims,outaims,col.names=T,quote=F,row.names=F,sep="\t")


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
