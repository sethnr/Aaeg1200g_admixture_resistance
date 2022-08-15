
library("tidyverse")
library("lostruct")
library("patchwork")
library("gridExtra")
library("grid")

library("getopt")

opttab <- matrix(c("dir","d","1","character",
                   "chr","c","1","character",
                   "out","o","1","character",
                   "vcf","v","1","character",
                   "meta","m","1","character",
                   "r2","r","2","numeric"
),byrow=T,ncol=4)
opt <- getopt(opttab)

indir <- opt$dir
chrom <- opt$chr
metafile <- opt$meta
vcffile <- opt$vcf
outprefix <- opt$out
minr2 <- 0.25 #minimum r2 for merging
blocksize<-5e05


if (!is.null(opt$r2)  ) {minr2 <- opt$r2}

# aimsfile <- "data/lostruct_aims/aims_chr1_Kenya.txt"
# invfile <- "data/lostruct_inversion_candidates/regions_chr1_Kenya.txt"
# assessfile <- "data/lostruct_inversions/inversions_chr1_Kenya.txt"
# metafile <- "resources/meta_Aaeg1kg_spp.txt"
# outprefix <- "merged_aims_chr1_Kenya"


#####
# get and cat blocks
#####

write(paste("getting positions from ",indir),stderr())
blockfiles <- list.files(indir,pattern = paste(".*chr",chrom,".*blocks.txt",sep=""),full.names = T)
countries <- gsub("invs_chr._","",gsub("_blocks.txt","",basename(blockfiles)))
names(blockfiles) <- countries

chromname <- c("NC_035107.1","NC_035108.1","NC_035109.1")
chromlen <- c(310827022,474425716,409777670)
names(chromlen) <- chromname

if(exists("allblocks")){rm(allblocks)}
for(C in countries) {
  if(file.size(blockfiles[C])>0) {
    cblocks = read.table(blockfiles[C],header=T)
    cblocks$region <- C
    if(exists("allblocks")) {
      allblocks <- rbind(allblocks,cblocks)
    } else{
      allblocks <- cblocks
    }
  }
}
allblocks$id <- paste(allblocks$chrom,allblocks$region,allblocks$inv,sep="_")

csizes <- as.data.frame(allblocks %>% group_by(id) %>% summarise("size"=n(),"mid"=mean(pos)))
cids <- csizes[order(csizes[,"mid"],decreasing = F),1]


#####
# get and cat aims
######

aimsfiles <- list.files(indir,pattern = paste(".*chr",chrom,".*aims.txt",sep=""),full.names = T)
countries <- gsub("invs_chr._","",gsub("_aims.txt","",basename(aimsfiles)))
names(aimsfiles) <- countries

write(aimsfiles,stderr())

if(exists("allaims")){rm(allaims)}
for(C in countries) {
  if(file.size(aimsfiles[C])>0) {
    write(paste(C,file.size(aimsfiles[C])),stderr())
    caims = read.table(aimsfiles[C],header=T)
    caims <- caims[,c("chrom","pos","i","inv","assoc")]
    caims$region <- C
    caims$id <- paste(chrom,C,caims$inv,sep="_")
    if(exists("allaims")) {
      allaims <- rbind(allaims,caims)
    } else{
      allaims <- caims
    }
  }
}
#samples <- colnames(allaims)[8:ncol(allaims)-1]

#allaims$id <- paste(match(allaims$chrom,chromname),allaims$region,allaims$inv,sep="_")
allaims$end <- allaims$pos
allaims$start <- allaims$pos

#write.table(head(allaims),stderr())
######
# get SNPs for all aims
######
#write("gathering meta",file=stderr())
metatab <- read.table(metafile,header=T, sep="\t")
samples <- metatab$sample


write("parsing all SNPs in inversion",file=stderr())
for(I in unique(allaims$id)) {
    write(paste("  ",I,sum(allaims$id==I)," SNPs"),file=stderr())
    invsnps <- vcf_query(vcffile,
                     samples=samples,
                     regions=allaims[allaims$id==I,c("chrom","start","end")])
    colnames(invsnps) <- samples
    if(!exists("allinvsnps")) {
        allinvsnps <- invsnps
    } else {
        allinvsnps <- rbind(allinvsnps,invsnps)
    }
}
write("binding all assoc SNPs",file=stderr())
write(dim(allaims),file=stderr())
write(dim(allinvsnps),file=stderr())
allaims <- cbind(allaims,allinvsnps)

samples <- samples[samples %in% colnames(allinvsnps)]


#######
# calc R2 for all inv pairs
#######
write("calculating linkage ",file=stderr())
if(file.exists(paste(outprefix,".r2.txt",sep=""))) {
  write(paste("reading from",paste(outprefix,".r2.txt",sep="")),stderr())
  r2df <- read.table(paste(outprefix,".r2.txt",sep=""),header=T,row.names=1)
  r2matrix <- as.matrix(r2df)
  colnames(r2matrix) <- row.names(r2df)
  colnames(r2df) <- row.names(r2df)
} else {

  r2matrix <- matrix(rep(NA,length(cids)^2),nrow=length(cids),dimnames=list(cids,cids))

  for (C1 in cids) {
    cat(paste("cf",C1),file=stderr())
    A1 <- as.matrix(allaims[allaims$id==C1,samples])
    if(nrow(A1)==0) {next}
    for (C2 in cids) {
      if(C2>=C1) {next}
      cat(".",file=stderr())
      A2 <- as.matrix(allaims[allaims$id==C2,samples])
      if(nrow(A2)==0) {next}

      r2s <- matrix(rep(-1,nrow(A1)*nrow(A2)),nrow=nrow(A1))
      for(si in c(1:nrow(A1))){
        for(sj in c(1:nrow(A2))){
          r2s[si,sj] = mean(cor(A1[si,],A2[sj,],use = "pairwise.complete.obs")^2)
        }
      }
      meanr2s <- mean(abs(r2s))
      #meanr2s <- mean(r2s)
      r2matrix[C1,C2] <- meanr2s
      r2matrix[C2,C1] <- meanr2s

    }
    write(paste(" "),stderr())
  }

  r2df <- as.data.frame(r2matrix)
  write.table(r2df,paste(outprefix,".r2.txt",sep=""),col.names=T,row.names=T,quote=F,sep="\t")
}


r2df$from <- rownames(r2matrix)
r2dfM <- pivot_longer(r2df,cols=any_of(cids),names_to = "to",values_to = "r2")
r2dfM$to <- factor(r2dfM$to,levels=cids,ordered=T)
r2dfM$from <- factor(r2dfM$from,levels=cids,ordered=T)
r2plot <- ggplot(r2dfM,aes(x=from,y=to,fill=r2)) + geom_raster() +
  coord_fixed() +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.position="bottom",axis.title=element_blank())



#merge inv calls based on LD of AIMs
rawnames <- colnames(r2matrix)
ldmerges <- list()

mergetab <- data.frame("cluster"=character(),
                       "inversion"=character(),
                       "r2"=character())

while(length(rawnames)>0) {
  C1 <- rawnames[1]
  ols=c(C1)
  for(C2 in rawnames) {
    if (C1==C2) {next}
    if (is.na(r2matrix[C1,C2])) {next}
    if(r2matrix[C1,C2] > minr2) {
      ols <- c(ols,C2)
      i <- nrow(mergetab)+1
      mergetab[i,c("cluster","inversion")] <- c(C1,C2)
      mergetab[i,"r2"] <- r2matrix[C1,C2]
    }
  }
  rawnames <- rawnames[!rawnames %in% ols]
  write(paste(" ",C1,"<-",paste(ols,collapse="/")),stderr())
  ldmerges[[C1]]<-ols
}


####
# write cluster IDs into file
#####

r2dfM$cluster = NA
allblocks$cluster = NA
allaims$cluster = NA
newnames<-c()

for(i in c(1:length(ldmerges))) {
  cname <- names(ldmerges)[i]
  ols <- ldmerges[[i]]

  allaims$cluster[allaims$id %in% ols] <- i
  allblocks$cluster[allblocks$id %in% ols] <- i

  #parse and create new cluster name
  olcountries = unique(gsub('[\\d\\_]*',"",ols,perl=T))
  olconts = unique(metatab$contgroup[metatab$country %in% olcountries)
  if(length(olcountries)==1) {
    cextent=olcountries[1]
  } else if (length(olconts)==1){
    cextent=olconts[1]
  } else {
    cextent="global"
  }
  ci <- length(grep(paste(chr,cextent,se="_"),newnames))+1
  newname <- paste(chr,cextent,ci,se="_")
  mergetab$newname[mergetab$cluster==cname] <- newname
  newnames <- c(newnames,newname)
}

write.table(allblocks,file=paste(outprefix,"blocks.txt",sep="_"),sep="\t",col.names=T,row.names=F,quote=F)

write.table(mergetab,file=paste(outprefix,"LDmerges.txt",sep="_"),col.names=T,quote=F,row.names=F,sep="\t")




clustersizes <- allblocks %>% group_by(cluster) %>%
                summarise(invs=n_distinct(inv),
                          blocks=n_distinct(block))

#bigclusts <- unique(na.omit(allblocks$cluster))
bigclusts <- clustersizes$cluster[clustersizes$invs>1]
clustcols <- sample(rainbow(length(bigclusts),s=0.6,v=0.9))
names(clustcols) <- bigclusts

allblocks$id <- factor(allblocks$id,levels=cids,ordered=T)
blockclustplot <- ggplot(allblocks,aes(x=pos,y=id,fill=as.factor(cluster))) +
  geom_raster() + scale_y_discrete(position="right") +
  theme(axis.title=element_blank(),legend.position="bottom",legend.title.align = 1) +
  scale_fill_manual(values = clustcols,name=(paste("cluster\nR^2>",minr2))) +
  xlim(0,chromlen[chrom])

r2plot | blockclustplot
ggsave(paste(outprefix,".png",sep=""),dpi = 300,width=400,height=220,units="mm")



#set aims for each cluster
clustaims <- allaims[!duplicated(allaims[,c("cluster","chrom","pos")]),] %>% select(!c(inv,assoc))
clustaims <- clustaims[order(clustaims$pos),] %>% rename("inv"="cluster")
for(i in unique(allaims$cluster)) {
  clustaims$i[clustaims$cluster==i] <- c(1:sum(clustaims$cluster==i))
}

write.table(clustaims,paste(outprefix,"aims.txt",sep="_"),sep="\t",quote=F,col.names=T,row.names=F)
