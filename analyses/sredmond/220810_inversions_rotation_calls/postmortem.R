library("tidyverse")

library("patchwork")
library("gridExtra")

library("zoo")
library("plyr")
library("getopt")



library("getopt")

opttab <- matrix(c("infile","i","1","character",
                   "dir","d","1","character",
                   "chr","c","1","character",
                   "out","o","1","character",
                   "vcf","v","1","character"
),byrow=T,ncol=4)
opt <- getopt(opttab)

indir <- opt$dir
ldmerges <- opt$infile
chrom <- opt$chr
#vcffile <- opt$vcf
outprefix <- opt$out


ldmergetab <- read.table(ldmerges,header=T) %>% dplyr::select(!"cluster") %>% dplyr::rename("poscluster"="inversion",
                                                              "ldcluster"="newname")

write(paste("getting merges from ",indir),stderr())
mergefiles <- list.files(indir,pattern=paste("invs_chr",chrom,"_.*merges.txt",sep=""),full.names = T)
mergefiles <- mergefiles[laply(mergefiles, file.size)>0]
posmergetab <- ldply(mergefiles, read.table,header=T) %>% 
                        dplyr::rename("poscluster"="cluster","hclust"="inversion")
posmergetab$country <- gsub("_.*","",gsub("^\\d_","",posmergetab$poscluster,perl=T))


assfiles <- list.files(indir,pattern=paste("assess_chr",chrom,"_.*txt",sep=""),full.names = T)
assfiles <- assfiles[!grepl("calls",assfiles)]
assfiles <- assfiles[laply(assfiles, file.size)>0]
names(assfiles) <- gsub(".txt","",gsub(".*chr._","",assfiles))
asstab <- ldply(assfiles, function(x) {read.table(x,header=T) %>% add_column(country=names(x))}) %>% 
  dplyr::rename("country"=".id","hclust"="cluster")


head(ldmergetab)
head(posmergetab)
head(asstab)

mergetab <- merge(ldmergetab,merge(posmergetab,asstab,by=c("hclust","country"),all=T),by="poscluster",all=T)



ggplot(mergetab,aes(x=d,y=mean_wss,color=admixed)) + 
  geom_point(color="grey") + 
  geom_point(data=subset(mergetab,hzgood | admixed))

ggplot(mergetab,aes(x=mean_wss,y=mean_bss,color=admixed)) + 
  geom_point(color="grey") + 
  geom_point(data=subset(mergetab,hzgood | admixed))


ggplot(mergetab,aes(x=hzaa,y=hzab,color=admixed)) + 
  geom_point(color="grey") + 
  geom_point(data=subset(mergetab,hzgood | admixed),alpha=0.5) +
  coord_fixed()


mergetab$hzabp <- mergetab$hzab/apply(mergetab[,c("hzaa","hzbb")],1,FUN=max)
ggplot(mergetab,aes(x=f3,y=hzabp,color=admixed)) + 
  geom_point(color="grey") + 
  geom_point(data=subset(mergetab,hzgood | admixed),alpha=0.5)
