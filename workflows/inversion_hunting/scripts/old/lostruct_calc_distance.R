# devtools::install_github("petrelharp/local_pca/lostruct")
library("lostruct")
library("tidyverse")

args = commandArgs(trailingOnly=TRUE)

inpcs <- args[1]
outdist <- args[2]
invs = args[3]
outpng <- args[4]

cores <- 4

intable <- read.table(inpcs,header=T)

allposns <- intable %>% select(chrom,start,end)
allpcs <- as.matrix(intable %>% select(-chrom,-start,-end))

write("calculating PC distances",stderr())
pcdist <- pc_dist(allpcs,npc=2,w=1,mc.cores=cores)


pcdistdf <- as.data.frame(pcdist)


blocks <- paste(allposns[,1],allposns[,2],sep=":")
allposns$block <- blocks
allposns$pos <- c(1:nrow(allposns))

colnames(pcdistdf) <- blocks
write("saving PC distances",stderr())
write.table(pcdistdf,file=outdist,col.names = T,quote=F,sep="\t")


pcdistdf$x <- blocks
pcdistflat <- pivot_longer(pcdistdf,cols=blocks,names_to = "y")
pcdistflat <- merge(merge(pcdistflat,allposns,by.x="x",by.y="block"),allposns,by.x="y",by.y="block",suffixes = c(".x",".y"))

#ggplot(pcdistflat,aes(x=x,y=y,fill=value)) + geom_raster() + coord_fixed()

ggplot(pcdistflat,aes(x=pos.x,y=pos.y,fill=value)) + geom_raster() + coord_fixed()
ggsave(outpng)
