# devtools::install_github("petrelharp/local_pca/lostruct")
library("lostruct")

args = commandArgs(trailingOnly=TRUE)

inpcs <- args[1]
outdist <- args[2]
outpng <- args[3]

intable <- read.table(inpcs,header=T)

allposns <- intable %>% select(chrom,start,end)
allpcs <- intable %>% select(-chrom,-start,-end)

pcdist <- pc_dist(allpcs,npc=2,w=1)


pcdistdf <- as.data.frame(pcdist)

blocks <- paste(allposns[,1],allposns[,2],sep=":")
allposns$block <- blocks
blocktab$pos <- as.numeric(blocktab$pos)

colnames(pcdistdf) <- blocks
pcdistdf$x <- blocks
pcdistflat <- pivot_longer(pcdistdf,cols=blocks,names_to = "y")
pcdistflat <- merge(merge(pcdistflat,blocktab,by.x="x",by.y="block"),blocktab,by.x="y",by.y="block",suffixes = c(".x",".y"))


ggplot(pcdistflat,aes(x=pos.x,y=pos.y,fill=value)) + geom_raster() + coord_fixed()
ggsave("test_snpable_AC3.png")
