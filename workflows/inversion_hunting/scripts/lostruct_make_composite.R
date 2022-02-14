
library("tidyverse")

args = commandArgs(trailingOnly=TRUE)

distsf <- args[1]
region <- args[2]
chrom <- as.numeric(args[3])
sppfile <- args[4]
invfile <- args[5]
outpng <- args[6]
outinvs <- args[7]

# distsf <- "redmond-lab-aaeg1000g/results_lostruct/lostruct_all_chr/"
# region <- "wafrica"
# sppfile <- "resources/meta_Aaeg1kg_spp.txt"
# outpng <- "~/desktop/wafrica_composite_test.png"

write(paste("dists folder: ",distsf),file=stderr())
write(paste("region: ",region),file=stderr())
write(paste("sppfile: ",sppfile),file=stderr())
write(paste("png: ",outpng),file=stderr())



chromname <- c("NC_035107.1","NC_035108.1","NC_035109.1")
chromlen <- c(310827022,474425716,409777670)
names(chromlen) <- chromname

pcblocksize <- 5e5

spptab <- read.table(sppfile,header=T, sep="\t")

countries <- unique(spptab$country[tolower(spptab$contgroup)==region])

distfolderlist <- list.files(distsf,full.names = T)
alldistfiles <- vector()

#for(country in c(region,countries)) {
for(country in c(countries)) {
    cdist <- distfolderlist[grep(paste("chr",chrom,"_",country,".txt",sep=""),distfolderlist)]
  if(length(cdist)==1) {
    write(paste(country,"->",cdist[1]),file=stderr())
    alldistfiles[country] <- cdist[1]
  } else {
    write(paste("found ",length(cdist),"files for ",country),file=stderr())
  }
}

rm("allpcdists")

for(country in names(alldistfiles)) {
  write(paste("parsing",country,"::",alldistfiles[country]),file=stderr())
  regionfile <- alldistfiles[country]
  regiondf <- read.table(regionfile,header=T)
  pcdistdf <- regiondf[,c(1:dim(regiondf)[1])]

  pos <- pcblocksize * 1:dim(regiondf)[1]
  blocks <- paste(chrom,pos,sep=":")

  regions <- data.frame("pos"=pos,
                      "block"=blocks)
  colnames(pcdistdf) <- blocks

  pcdistdf$x <- blocks
  pcdistflat <- pivot_longer(pcdistdf,cols=all_of(blocks),names_to = "y")
  pcdistflat <- merge(merge(pcdistflat,regions,by.x="x",by.y="block"),regions,by.x="y",by.y="block",suffixes = c(".x",".y"))
  pcdistflat$set <- country

  if(exists("allpcdists")) {
    allpcdists <- rbind(allpcdists,pcdistflat)
  } else {
    allpcdists <- pcdistflat
  }
}

nrows = round(sqrt(length(alldistfiles)))
axlab <- as_mapper(~ .x /1e6)


ggplot(allpcdists,aes(x=pos.x,y=pos.y,fill=value)) + geom_raster() + coord_fixed() +
  scale_x_continuous(expand = c(0,0),labels=axlab) + scale_y_continuous(expand = c(0,0),labels=axlab) +
  facet_wrap("set",nrow=nrows) + theme(legend.position="none")

ggsave(outpng)


#######
# get inversion calls
#######

invcalls <- read.table(invfile,header=T,sep="\t")
colnames(invcalls) <- tolower(colnames(invcalls))
invcalls <- invcalls[invcalls$chrom==chrom,]


ggplot(allpcdists,aes(x=pos.x,y=pos.y,fill=value)) + geom_raster() + coord_fixed() +
  geom_rect(aes(xmin=start,xmax=end,ymin=start,ymax=end),
          data=invcalls,
          inherit.aes=F,fill=NA,color="orange") +
  scale_x_continuous(expand = c(0,0),labels=axlab) + scale_y_continuous(expand = c(0,0),labels=axlab) +
  facet_wrap("set",nrow=nrows) + theme(legend.position="none")
ggsave(outinvs)

