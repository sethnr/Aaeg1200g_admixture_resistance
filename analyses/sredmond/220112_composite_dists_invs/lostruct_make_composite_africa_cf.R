
library("tidyverse")
library("patchwork")

args = commandArgs(trailingOnly=TRUE)

setwd("~/Gits/Aaeg1000g_analyses/analyses/sredmond/220112_composite_dists_invs/")
distsf <- "redmond-lab-aaeg1000g/results_lostruct/lostruct_all_chr/"
sppfile <- "resources/meta_Aaeg1kg_spp.txt"
invfile <- "resources/redmond_2020_inversion_calls.txt"


chromname <- c("NC_035107.1","NC_035108.1","NC_035109.1")
chromlen <- c(310827022,474425716,409777670)
names(chromlen) <- chromname

pcblocksize <- 5e5

spptab <- read.table(sppfile,header=T, sep="\t")


rm("allpcdists")
allplots = list();

for(cont in c("eafrica","wafrica","asia","americas")) {
  countries <- unique(spptab$country[tolower(spptab$contgroup)==cont])
  
  distfolderlist <- list.files(distsf,full.names = T)
  alldistfiles <- vector()
  for(chrom in c(1,2,3)) {
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
    
    for(country in names(alldistfiles)) {
      write(paste("parsing",country,"::",alldistfiles[country]),file=stderr())
      regionfile <- alldistfiles[country]
      regiondf <- read.table(regionfile,header=T)
      pcdistdf <- regiondf[,c(1:dim(regiondf)[1])]
    
      pos <- pcblocksize * 1:dim(regiondf)[1]
      blocks <- paste(chrom,pos,sep=":")
    
      regions <- data.frame(
                          "chrom"=rep(chrom,length(pos)),
                          "pos"=pos,
                          "block"=blocks)
      colnames(pcdistdf) <- blocks
    
      pcdistdf$x <- blocks
      pcdistflat <- pivot_longer(pcdistdf,cols=all_of(blocks),names_to = "y")
      pcdistflat <- merge(merge(pcdistflat,regions,by.x="x",by.y="block"),regions,by.x="y",by.y="block",suffixes = c(".x",".y"))
      pcdistflat$set <- country
      pcdistflat$continent <- cont
      
      if(exists("allpcdists")) {
        allpcdists <- rbind(allpcdists,pcdistflat)
      } else {
        allpcdists <- pcdistflat
      }
    }
  }
  
 }

for (cont in c("eafrica","wafrica","asia","americas")){

c1p <- ggplot(subset(allpcdists,continent==cont & chrom.x==1),aes(x=pos.x,y=pos.y,fill=value)) + geom_raster() + coord_fixed() +
  scale_x_continuous(expand = c(0,0),labels=axlab) + scale_y_continuous(expand = c(0,0),labels=axlab) +
  facet_grid("chrom.x ~ set") + theme(legend.position="none",axis.title = element_blank())
c2p <- ggplot(subset(allpcdists,continent==cont & chrom.x==2),aes(x=pos.x,y=pos.y,fill=value)) + geom_raster() + coord_fixed() +
  scale_x_continuous(expand = c(0,0),labels=axlab) + scale_y_continuous(expand = c(0,0),labels=axlab) +
  facet_grid("chrom.x ~ set") + theme(legend.position="none",axis.title = element_blank())
c3p <- ggplot(subset(allpcdists,continent==cont & chrom.x==3),aes(x=pos.x,y=pos.y,fill=value)) + geom_raster() + coord_fixed() +
  scale_x_continuous(expand = c(0,0),labels=axlab) + scale_y_continuous(expand = c(0,0),labels=axlab) +
  facet_grid("chrom.x ~ set") + theme(legend.position="none",axis.title = element_blank())

c1p / c2p / c3p
ggsave(paste("lostruct_allchrom_",cont,"_composite.png",sep=""),width = 250,height=250,units="mm")

}




