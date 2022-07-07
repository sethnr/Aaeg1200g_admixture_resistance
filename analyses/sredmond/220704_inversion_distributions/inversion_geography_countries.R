library("tidyverse")

library("patchwork")
library("gridExtra")
library("grid")

library("genetics")
library("getopt")
library("plyr")

# library(rgdal)
# library(rgeos)
library(maptools)
library(raster)
library(ggmap)


indir <- "."
outprefix <- "invs_ldcf"

write(paste("getting freqs from ",indir),stderr())
freqfiles <- list.files(indir,pattern = paste(".*inv_freqs.txt",sep=""),full.names = T)
chrs <- gsub("invs_ldcf_chr","",gsub("_inv_freqs.txt","",basename(freqfiles)))

names(freqfiles) <- chrs

#read in all inv freqs, add chr name as col
invfreqs <- ldply(chrs, function(C) {read.table(freqfiles[C],header=T) %>% add_column("chrom"=C)})

countryfile <- "country_geocoords.txt"
if(file.exists(countryfile)){
  locations <- read.table(countryfile,header=T) 
} else {
  register_google(key="AIzaSyBFbT5X6MjrepPAlzIj3zz1HS58hnIPtEM", write=T)
  countries <- unique(invfreqs$country)
  locations <- geocode(countries)
  locations$country <- countries
  write.table(locations,countryfile,sep="\t",row.names=F,col.names=T)
}

countryfreqs <- merge(subset(invfreqs,pop=="all"),locations)


world <- map_data("world")



for(C in chrs) {
  cinvs <- subset(countryfreqs,chrom==C)
  ncols <- floor(sqrt(length(unique(cinvs$inv))))
  ggplot() +
    geom_polygon(data = world, aes(x=long, y = lat, group=group), fill='grey',color="black",size=0.2) +
    geom_point(data=cinvs, aes(x=lon, y=lat, fill=maf,size=n), shape=21,inherit.aes=F) +
    #geom_text(data=sample_table, aes(x=lon, y=lat, fill=set,label=Population), shape=21,inherit.aes=F) +
    ggtitle(paste("inversion distributions c",C," invs",sep="")) +
    scale_fill_distiller(limits=c(0,1),palette="Reds",direction = 1)+
    coord_fixed(ylim=c(-48,48),xlim=c(-155,140))+
    theme(axis.text=element_blank(),
          panel.border=element_rect(fill=NA, color="black"),
          panel.background=element_rect(fill='light blue', color=NA),
          axis.title=element_blank(),
          legend.position="bottom") +
    facet_wrap(inv ~ .,ncol=ncols)

  
  ggsave(paste(outprefix,"_c",C,"_inv_dists.png",sep=""),
         dpi=300,width=350,height=200,units='mm')
}



#do subsets for chr1 based on correlations / positions

cinvs <- subset(countryfreqs,chrom==1 & inv %in% c("13","14","15","16"))
ncols <- floor(sqrt(length(unique(cinvs$inv))))
ggplot() +
  geom_polygon(data = world, aes(x=long, y = lat, group=group), fill='grey',color="black",size=0.2) +
  geom_point(data=cinvs, aes(x=lon, y=lat, fill=maf,size=n), shape=21,inherit.aes=F) +
  #geom_text(data=sample_table, aes(x=lon, y=lat, fill=set,label=Population), shape=21,inherit.aes=F) +
  ggtitle(paste("inversion distributions c1 Q cluster",sep="")) +
  scale_fill_distiller(limits=c(0,1),palette="Reds",direction = 1)+
  coord_fixed(ylim=c(-48,48),xlim=c(-155,140))+
  theme(axis.text=element_blank(),
        panel.border=element_rect(fill=NA, color="black"),
        panel.background=element_rect(fill='light blue', color=NA),
        axis.title=element_blank(),
        legend.position="bottom") +
  facet_wrap(inv ~ .,ncol=ncols)
ggsave(paste(outprefix,"_c",1,"_Qcluster_inv_dists.png",sep=""),
       dpi=300,width=350,height=200,units='mm')





cinvs <- subset(countryfreqs,chrom==1 & inv %in% c("4","5","6","7","10","11"))
ncols <- floor(sqrt(length(unique(cinvs$inv))))
ggplot() +
  geom_polygon(data = world, aes(x=long, y = lat, group=group), fill='grey',color="black",size=0.2) +
  geom_point(data=cinvs, aes(x=lon, y=lat, fill=maf,size=n), shape=21,inherit.aes=F) +
  #geom_text(data=sample_table, aes(x=lon, y=lat, fill=set,label=Population), shape=21,inherit.aes=F) +
  ggtitle(paste("inversion distributions c1 C cluster",sep="")) +
  scale_fill_distiller(limits=c(0,1),palette="Reds",direction = 1)+
  coord_fixed(ylim=c(-48,48),xlim=c(-155,140))+
  theme(axis.text=element_blank(),
        panel.border=element_rect(fill=NA, color="black"),
        panel.background=element_rect(fill='light blue', color=NA),
        axis.title=element_blank(),
        legend.position="bottom") +
  facet_wrap(inv ~ .,ncol=ncols)
ggsave(paste(outprefix,"_c",1,"_Ccluster_inv_dists.png",sep=""),
       dpi=300,width=350,height=200,units='mm')


cinvs <- subset(countryfreqs,chrom==1 & inv %in% c("1","2","3","12"))
ncols <- floor(sqrt(length(unique(cinvs$inv))))
ggplot() +
  geom_polygon(data = world, aes(x=long, y = lat, group=group), fill='grey',color="black",size=0.2) +
  geom_point(data=cinvs, aes(x=lon, y=lat, fill=maf,size=n), shape=21,inherit.aes=F) +
  #geom_text(data=sample_table, aes(x=lon, y=lat, fill=set,label=Population), shape=21,inherit.aes=F) +
  ggtitle(paste("inversion distributions c1 unclustered",sep="")) +
  scale_fill_distiller(limits=c(0,1),palette="Reds",direction = 1)+
  coord_fixed(ylim=c(-48,48),xlim=c(-155,140))+
  theme(axis.text=element_blank(),
        panel.border=element_rect(fill=NA, color="black"),
        panel.background=element_rect(fill='light blue', color=NA),
        axis.title=element_blank(),
        legend.position="bottom") +
  facet_wrap(inv ~ .,ncol=ncols)
ggsave(paste(outprefix,"_c",1,"_unclustered_inv_dists.png",sep=""),
       dpi=300,width=350,height=200,units='mm')




sites <- read.table("resources/aegy.wgs.pops.list.csv",sep=",",header=T)
#cityfreqs <- merge(subset(invfreqs,pop!="all"),locations)

cityfreqs <- subset(invfreqs,pop!="all")
cityfreqs$pop[cityfreqs$pop == "Lope_village"] <- "LopeVillage"
cityfreqs$pop[cityfreqs$pop == "Virhembe"] <- "Virembe"
cityfreqs <- merge(cityfreqs,sites,by.x="pop",by.y="Pop",all.x=T)

ggplot() +
  geom_polygon(data = world, aes(x=long, y = lat, group=group), fill='grey',color="black",size=0.2) +
  geom_point(data=cityfreqs, aes(x=Long, y=Lat, fill=maf,size=n), shape=21,inherit.aes=F) +
  #geom_text(data=sample_table, aes(x=lon, y=lat, fill=set,label=Population), shape=21,inherit.aes=F) +
  ggtitle(paste("inversion distributions c1 unclustered",sep="")) +
  scale_fill_distiller(limits=c(0,1),palette="Reds",direction = 1)+
  coord_fixed(ylim=c(-48,48),xlim=c(-155,140))+
  theme(axis.text=element_blank(),
        panel.border=element_rect(fill=NA, color="black"),
        panel.background=element_rect(fill='light blue', color=NA),
        axis.title=element_blank(),
        legend.position="bottom") +
  facet_wrap(inv ~ .,ncol=ncols)

