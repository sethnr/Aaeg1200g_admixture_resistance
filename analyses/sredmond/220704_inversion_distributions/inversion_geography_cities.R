library("tidyverse")

library("patchwork")
library("gridExtra")
library("grid")

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

invfreqs$inv <- paste(invfreqs$chrom,invfreqs$inv,sep="-")

hwecountries <- unique(invfreqs$country[invfreqs$HWE==0])


countryfile <- "country_geocoords.txt"
if(file.exists(countryfile)){
  countrylocs <- read.table(countryfile,header=T) 
} else {
  register_google(key="AIzaSyBFbT5X6MjrepPAlzIj3zz1HS58hnIPtEM", write=T)
  countries <- unique(invfreqs$country)
  countrylocs <- geocode(countries)
  countrylocs$country <- countries
  write.table(countrylocs,countryfile,sep="\t",row.names=F,col.names=T)
}


cmaps <- list()
bigcountries <- c("Brazil","SaudiArabia","Argentina")
for(C in hwecountries) {
  lon <- countrylocs$lon[countrylocs$country==C]
  lat <- countrylocs$lat[countrylocs$country==C]
  if(C %in% bigcountries) {
    cmap <- get_map(c(lon,lat),zoom=4)
  } else {
    cmap <- get_map(c(lon,lat),zoom=6)}
  cmaps[[C]] <- cmap
}


sites <- read.table("resources/aegy.wgs.pops.list.csv",sep=",",header=T)
cityfreqs <- subset(invfreqs,pop!="all")
cityfreqs$pop[cityfreqs$pop == "Lope_village"] <- "LopeVillage"
cityfreqs$pop[cityfreqs$pop == "Virhembe"] <- "Virembe"
cityfreqs <- merge(cityfreqs,sites,by.x="pop",by.y="Pop",all.x=T)


for(I in unique(invfreqs$inv[invfreqs$HWE==0])) {
  ggmaps <- list()
  for(C in unique(invfreqs$country[invfreqs$inv==I & invfreqs$HWE==0])) {
    #get inv freqs for _just_ one country
    ccfreqs <- cityfreqs[cityfreqs$country==C & cityfreqs$inv==I,]
    ggmaps[[C]] <- ggmap(cmaps[[C]]) +
      geom_point(data=ccfreqs, aes(x=Long, y=Lat, fill=maf,size=n), shape=21,inherit.aes=F) +
    #   #geom_text(data=sample_table, aes(x=lon, y=lat, fill=set,label=Population), shape=21,inherit.aes=F) +
      scale_fill_distiller(limits=c(0,1),palette="Reds",direction = 1)+
      scale_size_continuous(limits=c(0,50))+
      theme(axis.text=element_blank(),
            panel.border=element_rect(fill=NA, color="black"),
            panel.background=element_rect(fill='light blue', color=NA),
            axis.title=element_blank(),
            legend.position="none")
  }
  png(paste(outprefix,"_site_distribution_",I,".png",sep=""),res=300,width=200,height=200,units='mm')
  do.call(grid.arrange,ggmaps)
  dev.off()
  
}




for(C in c("Kenya")) {
  ggmaps <- list()
  for(I in unique(invfreqs$inv[invfreqs$HWE==0 & invfreqs$country==C])) {
    #get inv freqs for _just_ one country
    ccfreqs <- cityfreqs[cityfreqs$country==C & cityfreqs$inv==I,]
    ggmaps[[I]] <- ggmap(cmaps[[C]]) +
      geom_point(data=ccfreqs, aes(x=Long, y=Lat, fill=maf,size=n), shape=21,inherit.aes=F) +
      #   #geom_text(data=sample_table, aes(x=lon, y=lat, fill=set,label=Population), shape=21,inherit.aes=F) +
      scale_fill_distiller(limits=c(0,1),palette="Reds",direction = 1)+
      scale_size_continuous(limits=c(0,50))+
      ggtitle(I) +
      theme(axis.text=element_blank(),
            panel.border=element_rect(fill=NA, color="black"),
            panel.background=element_rect(fill='light blue', color=NA),
            axis.title=element_blank(),
            legend.position="none")
  }
  png(paste(outprefix,"site_distribution_Kenya.png",sep="_"),res=300,width=200,height=200,units='mm')
  do.call(grid.arrange,ggmaps)
  dev.off()
  
}

# ggmap(cmaps[["Kenya"]]) +
#   geom_point(data=unique(cityfreqs[cityfreqs$Country=="Kenya",c("Country","Search","Lat","Long")]), 
#              aes(x=Long, y=Lat), shape=21,inherit.aes=F) +
#   geom_text(data=unique(cityfreqs[cityfreqs$Country=="Kenya",c("Country","Search","Lat","Long")]), 
#              aes(x=Long, y=Lat,label=Search), inherit.aes=F) +
#   scale_size_continuous(limits=c(0,50))+
#   ggtitle(I) +
#   theme(axis.text=element_blank(),
#         panel.border=element_rect(fill=NA, color="black"),
#         panel.background=element_rect(fill='light blue', color=NA),
#         axis.title=element_blank(),
#         legend.position="none")
# 

