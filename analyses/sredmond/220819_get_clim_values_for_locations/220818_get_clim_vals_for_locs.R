library("tidyverse")

library("patchwork")
library("gridExtra")
library("grid")

library("genetics")

library(maptools)
library(ggmap)

library("raster")


poplocs <- read.table("../../../resources/aegy.wgs.pops.list.csv",sep=",",header=T) %>% rename_with(tolower)


rasterfile <- "resources/rasters/ppp_2020_1km_Aggregated.tif"
density <- raster(rasterfile)

poplocs$popdens <- extract(density,poplocs[,c("long","lat")],method="bilinear")


seasonality <- raster("resources/rasters/wc2.1_2.5m_bio_4.tif")
maxtemp <- raster("resources/rasters/wc2.1_2.5m_bio_5.tif")
mintemp <- raster("resources/rasters/wc2.1_2.5m_bio_6.tif")
precipitation <- raster("resources/rasters/wc2.1_2.5m_bio_12.tif")
tempseasonality <- raster("resources/rasters/wc2.1_2.5m_bio_4.tif")
precseasonality <- raster("resources/rasters/wc2.1_2.5m_bio_15.tif")

poplocs$maxtemp <- extract(maxtemp,poplocs[,c("long","lat")],method="bilinear")
poplocs$mintemp <- extract(mintemp,poplocs[,c("long","lat")],method="bilinear")
poplocs$precipitation <- extract(precipitation,poplocs[,c("long","lat")],method="bilinear")
poplocs$precseason <- extract(seasonality,poplocs[,c("long","lat")],method="bilinear")
poplocs$tempseason <- extract(seasonality,poplocs[,c("long","lat")],method="bilinear")


colnames(poplocs)

write.table(poplocs[,c("index","label","lat","long","popdens","maxtemp","mintemp","precipitation","precseason","tempseason")],
            "meta_loc_demo_clim_vals.txt",sep="\t",quote=F,row.names=F,col.names=T)
