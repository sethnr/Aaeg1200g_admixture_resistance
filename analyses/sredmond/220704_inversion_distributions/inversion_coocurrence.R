library("tidyverse")

library("patchwork")
library("gridExtra")
library("grid")

library("genetics")
library("getopt")

opttab <- matrix(c("calls","i","1","character",
                   "outfile","o","1","character"),byrow=T,ncol=4)
opt <- getopt(opttab)

callsfile <- opt$calls
outprefix <- opt$outfile

calldf <- read.table(callsfile,header=T,sep="\t")
invnames <- colnames(calldf)[!colnames(calldf) %in% c("sample","spp","pop","region","contgroup","country")]

countries <- unique(calldf$country)

invfreqs <- aggregate(. ~ country,data=calldf[,c("country",as.character(invnames))],FUN=function(x) {sum(x) / (length(x)*2)})
invfreqscont <- aggregate(. ~ contgroup,data=calldf[,c("contgroup",as.character(invnames))],FUN=function(x) {sum(x) / (length(x)*2)})

cormatrix <- matrix(rep(NA,length(invnames)^2),
                    dimnames = list(invnames,invnames),
                    nrow=length(invnames))

cortab <- data.frame(
  inv1=character(),
  inv2=character(),
  occurpops=logical(),
  p=numeric(),
  cor=numeric()
  )

write("getting correlations",stderr())
minmaf <- 0.05
for(I1 in invnames) {
  for(I2 in invnames) {
    if(match(I1,invnames) >= match(I2,invnames)) {next}
    incountries <- invfreqs$country[invfreqs[,I1] > minmaf & invfreqs[,I2] > minmaf]
    
    cortest <- cor.test(calldf[,I1],
                        calldf[,I2])
    
    i <- nrow(cortab)+1
    cortab[i,c("inv1","inv2")] <- c(I1,I2)
    cortab[i,c("p","cor")] <- c(cortest$p.value,cortest$estimate)
    cortab[i,c("occurpops")] <- FALSE
    
    if(length(incountries)>0) {
      cortest <- cor.test(calldf[calldf$country %in% incountries,I1],
                          calldf[calldf$country %in% incountries,I2])
      
      i <- nrow(cortab)+1
      cortab[i,c("inv1","inv2")] <- c(I2,I1)
      cortab[i,c("occurpops")] <- TRUE
      cortab[i,c("p","cor")] <- c(cortest$p.value,cortest$estimate)
    }
  }
}


occurtab <- data.frame(
  inv1=character(),
  inv2=character(),
  continent=logical(),
  country=logical(),
  exclusive=logical()
)
write("getting occurrences",stderr())

for(I1 in invnames) {
  for(I2 in invnames) {
    if(I1==I2) {next}
    incountries <- invfreqs$country[invfreqs[,I1] > minmaf & invfreqs[,I2] > minmaf]
    inconts <- invfreqscont$contgroup[invfreqscont[,I1] > minmaf & invfreqscont[,I2] > minmaf]
    
    excl <- (all(calldf[calldf[,I1]==2,I2]==0) &
               all(calldf[calldf[,I2]==2,I1]==0)) &
            (max(calldf[calldf[,I1]==1,I2])==1 &
               max(calldf[calldf[,I2]==1,I1])==1)
    
    occurcountry <- length(incountries)>0
    occurcont <- length(inconts)>0
    
    i <- nrow(occurtab)+1
    occurtab[i,c("inv1","inv2")] <- c(I1,I2)
    occurtab[i,c("continent","country","exclusive")] <- c(occurcont,occurcountry,excl)

  }    
}

# merge(subset(occurtab,exclusive),
#       subset(cortab,country == "all"),
#       by=c("inv1","inv2"))

cortab$inv1 <- factor(cortab$inv1,levels=invnames,ordered=T)
cortab$inv2 <- factor(cortab$inv2,levels=invnames,ordered=T)
ggplot(subset(cortab),aes(x=inv1,y=inv2,fill=cor)) + 
  geom_raster() + 
  geom_tile(data=subset(occurtab,exclusive),fill=NA,color="orange",size=1)+ 
  coord_fixed() + scale_fill_distiller(palette="RdBu",limits=c(-1,1)) +
  ylab("all samples") + xlab(paste("maf >",minmaf)) +
  theme(panel.background =element_rect("dark grey"),
        panel.grid = element_blank())

ggsave(paste(outprefix,"cor.png",sep="_"))
write.table(cortab,paste(outprefix,"cor.txt",sep="_"),sep="\t",quote=F,row.names=F,col.names=T)
write.table(occurtab,paste(outprefix,"occur.txt",sep="_"),sep="\t",quote=F,row.names=F,col.names=T)
