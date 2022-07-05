
library("tidyverse")
library("lostruct")

#library("patchwork")
library("gridExtra")
library("grid")

library("getopt")

opttab <- matrix(c("infile","i","1","character",
                   "samples","s","1","character",
                   "outfile","o","1","character"
),byrow=T,ncol=4)
opt <- getopt(opttab)

invaimsfile <- opt$infile
metafile <- opt$samples
outsnpspng <- opt$outfile


write("gathering meta",file=stderr())

metatab <- read.table(metafile,header=T, sep="\t")
metatab$contgroup <- factor(metatab$contgroup,levels=c("Wafrica","Safrica","Eafrica","Americas","Asia"),ordered=T)
metatab$region <- factor(metatab$region,levels=c("East Africa","West Africa","South America","Carribean","North America","Middle East","Asia","Pacific"),ordered=T)

countrysorttab <- unique(metatab[,c("country","contgroup","region")])
metatab$country <- factor(metatab$country,levels=unique(metatab$country[order(metatab$region)]),ordered=T)


# invcands <- read.table(invfile,header=T)
# invcands <- subset(invcands,as.logical(valid))

write("reading inv snps",file=stderr())

if (file.size(invaimsfile)>0) {
  allinvsnps <- read.table(invaimsfile,header=T,check.names = F)
  write.table(table(allinvsnps$inv),file=stderr(),row.names = F,quote=F,col.names = F)
  write(paste(length(unique(allinvsnps$inv)),"invs found in",invaimsfile),file=stderr())
  allinvnames <- unique(allinvsnps$inv)
} else {
  allinvnames <- c()
}

nsamples <- nrow(metatab)
metatab <- metatab[metatab$sample %in% colnames(allinvsnps),]
write(paste("found",nrow(metatab),"samples of",nsamples),stderr())
samples <- metatab$sample
#write(samples,stderr())

#order countries in metatable by region (w african, eafrican, american, asian)
cntorder <- unique(metatab$country[order(metatab$contgroup,metatab$region)])


# mlpreds <- read.table(mlcallsfile,header=F,col.names = unique(allinvsnps$inv))
# meanpreds <- read.table(meancallsfile,header=T)

aimplots=list()
for(invname in allinvnames) {
    write(paste(" ",invname),file=stderr())
    invnamesafe <- paste("X",gsub("\\D",".",invname,perl=T),sep="")
    invsnps <- subset(allinvsnps,inv==invname)

    #reset inv positions per set in case using merged aim blocks

    write(paste("making sort arrays (country, mean SNP value)"),stderr())
    meansnp <- apply(invsnps[,samples],2,FUN=function(x) {mean(na.omit(x))})
    write(paste("got mean snp for ",dim(invsnps)[1],"snps"),stderr())

    write("sorting SNPs",stderr())
    invsnps$i <- order(invsnps$pos)

    aimsM <- pivot_longer(invsnps,any_of(samples),names_to = "sample") %>% rename("invcountry"="country")
    aimsM <- merge(aimsM,metatab,by="sample")

    write("sorting countries",stderr())
    aimsM$country <- factor(aimsM$country,levels = cntorder,ordered=T)

    write("sorting samples",stderr())
    #order SNPs by country, then mean call in this inversion
    snporder <- metatab$sample[order(metatab$country,meansnp[metatab$sample])]

    #write(snporder,stderr())
    aimsM$sample <- factor(aimsM$sample,levels = snporder,ordered=T)
    aimsM$cncode <- aimsM$country
    levels(aimsM$cncode) <- substr(levels(aimsM$country),0,3)

    aimplot <- ggplot(aimsM,aes(x=i,y=as.numeric(sample),fill=as.factor(value))) + geom_raster() +
       ylab("samples") + xlab("SNPs")+ theme(legend.position="none")+
       scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) +
      scale_fill_manual(values=c("0"="blue","1"="purple","2"="red")) +
      facet_grid("cncode ~ .",scale="free_y",space="free_y") +
      ggtitle(paste(invname,"aim calls","( top",nrow(invsnps),")"))+
      theme(panel.spacing = unit(0.2, "mm"),
            axis.text.y=element_blank(),
            axis.title.y=element_blank(),
            axis.ticks.y=element_blank(),
            strip.background = element_blank(),
            strip.text = element_blank(),
            plot.margin = unit(c(0,0,0,0), "mm"))
    aimplots[[invnamesafe]] <- aimplot



}

if(length(aimplots)>0) {
  print(length(aimplots),stderr())
  png(outsnpspng,res=400,width=400,height=200,units='mm')
    do.call("grid.arrange", c(aimplots, nrow=1))
  dev.off()
} else {file.create(outsnpspng)}
