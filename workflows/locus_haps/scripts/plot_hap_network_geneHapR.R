# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.14")
# 
# 
#BiocManager::install("GenomicRanges",update=F)
#BiocManager::install("Biostrings",update=F)
#BiocManager::install("rtracklayer",update=F)
#BiocManager::install("trackViewer",update=F)
#install.packages("geneHapR")
# 

library("geneHapR")
library("vcfR")
library("pegas")
library("getopt")
library("tidyverse")

opttab <- matrix(c("vcf","i","1","character",
                   "out","o","1","character",
                   "meta","m","1","character"
),byrow=T,ncol=4)
opt <- getopt(opttab)

#f <-"locushap_shapeit/snpeff.kdr3.m0.01_noUTR.vcf.gz"
#prefix <- "kdr3"


vcffile <- opt$vcf
prefix <- opt$out
metafile <- opt$meta

meta <- read.table(metafile,sep="\t",header=T)

metahap <- rbind(meta %>% mutate(sample=paste(sample,"0",sep="_")),
                 meta %>% mutate(sample=paste(sample,"1",sep="_")))
rownames(metahap) <- metahap$sample

poplocs <- read.table("resources/aegy.wgs.pops.list.csv",sep=",",header=T) %>% 
                    select("Pop","Lat","Long") %>% 
                    rename("pop"="Pop","poplat"="Lat","poplong"="Long")
cntlocs <- read.table("resources/country_geocoords.txt",sep="\t",header=T) %>% 
                    rename("cntlat"="lat","cntlong"="lon")


metahaploc <- merge(merge(metahap,cntlocs,by="country"),
                    poplocs,by="pop") %>% rename("Hap"="sample")
rownames(metahaploc) <- metahaploc$Hap


# vcf <- import_vcf(vcffile)
# hapResult <- vcf2hap(vcf,
#                      hapPrefix = "H",
#                      hetero_remove = FALSE,
#                      na_drop = T)

#vcf <- pegas::read.vcf(vcffile)
loci <- VCFloci(vcffile)

vcfr <- read.vcfR(vcffile)
gts <- extract.gt(vcfr)
haps <- extract.haps(vcfr)

#only include samples in haps file (removed related, etc)
haps <- haps[,colnames(haps) %in% metahap$sample]

#paste(rep(VCFlabels(vcffile),each=2),c("0","1"),sep="_")
haptab <- cbind(loci[,c(1,2,4,5,8)],haps)

hapResult <- table2hap(haptab[sort(c(grep("MODERATE",haptab$INFO),
                                grep("HIGH",haptab$INFO),
                                grep("315939224",haptab$POS),
                                grep("315999297",haptab$POS) )),])

                                


hapSummary <- hap_summary(hapResult)

W=300
H=200
R=400
png(paste(prefix,"haptable.png",sep="_"),width=W,height=H,res=R,units="mm")
#plot haplotype table
plotHapTable(hapSummary)
dev.off()


#plot haplotype network
hapNet <- get_hapNet(hapSummary,
                     AccINFO = metahap,
                     groupName = "region")

# hapcols = c("West Africa"="dark blue","East Africa"="light blue",
#             "Carribean"="Purple","South America"="Dark Red","North America"="Light Red",
#             "Middle East"="Light Green","Asia"="Green","Pacific"="Dark Green")

png(paste(prefix,"hapnet.png",sep="_"),width=W,height=H,res=R,units="mm")
plotHapNet(hapNet,
           show.mutation = 2,
           legend=c(10,-15),
           show_size_legend=F,
           scale='log2',
           threshold=0) #+ scale_fill_manual(values=hapcols)
dev.off()

png(paste(prefix,"geodist_pop.png",sep="_"),width=W,height=H,res=R,units="mm")
hapDistribution(hapResult,metahaploc,"poplong","poplat",
                c("H001","H002","H003","H004","H005","H006"),
                legend=TRUE)
dev.off()

png(paste(prefix,"geodist_cnt.png",sep="_"),width=W,height=H,res=R,units="mm")
hapDistribution(hapResult,metahaploc,"cntlong","cntlat",
                c("H001","H002","H003","H004","H005",
                  "H006","H007","H008","H009","H010"),
                legend=TRUE)
dev.off()

