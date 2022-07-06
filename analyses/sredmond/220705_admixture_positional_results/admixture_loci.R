library("tidyverse")

library("patchwork")
library("gridExtra")
library("plyr")

library("getopt")



indir <- "admix_tests"
outprefix <- "admix_loci"

#make function to read/parse all vars and create table
parse_3pop <- function(x) {
  fnspl <- strsplit(basename(x),"_")[[1]]
  chr = fnspl[2]
  mixpop <- fnspl[5]
  pop1 <- fnspl[3]
  pop2 <- fnspl[4]
  read.table(x,col.names = c("start","mid","end","f3")) %>% 
    add_column("p1"=pop1,"p2"=pop2,"p3"=mixpop,"chr"=chr)
}


write(paste("getting blocks from ",indir),stderr())
pop3files <- list.files(indir,pattern="3pop.*_blocks.txt",full.names = T)
pop3tab <- ldply(pop3files, parse_3pop)
pop3tab[pop3tab$p2>pop3tab$p1,c("p1","p2")] <- pop3tab[pop3tab$p2>pop3tab$p1,c("p2","p1")]
pop3tab$parent <- paste(pop3tab$p1,pop3tab$p2,sep="/")
pop3tab$set <- paste(pop3tab$p3,pop3tab$p1,pop3tab$p2,sep="/")


pop3files <- list.files(indir,pattern="3pop.*_summary.txt",full.names = T)
pop3sum <- ldply(pop3files, function(x) {read.table(x,col.names=c("p3","p1","p2","chrom","f3","f3sd","f3z","F3sig"))})
pop3sum$set <- paste(pop3sum$p3,pop3sum$p1,pop3sum$p2,sep="/")
pop3sum$F3sig <- as.logical(pop3sum$F3sig)
sigpops <- pop3sum$set[pop3sum$F3sig]

ggplot(pop3tab,aes(x=mid,y=f3,group=set,color=parent)) + geom_line() + facet_grid(p3 ~ chr,scale="free",space="free_x")


ggplot(subset(pop3tab,set %in% sigpops),aes(x=mid,y=f3,group=set,color=parent)) + 
         geom_line() + 
         facet_grid(p3 ~ chr,scale="free",space="free_x")


pop3tabneg <- pop3tab
pop3tabneg$f3[pop3tabneg$f3>0] <- 0
ggplot(pop3tabneg,aes(x=mid,y=f3,group=set,color=parent)) + 
  geom_line() + facet_grid(p3 ~ chr,scale="free_x",space="free_x")
