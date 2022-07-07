library("tidyverse")

library("patchwork")
library("gridExtra")
library("plyr")

library("getopt")



indir <- "admix_tests"
outprefix <- "3pop_admix_loci"

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
pop3sum <- ldply(pop3files, function(x) {read.table(x,col.names=c("p3","p1","p2","chr","f3","f3sd","f3z","f3sig"))})
pop3sum$set <- paste(pop3sum$p3,pop3sum$p1,pop3sum$p2,sep="/")
pop3sum$f3sig <- as.logical(pop3sum$f3sig)

sigpops <- unique(pop3sum$p3[pop3sum$f3sig])
#sigsets <- pop3sum$set[pop3sum$f3sig]
#get those that show sig f3 but not from admixed parents
sigsets <- pop3sum$set[pop3sum$f3sig & ((!pop3sum$p1 %in% sigpops) & (!pop3sum$p2 %in% sigpops))]

ggplot(pop3tab,aes(x=mid,y=f3,group=set,color=parent)) + geom_line() + facet_grid(p3 ~ chr,scale="free",space="free_x")


ggplot(subset(pop3tab,set %in% sigsets),aes(x=mid,y=f3,group=set,color=parent)) + 
         geom_line() + 
         facet_grid(p3 ~ chr,scale="free",space="free_x")
ggsave(paste(outprefix,"f_sigpops.png",sep="_"),dpi = 300,width=250,height=175,units="mm")

pop3tabneg <- pop3tab
pop3tabneg$f3[pop3tabneg$f3>0] <- 0
ggplot(subset(pop3tabneg,set %in% sigsets),aes(x=mid,y=f3,group=set,color=parent)) + 
  geom_line() + facet_grid(p3 ~ chr,scale="free_x",space="free_x")
ggsave(paste(outprefix,"f_sigpops_neg.png",sep="_"),dpi = 300,width=250,height=175,units="mm")


#######
# make and plot z scores
####### 

pop3tab2 <- merge(pop3tab,pop3sum[,c("set","chr","f3sd","f3sig")],by=c("set","chr"))
pop3tab2$f3z <- pop3tab2$f3/pop3tab2$f3sd

pop3tabneg <- pop3tab2
pop3tabneg$f3z[pop3tabneg$f3z>0] <- 0
ggplot(subset(pop3tabneg,set %in% sigsets),aes(x=mid,y=f3z,group=set,color=parent)) + 
  geom_line() + facet_grid(p3 ~ chr,scale="free_x",space="free_x")
ggsave(paste(outprefix,"fz_sigpops_neg.png",sep="_"),dpi = 300,width=250,height=175,units="mm")


ggplot(subset(pop3tabneg,!p1 %in% sigpops & !p2 %in% sigpops & !set %in% sigsets),aes(x=mid,y=f3z,group=set,color=parent)) + 
  geom_line() + 
  facet_grid(p3 ~ chr,scale="free_x",space="free_x")
ggsave(paste(outprefix,"fz_nonsig_neg.png",sep="_"),dpi = 300,width=250,height=175,units="mm")


