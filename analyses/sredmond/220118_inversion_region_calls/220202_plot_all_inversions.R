library(tidyverse)
library(patchwork)

allinvs <- read.table("220203_all_inv_candidates.txt",header=F,
                      col.names = c("chrom","start","end","country","chromname","valid","invname","d","mean_wss","bss_tot"))

corder <- c("Philippines",
            "Vietnam",
            "SaudiArabia",
            "USA",
            "PuertoRico",
            "Trinidad",
            "Brazil",
            "Senegal",
            "BurkinaFaso",
            "Ghana",
            "Gabon",
            "Uganda",
            "Kenya")
allinvs$country <- factor(allinvs$country,levels=corder,ordered=T)

chroms <- data.frame("chrom"=c(1,2,3),
                     "start"=c(1,1,1),
                     "end"=c(310827022,474425716,409777670))
chroms <- merge(chroms,unique(allinvs$country)) %>% rename("country"="y")



allinvsp <- ggplot(data=subset(allinvs,!valid),aes(y=country,yend=country,x=start,xend=end,color=valid)) + 
  geom_segment(data=chroms,color="grey") +
  geom_curve() + 
  geom_curve(data=subset(allinvs,valid)) + 
  facet_grid("chrom ~ .") + theme(legend.position="bottom", panel.grid.major.y=element_blank(),axis.title=element_blank(),)
allinvsp
ggsave("220203_all_inv_candidates.png",width=200,height=150,units="mm")

group=0
allinvs$group=NA

OLMIN = 0.90
for(i in c(1:(nrow(allinvs)-1))) {
  for(j in c((i+1):nrow(allinvs))) {
    ist <- allinvs[i,"start"]
    ien <- allinvs[i,"end"]
    jst <- allinvs[j,"start"]
    jen <- allinvs[j,"end"]
    
    ival <- allinvs[i,"valid"]
    jval <- allinvs[j,"valid"]
    
    igrp <- allinvs[i,"group"]
    jgrp <- allinvs[j,"group"]
    
    
    ilen <- ien-ist
    jlen <- jen-jst
    ol <- min(ien,jen)-max(ist,jst)
    
    if((ol > ilen*OLMIN & ol > jlen*OLMIN) & (ival | jval)) {
      write(paste(i,j,nrow(allinvs)),file=stderr())
      if(is.na(igrp) & is.na(jgrp)) {
        group <- group+1
        allinvs[i,"group"] <- group
        allinvs[j,"group"] <- group
        write("\tassigning both",file=stderr())
      } else if (!is.na(igrp) & is.na(jgrp)) {
        allinvs[j,"group"] <- igrp
        write("\tassigning i > j",file=stderr())
      } else if (!is.na(igrp) & igrp==jgrp) {
        write("\tboth assigned",file=stderr())
      } else if (!is.na(igrp) & igrp!=jgrp) {
        write("\tconflict!",file=stderr())
      }
    }
    
  }
}

for(i in c(1:(nrow(allinvs)))) {
  ival <- allinvs[i,"valid"]
  igrp <- allinvs[i,"group"]
  if(ival & is.na(igrp)) {
    group<-group+1
    allinvs[i,"group"] <- group
  }
}

goodinvs <- subset(allinvs,valid | group>0)

gorder <- unique(goodinvs$group[order(goodinvs$chrom,goodinvs$start,goodinvs$end)])
allgroups <- unique(goodinvs$group[order(goodinvs$group)])
names(allgroups) <- gorder
goodinvs$group <- factor(allgroups[as.character(goodinvs$group)],levels=allgroups,ordered=T)

goodinvsp <- ggplot(goodinvs,aes(y=country,yend=country,x=start,xend=end,color=valid)) + 
  geom_segment(data=chroms,color="grey") +
  geom_curve() + 
  geom_curve(data=subset(allinvs,valid)) + scale_color_manual(values=c("green","dark green"),name="valid") +
  facet_grid("chrom ~ .") + theme(legend.position="bottom", panel.grid.major.y=element_blank(), 
                                  axis.title=element_blank())
goodinvsp
ggsave("220203_all_inv_candidates_valid.png",width=200,height=150,units="mm")

write.table(goodinvs,
            file="220203_all_inv_candidates_grouped.txt",sep="\t",row.names=F)

countinvsp <- ggplot(data=subset(allinvs,valid | group>0),aes(x=country,fill=as.factor(group))) + 
  geom_bar() +guides(fill=guide_legend("inv",nrow=2),) + 
  facet_grid("chrom ~ .") + theme(legend.position="bottom",
                                  axis.text.x = element_text(angle=45,hjust=1),
                                  axis.title=element_blank())
countinvsp
ggsave("220203_all_inv_candidates_counts.png",width=150,height=150,units="mm")


invgroupcols <- c("#AAAAAA","#777777","#555555")[as.numeric(levels(goodinvs$group))%%3+1]
names(invgroupcols) <- levels(goodinvs$group)
invgroupcols[which(table(goodinvs$group)>1)] <- rainbow(sum(table(goodinvs$group)>1))


countinvsp_flip <- ggplot(unique(goodinvs[c("country","group","chrom")]),aes(x=country,fill=forcats::fct_rev(as.factor(group)))) + 
  geom_bar() + coord_flip() + scale_x_discrete(drop=F) + scale_fill_manual(values=invgroupcols) + guides(fill=guide_legend("inv",nrow=2),) + 
  facet_grid("chrom ~ .") + theme(legend.position="bottom",
                                  axis.title=element_blank(),
                                  axis.text.y=element_blank())
countinvsp_flip
allinvsp + plot_spacer() + plot_layout(nrow=1,widths=c(10,3),guides = "collect") & theme(legend.position='bottom')
ggsave("220203_all_inv_candidates_comb.png",width=200,height=150,units="mm")

goodinvsp + countinvsp_flip + plot_layout(nrow=1,widths=c(10,3),guides = "collect") & theme(legend.position='bottom')
ggsave("220203_good_inv_candidates_comb.png",width=200,height=150,units="mm")
