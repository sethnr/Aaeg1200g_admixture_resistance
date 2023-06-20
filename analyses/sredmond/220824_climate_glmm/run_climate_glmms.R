library("tidyverse")

library("patchwork")
library("gridExtra")
library("grid")

library("plyr")

library(glmm)





callsfiles <- list.files("lostruct_merge",pattern="invs_chr.*_inv_calls.txt",full.names = T)
callsfiles <- callsfiles[laply(callsfiles, file.size)>0]
rm("calltab")

for(x in callsfiles) {
  icalls <- read.table(x,header=T,sep="\t")
  invnames <- colnames(icalls)[!colnames(icalls) %in% c("sample","spp","pop","region","contgroup","country")]
  sample <- icalls$sample
  icallmat <- icalls[,invnames]
  if(!exists("calltab")) {
    calltab <- cbind(sample,icallmat)
  } else {
    calltab <- cbind(calltab,icallmat)
  }
  
}

meta <- icalls[,c("sample","pop","region","country","contgroup")]
meta$region[meta$region == "Caribean"]<-"Caribbean"
calltab <- merge(calltab,meta)


countries <- unique(calltab$country)
invnames <- colnames(calltab)[!colnames(calltab) %in% c("sample","spp","pop","region","contgroup","country")]
maf <- function(x) {round(sum(x) / (length(x)*2),2)}

popfreqs <- aggregate(. ~ pop,data=calltab[,c("pop",as.character(invnames))],FUN=maf)
popN <- as.data.frame(table(calltab[,c("pop")]))
colnames(popN) <- c("pop","n")
popfreqs <- merge(popfreqs,popN)

poplocs <- read.table("../../resources/aegy.wgs.pops.list.csv",sep=",",header=T) %>% rename_with(tolower)
popvars <- read.table("meta_loc_demo_clim_vals.txt",header=T,sep="\t")
popvars <- merge(poplocs[,c("pop","country","region")],popvars)

popfreqs <- merge(popfreqs,popvars)
popfreqs$abslat <- abs(popfreqs$lat)

allvars <- c("abslat","popdens","maxtemp","mintemp","precipitation","precseason","tempseason")





MINMAF <- 0.05

cntfreqs <- aggregate(. ~ country,data=calltab[,c("country",as.character(invnames))],FUN=maf)
ctypresent = list()
for(I in invnames) {
  ctypresent[[I]] <- cntfreqs$country[cntfreqs[,I] > MINMAF]
  
}


glmmresults <- list()
for(I in invnames) {
  write(I,stderr())
  write(paste(ctypresent[[I]],collapse="/"),stderr())
  write(paste(popfreqs$pop[popfreqs$country %in% ctypresent[[I]]],collapse="/"),stderr())
  if(length(ctypresent[[I]])>1) {
    glmminv <- glmm(as.formula(paste(I,"~ abslat + popdens + maxtemp + mintemp + precipitation + tempseason + precseason")), 
                    random=list(~ 0 + country), varcomps.names=c("region"),
                    #data=subset(popfreqs,region %in% [[I]]),
                    data=popfreqs,
                    family.glmm=binomial.glmm,m=10^3,debug=T)
  } else {
    write("  skipping...",stderr())
  }
  
  glmmresults[[I]]  <- glmminv
}

glmmsummscty <- ldply(invnames,function(I) {s = summary(glmmresults[[I]]); 
as.data.frame(s$coefmat) %>% add_column("var"=rownames(s$coefmat),
                                        "inv"=I,
                                        "converged"=s$trust.converged)})
```




```{r  include=FALSE}

MINMAF <- 0.05
library(glmm)

ctyfreqs <- aggregate(. ~ country,data=calltab[,c("country",as.character(invnames))],FUN=maf)
ctypresent = list()
for(I in invnames) {
  presents <- ctyfreqs$country[ctyfreqs[,I] > MINMAF]
  if(any(presents %in% bigcountry)) {
    ctypresent[[I]] <- presents
  }
}

cities_in_country <- aggregate(pop ~ country,calltab,FUN=function(x) {length(unique(x))})
bigcountry <- cities_in_country$country[cities_in_country$pop>1]

glmmresults <- list()
for(I in invnames) {
  write(I,stderr())
  write(paste(cntpresent[[I]],collapse="/"),stderr())
  write(length(popfreqs$pop[popfreqs$region %in% cntpresent[[I]]]),stderr())
  if(length(cntpresent[[I]])>1) {
    vars <- allvars
    while(length(vars)>=1) {
      formstr <- paste(I,"~",paste(vars,collapse=" + "))
      write(formstr,stderr())
      glmminv <- glmm(as.formula(formstr), 
                      random=list(~ 0 + region), varcomps.names=c("region"),
                      #data=subset(popfreqs,region %in% cntpresent[[I]]),
                      #data=subset(popfreqs,country %in% bigcountry),
                      data=popfreqs,
                      family.glmm=binomial.glmm,m=10^5,debug=T)
      glmmresults[[formstr]]  <- glmminv
      
      #get worst performing var and remove
      coef <- as.data.frame(summary(glmminv)$coefmat)
      coef <- coef[2:nrow(coef),]
      vars <- rownames(coef)[coef$`Pr(>|z|)` < max(coef$`Pr(>|z|)`)]
    }
  }
}

glmmsummscty <- ldply(names(glmmresults),function(I) {s = summary(glmmresults[[I]]); 
as.data.frame(s$coefmat) %>% add_column("var"=rownames(s$coefmat),
                                        "inv"=as.character(s$fixedcall)[2],
                                        "vars"=as.character(s$fixedcall)[3],
                                        "converged"=s$trust.converged)})

glmmsummscty[glmmsummscty$`Pr(>|z|)`<0.05,]

```
