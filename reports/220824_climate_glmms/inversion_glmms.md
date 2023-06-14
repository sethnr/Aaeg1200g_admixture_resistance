---
output: html_document
editor_options: 
  chunk_output_type: console
---



```r
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
```






```r
glmmresults <- list()
for(I in invnames) {
  write(paste(I,"(",length(popfreqs$pop[popfreqs$country %in% ctypresent[[I]]]),")"),stderr())
  if(length(rgnpresent[[I]])>1) {
    vars <- allvars
    while(length(vars)>=1) {
      formstr <- paste(I,"~",paste(vars,collapse=" + "))
      write(formstr,stderr())
      ws <- c(popfreqs$n)
      glmminv <- glmm(as.formula(formstr), 
                      random=list(~ 0 + country), varcomps.names=c("country"),
                      data=subset(popfreqs,region %in% rgnpresent[[I]]),
                      #data=popfreqs,
                      family.glmm=binomial.glmm,m=ITERS,debug=T,
                      doPQL = T)
      
      glmmresults[[formstr]]  <- glmminv
      
      #get worst performing var and remove
      coef <- as.data.frame(summary(glmminv)$coefmat)
      coef <- coef[2:nrow(coef),]
      vars <- rownames(coef)[coef$`Pr(>|z|)` < max(coef$`Pr(>|z|)`)]
    }
  }
}
```

```
## Error in mcvcov(object): NA/NaN/Inf in foreign function call (arg 29)
```

```r
glmmsummscty <- ldply(names(glmmresults),function(G) {m = glmmresults[[G]]
                                          s = summary(m); 
                                          t <- as.data.frame(s$coefmat) %>% add_column("var"=rownames(s$coefmat),
                                          "inv"=as.character(s$fixedcall)[2],
                                          "vars"=as.character(s$fixedcall)[3],
                                          "converged"=s$trust.converged)
                                          t})
```

```
## Error in mcvcov(object): NA/NaN/Inf in foreign function call (arg 29)
```

```r
unique(glmmsummscty[glmmsummscty$`Pr(>|z|)`<0.05 & glmmsummscty$var != "(Intercept)",c("vars","inv")])
```

```
##                              vars        inv
## 30 popdens + mintemp + tempseason X1_Gabon_1
## 33           mintemp + tempseason X1_Gabon_1
```

```r
unique(glmmsummscty[glmmsummscty$`Pr(>|z|)`<0.05 & glmmsummscty$var != "(Intercept)",])
```

```
##    Estimate Std. Error   z value   Pr(>|z|)        var        inv
## 30  -0.0195     0.0088 -2.214725 0.02677897 tempseason X1_Gabon_1
## 33  -0.0204     0.0088 -2.320469 0.02031552 tempseason X1_Gabon_1
##                              vars converged
## 30 popdens + mintemp + tempseason      TRUE
## 33           mintemp + tempseason      TRUE
```

```r
write.table(glmmsummscty,"./glmm_pql_results_country_climate_variables.txt",sep="\t",quote=F,row.names=F)
```



```r
glmmresults <- list()
for(I in invnames) {
#for(I in likelies) {
  write(paste(I,"(",length(popfreqs$pop[popfreqs$country %in% ctypresent[[I]]]),")"),stderr())
  if(length(ctypresent[[I]])>1) {
    vars <- allvars[allvars!="abslat"]
    while(length(vars)>=1) {
      formstr <- paste(I,"~",paste(vars,collapse=" + "))
      write(formstr,stderr())
      ws <- c(popfreqs$n)
      glmminv <- glmm(as.formula(formstr), 
                      random=list(~ 0 + region), varcomps.names=c("region"),
                      data=subset(popfreqs,region %in% rgnpresent[[I]]),
                      #data=popfreqs,
                      family.glmm=binomial.glmm,m=ITERS,debug=T,
                      doPQL = T)
      
      glmmresults[[formstr]]  <- glmminv
      
      #get worst performing var and remove
      coef <- as.data.frame(summary(glmminv)$coefmat)
      coef <- coef[2:nrow(coef),]
      vars <- rownames(coef)[coef$`Pr(>|z|)` < max(coef$`Pr(>|z|)`)]
    }
  }
}
```

```
## Error in `contrasts<-`(`*tmp*`, value = contr.funs[1 + isOF[nn]]): contrasts can be applied only to factors with 2 or more levels
```

```r
glmmsummscty <- ldply(names(glmmresults),function(G) {m = glmmresults[[G]]
                                          s = summary(m); 
                                          t <- as.data.frame(s$coefmat) %>% add_column("var"=rownames(s$coefmat),
                                          "inv"=as.character(s$fixedcall)[2],
                                          "vars"=as.character(s$fixedcall)[3],
                                          "converged"=s$trust.converged)
                                          t})
unique(glmmsummscty[glmmsummscty$`Pr(>|z|)`<0.05 & glmmsummscty$var != "(Intercept)",c("vars","inv")])
```

```
##                              vars        inv
## 22 popdens + mintemp + tempseason X1_Gabon_1
```

```r
unique(glmmsummscty[glmmsummscty$`Pr(>|z|)`<0.05 & glmmsummscty$var != "(Intercept)",])
```

```
##    Estimate Std. Error   z value   Pr(>|z|)        var        inv
## 22   -0.021     0.0089 -2.361122 0.01821973 tempseason X1_Gabon_1
##                              vars converged
## 22 popdens + mintemp + tempseason      TRUE
```

```r
write.table(glmmsummscty,"./glmm_pql_results_region_climate_variables_nolat.txt",sep="\t",quote=F,row.names=F)
```
