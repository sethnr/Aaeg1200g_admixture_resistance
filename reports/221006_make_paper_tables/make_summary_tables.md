


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

meta <- icalls[,c("sample","pop","region","country")]
calltab <- merge(calltab,meta)
```



```r
countries <- unique(calltab$country)
invnames <- colnames(calltab)[!colnames(calltab) %in% c("sample","spp","pop","region","contgroup","country")]
maf <- function(x) {round(sum(x) / (length(x)*2),2)}
#order countries in metatable by region (w african, eafrican, american, asian)
regorder <- c("East Africa","West Africa","Caribean","South America","North America","Middle East","Asia","Pacific")
meta$region <- factor(meta$region,levels=regorder,ordered=T)
cntorder <- unique(meta$country[order(meta$region)])
```


#calculate HWE for pops

```r
hwetab <- data.frame("inversion"=factor(levels=invnames),
           "country"=factor(levels=unique(calltab$country)),
           "pop"=factor(levels=unique(calltab$pop)),
           "N"=numeric(),
           "aa"=numeric(),
           "ab"=numeric(),
           "bb"=numeric(),
           "MAF"=numeric(),
           "F"=numeric(),
           "HWE"=numeric()
           )

maf <- function(x) {round(sum(x) / (length(x)*2),2)}
f <- function(x) {
    Ho<-sum(na.omit(x==1))/length(na.omit(x));
    He<-2*(maf(x)*(1-maf(x)));
    f <- Ho/He;
    f}

for (C in unique(calltab$country)) {
      countrycalls <- calltab[calltab$country==C,]
      
      for(I in invnames) {
        ICcalls <- calltab[calltab$country==C,I]
        geno <- genotype(c("a/a","a/b","b/b")[calltab[calltab$country==C,I]+1])
        if(nallele(geno)==2) {
          test <- HWE.test(geno)
          Pval <- test$test$p.value
        } else{
          Pval <- 1
        }
        i <- nrow(hwetab)+1
        hwetab[i,c("aa","ab","bb")] <- c(sum(geno=="a/a"),sum(geno=="a/b"),sum(geno=="b/b"))
        hwetab[i,c("country","inversion")] <- c(C,I)
        hwetab[i,c("HWE")] <- Pval
        hwetab[i,c("MAF")] <- maf(ICcalls)
        hwetab[i,c("F")] <- f(ICcalls)
        
      }
}
```

```
## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call
```

```r
for (P in unique(calltab$pop)) {
      #countrycalls <- calltab[calltab$country==C,]
      C <- unique(calltab[calltab$pop==P,"country"])
      for(I in invnames) {
        ICcalls <- calltab[calltab$pop==P,I]
        geno <- genotype(c("a/a","a/b","b/b")[calltab[calltab$pop==P,I]+1])
        if(nallele(geno)==2) {
          test <- HWE.test(geno)
          Pval <- test$test$p.value
        } else{
          Pval <- 1
        }
        i <- nrow(hwetab)+1
        hwetab[i,c("aa","ab","bb")] <- c(sum(geno=="a/a"),sum(geno=="a/b"),sum(geno=="b/b"))
        hwetab[i,c("country","pop","inversion")] <- c(C,P,I)
        hwetab[i,c("HWE")] <- Pval
        hwetab[i,c("MAF")] <- maf(ICcalls)
        hwetab[i,c("F")] <- f(ICcalls)
        
      }
}
```

```
## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call
```

```
## Warning in diseq.ci(x, R = ci.B, conf = conf): One or more observed value outide of confidence interval. Check results.
```

```
## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call

## Warning in diseq.ci(x, R = ci.B, conf = conf): NAs returned from diseq call
```

```r
hwetabP <- subset(hwetab,HWE<0.01)
```


#make f/af cols with signif asterisks for paper

```r
hwetab$f <- round(hwetab$F,2)
hwetab$f[is.nan(hwetab$f)] <- "-"
hwetab$f[hwetab$HWE<0.05] <- paste(hwetab$f[hwetab$HWE<0.05],"*",sep="")
hwetab$f[hwetab$HWE<0.01] <- paste(hwetab$f[hwetab$HWE<0.01],"*",sep="")
hwetab$AF <- hwetab$MAF

#add popn sizes to table
popsizes <- table(calltab$pop)
cntsizes <- table(calltab$country)
hwetab$N <- popsizes[as.character(hwetab$pop)]
hwetab$N[is.na(hwetab$pop)] <- cntsizes[as.character(hwetab$country[is.na(hwetab$pop)])]


hwetabCP <- pivot_wider(hwetab[,c("inversion","country","pop","N","AF","f")],names_from = inversion,values_from = c(AF,f),names_glue = "{inversion}::{.value}",)
colorder <- c("country","pop","N",colnames(hwetabCP)[order(colnames(hwetabCP)[4:ncol(hwetabCP)])+3])
hwetabCP$country <- factor(hwetabCP$country,levels=cntorder,ordered=T)
hwetabCP <- hwetabCP[order(hwetabCP$country),colorder]
hwetabCP
```

```
## # A tibble: 102 × 23
##    country     pop         N       `X1_Gabon_1::AF` `X1_Gabon_1::f` `X1_global_1::AF` `X1_global_1::f` `X1_PuertoRico_1…` `X1_PuertoRico…`
##    <ord>       <fct>       <table>            <dbl> <chr>                       <dbl> <chr>                         <dbl> <chr>           
##  1 SouthAfrica <NA>         11                 0    -                            0.18 1.23                           0.23 1.28            
##  2 SouthAfrica Skukusa      11                 0    -                            0.18 1.23                           0.23 1.28            
##  3 Kenya       <NA>        132                 0.26 0.37**                       0.27 1                              0.03 1.17            
##  4 Kenya       Arabuko      18                 0    -                            0.36 1.33                           0.06 0.99            
##  5 Kenya       Ganda         8                 0.06 1.11                         0.19 1.22                           0    -               
##  6 Kenya       Kakamega     18                 0.83 0.79                         0.5  0.89                           0    -               
##  7 Kenya       KayaBomu     17                 0.12 0.56                         0.15 1.15                           0    -               
##  8 Kenya       Kwale        19                 0.11 1.08                         0.26 1.37                           0.05 1.11            
##  9 Kenya       Rabai_2017   17                 0.12 1.11                         0.35 1.03                           0.15 1.15            
## 10 Kenya       ShimbaHills   8                 0    -                            0    -                              0    -               
## # … with 92 more rows, and 14 more variables: `X1_PuertoRico_3::AF` <dbl>, `X1_PuertoRico_3::f` <chr>, `X1_Senegal_1::AF` <dbl>,
## #   `X1_Senegal_1::f` <chr>, `X1_Uganda_2::AF` <dbl>, `X1_Uganda_2::f` <chr>, `X1_Wafrica_1::AF` <dbl>, `X1_Wafrica_1::f` <chr>,
## #   `X2_BurkinaFaso_1::AF` <dbl>, `X2_BurkinaFaso_1::f` <chr>, `X2_SaudiArabia_1::AF` <dbl>, `X2_SaudiArabia_1::f` <chr>,
## #   `X2_Senegal_1::AF` <dbl>, `X2_Senegal_1::f` <chr>
```

```r
write.table(hwetabCP,file = "inversion_freqs_cnt_pops.txt",sep="\t",col.names = T,row.names = F,quote=F)
```


#get inversion locations

```r
blocks <- ldply(as.list(list.files("./lostruct_merge","*blocks.txt",full.names = T)),read_tsv)
```

```
## Rows: 233 Columns: 7
## ── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (5): inv, block, region, id, cluster
## dbl (2): chrom, pos
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
## Rows: 129 Columns: 7
## ── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (5): inv, block, region, id, cluster
## dbl (2): chrom, pos
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
## Rows: 0 Columns: 0
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
blocks$inversion <- paste("X",blocks$cluster,sep="")
blocks <- blocks[blocks$inversion %in% invnames,]

blockextents <- blocks %>% dplyr::group_by(inversion) %>% 
                           dplyr::summarise(chrom=min(chrom),
                                            start=min(pos),
                                            end=max(pos)+0.5e6,
                                            size=round((0.5e6+max(pos)-min(pos))/1e6,2)) %>% 
                            mutate(region=paste(chrom,":",as.integer(start),"-",as.integer(end),sep="")) %>%
                            dplyr::relocate(region,.before = size)

blockextents
```

```
## # A tibble: 10 × 6
##    inversion        chrom     start       end region                 size
##    <chr>            <dbl>     <dbl>     <dbl> <chr>                 <dbl>
##  1 X1_Gabon_1           1   1500000  17000000 1:1500000-17000000     15.5
##  2 X1_global_1          1 148000000 160500000 1:148000000-160500000  12.5
##  3 X1_PuertoRico_1      1 107000000 109500000 1:107000000-109500000   2.5
##  4 X1_PuertoRico_3      1 248500000 252000000 1:248500000-252000000   3.5
##  5 X1_Senegal_1         1 264500000 309000000 1:264500000-309000000  44.5
##  6 X1_Uganda_2          1 302500000 308000000 1:302500000-308000000   5.5
##  7 X1_Wafrica_1         1 162000000 170500000 1:162000000-170500000   8.5
##  8 X2_BurkinaFaso_1     2 464500000 473000000 2:464500000-473000000   8.5
##  9 X2_SaudiArabia_1     2 227000000 239000000 2:227000000-239000000  12  
## 10 X2_Senegal_1         2 306500000 317000000 2:306500000-317000000  10.5
```



```r
maf <- function(x) {round(sum(x) / (length(x)*2),3)}

calltab$african <- c("RoW","Africa")[+1]

regioncount <- data.frame("inversion"=factor(levels=invnames),
           "f"=numeric(),
           "f_Afr"=numeric(),
           "f_Row"=numeric()
           )

for(I in invnames) {
  ICcalls <- 
  i <- nrow(regioncount)+1
  regioncount[i,c("inversion")] <- I
  regioncount[i,c("f")] <- maf(calltab[,I])
  regioncount[i,c("f_Afr")] <- maf(calltab[calltab$region %in% c("West Africa","East Africa"),I])
  regioncount[i,c("f_Row")] <- maf(calltab[!calltab$region %in% c("West Africa","East Africa"),I])
  
}

regioncount <- merge(blockextents[,c("inversion","region","size")],regioncount)
write.table(regioncount,file = "inversion_extents_freqs.txt",sep="\t",col.names = T,row.names = F,quote=F)
```
