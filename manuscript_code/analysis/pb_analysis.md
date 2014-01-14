# Personal Microbial Cloud

## Beta-diversity Data Analysis

Read in prepared workspace and load packages.


```r
load("../data/pb_data.RData")
library(vegan)
```

```
## Loading required package: permute This is vegan 2.0-8
```

```r
library(labdsv)
```

```
## Loading required package: mgcv This is mgcv 1.7-22. For overview type
## 'help("mgcv-package")'. Loading required package: MASS
## 
## Attaching package: 'labdsv'
## 
## The following object is masked from 'package:stats':
## 
## density
```

```r
library(ape)
library(xtable)
source("../data/pb_functions_for_manuscript.R")
```



Unpack list of groups.


```r
for (i in 1:length(names(groups))) {
    assign(names(groups)[i], unlist(groups[[i]]))
}
rm(i)
```



How do samples fall into groups? I.e., do we have balance in treatments?


```r
table(pb.map$person, pb.map$sampleType)
```

```
##     
##      filter petri.dish
##   s1     48         24
##   s2     46         23
##   s3     46         24
```

```r
table(pb.map$location)
```

```
## 
##   occ unocc 
##   108   103
```

```r
table(pb.map$duration)
```

```
## 
##   2   4 
## 105 106
```

```r
table(pb.map$sampleType)
```

```
## 
##     filter petri.dish 
##        140         71
```




The canberra distance metric is used to calculate beta diversity. Then 2d NMDS objects are created for the whole experiment and for each treatment individually.


```r
pb.can <- vegdist(pb.1000, "canberra")
```



Since most abundant taxa were removed as contaminants, try Bray-Curtis. 


```r
pb.bc <- vegdist(pb.1000, "bray")
```







```r
# Filters 4h diff people
occ4f <- c(s24of, s34of, s14of, s24uf, s34uf, s14uf)
occ4f.table <- pb.1000[occ4f, ]
occ4f.map <- pb.map[occ4f, ]
occ4f.can <- as.dist(as.matrix(pb.can)[occ4f, occ4f])
occ4f.nmds <- bestnmds(occ4f.can)
occ <- which(occ4f.map$location == "occ")
unocc <- which(occ4f.map$location == "unocc")

## 2h diff people
occ2f <- c(s22of, s32of, s12of, s22uf, s32uf, s12uf)
occ2f.table <- pb.1000[occ2f, ]
occ2f.map <- pb.map[occ2f, ]
occ2f.can <- as.dist(as.matrix(pb.can)[occ2f, occ2f])
occ2f.nmds <- bestnmds(occ2f.can)
occ2 <- which(occ2f.map$location == "occ")
unocc2 <- which(occ2f.map$location == "unocc")

# Petri dishes 4h diff people
occ4p <- c(s24op, s34op, s14op, s24up, s34up, s14up)
occ4p.table <- pb.1000[occ4p, ]
occ4p.map <- pb.map[occ4p, ]
occ4p.can <- as.dist(as.matrix(pb.can)[occ4p, occ4p])
occ4p.nmds <- bestnmds(occ4p.can)
occ.p <- which(occ4p.map$location == "occ")
unocc.p <- which(occ4p.map$location == "unocc")

## 2h diff people
occ2p <- c(s22op, s32op, s12op, s22up, s32up, s12up)
occ2p.table <- pb.1000[occ2p, ]
occ2p.map <- pb.map[occ2p, ]
occ2p.can <- as.dist(as.matrix(pb.can)[occ2p, occ2p])
occ2p.nmds <- bestnmds(occ2p.can)
occ2.p <- which(occ2p.map$location == "occ")
unocc2.p <- which(occ2p.map$location == "unocc")
```



Pretty much all of these NMDS ordinations have high stress levels:


```r

nmds.stress <- data.frame(treatment = c("4hr filters", "2hr filters", "4hr petri dishes", 
    "2hr petri dishes"), stress = c(occ4f.nmds$stress, occ2f.nmds$stress, occ4p.nmds$stress, 
    occ2p.nmds$stress), n = c(length(occ4f), length(occ2f), length(occ4p), length(occ2p)))
print(xtable(nmds.stress), type = "html")
```

<!-- html table generated in R 3.0.1 by xtable 1.7-1 package -->
<!-- Sun Jan  5 15:30:25 2014 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> treatment </TH> <TH> stress </TH> <TH> n </TH>  </TR>
  <TR> <TD align="right"> 1 </TD> <TD> 4hr filters </TD> <TD align="right"> 27.91 </TD> <TD align="right">  71 </TD> </TR>
  <TR> <TD align="right"> 2 </TD> <TD> 2hr filters </TD> <TD align="right"> 30.94 </TD> <TD align="right">  69 </TD> </TR>
  <TR> <TD align="right"> 3 </TD> <TD> 4hr petri dishes </TD> <TD align="right"> 26.00 </TD> <TD align="right">  35 </TD> </TR>
  <TR> <TD align="right"> 4 </TD> <TD> 2hr petri dishes </TD> <TD align="right"> 25.73 </TD> <TD align="right">  36 </TD> </TR>
   </TABLE>


But since statistics back up the visible patterns, it doesn't present much of a problem. 


The functions file (`pb_functions_for_mansucript.R`) contains a simple function to calculate richness, Shannon-Weiner diversity, and Evenness for each sample. 


```r
pb.hrj <- Evenness(pb.1000)

occ4f.map$R <- pb.hrj[occ4f, "R"]
occ4f.map$H1 <- pb.hrj[occ4f, "H1"]

occ2f.map$R <- pb.hrj[occ2f, "R"]
occ2f.map$H1 <- pb.hrj[occ2f, "H1"]

occ4p.map$R <- pb.hrj[occ4p, "R"]
occ4p.map$H1 <- pb.hrj[occ4p, "H1"]

occ2p.map$R <- pb.hrj[occ2p, "R"]
occ2p.map$H1 <- pb.hrj[occ2p, "H1"]
```



One last piece to assess the spread of points with beta-dispersion tests. It appears that s1 and s3 points are tightly grouped but s2 is not - this gives us a metric to report. 


```r
# 4of
pb.onlyocc.can <- as.dist(as.matrix(pb.can)[c(s24of, s34of, s14of), c(s24of, 
    s34of, s14of)])
pb.onlyocc.map <- pb.map[c(s24of, s34of, s14of), ]
occ4f.bd <- betadisper(pb.onlyocc.can, pb.onlyocc.map$person)
occ4f.betad <- data.frame(dists = occ4f.bd$distances, person = pb.onlyocc.map$person)
occ4f.betad.aov <- anova(occ4f.bd)
# boxplot(occ4f.bd$distances~pb.onlyocc.map$person)

# 2of
pb.onlyocc2.can <- as.dist(as.matrix(pb.can)[c(s22of, s32of, s12of), c(s22of, 
    s32of, s12of)])
pb.onlyocc2.map <- pb.map[c(s22of, s32of, s12of), ]
occ2f.bd <- betadisper(pb.onlyocc2.can, pb.onlyocc2.map$person)
occ2f.betad <- data.frame(dists = occ2f.bd$distances, person = pb.onlyocc2.map$person)
occ2f.betad.aov <- anova(occ2f.bd)
# boxplot(occ2f.bd$distances~pb.onlyocc2.map$person)

# 4op
pb.onlyoccp.can <- as.dist(as.matrix(pb.can)[c(s24op, s34op, s14op), c(s24op, 
    s34op, s14op)])
pb.onlyoccp.map <- pb.map[c(s24op, s34op, s14op), ]
occ4p.bd <- betadisper(pb.onlyoccp.can, pb.onlyoccp.map$person)
occ4p.betad <- data.frame(dists = occ4p.bd$distances, person = pb.onlyoccp.map$person)
occ4p.betad.aov <- anova(occ4p.bd)
# boxplot(occ4p.bd$distances~pb.onlyoccp.map$person)

# 2of
pb.onlyocc2p.can <- as.dist(as.matrix(pb.can)[c(s22op, s32op, s12op), c(s22op, 
    s32op, s12op)])
pb.onlyocc2p.map <- pb.map[c(s22op, s32op, s12op), ]
occ2p.bd <- betadisper(pb.onlyocc2p.can, pb.onlyocc2p.map$person)
occ2p.betad <- data.frame(dists = occ2p.bd$distances, person = pb.onlyocc2p.map$person)
occ2p.betad.aov <- anova(occ2p.bd)
# boxplot(occ2p.bd$distances~pb.onlyocc2p.map$person)

betad.aov.list <- list(occ4f.betad.aov, occ2f.betad.aov, occ4p.betad.aov, occ2p.betad.aov)
betad.table.list <- list(occ4f.betad, occ2f.betad, occ4p.betad, occ2p.betad)

rm(pb.onlyocc.can, pb.onlyocc.map, occ4f.bd, pb.onlyocc2.can, pb.onlyocc2.map, 
    occ2f.bd, pb.onlyoccp.can, pb.onlyoccp.map, occ4p.bd, pb.onlyocc2p.can, 
    pb.onlyocc2p.map, occ2p.bd, occ4f.betad, occ4f.betad.aov, occ4p.betad, occ4p.betad.aov, 
    occ2f.betad, occ2f.betad.aov, occ2p.betad, occ2p.betad.aov)
```



After all of these data manipulations, some of these objects now go to a new script to make figures.  


```r
ls()
```

```
##  [1] "betad.aov.list"   "betad.table.list" "cons"            
##  [4] "Evenness"         "gr.name"          "gr2"             
##  [7] "groups"           "makeTaxo"         "nmds.stress"     
## [10] "occ"              "occ.p"            "occ2"            
## [13] "occ2.p"           "occ2f"            "occ2f.can"       
## [16] "occ2f.map"        "occ2f.nmds"       "occ2f.table"     
## [19] "occ2p"            "occ2p.can"        "occ2p.map"       
## [22] "occ2p.nmds"       "occ2p.table"      "occ4f"           
## [25] "occ4f.can"        "occ4f.map"        "occ4f.nmds"      
## [28] "occ4f.table"      "occ4p"            "occ4p.can"       
## [31] "occ4p.map"        "occ4p.nmds"       "occ4p.table"     
## [34] "pb.1000"          "pb.bc"            "pb.can"          
## [37] "pb.hrj"           "pb.map"           "rarefied.otus"   
## [40] "rarefied.samples" "rarefied.seqs"    "s12of"           
## [43] "s12op"            "s12uf"            "s12up"           
## [46] "s14of"            "s14op"            "s14uf"           
## [49] "s14up"            "s22of"            "s22op"           
## [52] "s22uf"            "s22up"            "s24of"           
## [55] "s24op"            "s24uf"            "s24up"           
## [58] "s32of"            "s32op"            "s32uf"           
## [61] "s32up"            "s34of"            "s34op"           
## [64] "s34uf"            "s34up"            "se"              
## [67] "taxo"             "total.otus"       "total.samples"   
## [70] "total.seqs"       "unocc"            "unocc.p"         
## [73] "unocc2"           "unocc2.p"
```

```r

# clean groupings.
rm(list = ls()[which(ls() %in% names(groups))])
rm(occ, unocc, occ.p, unocc.p, occ2, unocc2, occ2.p, unocc2.p, occ2f, occ2p, 
    occ4f, occ4p)
ls()
```

```
##  [1] "betad.aov.list"   "betad.table.list" "cons"            
##  [4] "Evenness"         "gr.name"          "gr2"             
##  [7] "groups"           "makeTaxo"         "nmds.stress"     
## [10] "occ2f.can"        "occ2f.map"        "occ2f.nmds"      
## [13] "occ2f.table"      "occ2p.can"        "occ2p.map"       
## [16] "occ2p.nmds"       "occ2p.table"      "occ4f.can"       
## [19] "occ4f.map"        "occ4f.nmds"       "occ4f.table"     
## [22] "occ4p.can"        "occ4p.map"        "occ4p.nmds"      
## [25] "occ4p.table"      "pb.1000"          "pb.bc"           
## [28] "pb.can"           "pb.hrj"           "pb.map"          
## [31] "rarefied.otus"    "rarefied.samples" "rarefied.seqs"   
## [34] "se"               "taxo"             "total.otus"      
## [37] "total.samples"    "total.seqs"
```

```r
# save the rest.
save.image(file = "../data/pb_analysis_to_figures.RData")
getwd()
```

```
## [1] "/Users/jfmeadow/Dropbox/pb_manuscript/manuscript_code/analysis"
```









