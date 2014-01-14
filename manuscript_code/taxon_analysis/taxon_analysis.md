# Personal Microbial Cloud

## Taxon Analysis

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





Quick refrerence vector of names for figures and whatnot. Some genus, families, etc are empty, so search for the finest taxonomic level and use that name. `cons` is in functions script sourced above. 


```r
consensus <- apply(taxo[, 1:6], 1, cons)
consensus[1:10]
```

```
##                     2                     6                    12 
##      "Halomonadaceae" "Solirubrobacterales"         "Rhizobiales" 
##                    21                    24                    25 
##            "Bacteria" "Gammaproteobacteria"            "Bacteria" 
##                    35                    73                   120 
##          "Bacillales"  "Betaproteobacteria"            "Bacteria" 
##                   122 
##           "Ralstonia"
```



Generate a list of the 10 most abundant OTUs overall, but also tack on OTUs that are among the 10 most abundant in any particular group of samples. Then a matching version of the mapping table. 


```r
tops <- names(rev(sort(colSums(pb.1000))))[1:10]
for (i in 1:length(gr2)) {
    tops.plus <- names(rev(sort(colSums(pb.1000[gr2[[i]], ]))))[1:10]
    tops <- unique(c(tops, tops.plus))
}
tops <- tops[rev(order(taxo[tops, 7]))]
taxo.tops <- taxo[tops, ]
length(tops)
```

```
## [1] 36
```

```r
consensus[tops]
```

```
##               234992                25506               121358 
##   "Methylobacterium"     "Staphylococcus"       "Tumebacillus" 
##               201443                94008                80526 
##   "Stenotrophomonas"    "Corynebacterium"    "Corynebacterium" 
##               151122                72518               136033 
## "Enterobacteriaceae"    "Corynebacterium"      "Acinetobacter" 
##               129988               172108               199152 
##    "Corynebacterium"        "Pseudomonas"      "Lactobacillus" 
##               129498               134501               229364 
##      "Streptococcus"          "Comamonas"        "Micrococcus" 
##               204422                57760                41465 
##       "Sphingomonas"  "Sphingomonadaceae"     "Dolosigranulum" 
##                63496               235029               205036 
##    "Actinomycetales"            "Kocuria"         "Finegoldia" 
##               164584               189711                92831 
##     "Alcaligenaceae"        "Cupriavidus"            "Pantoea" 
##                  264               133239               196388 
##   "Methylobacterium"            "Gemella"        "Leuconostoc" 
##               100717               187238                58876 
##          "Neisseria"           "Bacteria"          "Aeromonas" 
##                 4557               200842               102650 
##    "Cloacibacterium"  "Flavobacteriaceae"        "Rhizobiales" 
##                55224               169608               102857 
##        "Rhizobiales"         "Inquilinus"         "Pedobacter"
```



Create a reordered map and otu table. 


```r
pb.order <- data.frame(pb.1000[unlist(gr2), tops])
pb.order$gr.num <- c(rep(gr.name, c(unlist(lapply(gr2, length)))))
map.order <- pb.map[row.names(pb.order), ]
pb.order$location <- factor(map.order$location, levels = c("occ", "unocc"))
```



Make new trimmed OTU table and assess with boxplots. Then send the statistically interesting otus to figures and tables script. 


```r
pb.order <- pb.order[map.order$duration == "4", ]
# for (i in 1:(ncol(otus.to.barplot.tmp)-2)) {
# boxplot(otus.to.barplot.tmp[, i] ~ factor(otus.to.barplot.tmp$gr.num),
# las=2, main=consensus[gsub('X', '', names(otus.to.barplot.tmp)[i])],
# sub=i) readline('Press return for next plot') }
pb.order <- pb.order[, c(1, 2, 3, 4, 5, 6, 8, 10, 12, 18, 21, ncol(pb.order) - 
    1, ncol(pb.order))]
names(pb.order) <- gsub("X", "", names(pb.order))
consensus[names(pb.order)]
```

```
##             234992              25506             121358 
## "Methylobacterium"   "Staphylococcus"     "Tumebacillus" 
##             201443              94008              80526 
## "Stenotrophomonas"  "Corynebacterium"  "Corynebacterium" 
##              72518             129988             199152 
##  "Corynebacterium"  "Corynebacterium"    "Lactobacillus" 
##              41465             205036               <NA> 
##   "Dolosigranulum"       "Finegoldia"                 NA 
##               <NA> 
##                 NA
```

```r
rm(cons, Evenness, gr.name, gr2, tops.plus, makeTaxo)
```




Now for indicator analysis - this takes a long time and lots of memory! This uses the `indval()` function in the `labdsv` package. This gives output for all OTUs, the vast majority of which are unimportant to analysis. So there are a few steps to cut the output down to just the statistically salient OTUs. 


Unpack list of groups and pull out just 4 hour air filter dataset. 


```r
for (i in 1:length(names(groups))) {
    assign(names(groups)[i], unlist(groups[[i]]))
}

occ4f <- c(s14of, s24of, s34of, s14uf, s24uf, s34uf)
occ4f.table <- pb.1000[occ4f, ]
occ4f.table <- occ4f.table[, -which(colSums(occ4f.table) == 0)]
occ4f.map <- pb.map[occ4f, ]
```



Split this into two chunks - each takes a long time. 


```r
iv.location <- indval(occ4f.table, occ4f.map$location)
```



```r
iv.person <- indval(occ4f.table, occ4f.map$occ.person2)
```


levels(occ4f.map$occ.person2)


Keep those that are significant before p-value adjustment. 


```r
p <- 0.05

these <- which(iv.location$pval <= p)
iv.location.table <- data.frame(iv.location$maxcls[these], round(iv.location$indcls[these], 
    4), iv.location$pval[these])
names(iv.location.table) <- c("cluster", "indicator_value", "probability")
iv.location.table$cluster[iv.location.table$cluster == "1"] <- "occupied"
iv.location.table$cluster[iv.location.table$cluster == "2"] <- "unoccupied"

these <- which(iv.person$pval <= p)
iv.person.table <- data.frame(iv.person$maxcls[these], round(iv.person$indcls[these], 
    4), iv.person$pval[these])
names(iv.person.table) <- c("cluster", "indicator_value", "probability")
iv.person.table$cluster[iv.person.table$cluster == "1"] <- "Subject1"
iv.person.table$cluster[iv.person.table$cluster == "2"] <- "Subject2"
iv.person.table$cluster[iv.person.table$cluster == "3"] <- "Subject3"
iv.person.table$cluster[iv.person.table$cluster == "4"] <- "unoccupied"
```




Order by indicator value and tack on taxonomic information. 


```r
ivL <- iv.location.table[rev(order(iv.location.table$indicator_value)), ]
ivP <- iv.person.table[rev(order(iv.person.table$indicator_value)), ]

head(ivP)
```

```
##          cluster indicator_value probability
## X41465  Subject1          0.9620       0.001
## X199152 Subject3          0.8948       0.001
## X154708 Subject3          0.8193       0.001
## X63496  Subject1          0.7999       0.001
## X172036 Subject1          0.7554       0.001
## X186045 Subject3          0.7500       0.001
```

```r
ivL.names <- gsub("X", "", row.names(ivL))
ivL <- data.frame(cbind(ivL[, ], taxo[ivL.names, ]))

ivP.names <- gsub("X", "", row.names(ivP))
ivP <- cbind(ivP[, ], taxo[ivP.names, ])

ivL <- ivL[, c(2, 3, 10, 1, 9, 8, 7, 6, 5, 4)]
ivP <- ivP[, c(2, 3, 10, 1, 9, 8, 7, 6, 5, 4)]
```



Print out for later reference. 


```r
print(xtable(head(ivL, 20), digits = 3), type = "html")
```

<!-- html table generated in R 3.0.1 by xtable 1.7-1 package -->
<!-- Sun Jan  5 15:55:44 2014 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> indicator_value </TH> <TH> probability </TH> <TH> abundance </TH> <TH> cluster </TH> <TH> genus </TH> <TH> family </TH> <TH> order </TH> <TH> class </TH> <TH> phylum </TH> <TH> kingdom </TH>  </TR>
  <TR> <TD align="right"> X94008 </TD> <TD align="right"> 0.830 </TD> <TD align="right"> 0.001 </TD> <TD align="right"> 10134.000 </TD> <TD> occupied </TD> <TD> Corynebacterium </TD> <TD> Corynebacteriaceae </TD> <TD> Actinomycetales </TD> <TD> Actinobacteria </TD> <TD> Actinobacteria </TD> <TD> Bacteria </TD> </TR>
  <TR> <TD align="right"> X72518 </TD> <TD align="right"> 0.795 </TD> <TD align="right"> 0.001 </TD> <TD align="right"> 5093.000 </TD> <TD> occupied </TD> <TD> Corynebacterium </TD> <TD> Corynebacteriaceae </TD> <TD> Actinomycetales </TD> <TD> Actinobacteria </TD> <TD> Actinobacteria </TD> <TD> Bacteria </TD> </TR>
  <TR> <TD align="right"> X25506 </TD> <TD align="right"> 0.770 </TD> <TD align="right"> 0.001 </TD> <TD align="right"> 19463.000 </TD> <TD> occupied </TD> <TD> Staphylococcus </TD> <TD> Staphylococcaceae </TD> <TD> Bacillales </TD> <TD> Bacilli </TD> <TD> Firmicutes </TD> <TD> Bacteria </TD> </TR>
  <TR> <TD align="right"> X80526 </TD> <TD align="right"> 0.758 </TD> <TD align="right"> 0.001 </TD> <TD align="right"> 8760.000 </TD> <TD> occupied </TD> <TD> Corynebacterium </TD> <TD> Corynebacteriaceae </TD> <TD> Actinomycetales </TD> <TD> Actinobacteria </TD> <TD> Actinobacteria </TD> <TD> Bacteria </TD> </TR>
  <TR> <TD align="right"> X201443 </TD> <TD align="right"> 0.751 </TD> <TD align="right"> 0.001 </TD> <TD align="right"> 12640.000 </TD> <TD> unoccupied </TD> <TD> Stenotrophomonas </TD> <TD> Xanthomonadaceae </TD> <TD> Xanthomonadales </TD> <TD> Gammaproteobacteria </TD> <TD> Proteobacteria </TD> <TD> Bacteria </TD> </TR>
  <TR> <TD align="right"> X129498 </TD> <TD align="right"> 0.692 </TD> <TD align="right"> 0.001 </TD> <TD align="right"> 2393.000 </TD> <TD> occupied </TD> <TD> Streptococcus </TD> <TD> Streptococcaceae </TD> <TD> Lactobacillales </TD> <TD> Bacilli </TD> <TD> Firmicutes </TD> <TD> Bacteria </TD> </TR>
  <TR> <TD align="right"> X129988 </TD> <TD align="right"> 0.650 </TD> <TD align="right"> 0.010 </TD> <TD align="right"> 2893.000 </TD> <TD> occupied </TD> <TD> Corynebacterium </TD> <TD> Corynebacteriaceae </TD> <TD> Actinomycetales </TD> <TD> Actinobacteria </TD> <TD> Actinobacteria </TD> <TD> Bacteria </TD> </TR>
  <TR> <TD align="right"> X151122 </TD> <TD align="right"> 0.619 </TD> <TD align="right"> 0.006 </TD> <TD align="right"> 7579.000 </TD> <TD> unoccupied </TD> <TD>  </TD> <TD> Enterobacteriaceae </TD> <TD> Enterobacteriales </TD> <TD> Gammaproteobacteria </TD> <TD> Proteobacteria </TD> <TD> Bacteria </TD> </TR>
  <TR> <TD align="right"> X121358 </TD> <TD align="right"> 0.604 </TD> <TD align="right"> 0.035 </TD> <TD align="right"> 16283.000 </TD> <TD> unoccupied </TD> <TD> Tumebacillus </TD> <TD> Bacillaceae </TD> <TD> Bacillales </TD> <TD> Bacilli </TD> <TD> Firmicutes </TD> <TD> Bacteria </TD> </TR>
  <TR> <TD align="right"> X13338 </TD> <TD align="right"> 0.556 </TD> <TD align="right"> 0.001 </TD> <TD align="right"> 266.000 </TD> <TD> occupied </TD> <TD> Anaerococcus </TD> <TD> Incertae Sedis XI </TD> <TD> Clostridiales </TD> <TD> Clostridia </TD> <TD> Firmicutes </TD> <TD> Bacteria </TD> </TR>
  <TR> <TD align="right"> X84956 </TD> <TD align="right"> 0.551 </TD> <TD align="right"> 0.001 </TD> <TD align="right"> 251.000 </TD> <TD> occupied </TD> <TD> Peptoniphilus </TD> <TD> Incertae Sedis XI </TD> <TD> Clostridiales </TD> <TD> Clostridia </TD> <TD> Firmicutes </TD> <TD> Bacteria </TD> </TR>
  <TR> <TD align="right"> X172949 </TD> <TD align="right"> 0.546 </TD> <TD align="right"> 0.001 </TD> <TD align="right"> 242.000 </TD> <TD> occupied </TD> <TD> Corynebacterium </TD> <TD> Corynebacteriaceae </TD> <TD> Actinomycetales </TD> <TD> Actinobacteria </TD> <TD> Actinobacteria </TD> <TD> Bacteria </TD> </TR>
  <TR> <TD align="right"> X234992 </TD> <TD align="right"> 0.544 </TD> <TD align="right"> 0.034 </TD> <TD align="right"> 28147.000 </TD> <TD> unoccupied </TD> <TD> Methylobacterium </TD> <TD> Methylobacteriaceae </TD> <TD> Rhizobiales </TD> <TD> Alphaproteobacteria </TD> <TD> Proteobacteria </TD> <TD> Bacteria </TD> </TR>
  <TR> <TD align="right"> X205036 </TD> <TD align="right"> 0.540 </TD> <TD align="right"> 0.001 </TD> <TD align="right"> 793.000 </TD> <TD> occupied </TD> <TD> Finegoldia </TD> <TD> Incertae Sedis XI </TD> <TD> Clostridiales </TD> <TD> Clostridia </TD> <TD> Firmicutes </TD> <TD> Bacteria </TD> </TR>
  <TR> <TD align="right"> X134501 </TD> <TD align="right"> 0.538 </TD> <TD align="right"> 0.003 </TD> <TD align="right"> 1560.000 </TD> <TD> unoccupied </TD> <TD> Comamonas </TD> <TD> Comamonadaceae </TD> <TD> Burkholderiales </TD> <TD> Betaproteobacteria </TD> <TD> Proteobacteria </TD> <TD> Bacteria </TD> </TR>
  <TR> <TD align="right"> X233136 </TD> <TD align="right"> 0.532 </TD> <TD align="right"> 0.001 </TD> <TD align="right"> 277.000 </TD> <TD> unoccupied </TD> <TD> Halomonas </TD> <TD> Halomonadaceae </TD> <TD> Oceanospirillales </TD> <TD> Gammaproteobacteria </TD> <TD> Proteobacteria </TD> <TD> Bacteria </TD> </TR>
  <TR> <TD align="right"> X230275 </TD> <TD align="right"> 0.519 </TD> <TD align="right"> 0.001 </TD> <TD align="right"> 322.000 </TD> <TD> occupied </TD> <TD> Anaerococcus </TD> <TD> Incertae Sedis XI </TD> <TD> Clostridiales </TD> <TD> Clostridia </TD> <TD> Firmicutes </TD> <TD> Bacteria </TD> </TR>
  <TR> <TD align="right"> X8516 </TD> <TD align="right"> 0.505 </TD> <TD align="right"> 0.001 </TD> <TD align="right"> 337.000 </TD> <TD> occupied </TD> <TD> Anaerococcus </TD> <TD> Incertae Sedis XI </TD> <TD> Clostridiales </TD> <TD> Clostridia </TD> <TD> Firmicutes </TD> <TD> Bacteria </TD> </TR>
  <TR> <TD align="right"> X63496 </TD> <TD align="right"> 0.500 </TD> <TD align="right"> 0.001 </TD> <TD align="right"> 1061.000 </TD> <TD> occupied </TD> <TD>  </TD> <TD>  </TD> <TD> Actinomycetales </TD> <TD> Actinobacteria </TD> <TD> Actinobacteria </TD> <TD> Bacteria </TD> </TR>
  <TR> <TD align="right"> X199152 </TD> <TD align="right"> 0.489 </TD> <TD align="right"> 0.001 </TD> <TD align="right"> 2554.000 </TD> <TD> occupied </TD> <TD> Lactobacillus </TD> <TD> Lactobacillaceae </TD> <TD> Lactobacillales </TD> <TD> Bacilli </TD> <TD> Firmicutes </TD> <TD> Bacteria </TD> </TR>
   </TABLE>

```r
print(xtable(head(ivP, 20), digits = 3), type = "html")
```

<!-- html table generated in R 3.0.1 by xtable 1.7-1 package -->
<!-- Sun Jan  5 15:55:44 2014 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> indicator_value </TH> <TH> probability </TH> <TH> abundance </TH> <TH> cluster </TH> <TH> genus </TH> <TH> family </TH> <TH> order </TH> <TH> class </TH> <TH> phylum </TH> <TH> kingdom </TH>  </TR>
  <TR> <TD align="right"> X41465 </TD> <TD align="right"> 0.962 </TD> <TD align="right"> 0.001 </TD> <TD align="right"> 1089.000 </TD> <TD> Subject1 </TD> <TD> Dolosigranulum </TD> <TD> Carnobacteriaceae </TD> <TD> Lactobacillales </TD> <TD> Bacilli </TD> <TD> Firmicutes </TD> <TD> Bacteria </TD> </TR>
  <TR> <TD align="right"> X199152 </TD> <TD align="right"> 0.895 </TD> <TD align="right"> 0.001 </TD> <TD align="right"> 2554.000 </TD> <TD> Subject3 </TD> <TD> Lactobacillus </TD> <TD> Lactobacillaceae </TD> <TD> Lactobacillales </TD> <TD> Bacilli </TD> <TD> Firmicutes </TD> <TD> Bacteria </TD> </TR>
  <TR> <TD align="right"> X154708 </TD> <TD align="right"> 0.819 </TD> <TD align="right"> 0.001 </TD> <TD align="right"> 64.000 </TD> <TD> Subject3 </TD> <TD> Corynebacterium </TD> <TD> Corynebacteriaceae </TD> <TD> Actinomycetales </TD> <TD> Actinobacteria </TD> <TD> Actinobacteria </TD> <TD> Bacteria </TD> </TR>
  <TR> <TD align="right"> X63496 </TD> <TD align="right"> 0.800 </TD> <TD align="right"> 0.001 </TD> <TD align="right"> 1061.000 </TD> <TD> Subject1 </TD> <TD>  </TD> <TD>  </TD> <TD> Actinomycetales </TD> <TD> Actinobacteria </TD> <TD> Actinobacteria </TD> <TD> Bacteria </TD> </TR>
  <TR> <TD align="right"> X172036 </TD> <TD align="right"> 0.755 </TD> <TD align="right"> 0.001 </TD> <TD align="right"> 362.000 </TD> <TD> Subject1 </TD> <TD> Corynebacterium </TD> <TD> Corynebacteriaceae </TD> <TD> Actinomycetales </TD> <TD> Actinobacteria </TD> <TD> Actinobacteria </TD> <TD> Bacteria </TD> </TR>
  <TR> <TD align="right"> X186045 </TD> <TD align="right"> 0.750 </TD> <TD align="right"> 0.001 </TD> <TD align="right"> 154.000 </TD> <TD> Subject3 </TD> <TD> Facklamia </TD> <TD> Aerococcaceae </TD> <TD> Lactobacillales </TD> <TD> Bacilli </TD> <TD> Firmicutes </TD> <TD> Bacteria </TD> </TR>
  <TR> <TD align="right"> X8516 </TD> <TD align="right"> 0.744 </TD> <TD align="right"> 0.001 </TD> <TD align="right"> 337.000 </TD> <TD> Subject3 </TD> <TD> Anaerococcus </TD> <TD> Incertae Sedis XI </TD> <TD> Clostridiales </TD> <TD> Clostridia </TD> <TD> Firmicutes </TD> <TD> Bacteria </TD> </TR>
  <TR> <TD align="right"> X9535 </TD> <TD align="right"> 0.667 </TD> <TD align="right"> 0.001 </TD> <TD align="right"> 133.000 </TD> <TD> Subject1 </TD> <TD> Corynebacterium </TD> <TD> Corynebacteriaceae </TD> <TD> Actinomycetales </TD> <TD> Actinobacteria </TD> <TD> Actinobacteria </TD> <TD> Bacteria </TD> </TR>
  <TR> <TD align="right"> X175349 </TD> <TD align="right"> 0.665 </TD> <TD align="right"> 0.001 </TD> <TD align="right"> 99.000 </TD> <TD> Subject3 </TD> <TD> Peptoniphilus </TD> <TD> Incertae Sedis XI </TD> <TD> Clostridiales </TD> <TD> Clostridia </TD> <TD> Firmicutes </TD> <TD> Bacteria </TD> </TR>
  <TR> <TD align="right"> X108686 </TD> <TD align="right"> 0.652 </TD> <TD align="right"> 0.001 </TD> <TD align="right"> 205.000 </TD> <TD> Subject1 </TD> <TD> Corynebacterium </TD> <TD> Corynebacteriaceae </TD> <TD> Actinomycetales </TD> <TD> Actinobacteria </TD> <TD> Actinobacteria </TD> <TD> Bacteria </TD> </TR>
  <TR> <TD align="right"> X129988 </TD> <TD align="right"> 0.640 </TD> <TD align="right"> 0.001 </TD> <TD align="right"> 2893.000 </TD> <TD> Subject1 </TD> <TD> Corynebacterium </TD> <TD> Corynebacteriaceae </TD> <TD> Actinomycetales </TD> <TD> Actinobacteria </TD> <TD> Actinobacteria </TD> <TD> Bacteria </TD> </TR>
  <TR> <TD align="right"> X196388 </TD> <TD align="right"> 0.589 </TD> <TD align="right"> 0.001 </TD> <TD align="right"> 253.000 </TD> <TD> Subject2 </TD> <TD> Leuconostoc </TD> <TD> Leuconostocaceae </TD> <TD> Lactobacillales </TD> <TD> Bacilli </TD> <TD> Firmicutes </TD> <TD> Bacteria </TD> </TR>
  <TR> <TD align="right"> X205036 </TD> <TD align="right"> 0.583 </TD> <TD align="right"> 0.002 </TD> <TD align="right"> 793.000 </TD> <TD> Subject1 </TD> <TD> Finegoldia </TD> <TD> Incertae Sedis XI </TD> <TD> Clostridiales </TD> <TD> Clostridia </TD> <TD> Firmicutes </TD> <TD> Bacteria </TD> </TR>
  <TR> <TD align="right"> X217189 </TD> <TD align="right"> 0.576 </TD> <TD align="right"> 0.001 </TD> <TD align="right"> 37.000 </TD> <TD> Subject3 </TD> <TD>  </TD> <TD> Incertae Sedis XI </TD> <TD> Clostridiales </TD> <TD> Clostridia </TD> <TD> Firmicutes </TD> <TD> Bacteria </TD> </TR>
  <TR> <TD align="right"> X96641 </TD> <TD align="right"> 0.564 </TD> <TD align="right"> 0.001 </TD> <TD align="right"> 283.000 </TD> <TD> Subject1 </TD> <TD> Anaerococcus </TD> <TD> Incertae Sedis XI </TD> <TD> Clostridiales </TD> <TD> Clostridia </TD> <TD> Firmicutes </TD> <TD> Bacteria </TD> </TR>
  <TR> <TD align="right"> X63291 </TD> <TD align="right"> 0.557 </TD> <TD align="right"> 0.001 </TD> <TD align="right"> 70.000 </TD> <TD> Subject1 </TD> <TD> Corynebacterium </TD> <TD> Corynebacteriaceae </TD> <TD> Actinomycetales </TD> <TD> Actinobacteria </TD> <TD> Actinobacteria </TD> <TD> Bacteria </TD> </TR>
  <TR> <TD align="right"> X201443 </TD> <TD align="right"> 0.554 </TD> <TD align="right"> 0.002 </TD> <TD align="right"> 12640.000 </TD> <TD> unoccupied </TD> <TD> Stenotrophomonas </TD> <TD> Xanthomonadaceae </TD> <TD> Xanthomonadales </TD> <TD> Gammaproteobacteria </TD> <TD> Proteobacteria </TD> <TD> Bacteria </TD> </TR>
  <TR> <TD align="right"> X9604 </TD> <TD align="right"> 0.540 </TD> <TD align="right"> 0.001 </TD> <TD align="right"> 260.000 </TD> <TD> Subject1 </TD> <TD>  </TD> <TD> Aerococcaceae </TD> <TD> Lactobacillales </TD> <TD> Bacilli </TD> <TD> Firmicutes </TD> <TD> Bacteria </TD> </TR>
  <TR> <TD align="right"> X94008 </TD> <TD align="right"> 0.538 </TD> <TD align="right"> 0.001 </TD> <TD align="right"> 10134.000 </TD> <TD> Subject1 </TD> <TD> Corynebacterium </TD> <TD> Corynebacteriaceae </TD> <TD> Actinomycetales </TD> <TD> Actinobacteria </TD> <TD> Actinobacteria </TD> <TD> Bacteria </TD> </TR>
  <TR> <TD align="right"> X72518 </TD> <TD align="right"> 0.532 </TD> <TD align="right"> 0.001 </TD> <TD align="right"> 5093.000 </TD> <TD> Subject3 </TD> <TD> Corynebacterium </TD> <TD> Corynebacteriaceae </TD> <TD> Actinomycetales </TD> <TD> Actinobacteria </TD> <TD> Actinobacteria </TD> <TD> Bacteria </TD> </TR>
   </TABLE>



Clean up and save important bits to go to figures and tables script. 


```r
rm(list = ls()[which(ls() %in% names(groups))])
rm(occ4f, occ4f.table, occ4f.map, i, iv.location, iv.person, p, these, iv.location.table, 
    iv.person.table, ivL.names, ivP.names)
save(ivL, ivP, pb.order, taxo.tops, tops, consensus, file = "../data/taxon_analysis.RData")
```




