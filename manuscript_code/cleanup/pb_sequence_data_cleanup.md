# Personal Microbial Cloud

## Sequence Data Cleanup

* load ecology packages
* read in a few functions that will be used later


```r
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
source("../data/pb_functions_for_manuscript.R")
```



Load data: 

* Classic OTU table. Manually remove hash character from first line. This also has taxonomy from GreenGenes database as last field. 
* Mapping file used in QIIME. Also removed hash from first line and removed comment line. 
* And then jive the order of samples in the OTU table and mapping file.



```r
pb <- read.delim("../data/otu_table_r.txt", head = TRUE, row.names = 1)
pb.tax <- pb[, ncol(pb)]
pb <- pb[, -ncol(pb)]
pb <- t(pb)

pb.map <- read.delim("../data/pb_map.txt", head = TRUE, row.names = 1)

all(row.names(pb.map) %in% row.names(pb))
```

```
## [1] TRUE
```

```r
dim(pb.map)
```

```
## [1] 300  12
```

```r
dim(pb)
```

```
## [1]    300 235688
```

```r

pb.map <- pb.map[row.names(pb), ]
dim(pb)
```

```
## [1]    300 235688
```

```r
dim(pb.map)
```

```
## [1] 300  12
```

```r
length(pb.tax)
```

```
## [1] 235688
```

```r
identical(row.names(pb), row.names(pb.map))
```

```
## [1] TRUE
```



Keep original copy before changing things. 


```r
pb.original <- pb
pb.tax.original <- pb.tax
pb.map.original <- pb.map
```




```
# saved just in case of screwup. 
pb <- pb.original
pb.tax <- pb.tax.original
pb.map <- pb.map.original
```

Take out plant sequences


```r
plants <- grep("Streptophyt", pb.tax)
pb <- pb[, -plants]
pb.tax <- pb.tax[-plants]
rm(plants)
```




Set up taxonomy table


```r
taxo <- makeTaxo(taxo.in = pb.tax, otu.table = pb)
```

```
## Warning: no non-missing arguments to max; returning -Inf
```

```r
dim(taxo)
```

```
## [1] 234213      7
```

```r
dim(pb)
```

```
## [1]    300 234213
```


The big dataset started with 1.0028 &times; 10<sup>7</sup> total sequences from 234213 OTUs and 300 samples.


Identify those samples not used in the study. 


```r
nonexp <- which(pb.map$location == "east" | pb.map$location == "west" | pb.map$location == 
    "out")
pb <- pb[-nonexp, ]
pb.map <- pb.map[-nonexp, ]
rm(nonexp)
dim(pb)
```

```
## [1]    243 234213
```

```r
dim(pb.map)
```

```
## [1] 243  12
```




Index for controls, and clean up mapping table. 


```r
pb.map <- pb.map[, c("person", "location", "duration", "sampleType")]
control <- which(pb.map$location == "control")
pb <- pb[-c(control), ]
pb.map <- pb.map[-c(control), ]
dim(pb)
```

```
## [1]    216 234213
```

```r
dim(pb.map)
```

```
## [1] 216   4
```

```r
identical(row.names(pb), row.names(pb.map))
```

```
## [1] TRUE
```






Remove 4 contaminants that were introduced by PCR reagents. 
And remove control samples. 
This is detailed in `pb_analysis_Redo.Rmd`. 


```r
contaminants <- c("190873", "64285", "96443", "129992")
contaminants.index <- which(colnames(pb) %in% contaminants)
dim(pb)
```

```
## [1]    216 234213
```

```r
dim(pb.map)
```

```
## [1] 216   4
```

```r
pb <- pb[, -contaminants.index]
taxo <- taxo[-contaminants.index, ]  # pb
rm(control, contaminants, contaminants.index)
identical(row.names(pb.map), row.names(pb))
```

```
## [1] TRUE
```

```r
dim(pb.map)
```

```
## [1] 216   4
```

```r
dim(pb)
```

```
## [1]    216 234209
```

```r
dim(taxo)
```

```
## [1] 234209      7
```




Rarefy to 1000 seqs per sample. 


```r
sort(rowSums(pb))
```

```
## uf04.2h.s3 uf05.2h.s3 uf09.2h.s2 uf05.4h.s2 up04.4h.s2 of04.4h.s2 
##         14         30        319        445        622       1035 
## uf07.2h.s1 uf08.2h.s2 up05.2h.s3 uf09.4h.s2 uf11.4h.s1 uf12.2h.s3 
##       1168       1209       1264       1298       1509       1832 
## of01.2h.s2 up05.2h.s2 uf05.2h.s1 of11.2h.s3 uf10.4h.s3 uf12.4h.s2 
##       1981       2100       2132       2313       2361       2411 
## up02.2h.s2 of04.4h.s1 uf03.4h.s3 uf03.2h.s1 up02.4h.s3 uf03.2h.s2 
##       2414       2563       2577       2690       2725       2755 
## op03.2h.s3 up06.4h.s3 of02.2h.s2 of03.4h.s2 up05.4h.s2 of11.2h.s1 
##       2780       2896       3276       3339       3429       3544 
## up04.2h.s1 uf02.2h.s3 of06.2h.s3 uf06.4h.s1 of05.2h.s2 uf09.2h.s1 
##       3688       3716       3770       3813       3929       3956 
## uf09.2h.s3 up04.4h.s3 uf08.4h.s2 up03.2h.s1 uf06.2h.s3 uf01.4h.s1 
##       3974       4049       4254       4255       4340       4379 
## of09.4h.s2 uf08.2h.s1 of09.2h.s2 up04.2h.s2 op06.4h.s1 op05.2h.s1 
##       4398       4431       4633       4664       4743       4774 
## uf03.4h.s2 uf07.2h.s3 uf12.4h.s3 uf06.4h.s2 uf10.4h.s2 uf02.4h.s1 
##       4858       5012       5403       5463       5465       5507 
## uf12.2h.s1 uf04.4h.s3 of07.4h.s2 op06.2h.s3 of02.2h.s3 up06.4h.s1 
##       5553       5567       5695       5766       5859       6039 
## uf01.2h.s3 uf12.4h.s1 of03.2h.s2 uf03.4h.s1 uf08.2h.s3 op05.2h.s2 
##       6122       6242       6471       6544       6603       6699 
## uf06.2h.s1 of12.4h.s2 uf03.2h.s3 up01.4h.s3 op04.4h.s2 uf10.2h.s3 
##       6853       6966       7024       7045       7058       7083 
## uf04.2h.s1 up03.4h.s3 uf11.4h.s2 up06.2h.s1 of03.2h.s3 uf01.2h.s1 
##       7331       7417       7708       8239       8321       8341 
## of12.2h.s3 of06.2h.s2 op05.2h.s3 uf09.4h.s1 up01.2h.s1 up06.4h.s2 
##       8599       8667       8832       8924       8944       9123 
## op05.4h.s3 op06.2h.s1 of10.2h.s3 uf11.2h.s1 of03.2h.s1 up05.4h.s3 
##       9728       9836      10098      10105      10287      10455 
## up03.2h.s2 of04.2h.s2 up03.4h.s1 uf07.4h.s2 of04.4h.s3 of03.4h.s3 
##      10719      10737      10874      10982      11224      11229 
## uf02.4h.s3 of03.4h.s1 up03.2h.s3 of10.4h.s3 of01.2h.s3 uf08.4h.s1 
##      11361      11362      11409      11440      11496      11595 
## uf07.4h.s1 up01.4h.s2 op04.4h.s3 uf01.4h.s2 uf04.2h.s2 of10.2h.s1 
##      11680      11794      11943      12079      12268      12565 
## of08.2h.s2 up02.4h.s2 uf04.4h.s1 op03.4h.s2 uf10.2h.s1 uf11.2h.s2 
##      12570      12576      12827      12839      12993      13104 
## op04.2h.s2 uf01.2h.s2 uf11.2h.s3 uf07.2h.s2 uf04.4h.s2 uf10.2h.s2 
##      14264      14296      14360      14393      15082      15325 
## of05.2h.s3 uf09.4h.s3 of04.2h.s3 uf02.2h.s2 uf06.4h.s3 of08.4h.s2 
##      15617      15760      15851      15868      16233      16497 
## up02.2h.s1 up01.4h.s1 op03.4h.s3 of11.4h.s2 of06.2h.s1 of09.2h.s3 
##      16667      17823      17842      18619      18804      19083 
## op04.2h.s1 of07.2h.s1 op03.2h.s1 of02.4h.s2 up04.4h.s1 uf02.2h.s1 
##      19087      19664      19949      20004      20114      20435 
## of05.4h.s1 of08.2h.s3 up06.2h.s3 uf05.4h.s1 op01.2h.s2 op02.2h.s1 
##      20546      21763      21847      22022      22716      22824 
## op01.4h.s3 of04.2h.s1 op01.2h.s3 uf05.2h.s2 uf10.4h.s1 of06.4h.s2 
##      23185      23190      23312      23526      23916      24057 
## of05.2h.s1 uf02.4h.s2 of09.4h.s3 op02.4h.s1 of05.4h.s2 of07.2h.s3 
##      24064      24119      26054      26422      26816      27521 
## of12.4h.s1 uf05.4h.s3 of11.4h.s3 up02.2h.s3 op01.4h.s1 up01.2h.s3 
##      27893      27964      28939      29254      29271      29317 
## of12.2h.s2 uf12.2h.s2 op04.4h.s1 of05.4h.s3 uf06.2h.s2 uf08.4h.s3 
##      30192      30442      30820      30826      31470      33236 
## of02.4h.s1 of12.2h.s1 of02.2h.s1 up02.4h.s1 of10.4h.s1 up05.4h.s1 
##      34169      34801      36046      37192      37617      37713 
## op02.2h.s2 of01.4h.s1 uf01.4h.s3 up05.2h.s1 of08.2h.s1 of11.4h.s1 
##      38675      38721      39095      39258      40386      41980 
## up04.2h.s3 of01.2h.s1 uf11.4h.s3 op06.2h.s2 of07.2h.s2 of10.4h.s2 
##      43128      43367      44531      44663      44820      45242 
## up01.2h.s2 up06.2h.s2 op02.4h.s2 up03.4h.s2 of06.4h.s1 of08.4h.s3 
##      45419      45429      45633      45687      45882      47273 
## of07.4h.s1 of09.2h.s1 op04.2h.s3 of11.2h.s2 of10.2h.s2 of02.4h.s3 
##      47494      48523      48552      48748      49228      55685 
## of12.4h.s3 op05.4h.s1 of09.4h.s1 op06.4h.s2 of08.4h.s1 op01.2h.s1 
##      56417      58353      62013      64385      64649      69957 
## uf07.4h.s3 of07.4h.s3 of01.4h.s2 op03.4h.s1 op03.2h.s2 op01.4h.s2 
##      76818      78604      82810      87558      87843      90912 
## op05.4h.s2 op06.4h.s3 op02.2h.s3 of06.4h.s3 of01.4h.s3 op02.4h.s3 
##     100535     104632     113890     144827     173563     298953
```

```r
pb <- pb[-which(rowSums(pb) < 1000), ]
pb.1000 <- rrarefy(pb, 1000)
pb.1000 <- pb.1000[, -which(colSums(pb.1000) == 0)]
pb.map <- pb.map[row.names(pb.1000), ]
taxo <- taxo[colnames(pb.1000), ]
dim(pb.1000)
```

```
## [1]   211 19226
```

```r
dim(pb.map)
```

```
## [1] 211   4
```

```r
dim(taxo)
```

```
## [1] 19226     7
```




Put them all in alphabetical order, since that gives nice groupings. 
Yay for labelling schemes!


```r
pb.map <- pb.map[order(row.names(pb.map)), ]
pb.1000 <- pb.1000[row.names(pb.map), ]
```



Fix the abundance column - it still has abundance from full dataset, but abundance from rarefied dataset more meaningful. 


```r
identical(row.names(taxo), colnames(pb.1000))  # just to make sure
```

```
## [1] TRUE
```

```r
taxo$abundance <- colSums(pb.1000)
```



A few metrics to reference later.


```r
total.seqs <- sum(pb)
total.samples <- nrow(pb)
total.otus <- ncol(pb)

rarefied.seqs <- sum(pb.1000)
rarefied.samples <- nrow(pb.1000)
rarefied.otus <- ncol(pb.1000)
```





Clean up factor levels. 


```r
pb.map$person <- factor(pb.map$person)
pb.map$location <- factor(pb.map$location)
pb.map$duration <- factor(pb.map$duration)
pb.map$sampleType <- factor(pb.map$sampleType)
pb.map$occ.person <- factor(paste(pb.map$location, pb.map$person, sep = ""))

pb.map$occ.person2 <- as.character(pb.map$occ.person)
pb.map$occ.person2[pb.map$location == "unocc"] <- "unocc"
pb.map$occ.person2 <- factor(pb.map$occ.person2)
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



Indexing lists for quick parsing. 


```r
s22of <- c(1:nrow(pb.map))[pb.map$person == "s2" & pb.map$duration == 2 & pb.map$location == 
    "occ" & pb.map$sampleType == "filter"]
s24of <- c(1:nrow(pb.map))[pb.map$person == "s2" & pb.map$duration == 4 & pb.map$location == 
    "occ" & pb.map$sampleType == "filter"]
s22op <- c(1:nrow(pb.map))[pb.map$person == "s2" & pb.map$duration == 2 & pb.map$location == 
    "occ" & pb.map$sampleType == "petri.dish"]
s24op <- c(1:nrow(pb.map))[pb.map$person == "s2" & pb.map$duration == 4 & pb.map$location == 
    "occ" & pb.map$sampleType == "petri.dish"]
s22uf <- c(1:nrow(pb.map))[pb.map$person == "s2" & pb.map$duration == 2 & pb.map$location == 
    "unocc" & pb.map$sampleType == "filter"]
s24uf <- c(1:nrow(pb.map))[pb.map$person == "s2" & pb.map$duration == 4 & pb.map$location == 
    "unocc" & pb.map$sampleType == "filter"]
s22up <- c(1:nrow(pb.map))[pb.map$person == "s2" & pb.map$duration == 2 & pb.map$location == 
    "unocc" & pb.map$sampleType == "petri.dish"]
s24up <- c(1:nrow(pb.map))[pb.map$person == "s2" & pb.map$duration == 4 & pb.map$location == 
    "unocc" & pb.map$sampleType == "petri.dish"]

##
s32of <- c(1:nrow(pb.map))[pb.map$person == "s3" & pb.map$duration == 2 & pb.map$location == 
    "occ" & pb.map$sampleType == "filter"]
s34of <- c(1:nrow(pb.map))[pb.map$person == "s3" & pb.map$duration == 4 & pb.map$location == 
    "occ" & pb.map$sampleType == "filter"]
s32op <- c(1:nrow(pb.map))[pb.map$person == "s3" & pb.map$duration == 2 & pb.map$location == 
    "occ" & pb.map$sampleType == "petri.dish"]
s34op <- c(1:nrow(pb.map))[pb.map$person == "s3" & pb.map$duration == 4 & pb.map$location == 
    "occ" & pb.map$sampleType == "petri.dish"]
s32uf <- c(1:nrow(pb.map))[pb.map$person == "s3" & pb.map$duration == 2 & pb.map$location == 
    "unocc" & pb.map$sampleType == "filter"]
s34uf <- c(1:nrow(pb.map))[pb.map$person == "s3" & pb.map$duration == 4 & pb.map$location == 
    "unocc" & pb.map$sampleType == "filter"]
s32up <- c(1:nrow(pb.map))[pb.map$person == "s3" & pb.map$duration == 2 & pb.map$location == 
    "unocc" & pb.map$sampleType == "petri.dish"]
s34up <- c(1:nrow(pb.map))[pb.map$person == "s3" & pb.map$duration == 4 & pb.map$location == 
    "unocc" & pb.map$sampleType == "petri.dish"]

##
s12of <- c(1:nrow(pb.map))[pb.map$person == "s1" & pb.map$duration == 2 & pb.map$location == 
    "occ" & pb.map$sampleType == "filter"]
s14of <- c(1:nrow(pb.map))[pb.map$person == "s1" & pb.map$duration == 4 & pb.map$location == 
    "occ" & pb.map$sampleType == "filter"]
s12op <- c(1:nrow(pb.map))[pb.map$person == "s1" & pb.map$duration == 2 & pb.map$location == 
    "occ" & pb.map$sampleType == "petri.dish"]
s14op <- c(1:nrow(pb.map))[pb.map$person == "s1" & pb.map$duration == 4 & pb.map$location == 
    "occ" & pb.map$sampleType == "petri.dish"]
s12uf <- c(1:nrow(pb.map))[pb.map$person == "s1" & pb.map$duration == 2 & pb.map$location == 
    "unocc" & pb.map$sampleType == "filter"]
s14uf <- c(1:nrow(pb.map))[pb.map$person == "s1" & pb.map$duration == 4 & pb.map$location == 
    "unocc" & pb.map$sampleType == "filter"]
s12up <- c(1:nrow(pb.map))[pb.map$person == "s1" & pb.map$duration == 2 & pb.map$location == 
    "unocc" & pb.map$sampleType == "petri.dish"]
s14up <- c(1:nrow(pb.map))[pb.map$person == "s1" & pb.map$duration == 4 & pb.map$location == 
    "unocc" & pb.map$sampleType == "petri.dish"]


## A few groupings to call up later

groups <- list(s22of, s32of, s12of, s22op, s32op, s12op, s22uf, s32uf, s12uf, 
    s22up, s32up, s12up, s24of, s34of, s14of, s24op, s34op, s14op, s24uf, s34uf, 
    s14uf, s24up, s34up, s14up)

names(groups) <- c("s22of", "s32of", "s12of", "s22op", "s32op", "s12op", "s22uf", 
    "s32uf", "s12uf", "s22up", "s32up", "s12up", "s24of", "s34of", "s14of", 
    "s24op", "s34op", "s14op", "s24uf", "s34uf", "s14uf", "s24up", "s34up", 
    "s14up")

gr2 <- list(s14of, s12of, s14op, s12op, s14uf, s12uf, s14up, s12up, s34of, s32of, 
    s34op, s32op, s34uf, s32uf, s34up, s32up, s24of, s22of, s24op, s22op, s24uf, 
    s22uf, s24up, s22up)

gr.name <- c("s14of", "s12of", "s14op", "s12op", "s14uf", "s12uf", "s14up", 
    "s12up", "s34of", "s32of", "s34op", "s32op", "s34uf", "s32uf", "s34up", 
    "s32up", "s24of", "s22of", "s24op", "s22op", "s24uf", "s22uf", "s24up", 
    "s22up")

# clean up group names
rm(list = ls()[which(ls() %in% names(groups))])
```





Now save the image that can be used for analysis and figures. 


```r
rm(pb.original, pb.map.original, pb.tax.original, pb, pb.tax)
save.image("../data/pb_data.RData")
ls()
```

```
##  [1] "cons"             "Evenness"         "gr.name"         
##  [4] "gr2"              "groups"           "makeTaxo"        
##  [7] "pb.1000"          "pb.map"           "rarefied.otus"   
## [10] "rarefied.samples" "rarefied.seqs"    "se"              
## [13] "taxo"             "total.otus"       "total.samples"   
## [16] "total.seqs"
```





