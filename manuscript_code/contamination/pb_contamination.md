### Picklebox - look for contamination


Load data. 


```r
setwd("~/Dropbox/pb_shared_markdown/manuscript_code/contamination/")

library(vegan)
library(labdsv)

source("../data/pb_functions_for_manuscript.R")
```





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





Identify those not used in the study. 


```r
nonexp <- which(pb.map$location == "east" | pb.map$location == "west" | pb.map$location == 
    "out")
```



Keep original copy


```r
pb.original <- pb
pb.tax.original <- pb.tax
```



Take out plant sequences


```r
plants <- grep("Streptophyt", pb.tax)
pb <- pb[, -plants]
pb.tax <- pb.tax[-plants]  # remove from taxonomy to line up
```



Make a full rarefied dataset including controls and contaminants. 


```r
sort(rowSums(pb.original))
```

```
##  uf04.2h.ab  uf05.2h.ab  uf09.2h.aa  uf05.4h.aa    cef08.4h  uf07.2h.jm 
##          14          41         723         781        1522        1825 
##  uf11.4h.jm  up04.4h.aa  up05.2h.aa     extr.c3  of01.2h.aa  uf09.4h.aa 
##        1976        2520        2954        2980        3563        3624 
##  uf12.2h.ab  of11.2h.jm     extr.c4  of03.4h.aa    cep06.4h  of04.4h.jm 
##        4248        4310        4314        4343        4436        4544 
##     extr.c7  uf06.4h.jm  uf03.4h.ab    cwp02.4h  uf08.2h.aa  of11.2h.ab 
##        4805        4831        4932        5037        5191        5197 
##  neg.dish.4  uf08.2h.jm    extr.c11    cef03.4h  uf05.2h.jm  of09.2h.aa 
##        5260        5445        5452        5654        5697        6003 
## out02.2h.ab  uf08.4h.aa    extr.c10  of02.2h.aa  uf09.2h.ab  up04.2h.aa 
##        6241        6246        6253        6266        6389        6450 
##  uf03.2h.aa     extr.c1    cef11.4h    cwf01.4h    cwp01.4h     extr.c6 
##        6487        6519        6589        6717        6804        6920 
##  op05.2h.jm  uf03.2h.jm    cef09.4h  uf06.4h.aa  op06.4h.jm  up05.2h.ab 
##        6931        7057        7065        7139        7214        7247 
##  uf09.2h.jm  uf12.2h.jm  uf03.4h.aa  neg.swab.2  of07.4h.aa  uf12.4h.aa 
##        7420        7594        7691        7873        7885        7897 
##  uf01.2h.ab    cwf12.4h    cwf10.4h    extr.c12  uf03.4h.jm  neg.swab.3 
##        7945        7995        8132        8474        8581        8730 
##    cwf02.4h  neg.swab.4  up03.4h.ab  of05.2h.aa  neg.filt.4    cef12.4h 
##        8747        8772        8820        8980        9038        9241 
##    cep01.4h  up04.4h.ab  of09.4h.aa  of06.2h.ab  uf10.4h.ab  up02.4h.ab 
##        9539        9593        9798        9836       10006       10009 
## out01.2h.aa  up02.2h.aa    cep02.4h  op03.2h.ab  neg.filt.3  op06.2h.ab 
##       10018       10025       10199       10245       10533       10967 
##  up06.4h.jm    cwf07.4h  uf07.2h.ab  uf01.4h.jm  op04.4h.aa  of04.4h.aa 
##       11058       11096       11201       11275       11413       11651 
##  of03.2h.aa  uf11.2h.jm  uf03.2h.ab  neg.dish.3  uf10.4h.aa    cwp04.4h 
##       11746       11792       11877       12029       12067       12342 
##     extr.c5 out03.2h.jm  uf02.2h.ab  up05.4h.aa  uf10.2h.ab  up04.2h.jm 
##       12369       12409       12525       12736       12855       12920 
##    extr.c15  op05.4h.ab  up06.2h.jm  of03.4h.jm    extr.c13  uf02.4h.ab 
##       13047       13199       13262       13266       13319       13667 
##  uf02.4h.jm  uf06.2h.jm  up06.4h.ab  uf12.4h.jm  of10.2h.jm  up01.4h.ab 
##       14003       14050       14219       14321       14336       14460 
##  of10.2h.ab  of04.4h.ab  uf09.4h.jm  uf08.4h.jm    cep04.4h  uf04.4h.jm 
##       14484       14562       14784       14875       15011       15186 
##    cwf09.4h  of12.4h.aa  op05.2h.ab  of03.2h.jm out01.4h.jm out01.2h.jm 
##       15187       15593       15660       15663       15995       16003 
##  uf11.4h.aa  op04.2h.aa  uf01.4h.aa out02.2h.jm  of08.2h.aa  uf08.2h.ab 
##       16170       16334       16349       16354       16403       16425 
##  of02.2h.ab  of04.2h.aa  op06.2h.jm  uf01.2h.jm  of12.2h.ab  of10.4h.ab 
##       16627       16801       16933       16960       17098       17370 
##    cwf11.4h    cef02.4h  up03.4h.jm  up01.2h.jm  up06.4h.aa    cwp06.4h 
##       17528       17584       17988       18043       18173       18405 
##  neg.filt.2  of03.2h.ab out01.2h.ab  uf09.4h.ab  uf04.4h.ab  neg.dish.2 
##       18475       19093       19358       19378       19448       19501 
##  of03.4h.ab  op03.4h.aa  up03.2h.jm  op05.2h.aa  neg.dish.1  up01.4h.aa 
##       19629       19834       19851       19923       20392       20597 
##  uf10.2h.jm  up03.2h.ab  uf07.4h.jm    extr.c14  op03.4h.ab  uf10.2h.aa 
##       20620       20753       20805       21268       21729       21865 
##  op04.2h.jm  uf07.2h.aa  uf12.4h.ab  of06.2h.jm  up05.4h.ab  of01.2h.ab 
##       22021       22056       22110       22347       22419       22473 
##  of09.2h.ab  of08.4h.aa  uf11.2h.aa  op04.4h.ab  up03.2h.aa  uf11.2h.ab 
##       22768       22833       23313       23317       23579       23582 
##     extr.c9    cwf05.4h  of07.2h.jm  uf06.2h.ab  neg.filt.1  uf04.2h.jm 
##       23755       24256       24995       25019       25033       25113 
##  of06.2h.aa out01.4h.ab  op01.4h.ab  up02.2h.jm  of11.4h.aa  uf04.4h.aa 
##       25960       26001       26460       27156       27305       27412 
##  of05.2h.ab out03.2h.ab  of04.2h.jm  of09.4h.ab  of02.4h.aa  of05.4h.aa 
##       28213       28969       29019       29404       29841       30280 
##  up01.4h.jm  neg.swab.1 out02.4h.jm     extr.c8  of11.4h.ab  op01.2h.ab 
##       30325       30427       30760       31143       31996       32111 
##    cep05.4h  uf01.2h.aa  op02.2h.jm  of05.4h.jm  of12.4h.jm    cef01.4h 
##       32448       33241       34365       34840       35015       35256 
##    cwf04.4h    cef05.4h  uf07.4h.aa  of06.4h.aa  op03.2h.jm  of05.2h.jm 
##       35369       35792       36484       36688       36832       36981 
##  uf02.4h.aa  op01.2h.aa  of07.2h.ab  op02.4h.jm  of04.2h.ab  uf10.4h.jm 
##       37003       37791       38035       38617       39740       40604 
##  op01.4h.jm c.out.03.4h  of05.4h.ab    cef07.4h c.out.01.4h  of02.4h.jm 
##       40741       41508       41515       42017       42363       42711 
##    cep03.4h    cef04.4h  up05.2h.jm  op04.4h.jm  up02.2h.ab  of12.2h.jm 
##       42888       43308       43609       43673       43739       43743 
##  of12.2h.aa  uf04.2h.aa  op02.4h.aa  up06.2h.ab  of10.4h.jm  up04.4h.jm 
##       46188       47668       47853       47932       48189       48228 
##  up02.4h.aa  of11.4h.jm  of08.2h.jm out03.4h.aa out01.4h.aa  uf02.2h.jm 
##       48508       48928       49057       49166       50494       50814 
##  uf05.4h.jm  uf12.2h.aa    cwp03.4h  uf06.2h.aa    cwf06.4h  op02.2h.aa 
##       50884       51119       51157       51581       51733       51818 
##  of08.2h.ab out02.4h.aa  uf05.2h.aa  uf02.2h.aa  up06.2h.aa  of02.2h.jm 
##       52636       52806       53026       53420       53830       54590 
##  up01.2h.ab out03.4h.ab  of01.2h.jm  of10.4h.aa  of08.4h.ab    cef10.4h 
##       55121       56066       56205       58669       60012       60368 
##  up05.4h.jm  of07.4h.jm  uf05.4h.ab  of10.2h.aa  of02.4h.ab    cwf08.4h 
##       60551       61009       61050       61616       62069       62276 
##    cwp05.4h  of07.2h.aa  of11.2h.aa  uf06.4h.ab  op05.4h.jm  of01.4h.jm 
##       62606       63211       63559       64060       64454       65531 
##  op04.2h.ab  of09.2h.jm  of12.4h.ab  op06.2h.aa  uf01.4h.ab  up02.4h.jm 
##       68925       69222       71576       71693       73705       73726 
##  up03.4h.aa  of06.4h.jm  uf08.4h.ab c.out.02.4h  up04.2h.ab  uf11.4h.ab 
##       73920       73954       74676       76285       78288       80533 
## out03.4h.jm  op06.4h.aa out02.2h.aa  of07.4h.ab  up01.2h.aa out02.4h.ab 
##       80593       85142       85418       85503       87362       89350 
##  of09.4h.jm  op01.4h.aa     extr.c2    cef06.4h  op01.2h.jm  of08.4h.jm 
##       89773       95472      100420      103197      104357      108328 
##  op03.4h.jm  of01.4h.aa  op03.2h.aa  uf07.4h.ab  op02.2h.ab out03.2h.aa 
##      108847      113924      116025      117705      119490      121248 
##  op06.4h.ab  op05.4h.aa  of06.4h.ab  of01.4h.ab    cwf03.4h  op02.4h.ab 
##      131677      142861      186832      205424      235568      311354
```

```r
sort(rowSums(pb))
```

```
##  uf04.2h.ab  uf05.2h.ab  uf09.2h.aa  uf05.4h.aa    cef08.4h  uf07.2h.jm 
##          14          41         723         781        1522        1825 
##  uf11.4h.jm  up04.4h.aa  up05.2h.aa     extr.c3  of01.2h.aa  uf09.4h.aa 
##        1976        2520        2951        2954        3563        3624 
##  uf12.2h.ab  of11.2h.jm     extr.c4  of03.4h.aa    cep06.4h  of04.4h.jm 
##        4248        4310        4314        4332        4436        4516 
##     extr.c7  uf06.4h.jm  uf03.4h.ab    cwp02.4h  of11.2h.ab  uf08.2h.aa 
##        4805        4831        4932        5015        5188        5191 
##  neg.dish.4  uf08.2h.jm    extr.c11    cef03.4h  uf05.2h.jm  of09.2h.aa 
##        5256        5379        5452        5654        5692        5981 
##  of02.2h.aa  uf08.4h.aa out02.2h.ab    extr.c10  uf09.2h.ab  uf03.2h.aa 
##        6217        6238        6241        6253        6350        6391 
##  up04.2h.aa    cef11.4h     extr.c1    cwf01.4h  op05.2h.jm    cwp01.4h 
##        6450        6471        6518        6717        6751        6804 
##  op06.4h.jm     extr.c6  uf03.2h.jm    cef09.4h  uf06.4h.aa  up05.2h.ab 
##        6906        6920        7057        7057        7087        7247 
##  uf09.2h.jm  uf12.2h.jm  uf03.4h.aa  neg.swab.2  of07.4h.aa  uf12.4h.aa 
##        7387        7508        7634        7873        7873        7897 
##  uf01.2h.ab    cwf12.4h    cwf10.4h    extr.c12  neg.swab.4  uf03.4h.jm 
##        7914        7995        8131        8474        8478        8518 
##  neg.swab.3    cwf02.4h  up03.4h.ab  of05.2h.aa  neg.filt.4    cef12.4h 
##        8730        8747        8819        8880        9038        9215 
##    cep01.4h  up04.4h.ab  of09.4h.aa  of06.2h.ab out01.2h.aa  uf10.4h.ab 
##        9539        9576        9798        9834        9951       10005 
##  up02.4h.ab  up02.2h.aa    cep02.4h  op03.2h.ab  neg.filt.3  op06.2h.ab 
##       10008       10025       10186       10227       10517       10921 
##  up06.4h.jm    cwf07.4h  uf07.2h.ab  uf01.4h.jm  op04.4h.aa  of03.2h.aa 
##       11058       11095       11160       11260       11413       11538 
##  uf03.2h.ab  of04.4h.aa  uf11.2h.jm  neg.dish.3  uf10.4h.aa    cwp04.4h 
##       11545       11564       11792       12028       12067       12342 
##     extr.c5 out03.2h.jm  uf02.2h.ab  up01.4h.ab  up05.4h.aa  uf10.2h.ab 
##       12369       12390       12525       12640       12675       12836 
##  up04.2h.jm    extr.c15  of03.4h.jm  op05.4h.ab  up06.2h.jm    extr.c13 
##       12920       12988       13188       13198       13259       13294 
##  uf02.4h.ab  uf02.4h.jm  uf06.2h.jm  up06.4h.ab  of10.2h.jm  uf12.4h.jm 
##       13666       13937       14050       14219       14301       14321 
##  of10.2h.ab  of04.4h.ab  uf09.4h.jm  uf08.4h.jm    cep04.4h  uf04.4h.jm 
##       14456       14538       14640       14784       15010       15116 
##    cwf09.4h out01.4h.jm  of12.4h.aa  op05.2h.ab  of03.2h.jm out01.2h.jm 
##       15186       15417       15568       15660       15663       16003 
##  uf11.4h.aa out02.2h.jm  op04.2h.aa  uf01.4h.aa  of08.2h.aa  uf08.2h.ab 
##       16104       16221       16316       16349       16363       16425 
##  of02.2h.ab  of04.2h.aa  uf01.2h.jm    cwf11.4h  op06.2h.jm  of12.2h.ab 
##       16601       16743       16834       16914       16932       17060 
##  of10.4h.ab    cef02.4h  up03.2h.jm  up03.4h.jm  up01.2h.jm  up06.4h.aa 
##       17255       17583       17865       17951       18005       18090 
##    cwp06.4h  neg.filt.2  of03.2h.ab out01.2h.ab  op05.2h.aa  uf09.4h.ab 
##       18403       18472       19093       19245       19319       19378 
##  uf04.4h.ab  neg.dish.2  of03.4h.ab  op03.4h.aa  neg.dish.1  uf10.2h.jm 
##       19448       19501       19628       19834       20297       20469 
##  up01.4h.aa  up03.2h.ab  uf07.4h.jm    extr.c14  op03.4h.ab  uf10.2h.aa 
##       20531       20752       20805       21268       21727       21865 
##  uf12.4h.ab  op04.2h.jm  uf07.2h.aa  of06.2h.jm  of01.2h.ab  up05.4h.ab 
##       21908       21964       22021       22346       22366       22376 
##  of09.2h.ab  of08.4h.aa  uf11.2h.aa  op04.4h.ab  up03.2h.aa  uf11.2h.ab 
##       22706       22832       23312       23317       23365       23582 
##     extr.c9    cwf05.4h  of07.2h.jm  uf06.2h.ab  neg.filt.1  uf04.2h.jm 
##       23754       24255       24993       25018       25032       25112 
##  of06.2h.aa out01.4h.ab  op01.4h.ab  of11.4h.aa  up02.2h.jm  uf04.4h.aa 
##       25959       26000       26460       27047       27156       27232 
##  of05.2h.ab out03.2h.ab  of04.2h.jm  of09.4h.ab  of02.4h.aa out02.4h.jm 
##       28157       28771       29018       29385       29825       30148 
##  of05.4h.aa  up01.4h.jm  neg.swab.1     extr.c8  of11.4h.ab  op01.2h.ab 
##       30228       30325       30427       31143       31996       32111 
##    cep05.4h  of05.2h.jm  uf01.2h.aa  op02.2h.jm  of05.4h.jm  of12.4h.jm 
##       32419       32433       33129       34332       34700       34943 
##    cef01.4h    cwf04.4h    cef05.4h  uf02.4h.aa  uf07.4h.aa  of06.4h.aa 
##       35256       35367       35521       36051       36484       36516 
##  op03.2h.jm  op01.2h.aa  of07.2h.ab  op02.4h.jm  of04.2h.ab  uf10.4h.jm 
##       36521       37789       37853       38542       39565       40457 
##  op01.4h.jm c.out.03.4h  of05.4h.ab    cef07.4h c.out.01.4h  of02.4h.jm 
##       40679       41086       41274       42017       42229       42694 
##    cep03.4h    cef04.4h  up05.2h.jm  op04.4h.jm  of12.2h.jm  up02.2h.ab 
##       42885       43112       43608       43672       43735       43739 
##  of12.2h.aa  uf04.2h.aa  op02.4h.aa  of10.4h.jm  up06.2h.ab  up04.4h.jm 
##       46181       47666       47810       47855       47931       48227 
##  up02.4h.aa  of11.4h.jm  of08.2h.jm out03.4h.aa out01.4h.aa  uf05.4h.jm 
##       48339       48927       49056       49157       50320       50777 
##  uf02.2h.jm  uf12.2h.aa    cwp03.4h  uf06.2h.aa    cwf06.4h  op02.2h.aa 
##       50813       51118       51128       51580       51608       51663 
##  of08.2h.ab out02.4h.aa  uf05.2h.aa  uf02.2h.aa  up06.2h.aa  of02.2h.jm 
##       52635       52806       53026       53181       53829       54514 
##  up01.2h.ab out03.4h.ab  of01.2h.jm  of10.4h.aa  up05.4h.jm  of08.4h.ab 
##       55121       55624       56067       58667       59347       59620 
##    cef10.4h  of07.4h.jm  uf05.4h.ab  of10.2h.aa  of02.4h.ab  of07.2h.aa 
##       60368       61008       61050       61256       61834       62152 
##    cwf08.4h    cwp05.4h  of11.2h.aa  of01.4h.jm  uf06.4h.ab  op05.4h.jm 
##       62274       62603       63384       64005       64060       64454 
##  op06.2h.aa  op04.2h.ab  of09.2h.jm  of12.4h.ab  uf01.4h.ab  up02.4h.jm 
##       68753       68922       69219       71576       73548       73725 
##  up03.4h.aa  of06.4h.jm  uf08.4h.ab c.out.02.4h  up04.2h.ab out03.4h.jm 
##       73919       73954       74232       76112       77925       79188 
##  uf11.4h.ab  op06.4h.aa out02.2h.aa  of07.4h.ab  up01.2h.aa out02.4h.ab 
##       80253       80564       84876       85444       87130       88968 
##  of09.4h.jm  op01.4h.aa     extr.c2    cef06.4h  op01.2h.jm  of08.4h.jm 
##       89578       94784      100418      103076      104283      108327 
##  op03.4h.jm  of01.4h.aa  op03.2h.aa  uf07.4h.ab  op02.2h.ab out03.2h.aa 
##      108694      113922      115487      117704      119462      120473 
##  op06.4h.ab  op05.4h.aa  of06.4h.ab  of01.4h.ab    cwf03.4h  op02.4h.ab 
##      131675      142713      186451      205422      235565      311310
```

```r
pb.tmp <- pb[-nonexp, ]
pb.tmp <- pb.tmp[-which(rowSums(pb.tmp) < 3000), ]
pb.3500 <- rrarefy(pb.tmp, 3500)
pb.map <- pb.map[row.names(pb.3500), c("person", "location", "duration", "sampleType")]
dim(pb.3500)
```

```
## [1]    234 234213
```

```r
dim(pb.map)
```

```
## [1] 234   4
```

```r
identical(row.names(pb.3500), row.names(pb.map))
```

```
## [1] TRUE
```



Make taxonomy data frame for indexing. 


```r
taxo <- makeTaxo(pb.tax, pb.3500)
```

```
## Warning: no non-missing arguments to max; returning -Inf
```

```r
head(taxo)
```

```
##    kingdom         phylum               class              order
## 0 Bacteria  Bacteroidetes     Sphingobacteria Sphingobacteriales
## 1 Bacteria Proteobacteria Gammaproteobacteria                   
## 2 Bacteria Proteobacteria Gammaproteobacteria  Oceanospirillales
## 3 Bacteria Actinobacteria      Actinobacteria    Actinomycetales
## 4 Bacteria Actinobacteria      Actinobacteria    Actinomycetales
## 5 Bacteria                                                      
##             family genus abundance
## 0 Chitinophagaceae               0
## 1                                0
## 2   Halomonadaceae               1
## 3                                0
## 4                                1
## 5                                0
```



Get rid of empty OTUs to reduce computing demand. 


```r
pb.3500 <- pb.3500[, -which(colSums(pb.3500) == 0)]
taxo <- taxo[colnames(pb.3500), ]
head(taxo)
```

```
##     kingdom         phylum               class             order
## 2  Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 4  Bacteria Actinobacteria      Actinobacteria   Actinomycetales
## 10 Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
## 15 Bacteria Proteobacteria  Betaproteobacteria                  
## 16 Bacteria     Firmicutes             Bacilli        Bacillales
## 19 Bacteria                                                     
##            family        genus abundance
## 2  Halomonadaceae                      1
## 4                                      1
## 10 Halomonadaceae                      2
## 15                                     1
## 16    Bacillaceae Tumebacillus         1
## 19                                     1
```

```r
dim(taxo)
```

```
## [1] 33504     7
```

```r
dim(pb.3500)
```

```
## [1]   234 33504
```

```r
identical(row.names(taxo), colnames(pb.3500))
```

```
## [1] TRUE
```

```r
pb.3500.ra <- pb.3500/35
```



Quick refrerence vector of names for figures and whatnot. Some genus, families, etc are empty, so search for the finest taxonomic level and use that name. 


```r
cons <- function(x) {
    l <- length(x)
    while (x[l] == "") {
        l = l - 1
    }
    name <- x[l]
}
consensus <- apply(taxo[, 1:6], 1, cons)
consensus[1:10]
```

```
##                    2                    4                   10 
##     "Halomonadaceae"    "Actinomycetales"     "Halomonadaceae" 
##                   15                   16                   19 
## "Betaproteobacteria"       "Tumebacillus"           "Bacteria" 
##                   21                   34                   70 
##           "Bacteria"        "Bacillaceae"     "Proteobacteria" 
##                   79 
##     "Comamonadaceae"
```



Some indexing and metadata stuff for figures and such. 


```r
occ <- which(pb.map$location == "occ")
unocc <- which(pb.map$location == "unocc")
control <- which(pb.map$location == "control")

pb.map$pch <- 21
table(pb.map$sampleType)
```

```
## 
## extraction.control             filter         petri.dish 
##                 12                142                 74 
##    reagent.control               swab 
##                  2                  4
```

```r
pb.map$pch2 <- pb.map$pch
pb.map$pch2[pb.map$sampleType == "extraction.control"] <- 24
pb.map$pch2[pb.map$sampleType == "reagent.control"] <- 23
pb.map$pch2[pb.map$sampleType == "petri.dish" & pb.map$location == "control"] <- 22

pb.map$bg <- "gray60"
pb.map$bg[occ] <- "cornflowerblue"
pb.map$bg[unocc] <- "darkorange"
```



Make dissimilarity matrix for big dataset - this takes a while. 
Make with both Canberra and Bray-Curtis, and then make NMDS ordinations. 
The Bray-Curtis will be more useful for identifying the influence of the most abundant taxa, 
while Canberra is used for analysis. 
Here, Canberra has high stress, so is not used much. 
Also quick plot to make sure that BC captures differences between groups. 


```r
can <- vegdist(pb.3500, "canberra")
bc <- vegdist(pb.3500)

nmds.can <- nmds(can)
```

```
## initial  value 40.307806 
## iter   5 value 32.361958
## iter  10 value 31.603218
## final  value 31.348427 
## converged
```

```r
nmds.bc <- nmds(bc)
```

```
## initial  value 24.728368 
## iter   5 value 18.718472
## final  value 18.288598 
## converged
```




The BC ordination clearly shows that controls are tightly grouped (mostly), while occupied and unoccupied are defintely statistically different. Even with contaminants. Strong gradient along x-axis >> will be used for regression below. 


```r
# pdf('~/Desktop/contamination_ordination.pdf', useDingbats=TRUE)
plot(nmds.bc$points, col = pb.map$location, type = "n")
points(nmds.bc$points[occ, ], pch = 16, col = pb.map$bg[occ])
points(nmds.bc$points[unocc, ], pch = 16, col = pb.map$bg[unocc])
points(nmds.bc$points[control, ], pch = 21, bg = pb.map$bg[control])
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-121.png) 

```r
# dev.off()

plot(nmds.bc$points, type = "n", xlim = c(-0.5, -0.3), ylim = c(-0.2, 0.1))
points(nmds.bc$points[occ, ], pch = 16, col = pb.map$bg[occ])
points(nmds.bc$points[unocc, ], pch = 16, col = pb.map$bg[unocc])
points(nmds.bc$points[control, ], pch = 21, bg = pb.map$bg[control])
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-122.png) 



Find the most abundant taxa. This is done for 10, 100, 1000. Use 10 for plotting, and take out those contaminants, but use 1000 since some contaminants are really rare in actual samples. 


```r
control.taxa.10 <- rev(sort(colSums(pb.3500[control, ])))[1:10]/sum(pb.3500[control, 
    ])
occ.taxa.10 <- rev(sort(colSums(pb.3500[occ, ])))[1:10]/sum(pb.3500[occ, ])
unocc.taxa.10 <- rev(sort(colSums(pb.3500[unocc, ])))[1:10]/sum(pb.3500[unocc, 
    ])

control.taxa.10["other"] <- 1 - sum(control.taxa.10)
occ.taxa.10["other"] <- 1 - sum(occ.taxa.10)
unocc.taxa.10["other"] <- 1 - sum(unocc.taxa.10)

top.10 <- cbind(rev(control.taxa.10), rev(unocc.taxa.10), rev(occ.taxa.10))
# mids <- barplot(top.10, col=c('gray30', rep('gray70', 10)),
# border='gray30', space=1)

con.cum <- cumsum(rev(control.taxa.10))
occ.cum <- cumsum(rev(occ.taxa.10))
unocc.cum <- cumsum(rev(unocc.taxa.10))

control.taxa.100 <- rev(sort(colSums(pb.3500[control, ])))[1:100]/sum(pb.3500[control, 
    ])
occ.taxa.100 <- rev(sort(colSums(pb.3500[occ, ])))[1:100]/sum(pb.3500[occ, ])
unocc.taxa.100 <- rev(sort(colSums(pb.3500[unocc, ])))[1:100]/sum(pb.3500[unocc, 
    ])

control.taxa.100["other"] <- 1 - sum(control.taxa.100)
occ.taxa.100["other"] <- 1 - sum(occ.taxa.100)
unocc.taxa.100["other"] <- 1 - sum(unocc.taxa.100)

con.cum.100 <- cumsum(rev(control.taxa.100))
occ.cum.100 <- cumsum(rev(occ.taxa.100))
unocc.cum.100 <- cumsum(rev(unocc.taxa.100))

occ.taxa.1000 <- rev(sort(colSums(pb.3500[occ, ])))[1:1000]/sum(pb.3500[occ, 
    ])
unocc.taxa.1000 <- rev(sort(colSums(pb.3500[unocc, ])))[1:1000]/sum(pb.3500[unocc, 
    ])

occ.taxa.1000["other"] <- 1 - sum(occ.taxa.1000)
unocc.taxa.1000["other"] <- 1 - sum(unocc.taxa.1000)

occ.cum.1000 <- cumsum(rev(occ.taxa.1000))
unocc.cum.1000 <- cumsum(rev(unocc.taxa.1000))
```






Halomonas is the most abundant across the board. It drives the BC NMDS. 


```r
rev(sort(colSums(pb.3500)))[1:10]
```

```
## 190873 234992  64285  25506 121358 201443  94008  80526 151122  72518 
## 279828  66921  60454  45640  37840  28903  24596  20070  17032  12480
```

```r
rev(sort(colSums(pb.3500)))[1]/sum(pb.3500)  # 34%
```

```
## 190873 
## 0.3417
```

```r
halo <- "190873"
taxo["190873", ]
```

```
##         kingdom         phylum               class             order
## 190873 Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
##                family     genus abundance
## 190873 Halomonadaceae Halomonas    279828
```

```r
pb.3500[, "190873"]
```

```
## op01.4h.jm   extr.c10 uf10.4h.aa neg.dish.3 op01.4h.ab op06.2h.aa 
##        867       2151       1609       2134        350        867 
## of09.2h.ab of06.4h.ab op02.4h.ab uf07.4h.ab up05.4h.jm uf03.4h.ab 
##        440        591         89        871       1028       1494 
## op06.4h.ab uf05.2h.aa uf02.2h.jm op06.4h.aa of12.2h.aa up05.2h.ab 
##        564       1621       1609        488       1040       2310 
## of06.4h.jm    extr.c2 of01.2h.ab of09.2h.aa uf04.4h.aa op02.2h.ab 
##       1056        183       1358        619       1245        113 
## of01.2h.jm uf01.2h.aa of05.2h.ab uf11.4h.ab uf08.2h.aa op04.2h.ab 
##        638       1497       1281       1240       2236        889 
## of03.2h.ab uf02.2h.aa neg.filt.1 of12.4h.ab of07.2h.ab uf10.4h.jm 
##       1652       2016       1866        611        771       1270 
## up02.4h.jm of08.4h.jm up06.2h.aa of05.4h.jm op05.4h.aa of11.4h.jm 
##       1327       1097        378       1111        741        359 
## of08.2h.jm up04.2h.ab    extr.c8 uf04.2h.aa uf06.4h.aa of08.2h.aa 
##        503       1163       2406       2041        685        744 
## of01.4h.aa up01.2h.aa of07.4h.jm uf01.2h.jm of02.4h.jm of04.2h.jm 
##        798       1263        681       1552        545        587 
## up02.4h.aa of10.4h.aa uf03.4h.jm of01.4h.jm neg.dish.2 of06.2h.jm 
##       2188        638        649       1113       2284        490 
## op05.4h.jm of12.2h.jm uf05.4h.ab uf07.2h.ab op04.4h.jm op03.4h.jm 
##        220        558       1426       1721        799        541 
## of08.2h.ab of09.4h.jm uf07.4h.aa op01.4h.aa of07.4h.ab of01.4h.ab 
##       1543        923       1976        136        259        351 
## uf04.4h.ab uf06.2h.aa of02.4h.ab uf03.4h.aa of10.2h.aa of05.2h.jm 
##       2189       1115        302        910        542        685 
## up05.2h.jm up03.4h.aa uf03.2h.ab uf04.2h.jm of03.4h.ab of11.4h.aa 
##        292       1069       1224       1955       1223        923 
## op04.4h.aa op05.2h.ab of10.2h.ab op02.2h.jm of04.2h.ab of09.4h.aa 
##       1047       1263        933        909       1780       1642 
## op01.2h.jm op02.4h.aa uf05.2h.jm op03.2h.jm uf12.2h.ab uf01.4h.aa 
##        864        118       1759       1227       1587        749 
## op05.2h.aa op03.2h.aa of03.4h.jm uf02.4h.aa uf06.2h.ab of02.2h.jm 
##       1647        655        369        974       2303       1031 
## of09.2h.jm of10.4h.jm uf07.2h.aa uf01.4h.ab op03.4h.ab of08.4h.ab 
##        953        623        967       1373        530        566 
##   extr.c15 op03.4h.aa op01.2h.aa of11.2h.aa up04.4h.jm of09.4h.ab 
##       2030        905       1232        599       1658        321 
## uf09.4h.jm uf08.4h.ab uf05.4h.jm uf08.4h.jm uf06.4h.ab up06.2h.ab 
##       1062       1423       1372        636       2099       1605 
## up04.4h.ab uf11.4h.aa neg.filt.2 op02.4h.jm uf03.2h.aa op02.2h.aa 
##       1601       1582       2066        906       1503        664 
## op03.2h.ab uf09.4h.ab neg.swab.3 up02.2h.jm of07.2h.jm up01.2h.ab 
##       2039        550       2221       1175        553       1295 
## up05.4h.ab of03.2h.aa    extr.c9 op04.4h.ab of05.4h.aa up06.4h.aa 
##       1529       1303       2098       1191        335       1279 
## of06.2h.aa uf09.2h.ab   extr.c14 up01.2h.jm up03.2h.jm up04.2h.aa 
##       2018       1002       1690       1289       2181        828 
## of04.2h.aa uf03.2h.jm op04.2h.jm uf01.4h.jm of05.4h.ab of04.4h.ab 
##        976       1781        308       1708        670        615 
## of12.4h.jm uf10.2h.ab uf02.4h.ab of11.4h.ab uf10.2h.aa of02.2h.aa 
##        606       1365        432        281        914       1389 
## up06.2h.jm up04.2h.jm up03.4h.jm uf12.2h.jm uf10.2h.jm of10.2h.jm 
##       1003       1975       1191        864       1085        334 
## op06.2h.jm uf07.4h.jm of02.4h.aa uf09.2h.jm uf02.4h.jm uf04.4h.jm 
##       1052       1273        812       1222       1888        393 
## uf08.4h.aa of08.4h.aa op05.4h.ab of07.2h.aa uf08.2h.jm uf11.2h.ab 
##       1031        774        744        832        495       1037 
## up03.4h.ab uf10.4h.ab up01.4h.jm of06.2h.ab uf12.4h.jm neg.swab.1 
##        497       2067        997       1812       1559       1698 
## up03.2h.aa up02.2h.ab of06.4h.aa op06.2h.ab up01.4h.aa neg.dish.1 
##       1355        870        938       1259       1227       2586 
## uf12.2h.aa of03.2h.jm    extr.c6 up06.4h.jm op05.2h.jm op01.2h.ab 
##       1175        960       1490       1218        846        659 
## op04.2h.aa up02.2h.aa of04.4h.aa neg.filt.4 of01.2h.aa of11.2h.ab 
##        395       2075       2804       2317       1148       1561 
## up05.4h.aa uf08.2h.ab    extr.c5   extr.c11 of10.4h.ab uf02.2h.ab 
##       2135       1908        505       2191        937       2051 
##   extr.c12 up03.2h.ab    extr.c1 uf12.4h.aa uf01.2h.ab of02.2h.ab 
##       2092       1180       1315       2037        726       1684 
## of03.4h.aa up01.4h.ab   extr.c13 neg.swab.4 of04.4h.jm of11.2h.jm 
##        739       1159       1649       2430       1220        404 
## neg.swab.2 uf11.2h.jm up02.4h.ab of07.4h.aa uf09.4h.aa of12.2h.ab 
##       2385        431       2036        859       1466       1401 
## uf06.2h.jm op06.4h.jm neg.filt.3 of05.2h.aa uf12.4h.ab    extr.c7 
##       1440        790       2246       1583       2223       2306 
## of12.4h.aa uf11.2h.aa up06.4h.ab uf06.4h.jm neg.dish.4    extr.c4 
##       1490       1355       2191        552       2310       1241
```

```r
haloCol <- pb.3500[, halo]

plot(haloCol ~ nmds.bc$points[, 1], type = "n")
points(haloCol[occ] ~ nmds.bc$points[occ, 1], pch = 16, col = pb.map$bg[occ])
points(haloCol[unocc] ~ nmds.bc$points[unocc, 1], pch = 16, col = pb.map$bg[unocc])
points(haloCol[control] ~ nmds.bc$points[control, 1], pch = 21, bg = pb.map$bg[control], 
    )
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14.png) 

```r
# identify(haloCol[control] ~ nmds.bc$points[control, 1],
# labels=pb.map$sampleType[control]) the mixed group of controls are all
# extraction blanks
```




Big panel figure showing each top contaminant in the dataset, and its influence on NMDS1


```r

# pdf('~/Desktop/contaminationTop10.pdf', height=10, width=10,
# useDingbats=TRUE)
layout(matrix(c(11:15, 1:10, 16:20), 5, 4), heights = c(1, 1, 1, 1, 1.5), widths = c(1, 
    2, 2, 1))
par(las = 1, mar = c(0, 4, 0, 0), fg = "gray40", col.axis = "gray40", col.lab = "gray40")
for (i in 1:10) {
    id <- names(control.taxa.10)[i]
    print(taxo[id, ])
    idCol <- pb.3500.ra[, id]
    print(sum(idCol)/sum(pb.3500.ra))
    y.lim <- round(range(idCol), 2)
    ybuff <- y.lim[2]/8
    x.lim <- round(range(nmds.bc$points[, 1]), 2)
    xbuff <- sum(abs(x.lim))/8
    
    if (i %in% c(1:4)) {
        par(las = 1, mar = c(0, 4, 0, 0))
    }
    if (i %in% c(5:8)) {
        par(las = 1, mar = c(0, 0, 0, 4))
    }
    if (i == 5) {
        par(las = 1, mar = c(4, 4, 0, 0))
    }
    if (i == 10) {
        par(las = 1, mar = c(4, 0, 0, 4))
    }
    plot(idCol ~ nmds.bc$points[, 1], ylim = c(y.lim + c(-ybuff, 2 * ybuff)), 
        xlim = c(x.lim + c(-xbuff, xbuff)), type = "n", xaxt = "n", yaxt = "n", 
        xlab = "", ylab = "")
    if (i %in% c(1:5)) {
        yax <- 2
    } else {
        yax <- 4
    }
    axis(side = yax, at = y.lim, labels = TRUE)
    if (i %in% c(5, 10)) {
        axis(side = 1, at = x.lim)
        mtext("NMDS 1", side = 1, line = 1, cex = 0.7)
    }
    if (i == 1) {
        mtext("relative\nabundance", side = 2, cex = 0.7, line = 0.5)
    }
    if (i == 6) {
        mtext("relative\nabundance", side = 4, cex = 0.7, line = 0.5)
        legend("right", legend = c("extraction control", "reagent control", 
            "petri control", "dish control", "occupied", "unoccupied"), pch = c(24, 
            23, 22, 21, 16, 16), col = c(1, 1, 1, 1, "cornflowerblue", "darkorange"), 
            pt.bg = rgb(0, 0, 0, 0.3))
    }
    points(idCol[occ] ~ nmds.bc$points[occ, 1], pch = 16, col = pb.map$bg[occ])
    points(idCol[unocc] ~ nmds.bc$points[unocc, 1], pch = 16, col = pb.map$bg[unocc])
    points(idCol[control] ~ nmds.bc$points[control, 1], pch = pb.map$pch2[control], 
        col = "black", bg = rgb(0, 0, 0, 0.3), cex = 2)
    mtext(paste(consensus[id], " (", id, ")", sep = ""), line = -1.5)
}
```

```
##         kingdom         phylum               class             order
## 190873 Bacteria Proteobacteria Gammaproteobacteria Oceanospirillales
##                family     genus abundance
## 190873 Halomonadaceae Halomonas    279828
## [1] 0.3417
```

```
##        kingdom         phylum               class           order
## 64285 Bacteria Proteobacteria Gammaproteobacteria Alteromonadales
##               family      genus abundance
## 64285 Shewanellaceae Shewanella     60454
## [1] 0.07381
```

```
##         kingdom         phylum               class       order
## 234992 Bacteria Proteobacteria Alphaproteobacteria Rhizobiales
##                     family            genus abundance
## 234992 Methylobacteriaceae Methylobacterium     66921
## [1] 0.08171
```

```
##         kingdom         phylum               class           order
## 201443 Bacteria Proteobacteria Gammaproteobacteria Xanthomonadales
##                  family            genus abundance
## 201443 Xanthomonadaceae Stenotrophomonas     28903
## [1] 0.03529
```

```
##         kingdom     phylum   class      order      family        genus
## 121358 Bacteria Firmicutes Bacilli Bacillales Bacillaceae Tumebacillus
##        abundance
## 121358     37840
## [1] 0.0462
```

```
##        kingdom         phylum              class           order
## 96443 Bacteria Proteobacteria Betaproteobacteria Burkholderiales
##                 family        genus abundance
## 96443 Burkholderiaceae Burkholderia      3582
## [1] 0.004374
```

```
##         kingdom         phylum              class           order
## 129992 Bacteria Proteobacteria Betaproteobacteria Burkholderiales
##                family genus abundance
## 129992 Comamonadaceae            5273
## [1] 0.006438
```

```
##         kingdom         phylum               class             order
## 151122 Bacteria Proteobacteria Gammaproteobacteria Enterobacteriales
##                    family genus abundance
## 151122 Enterobacteriaceae           17032
## [1] 0.0208
```

```
##         kingdom        phylum           class              order
## 180093 Bacteria Bacteroidetes Sphingobacteria Sphingobacteriales
##               family        genus abundance
## 180093 Cytophagaceae Hymenobacter       782
## [1] 0.0009548
```

```
##        kingdom         phylum          class           order
## 94008 Bacteria Actinobacteria Actinobacteria Actinomycetales
##                   family           genus abundance
## 94008 Corynebacteriaceae Corynebacterium     24596
## [1] 0.03003
```

```r

for (i in 11:2) {
    xup <- xdown <- yup <- ydown <- NULL
    
    yup <- rep(con.cum[i], 2)
    ydown <- rep(con.cum[i - 1], 2)
    xup <- xdown <- c(1, 2)
    
    if (names(con.cum)[i] %in% names(unocc.cum.1000)) {
        ucu <- which(names(unocc.cum.1000) == names(con.cum)[i])
        yup <- c(yup, rep(unocc.cum.1000[ucu], 2))
        ydown <- c(ydown, rep(unocc.cum.1000[ucu - 1], 2))
        xup <- xdown <- c(xup, 3, 4)
    }
    if (names(con.cum)[i] %in% names(occ.cum.1000)) {
        ocu <- which(names(occ.cum.1000) == names(con.cum)[i])
        yup <- c(yup, rep(occ.cum.1000[ocu], 2))
        ydown <- c(ydown, rep(occ.cum.1000[ocu - 1], 2))
        xup <- xdown <- c(xup, 5, 6)
    }
    par(mar = c(3, 3, 3, 3))
    if (i %in% c(7, 2)) {
        par(mar = c(8, 3, 3, 3))
    }
    barplot(top.10, col = c("gray30", rep("gray93", 10)), border = "gray30", 
        space = 1, yaxt = "n")
    if (i %in% c(7, 2)) {
        par(las = 2)
        mtext(c("control", "unoccupied", "occupied"), side = 1, at = c(1.5, 
            3.5, 5.5), cex = 0.8, line = 0.2)
    }
    par(las = 1)
    if (!i %in% c(7, 2)) {
        par(las = 2)
        mtext(c("c", "u", "o"), side = 1, at = c(1.5, 3.5, 5.5), cex = 0.8, 
            line = 0.2)
    }
    par(las = 1)
    if (i %in% 11:7) {
        axis(2, at = c(0, 0.5, 1), labels = c(0, 50, 100))
    }
    if (i %in% 2:6) {
        axis(4, at = c(0, 0.5, 1), labels = c(0, 50, 100))
    }
    polygon(c(xup, rev(xdown)), c(yup, rev(ydown)), col = rgb(0, 0, 0, 0.3))
}
```

![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-15.png) 

```r

# dev.off()
```




Now take out those most abundant in pcr controls.


```r
dim(pb.tmp)
```

```
## [1]    234 234213
```

```r
pb.c <- pb.tmp[-c(control), ]
pb.nc <- pb.tmp[-c(control), -which(colnames(pb.tmp) %in% names(control.taxa.10[c(1, 
    2, 6, 7)]))]
pb.map.nc <- pb.map[-control, ]

sort(rowSums(pb.nc))
```

```
## of04.4h.aa uf08.2h.aa up05.2h.ab uf09.4h.aa uf12.2h.ab of01.2h.aa 
##       1035       1209       1264       1298       1832       1981 
## uf05.2h.jm of11.2h.ab uf10.4h.ab uf12.4h.aa up02.2h.aa of04.4h.jm 
##       2132       2313       2361       2411       2414       2563 
## uf03.4h.ab uf03.2h.jm up02.4h.ab uf03.2h.aa op03.2h.ab up06.4h.ab 
##       2577       2690       2725       2755       2780       2896 
## of02.2h.aa of03.4h.aa up05.4h.aa of11.2h.jm up04.2h.jm uf02.2h.ab 
##       3276       3339       3429       3544       3688       3716 
## of06.2h.ab uf06.4h.jm of05.2h.aa uf09.2h.jm uf09.2h.ab up04.4h.ab 
##       3770       3813       3929       3956       3974       4049 
## uf08.4h.aa up03.2h.jm uf06.2h.ab uf01.4h.jm of09.4h.aa uf08.2h.jm 
##       4254       4255       4340       4379       4398       4431 
## of09.2h.aa up04.2h.aa op06.4h.jm op05.2h.jm uf03.4h.aa uf07.2h.ab 
##       4633       4664       4743       4774       4858       5012 
## uf12.4h.ab uf06.4h.aa uf10.4h.aa uf02.4h.jm uf12.2h.jm uf04.4h.ab 
##       5403       5463       5465       5507       5553       5567 
## of07.4h.aa op06.2h.ab of02.2h.ab up06.4h.jm uf01.2h.ab uf12.4h.jm 
##       5695       5766       5859       6039       6122       6242 
## of03.2h.aa uf03.4h.jm uf08.2h.ab op05.2h.aa uf06.2h.jm of12.4h.aa 
##       6471       6544       6603       6699       6853       6966 
## uf03.2h.ab up01.4h.ab op04.4h.aa uf10.2h.ab uf04.2h.jm up03.4h.ab 
##       7024       7045       7058       7083       7331       7417 
## uf11.4h.aa up06.2h.jm of03.2h.ab uf01.2h.jm of12.2h.ab of06.2h.aa 
##       7708       8239       8321       8341       8599       8667 
## op05.2h.ab uf09.4h.jm up01.2h.jm up06.4h.aa op05.4h.ab op06.2h.jm 
##       8832       8924       8944       9123       9728       9836 
## of10.2h.ab uf11.2h.jm of03.2h.jm up05.4h.ab up03.2h.aa of04.2h.aa 
##      10098      10105      10287      10455      10719      10737 
## up03.4h.jm uf07.4h.aa of04.4h.ab of03.4h.ab uf02.4h.ab of03.4h.jm 
##      10874      10982      11224      11229      11361      11362 
## up03.2h.ab of10.4h.ab of01.2h.ab uf08.4h.jm uf07.4h.jm up01.4h.aa 
##      11409      11440      11496      11595      11680      11794 
## op04.4h.ab uf01.4h.aa uf04.2h.aa of10.2h.jm of08.2h.aa up02.4h.aa 
##      11943      12079      12268      12565      12570      12576 
## uf04.4h.jm op03.4h.aa uf10.2h.jm uf11.2h.aa op04.2h.aa uf01.2h.aa 
##      12827      12839      12993      13104      14264      14296 
## uf11.2h.ab uf07.2h.aa uf04.4h.aa uf10.2h.aa of05.2h.ab uf09.4h.ab 
##      14360      14393      15082      15325      15617      15760 
## of04.2h.ab uf02.2h.aa uf06.4h.ab of08.4h.aa up02.2h.jm up01.4h.jm 
##      15851      15868      16233      16497      16667      17823 
## op03.4h.ab of11.4h.aa of06.2h.jm of09.2h.ab op04.2h.jm of07.2h.jm 
##      17842      18619      18804      19083      19087      19664 
## op03.2h.jm of02.4h.aa up04.4h.jm uf02.2h.jm of05.4h.jm of08.2h.ab 
##      19949      20004      20114      20435      20546      21763 
## up06.2h.ab uf05.4h.jm op01.2h.aa op02.2h.jm op01.4h.ab of04.2h.jm 
##      21847      22022      22716      22824      23185      23190 
## op01.2h.ab uf05.2h.aa uf10.4h.jm of06.4h.aa of05.2h.jm uf02.4h.aa 
##      23312      23526      23916      24057      24064      24119 
## of09.4h.ab op02.4h.jm of05.4h.aa of07.2h.ab of12.4h.jm uf05.4h.ab 
##      26054      26422      26816      27521      27893      27964 
## of11.4h.ab up02.2h.ab op01.4h.jm up01.2h.ab of12.2h.aa uf12.2h.aa 
##      28939      29254      29271      29317      30192      30442 
## op04.4h.jm of05.4h.ab uf06.2h.aa uf08.4h.ab of02.4h.jm of12.2h.jm 
##      30820      30826      31470      33236      34169      34801 
## of02.2h.jm up02.4h.jm of10.4h.jm up05.4h.jm op02.2h.aa of01.4h.jm 
##      36046      37192      37617      37713      38675      38721 
## uf01.4h.ab up05.2h.jm of08.2h.jm of11.4h.jm up04.2h.ab of01.2h.jm 
##      39095      39258      40386      41980      43128      43367 
## uf11.4h.ab op06.2h.aa of07.2h.aa of10.4h.aa up01.2h.aa up06.2h.aa 
##      44531      44663      44820      45242      45419      45429 
## op02.4h.aa up03.4h.aa of06.4h.jm of08.4h.ab of07.4h.jm of09.2h.jm 
##      45633      45687      45882      47273      47494      48523 
## op04.2h.ab of11.2h.aa of10.2h.aa of02.4h.ab of12.4h.ab op05.4h.jm 
##      48552      48748      49228      55685      56417      58353 
## of09.4h.jm op06.4h.aa of08.4h.jm op01.2h.jm uf07.4h.ab of07.4h.ab 
##      62013      64385      64649      69957      76818      78604 
## of01.4h.aa op03.4h.jm op03.2h.aa op01.4h.aa op05.4h.aa op06.4h.ab 
##      82810      87558      87843      90912     100535     104632 
## op02.2h.ab of06.4h.ab of01.4h.ab op02.4h.ab 
##     113890     144827     173563     298953
```

```r
pb.nc.1000 <- rrarefy(pb.nc, 1000)
pb.nc.1000 <- pb.nc.1000[, -which(colSums(pb.nc.1000) == 0)]
pb.1000 <- rrarefy(pb.c, 1000)
pb.1000 <- pb.1000[, -which(colSums(pb.1000) == 0)]
```





Make same NMDS objects to compare.


```r
can.nc <- vegdist(pb.nc.1000, "canberra")
bc.nc <- vegdist(pb.nc.1000)
can.c <- vegdist(pb.1000, "canberra")
bc.c <- vegdist(pb.1000)

nmds.can.nc <- nmds(can.nc)
```

```
## initial  value 39.616151 
## iter   5 value 30.761031
## iter  10 value 30.282891
## iter  10 value 30.278175
## iter  10 value 30.249207
## final  value 30.249207 
## converged
```

```r
nmds.bc.nc <- nmds(bc.nc)
```

```
## initial  value 33.884719 
## iter   5 value 23.374712
## final  value 21.536789 
## converged
```

```r
nmds.can.c <- nmds(can.c)
```

```
## initial  value 42.609751 
## iter   5 value 31.960807
## iter  10 value 31.226954
## final  value 31.151002 
## converged
```

```r
nmds.bc.c <- nmds(bc.c)
```

```
## initial  value 26.232957 
## iter   5 value 20.403711
## final  value 19.973245 
## converged
```



Still see really strong pattern, but dissimilarities are pretty different. So this says probably redo all analyses to see how big of a difference. 


```r

levels(pb.map$location)
```

```
## [1] "control" "east"    "occ"     "out"     "unocc"   "west"
```

```r
plot(nmds.can.nc$points, col = pb.map.nc$location)
```

![plot of chunk unnamed-chunk-18](figure/unnamed-chunk-181.png) 

```r
plot(nmds.bc.nc$points, col = pb.map.nc$location)
```

![plot of chunk unnamed-chunk-18](figure/unnamed-chunk-182.png) 

```r
plot(nmds.can.c$points, col = pb.map.nc$location)
```

![plot of chunk unnamed-chunk-18](figure/unnamed-chunk-183.png) 

```r
plot(nmds.bc.c$points, col = pb.map.nc$location)
```

![plot of chunk unnamed-chunk-18](figure/unnamed-chunk-184.png) 

```r

plot(can.nc, can.c)
```

![plot of chunk unnamed-chunk-18](figure/unnamed-chunk-185.png) 

```r
plot(bc.nc, bc.c)
```

![plot of chunk unnamed-chunk-18](figure/unnamed-chunk-186.png) 




Still a very big statistical difference. 


```r

adonis(can.nc ~ pb.map.nc$location)
```

```
## 
## Call:
## adonis(formula = can.nc ~ pb.map.nc$location) 
## 
## Terms added sequentially (first to last)
## 
##                     Df SumsOfSqs MeanSqs F.Model    R2 Pr(>F)    
## pb.map.nc$location   1       1.1   1.053    2.23 0.011  0.001 ***
## Residuals          206      97.0   0.471         0.989           
## Total              207      98.1                 1.000           
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
adonis(bc.nc ~ pb.map.nc$location)
```

```
## 
## Call:
## adonis(formula = bc.nc ~ pb.map.nc$location) 
## 
## Terms added sequentially (first to last)
## 
##                     Df SumsOfSqs MeanSqs F.Model    R2 Pr(>F)    
## pb.map.nc$location   1       5.5    5.46    18.6 0.083  0.001 ***
## Residuals          206      60.5    0.29         0.917           
## Total              207      66.0                 1.000           
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```



```r

plot(nmds.bc$points, col = pb.map$location, type = "n")
points(nmds.bc$points[occ, ], pch = 16, col = pb.map$bg[occ])
points(nmds.bc$points[unocc, ], pch = 16, col = pb.map$bg[unocc])
points(nmds.bc$points[control, ], pch = 21, bg = pb.map$bg[control])
```

![plot of chunk unnamed-chunk-20](figure/unnamed-chunk-201.png) 

```r

plot(nmds.bc$points, type = "n", xlim = c(-0.5, -0.3), ylim = c(-0.2, 0.1))
points(nmds.bc$points[occ, ], pch = 16, col = pb.map$bg[occ])
points(nmds.bc$points[unocc, ], pch = 16, col = pb.map$bg[unocc])
points(nmds.bc$points[control, ], pch = 21, bg = pb.map$bg[control])
```

![plot of chunk unnamed-chunk-20](figure/unnamed-chunk-202.png) 










