

# Quality checking for PickleBox2014 sequences. 

_James F Meadow_ (`jfmeadow at gmail dot com`)

This script looks at the number of sequences returned for each sample. 




```r
# setwd('~/Dropbox/Pickle2014/R/quality')
set.seed(42)
require(xtable)
require(wesanderson)
```


Read data


```r
qu <- read.delim('split_library_log.txt', sep='\t', head=TRUE, row.names=1)
```



Assign colors


```r
pal <- wes.palette(name='Zissou', type='continuous')
```


Setup plotting parameters. 


```r
cutoffs <- c(100, 500, 1000, 2000, 5000)
bottom <- c(0, 100, 500, 1000, 2000)
nr <- nrow(qu)
col5 <- rev(pal(5))
xcut <- rep(0, 5)
for(i in 1:5) {
  xcut[i] <- length(which(qu$seqs < cutoffs[i]))
  }
```



Make plots. 


```r
par(mfrow=c(2, 1))
par(xpd=FALSE, las=1, mar=c(5, 5, 4, 1))
plot(qu$seqs, las=1, type='n', bty='n', xaxs='i', yaxs='i', las=1, 
     ann=FALSE, xaxt='n')
for(i in 1:5) {
  polygon(c(1, nr, nr, 1), c(bottom[i], bottom[i], cutoffs[i], cutoffs[i]), 
          border='transparent', col=col5[i])
  }
segments(c(nr-xcut), cutoffs, c(nr-xcut), -100, 
         lty=3, lwd=2, col='gray40')
par(xpd=TRUE)
points(1:nrow(qu), qu$seqs, pch=21, bg=rgb(0,0,0,.3), col=rgb(0,0,0,.7))
text(c(0, nr-xcut), 0, c(nr, xcut), pos=1)
mtext('Seqs per sample', side=2, line=3.5, las=0)
mtext('How many samples lost?', side=1, line=2)
text(1, max(qu$seqs), max(qu$seqs), pos=4)

par(xpd=FALSE, las=1)
plot(qu$seqs, las=1, type='n', bty='n', xaxs='i', yaxs='i', 
     ylim=c(0, 5200), las=1, ann=FALSE, , yaxt='n',
     ylab='Number of sequences per sample', xaxt='n')
axis(side=2, at=c(cutoffs))
abline(h=cutoffs, col=col5, las=1)
for(i in 1:5) {
  polygon(c(1, nr, nr, 1), c(bottom[i], bottom[i], cutoffs[i], cutoffs[i]), 
          border='transparent', col=col5[i])
  }
par(xpd=TRUE)
segments(c(nr-xcut), cutoffs, c(nr-xcut), -20, 
         lty=3, lwd=2, col='gray40')
points(1:nrow(qu), qu$seqs, pch=21, bg=rgb(0,0,0,.3), col=rgb(0,0,0,.7))
text(c(0, nr-xcut), 0, c(nr, xcut), pos=1)
mtext('Seqs per sample', side=2, line=3.5, las=0)
mtext('How many samples lost?', side=1, line=2)
```

![plot of chunk plotSeqs](figure/plotSeqs.png) 







