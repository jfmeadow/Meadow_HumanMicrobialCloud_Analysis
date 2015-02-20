# cutoffs =    0.3   0.5     1   2.5     5    10+
rm(list=ls())

part <- read.delim('particlesCombined.txt', head=TRUE, sep='\t')
names(part)
part[, c(4:15, 17:28)] <- part[, c(4:15, 17:28)]/2.8
part$subTreat <- factor(paste0(part$subject, part$treatment))
part$diff10 <- part$Ch6DiffInside - part$Ch6DiffSupply
part$diff510 <- part$Ch5DiffInside - part$Ch5DiffSupply
part$diff35 <- part$Ch4DiffInside - part$Ch4DiffSupply
part <- part[, c('subTreat', 'diff10', 'diff510', 'diff35')]

#### 3 minute delay?
# plot(c(0, part$Ch6DiffInside) ~ c(part$Ch6DiffSupply, 0))
# test <- part[part$subTreat == 's1high', ]
# cor.test(test$Ch6DiffInside, test$Ch6DiffSupply)
# cor.test(test$Ch6DiffInside[2:90], test$Ch6DiffSupply[1:89])
# cor.test(test$Ch6DiffInside[3:90], test$Ch6DiffSupply[1:88])
# cor.test(test$Ch6DiffInside[4:90], test$Ch6DiffSupply[1:87])
# cor.test(test$Ch6DiffInside[5:90], test$Ch6DiffSupply[1:86])


dim(part)




cols <- rep('tomato', 10)
cols[c(1, 5, 6, 8)] <- 'cornflowerblue'
#males <- c(1, 5, 6, 8)
#females <- c(2, 3, 4, 7)

for (i in 1:nlevels(part$subTreat)) {
	lev <- levels(part$subTreat)[i]
	assign(lev, part[part$subTreat == lev, ])
	}
	
list.Low <- ls(pattern = 'low')
list.High <- ls(pattern = 'high')

for (i in 1:length(list.Low)) {
	print(dim(get(list.Low[i])))
	}
for (i in 1:length(list.High)) {
	print(dim(get(list.High[i])))
	}

plot(s1low$diff10, type='l')
plot(s1high$diff10, type='l')

## Low
Agg18 <- rep(letters[1:18], each = 5)
for (i in 1:length(list.Low)) {
	name <- paste("agg.", list.Low[i], sep = "")
	assign(name, aggregate(get(list.Low[i])[, -1], by = list(Agg18), FUN = "sum"))
}

Agg.df <- data.frame(matrix(0, 0, 4))
Aggs <- ls(pattern = "agg.")
for (i in 1:length(Aggs)) {
	Agg.tmp <- get(Aggs[i])
	print(get(Aggs[i]))
	Agg.df <- data.frame(rbind(Agg.df, 
		data.frame(rep(letters[i], nrow(Agg.tmp)), 
			Agg.tmp[, 2:4])))

}
names(Agg.df)[1] <- 'group'
Agg.df

rm(list=ls(pattern='agg.'))

## High
Agg18 <- rep(letters[1:18], each = 5)
for (i in 1:length(list.High)) {
	name <- paste("agg.", list.High[i], sep = "")
	assign(name, aggregate(get(list.High[i])[, -1], by = list(Agg18), FUN = "sum"))
}

Agg.df.High <- data.frame(matrix(0, 0, 4))
Aggs <- ls(pattern = "agg.")
for (i in 1:length(Aggs)) {
	Agg.tmp <- get(Aggs[i])
	print(get(Aggs[i]))
	Agg.df.High <- data.frame(rbind(Agg.df.High, 
		data.frame(rep(letters[i], nrow(Agg.tmp)), 
			Agg.tmp[, 2:4])))

}
names(Agg.df.High)[1] <- 'group'
Agg.df.High


dev.new()
boxplot(Agg.df$diff10 ~ Agg.df$group, cex=.6, col=cols, pch=16, 
	las=1, xaxt='n', 
	ylab='')
mtext(expression(paste('Particles (>10', mu, 'm) L' ^{-1},  ' min' ^{-1})), 
	side=2, line=2)
abline(h=0, col='gray', lty=3, lwd=2)
mtext('People', side=1)
mtext(c('Men', 'Women'), side=1, at=c(3.5, 5.5), line=1.5, font=2, cex=1.4,
	col=c('cornflowerblue', 'tomato'))


dev.new()
boxplot(Agg.df.High$diff10 ~ Agg.df.High$group, cex=.6, col=cols, pch=16, 
	las=1, xaxt='n', 
	ylab='')
mtext(expression(paste('Particles (>10', mu, 'm) L' ^{-1},  ' min' ^{-1})), 
	side=2, line=2)
abline(h=0, col='gray', lty=3, lwd=2)
mtext('People', side=1)
mtext(c('Men', 'Women'), side=1, at=c(3.5, 5.5), line=1.5, font=2, cex=1.4,
	col=c('cornflowerblue', 'tomato'))

dev.new()
boxplot(Agg.df$diff510 ~ Agg.df$group, cex=.6, col=cols, pch=16, 
	las=1, xaxt='n', 
	ylab='')
mtext(expression(paste('Particles (>10', mu, 'm) L' ^{-1},  ' min' ^{-1})), 
	side=2, line=2)
abline(h=0, col='gray', lty=3, lwd=2)
mtext('People', side=1)
mtext(c('Men', 'Women'), side=1, at=c(3.5, 5.5), line=1.5, font=2, cex=1.4,
	col=c('cornflowerblue', 'tomato'))


dev.new()
boxplot(Agg.df.High$diff510 ~ Agg.df.High$group, cex=.6, col=cols, pch=16, 
	las=1, xaxt='n', 
	ylab='')
mtext(expression(paste('Particles (>10', mu, 'm) L' ^{-1},  ' min' ^{-1})), 
	side=2, line=2)
abline(h=0, col='gray', lty=3, lwd=2)
mtext('People', side=1)
mtext(c('Men', 'Women'), side=1, at=c(3.5, 5.5), line=1.5, font=2, cex=1.4,
	col=c('cornflowerblue', 'tomato'))







