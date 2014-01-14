
########### richness, diversity and evenness in this function
Evenness <- function(mat) {
  require(vegan)
  H1 <- diversity(mat)
  R <- specnumber(mat)
  J <- H1 / log(R)
  hrj <- data.frame( H1, R, J )
  invisible(hrj)
}


################## makes taxonomy data.frame from taxon vector, sep='; '
makeTaxo <- function(taxo.in=tax.np, otu.table=pb.3500) {
	taxo.in <- as.character(taxo.in)
	tax.tmp.ls <- strsplit(taxo.in, split='; ')
	tax.lengths <- unlist(lapply(tax.tmp.ls, length))
	max(tax.lengths)
	tax.tmp.ls[[1]][1]
	
	# test
	# x <- c('a; b; c; d', 'e; f; g', 'i; j')
	# x2 <- strsplit(x, '; ')
	# x3 <- data.frame(one=sapply(x2, function(x){x[1]}),
					 # two=sapply(x2, function(x){x[2]}),
					 # three=sapply(x2, function(x){x[3]}),
					 # four=sapply(x2, function(x){x[4]}))
	# x3
	# x3$four <- as.character(x3$four)
	# x3$four[which(is.na(x3$four))] <- 'h'
	
	taxo <- data.frame(kingdom=sapply(tax.tmp.ls, function(x){x[1]}),
					   phylum=sapply(tax.tmp.ls, function(x){x[2]}),
					   class=sapply(tax.tmp.ls, function(x){x[3]}),
					   order=sapply(tax.tmp.ls, function(x){x[4]}),
					   family=sapply(tax.tmp.ls, function(x){x[5]}),
					   genus=sapply(tax.tmp.ls, function(x){x[6]}))
	
	taxo$kingdom <- as.character(taxo$kingdom)
	taxo$phylum <- as.character(taxo$phylum)
	taxo$class <- as.character(taxo$class)
	taxo$order <- as.character(taxo$order)
	taxo$family <- as.character(taxo$family)
	taxo$genus <- as.character(taxo$genus)
	
	for (i in 1:ncol(taxo)){
		taxo[which(is.na(taxo[, i])), i] <- '' 
		}
	
	# taxo.all <- taxo # save big one
	taxo$abundance <- colSums(otu.table)
	row.names(taxo) <- colnames(otu.table)
	
	invisible(taxo)
}

#############
# get a vector of names from taxo 

cons <- function(x) {
  l <- length(x)
  while(x[l] == '') {l = l-1}
  name <- x[l]
}
# consensus <- apply(taxo[, 1:6], 1, cons)



########  standard error
se <- function(x){(sd(x)/sqrt(length(x)))}













