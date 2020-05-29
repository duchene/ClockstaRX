collect.clocks <- function(loctrs, sptr, branch.support.threshold = 0, verbose = T){
	require(phangorn)
	if(is.rooted(sptr)) sptr <- unroot(sptr)

	# Collect rates
	ratesmat <- do.call(rbind, lapply(loctrs, function(x) get.locus.rates(x, sptr = sptr, branch.support.threshold = branch.support.threshold)))
	ratesmat <- apply(as.data.frame(ratesmat), 2, as.numeric)
	
	# Collect number of samples per branch
	locN <- apply(ratesmat, 2, function(x) sum(!is.na(x)))
	
	# Return rates and N for each species tree branch
	if(verbose) cat("Output includes:", fill = T)
	# Export tree with missing data as annotations
	trN <- sptr
	trN$edge.length <- locN
	trlen <- sptr
	trlen$edge.length <- apply(ratesmat, 2, function(x) if(all(is.na(x))) 0 else median(x, na.rm = T))
	reslist <- list(raw.rates.matrix = ratesmat, N.loci.per.branch = locN, N.samples.tree = trN, median.clock.tree = trlen)
	if(verbose) for(i in 1:length(reslist)) cat(paste0(i, ". ", names(reslist)[i]), fill = T)
	return(reslist)
}