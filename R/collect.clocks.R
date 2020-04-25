# Collect branch lengths from gene trees if present in species tree.

collect.clocks <- function(loctrs, sptr, make.plot = T, branch.support.threshold = 0){
	require(phangorn)
	if(is.rooted(sptr)) sptr <- unroot(sptr)

	# Collect rates
	ratesmat <- do.call(rbind, lapply(loctrs, function(x) get.locus.rates(x, sptr = sptr, branch.support.threshold = branch.support.threshold)))
	ratesmat <- apply(as.data.frame(ratesmat), 2, as.numeric)
	
	# Collect number of samples per branch
	locN <- apply(ratesmat, 2, function(x) sum(!is.na(x)))
	
	# Consider making plot showing rates distributions across branches.
	#if(make.plot){
	#	dev.new()
	#}

	# Return rates and N for each species tree branch
	reslist <- list(raw.rates.matrix = ratesmat, N.loci.per.branch = locN)
	return(reslist)
}