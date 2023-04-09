collect.clocks <- function(loctrs, sptr, sp.time.tree = T, branch.support.threshold = 0, branch.length.threshold = 5e-6, verbose = T){
	require(phangorn)
	if(is.rooted(sptr)) sptr <- unroot(sptr)

	# Collect rates
	ratesmat <- do.call(rbind, lapply(loctrs, function(x) get.locus.rates(x, sptr = sptr, branch.support.threshold = branch.support.threshold, branch.length.threshold = branch.length.threshold)))
	ratesmat <- apply(as.data.frame(ratesmat), 2, as.numeric)
	if(sp.time.tree) for(i in 1:Nedge(sptr)) ratesmat[,i] <- ratesmat[,i] / sptr$edge.length[i]
	
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
	reslist <- new("clockstarx", crx_list = reslist)
	return(reslist)
}
