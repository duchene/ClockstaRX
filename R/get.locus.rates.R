get.locus.rates <- function(locustr, sptr, branch.support.threshold = 0, branch.length.threshold = 5e-6){
	if(is.rooted(locustr)) locustr <- unroot(locustr)
	rates <- sapply(1:Nedge(sptr), function(x) branch.present(locus.tree = locustr, sp.tree = sptr, sptredge = x, branch.support.threshold = branch.support.threshold, branch.length.threshold = branch.length.threshold))
	return(rates)
}