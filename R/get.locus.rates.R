get.locus.rates <- function(locustr, sptr, branch.support.threshold = 0){
	if(is.rooted(locustr)) locustr <- unroot(locustr)
	rates <- sapply(1:Nedge(sptr), function(x) branch.present(locus.tree = locustr, sp.tree = sptr, sptredge = x, branch.support.threshold = branch.support.threshold))
	return(rates)
}