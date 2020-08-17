# Replace all branches in gene tree with those in species tree where present.

get.sptree.rates <- function(locustr, sptr, branch.support.threshold = 0, branch.length.threshold = 5e-6){
	if(is.rooted(locustr)) locustr <- unroot(locustr)
	if(is.rooted(sptr)) sptr <- unroot(sptr)

	rates <- sapply(1:nrow(locustr$edge), function(x) branch.present(locus.tree = sptr, sp.tree = locustr, sptredge = x, branch.support.threshold = branch.support.threshold, branch.length.threshold = branch.length.threshold))

	locustr$edge.length[!is.na(rates)] <- rates[!is.na(rates)]	
	
	# Return gene tree
	return(locustr)
}