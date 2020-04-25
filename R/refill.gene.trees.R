# This function refills gene trees with lengths of branches present in the species tree. The species tree should have branch lengths in terms of rates.

refill.gene.trees <- function(loctrs, sp.phylogram){
	if(is.rooted(sp.phylogram)) sp.phylogram <- unroot(sp.phylogram)
	if(!is.null(names(loctrs))) locnames <- names(loctrs)
	loctrs <- lapply(loctrs, function(x) get.sptree.rates(x, sptr = sp.phylogram))
	if(exists("locnames")) names(loctrs) <- locnames

	# return list of gene.trees
	class(loctrs) <- "multiPhylo"
	return(loctrs)
}