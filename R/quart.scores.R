quart.scores <- function(loctrs, sptr){
	     blenstab <- collect.clocks(loctrs, sptr, sp.time.tree = F, branch.support.threshold = 0, branch.length.threshold = 0, verbose = F)
	     qscores <- apply(blenstab[[1]], 1, function(x) sum(!is.na(x))) / Nedge(loctrs)
	     return(qscores)
}