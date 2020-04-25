weigh.clocks <- function(ratesmat){
	     bsdmatdim <- nrow(ratesmat)
	     Nbranches <- ncol(ratesmat)
	     meanlengths <- colMeans(ratesmat)
	     
	     wghts <- (meanlengths * Nbranches) / sum(meanlengths)
	     locus.brlens <- do.call(rbind, lapply(1:bsdmatdim, function(x) ratesmat[x,] / wghts))
	     if(!is.null(rownames(ratesmat))) rownames(locus.brlens) <- rownames(ratesmat)
	     return(locus.brlens)
}