# Calculates a distance matrix. The score ranges from 0 to 1, and takes as reference the smallest tree (i.e., it is the proportion of quartets in the smallest tree that occur in the larger tree). All trees require branch lengths (can be dummy, such as all 1).
quart.scores <- function(trs){
	     qscores <- matrix(NA, length(trs), length(trs))
	     nedges <- sapply(trs, Nedge)
	     for(i in 1:length(trs)){
	     	   for(j in i:length(trs)){
		   	 trstemp <- trs[c(i,j)][order(nedges[c(i,j)])]
			 brsavail <- get.locus.rates(trstemp[[1]], sptr = trstemp[[2]], branch.support.threshold = 0, branch.length.threshold = 0)
			 if(all(is.na(brsavail))) qscores[j,i] <- 0 else qscores[j,i] <- sum(!is.na(brsavail)) / min(nedges[c(i,j)])
		   }
	     }
	     qscores <- as.dist(qscores)
	     return(qscores)
}