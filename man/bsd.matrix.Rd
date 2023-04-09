\name{bsd.matrix}
\alias{bsd.matrix}
\title{
bsd.matrix
}
\description{
This function calculates pairwise matrix of bsd distances, assuming that all locus trees have the same topology as the species tree, or that missing branches have been imputed.
}
\usage{
bsd.matrix(locus.brlens, species.tree, mean.brlen = 0.05, ncore = 1)
}
\arguments{
  \item{locus.brlens}{
A matrix as produced by 'collect.clocks' of locus rates where rows are loci and columns are branches.
}
  \item{species.tree}{
Species tree topology structure of class 'phylo'.
}
  \item{mean.brlen}{
Scaling value to be passed on to bsd.dist.
}
  \item{ncore}{
Numeric indicating how many computer cores to use during the calculation of the matrix.
}
}
\details{
The function produces pairwise distances between loci in terms of their branch lengths. This requires an optimisation of their minimal divergence that can be computationally costly, and parallel computing is therefore advisable. The function is used internally by 'clock.space'.
}
\value{
A matrix of pairwise BSDmin distances among loci.
}
\author{
David Duchene
}

\seealso{
bsd.dist
}

## The function is currently defined as
function(locus.brlens, species.tree, mean.brlen = 0.05, ncore = 1){

        if(!is.null(rownames(locus.brlens))) locnames <- rownames(locus.brlens)

        bsdmatdim <- nrow(locus.brlens)

        trs <- apply(locus.brlens, 1, function(x){
            species.tree$edge.length <- x
            return(species.tree)
        })
        class(trs) <- "multiPhylo"

        if(ncore == 1){
                 bsd.mat <- do.call(rbind, lapply(2:bsdmatdim, function(x) c(sapply(1:(x-1), function(y) bsd.dist(trs[[x]], trs[[y]], mean.brlen = mean.brlen)), rep(NA, length(x:bsdmatdim)))))
        } else if(ncore > 1){
                 # Parallelisaton to be tested
                 require(foreach)
                 require(doParallel)
                 #sub.trees <- lapply(2:bsdmatdim, function(x) trs[1:(x-1)])
                 sub.trees <- list()
                 for(k in 2:length(trs)){
                       sub.trees[[k]] <- trs[1:k-1]
                 }
                 compute.tree.dists <- function(tree.sub.list, fix.tree){
                        res <- sapply(tree.sub.list, function(a) bsd.dist(fix.tree, a))
                        return(res)
                 }
                 cl <- makeCluster(ncore)
                 registerDoParallel(cl)
                 cat(paste("Starting parallel computing with", ncore, "cores"), fill = T)
                 res.par <- foreach(s.trees = sub.trees, j = 1:length(trs), .packages = c("phangorn", "ape"), .export = c("bsd.dist")) %dopar% c(compute.tree.dists(tree.sub.list = s.trees, fix.tree = trs[[j]]), rep(NA, bsdmatdim-j+1))
                 stopCluster(cl)
                 cat("Finished parallel computing", fill = T)
                 bsd.mat <- do.call(rbind, res.par)
        }
		#bsd.mat <- rbind(rep(NA, bsdmatdim), bsd.mat)

        if(!is.null(rownames(locus.brlens))) rownames(bsd.mat) <- colnames(bsd.mat) <- locnames

        #bsd.mat <- as.dist(bsd.mat)

        return(bsd.mat)
}
\keyword{ BSDmin comparisons }