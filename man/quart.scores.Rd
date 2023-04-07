\name{quart.scores}
\alias{quart.scores}
\title{
quart.scores
}
\description{
This function takes a list of trees and calcuates quartet scores among them, the output being a pairwise distance matrix.
}
\usage{
quart.scores(trs)
}
\arguments{
  \item{trs}{
A list of phylogenetic trees, each of class 'phylo'. The class of the list can be 'phylo' or 'multiPhylo'.
}
}
\details{
The function uses the function get.locus.rates to identify the shared quartets, such that all trees should contain branch lengths.
}
\value{
A pairwise distance matrix of class dist.
}
\author{
David Duchene
}
\seealso{
get.locus.rates
}
\examples{
set.seed(12345)
myTrees <- lapply(1:20, function(x) rtree(20))
mytrQs <- quart.scores(myTrees)
mytrQs


## The function is currently defined as
function (trs)
{
    qscores <- matrix(NA, length(trs), length(trs))
    nedges <- sapply(trs, Nedge)
    for (i in 1:length(trs)) {
        for (j in i:length(trs)) {
            trstemp <- trs[c(i, j)][order(nedges[c(i, j)])]
            brsavail <- get.locus.rates(trstemp[[1]], sptr = trstemp[[2]],
                branch.support.threshold = 0, branch.length.threshold = 0)
            if (all(is.na(brsavail)))
                qscores[j, i] <- 0
            else qscores[j, i] <- sum(!is.na(brsavail))/min(nedges[c(i,
                j)])
        }
    }
    qscores <- as.dist(qscores)
    return(qscores)
}
}
\keyword{ quartet distances }
