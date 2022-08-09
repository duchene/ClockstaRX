\name{collect.clocks}
\alias{collect.clocks}
\title{
collect.clocks
}
\description{
collect.clocks takes locus trees with branch lengths in substitutions per site and a species tree topology. It collects the lengths of quartet branches in locus trees that are present in the species tree.
}
\usage{
collect.clocks(loctrs, sptr, sp.time.tree = T, branch.support.threshold = 0, branch.length.threshold = 5e-6, verbose = T)
}
\arguments{
  \item{loctrs}{
A list of trees of class 'multiPhylo'
}
  \item{sptr}{
A phylogenetic tree of class 'phylo'.
}
  \item{make.plot}{
A logical (T/F) indicating whether the function should make a plot summarising the branch lengths data collected (this option has not been finished).
}
  \item{branch.support.threshold}{
A numeric indicating the branch support below which branches in locus tree should be considered as absent, even if present in the species tree topology.
}
}
\details{
The function plots...
}
\value{
A list with two items:
  \item{raw.rates.matrix}{A matrix of each of the loci (rows) including the lengths of branches present in the species tree (columns). The matrix might include missing data as NAs in cases where brnaches of the species tree are absent in the locus tree.}
  \item{N.loci.per.branch}{A vector with the number of loci that contained each of the brnaches of the species tree.}
}
\references{
To be added.

}
\author{
David Duchene
}
\note{
Please see the references for more detailed information.
}

\seealso{
Pending.
}
\examples{
set.seed(12345)
myTree <- rtree(10)

par(mfrow = c(1, 2))
plot(myTree)
nde <- pathnode(myTree)
nde



## The function is currently defined as
function (phylo, tipsonly = T) 
{
    require(phangorn)
    di.tr <- dist.nodes(phylo)
    root.tr <- phylo$edge[, 1][!(phylo$edge[, 1] \%in\% phylo$edge[, 
        2])][1]
    tr.depth <- max(di.tr[as.numeric(colnames(di.tr)) == root.tr, 
        ])
    if (tipsonly == TRUE) {
        roottotippath <- di.tr[as.numeric(rownames(di.tr)) == 
            root.tr, 1:length(phylo$tip.label)]
        nodesinpath <- sapply(1:length(phylo$tip.label), function(x) length(Ancestors(phylo, 
            x)))
    }
    else {
        roottotippath <- di.tr[as.numeric(rownames(di.tr)) == 
            root.tr, ]
	nodesinpath <- sapply(1:(length(phylo$tip.label)+phylo$Nnode), function(x) length(Ancestors(phylo, x)))
    }
    plot(roottotippath, nodesinpath, xlab = "Root-to-tip path length", 
        ylab = "Number of parent nodes", pch = 20)
    return(list(roottotippath = roottotippath, nodesinpath = nodesinpath))
  }
}
\keyword{ Node-density-effect }
