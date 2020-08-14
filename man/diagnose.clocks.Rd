\name{diganose.clocks}
\alias{diganose.clocks}
\title{
diganose.clocks
}
\description{
Perform a description of clocks signals across loci and lineages.
}
\usage{
diagnose.clocks(loctrs, sptr, sp.time.tree = T, branch.support.threshold = 0, log.branches = F, pca = T, mds = F, pca.permutations = 100, mean.scaling.brlen = 0.05, ncore = 1, make.plots = F, sammon.correction = F, clustering.boot.samps = 50, kmax = 10, other.data = NULL, pdf.file = "clock.diagnosis")
}
\arguments{
  \item{loctrs}{
A 'list' or 'multiPhylo' with the locus trees (class 'phylo') with branch lengths measured in expected number of substitutions per site
}
  \item{sptr}{
A tree of class 'phylo' representing the species tree of all taxa included in gene trees
}
  \item{sp.time.tree}{
A logical value (T / F) specifying whether the branches from the species tree should be taken to be calculated in time, in which case locus branch lengths will divided by species tree branch lengths and taken to be rates
}
  \item{branch.support.threshold}{
A numeric specifying the minimum branch support in locus trees required for a branch to be included in clocks analyses
}
  \item{log.branches}{
A logical value (T / F) specifying whether rates (or locus branch lengths) should be log-transformed
}
  \item{pca}{
A logical value (T / F) specifying whether to describe the data of rates by performing dimensionality reduction using principal components analysis
}
  \item{mds}{
A logical value (T / F) specifying whether to describe the data of rates by performing dimensionality reduction using multi-dimensional scaling
}
  \item{pca.permutations}{
A numeric value specifying the number of permutations of rates data for each branch for testing principal components for their power to explain the data, and the contribution of branches to each principal component (ignored if pca = F)
}
  \item{mean.scaling.brlen}{
A numeric value specifying the mean branch length to be used when performing dimensionality reduction using multi-dimensional scaling (ignored if mds = F)
}
  \item{ncore}{
A numeric value specifying the number of computer cores to use when calculating the matrix used for multi-dimensional scaling (ignored if mds = F)
}
  \item{make.plots}{
A logical value (T / F) specifying whether to plot results, and save PDFs with names according to the 'pdf.file' argument
}
  \item{sammon.correction}{
A logical value (T / F) indicating whether to perform Sammon's non-linear mapping correction to multi-dimensional scaling analysis (ignored if mds = F)
}
  \item{clustering.boot.samps}{
A numeric value specifying the number of bootstrap replicates to be done for selecting the number of clusters in the data
}
  \item{kmax}{
A numeric value specifying the maximum number of clusters to be considered in the data
}
  \item{other.data}{
A 'data.frame' where each column contains a 'numeric' or 'factor' variable for each locus, which will then be added to PDFs of clock space (ignored if make.plots = F)
}
  \item{pdf.file}{
A 'character' string specifying the prefix of the PDF files saved, including any file paths if desired
}


}
\details{
The function plots...
}
\value{
A list with two items:
  \item{roottotippath }{A vector with the distance from the root to the tips. If tipsonly==F, it will also include the distances for internal nodes.}
  \item{nodesinpath }{A vector with the number of nodes along a lineage. If tipsonly==F, it will also include the number of nodes for internal nodes. These are counted from the tips to the root.}
}
\references{
Venditti, Chris, Andrew Meade, and Mark Pagel. "Detecting the node-density artifact in phylogeny reconstruction." Systematic biology 55.4 (2006): 637-643.
 
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