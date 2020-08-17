\name{diganose.clocks}
\alias{diganose.clocks}
\title{
diganose.clocks
}
\description{
Perform a description of clocks signals across loci and lineages.
}
\usage{
diagnose.clocks(loctrs, sptr, sp.time.tree = T, branch.support.threshold = 0, branch.length.threshold = 1e-6, log.branches = F, pca = T, mds = F, pca.permutations = 100, mean.scaling.brlen = 0.05, ncore = 1, make.plots = F, sammon.correction = F, clustering.boot.samps = 50, kmax = 10, other.data = NULL, pdf.file = "clock.diagnosis")
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
  \item{branch.length.threshold}{
A numeric specifying the minimum branch length in locus trees required for a branch to be included in clocks analyses
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
A NULL or 'character' string specifying if NULL not to save summary PDFs, or if character specifying the prefix of the PDF files saved including any file paths if desired
}


}
\details{
The function the full set of steps in ClockstaRX, serving as a wrapper for the other primary functions. In brief, the steps are to (1) collect branch lengths from loci if they appear in the species tree and divide them by time if an independently-estimated species time-tree is available (collect.clocks function); (2) reduce the dimensionality of clock space using PCA and MDS, making test for the significance of principal components and the contributions of branches to each in the case of PCA (clock.space); (3) perform clustering of loci in the first two dimensions of dimensionality-reduced data (group.clocks); write PDFs of several summary results, including user input data for visual exploration of clock space (write.clocks.plots).
}
\value{
A list with the following items (examples of PCA included but not given if pca = F, some other possible results not shown):
  \item{raw.rates.matrix }{A matrix of the branch lengths collected from gene tree, including NAs as missing data in the case of species tree branches that were not found in the gene tree.}
  \item{N.loci.per.branch }{A vector with the number of gene trees found to contain each of the branches in the species tree.}
  \item{N.samples.tree }{A phylo tree object identical to the species tree where branch lengths are the number of gene trees found to contain each of the branches in the species tree.}
  \item{median.clock.tree }{A phylo tree object identical to the species tree where branch  lengths are the mean length found in gene trees.}
  \item{imputed.clocks }{A matrix of the branch lengths or rates collected with NA values filled using the mean across the branch, or across the whole data set if <2 samples are available in the branch.}
  \item{weighted.imputed.clocks }{A matrix identical to imputed.clocks, but with data weighted by the mean of each branch length as proposed by Bedford and Hartl (2008).}
  \item{pca.clock.space }{A list including PCA results of the matrix given in imputed.clocks, and results of permutation tests on the PCA including p-values for the full analysis (pPhi and pPsi), each PC (pPC), and each branch within each PC (pIL), as described in Björklund (2019) and Vieira (2012).}
  \item{weighted.pca.clock.space }{As pca.clock.space but applied to the matrix in weighted.imputed.clocks.}
  \item{pca.cluster.support }{Results of the clusGap function in the package cluster, showing the support for each number of clusters as tested for the data from pca.clock.space.}
  \item{weighted.pca.cluster.support }{As pca.cluster.support but applied to weighted.pca.clock.space.}
  \item{pca.clustering }{A matrix showing the cluster for each locus (rows) for each different number of clusters tested (columns) as applied to the data in pca.clock.space.}
  \item{weighted.pca.clustering }{As pca.clustering but applied to weighted.pca.clock.space.}
  \item{pca.best.k }{A numeric showing the best-supported number of clusters, where the criterion is the following: If there is a significant decline in gap from one cluster to two, then there is a single cluster; otherwise, the largest increase in gap leads to the optimum number of clusters.}
  \item{weighted.pca.best.k }{As pca.best.k but applied to weighted.pca.clock.space.}
  \item{pca.influential.branches }{A matrix (0/1) showing which branches (rows) significantly contribute to the variation in each principal component (columns), where 1 indicates a significant contribution.}
  \item{weighted.pca.influential.branches }{As pca.influential.branches but applied ot weighted.pca.clock.space.}
  \item{variable.colours }{A matrix with the colours defined for each locus (rows) and each of the default variables explored (number of taxa, branch support, tree length, etc.), and any of the variables given by the user with the other.data argument.}
}
\references{
Bedford, T., & Hartl, D. L. (2008). Overdispersion of the molecular clock: temporal variation of gene-specific substitution rates in Drosophila. Molecular biology and evolution, 25(8), 1631-1638.

Björklund, M. (2019). Be careful with your principal components. Evolution, 73(10), 2151-2158.

Vieira, V. N. C. S. (2012). Permutation tests to estimate significances on Principal Components Analysis. Computational Ecology and Software, 2(2), 103. 
}
\author{
David Duchene
}
\note{
Please see the references for more detailed information.
}

\seealso{
\code{\link[ClockstaRX]{collect.clocks}}, \code{\link[ClockstaRX]{clock.space}}, \code{\link[ClockstaRX]{group.clocks}}.
}
\examples{
set.seed(12345)
myTree <- rtree(10)

par(mfrow = c(1, 2))
plot(myTree)
nde <- pathnode(myTree)
nde

\keyword{ Node-density-effect }