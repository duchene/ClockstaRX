\name{bsd.dist}
\alias{bsd.dist}
\title{
bsd.dist
}
\description{
This function calculates the distance between two trees in terms of their branch lengths, assuming identical topologies.
}
\usage{
bsd.dist(tree1, tree2, mean.brlen = 0.05)
}
\arguments{
  \item{tree1}{
A phylogenetic tree of class 'phylo'.
}
  \item{tree2}{
A phylogenetic tree of class 'phylo'.
}
  \item{mean.brlen}{
Value used for scaling the trees to make them comparable, which can be arbitrary.
}
}
\details{
The function uses the branch score distance (Kuhner & Felsenstein 1994) between two trees, but first rescaling them to be as similar to each other as possible.
}
\value{
A numeric value that is the BSDmin distance.
}
\references{
Kuhner, M. K. and Felsenstein, J. (1994) A simulation comparison of phylogeny algorithms under equal and unequal evolutionary rates, Molecular Biology and Evolution, *11(3)*, 459-468
}
\author{
David Duchene
}

\seealso{
quart.scores
}
\examples{
set.seed(12345)
myTree1 <- rtree(20)
myTree2 <- myTree1
myTree2$edge.length <- runif(Nedge(myTree2), min = 0, max = 10)
branchlengthdist <- bsd.dist(myTree1, myTree2)
branchlengthdist



## The function is currently defined as
function (tree1, tree2, mean.brlen = 0.05)
{
    if (sum(tree1$edge.length) < sum(tree2$edge.length)) {
        tree3 <- tree1
        tree1 <- tree2
        tree2 <- tree3
    }
    tree.dist.opt <- function(x) {
        tree3 <- tree2
        tree3$edge.length <- tree2$edge.length * x
        return(KF.dist(tree1, tree3))
    }
    opt.dist <- optim(0, fn = tree.dist.opt, method = "Brent", lower = 0, upper = 50)
    min.bdi <- opt.dist$value
    scaling <- opt.dist$par
    tree2.scaled <- tree2
    tree2.scaled$edge.length <- tree2$edge.length * scaling
    root.scaling <- mean.brlen/mean(c(tree1$edge.length[tree1$edge.length > 1e-05], tree2.scaled$edge.length[tree2.scaled$edge.length > 1e-05]))
    tree1.root.scaled <- tree1
    tree2.root.scaled <- tree2.scaled
    tree1.root.scaled$edge.length <- tree1$edge.length * root.scaling
    tree2.root.scaled$edge.length <- tree2.scaled$edge.length * root.scaling
    min.bdi.root.scaled <- KF.dist(tree1.root.scaled, tree2.root.scaled)
    return(min.bdi.root.scaled)
}
}
\keyword{ branch length distance }