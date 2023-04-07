\name{branch.present}
\alias{branch.present}
\title{
branch.present
}
\description{
This function tests whether a branch in one tree (species tree) is present in the other (locus tree).
}
\usage{
branch.present(locus.tree, sp.tree, sptredge = 1, branch.support.threshold = 0, branch.length.threshold = 5e-06)
}
\arguments{
  \item{locus.tree}{
A phylogenetic tree of class 'phylo'.
}
  \item{sp.tree}{
A phylogenetic tree of class 'phylo'.
}
  \item{sptredge}{
The edge in the species tree to be assessed for its presence in the locus tree.
}
  \item{branch.support.threshold}{
Numeric value of branch support as branch annotation below which a branch will be considered absent in the locus tree even if present in the topology structure.
}
  \item{branch.length.threshold}{
Numeric value of branch support as branch annotation below which a branch will be considered absent in the locus tree even if present in the topology structure.
}
}
\details{
The function searches for a branch (quartet) exploring the descendants on both sides of the quartet, allowing for missing taxa.
}
\value{
If the branch is present, the function returns the length of the branch in the locus.tree, and NA otherwise.
}
\author{
David Duchene
}

\seealso{
get.locus.rates
}
\examples{
set.seed(12345)
locus.tree <- rtree(20)
species.tree <- rtree(20)
presence1 <- branch.present(locus.tree, species.tree, sptredge = 5)
presence1
presence2 <- branch.present(locus.tree, species.tree, sptredge = 6)
presence2



## The function is currently defined as
function (locus.tree, sp.tree, sptredge = 1, branch.support.threshold = 0,
    branch.length.threshold = 5e-06)
{
    subrootchildren <- Descendants(sp.tree, sp.tree$edge[sptredge,
        1], type = "children")
    all.sisters <- sp.tree$tip.label[Descendants(sp.tree, sp.tree$edge[sptredge,
        1], type = "tips")[[1]]]
    other.sister <- all.sisters[-which(all.sisters %in% sp.tree$tip.label[Descendants(sp.tree,
        sp.tree$edge[sptredge, 2], type = "tips")[[1]]])]
    if (!any(other.sister %in% locus.tree$tip.label))
        return(NA)
    if (!is.monophyletic(locus.tree, locus.tree$tip.label[which(locus.tree$tip.label %in%
        all.sisters)]))
        return(NA)
    if (!is.monophyletic(locus.tree, locus.tree$tip.label[which(locus.tree$tip.label %in%
        other.sister)]))
        return(NA)
    if (length(subrootchildren) == 3) {
        otherchildren <- subrootchildren[-which(subrootchildren ==
            sp.tree$edge[sptredge, 2])]
        other.sister1 <- sp.tree$tip.label[Descendants(sp.tree,
            otherchildren[1], type = "tips")[[1]]]
        if (!any(other.sister1 %in% locus.tree$tip.label))
            return(NA)
        other.sister2 <- sp.tree$tip.label[Descendants(sp.tree,
            otherchildren[2], type = "tips")[[1]]]
        if (!any(other.sister2 %in% locus.tree$tip.label))
            return(NA)
        if (!is.monophyletic(locus.tree, locus.tree$tip.label[which(locus.tree$tip.label %in%
            other.sister1)]))
            return(NA)
        if (!is.monophyletic(locus.tree, locus.tree$tip.label[which(locus.tree$tip.label %in%
            other.sister2)]))
            return(NA)
    }
    if (sp.tree$edge[sptredge, 2] <= Ntip(sp.tree)) {
        if (!sp.tree$tip.label[sp.tree$edge[sptredge, 2]] %in%
            locus.tree$tip.label)
            return(NA)
        tipname <- sp.tree$tip.label[sp.tree$edge[sptredge, 2]]
        edgelen <- locus.tree$edge.length[which(locus.tree$edge[,
            2] == which(locus.tree$tip.label == tipname))]
        if (edgelen > branch.length.threshold)
            return(edgelen)
        else return(NA)
    }
    sisteredges <- Descendants(sp.tree, sp.tree$edge[sptredge,
        2], type = "children")
    spedgetaxa1 <- sp.tree$tip.label[Descendants(sp.tree, sisteredges[1],
        type = "tips")[[1]]]
    spedgetaxa2 <- sp.tree$tip.label[Descendants(sp.tree, sisteredges[2],
        type = "tips")[[1]]]
    if (any(spedgetaxa1 %in% locus.tree$tip.label) && any(spedgetaxa2 %in%
        locus.tree$tip.label)) {
        loctax1 <- locus.tree$tip.label[which(locus.tree$tip.label %in%
            spedgetaxa1)]
        loctax2 <- locus.tree$tip.label[which(locus.tree$tip.label %in%
            spedgetaxa2)]
        if (is.monophyletic(locus.tree, locus.tree$tip.label[which(locus.tree$tip.label %in%
            c(loctax1, loctax2))])) {
            if (!is.monophyletic(locus.tree, loctax1) || !is.monophyletic(locus.tree,
                loctax2))
                return(NA)
            locus.tree <- root(locus.tree, outgroup = locus.tree$tip.label[which(!locus.tree$tip.label %in%
                c(loctax1, loctax2))], resolve.root = T)
            mrcanode <- getMRCA(locus.tree, c(loctax1, loctax2))
            if (!is.null(locus.tree$node.label))
                if (!is.na(suppressWarnings(as.numeric(locus.tree$node.label[mrcanode -
                  Ntip(locus.tree)]))))
                  if (as.numeric(locus.tree$node.label[mrcanode -
                    Ntip(locus.tree)]) < branch.support.threshold)
                    return(NA)
            edgelen <- locus.tree$edge.length[which(locus.tree$edge[,
                2] == mrcanode)]
            if (length(edgelen) == 0) {
                alldaughters <- Descendants(locus.tree, mrcanode,
                  "children")
                mrcanode <- alldaughters[which(alldaughters !=
                  getMRCA(locus.tree, loctax1) && alldaughters !=
                  getMRCA(locus.tree, loctax2))]
                edgelen <- locus.tree$edge.length[which(locus.tree$edge[,
                  2] == mrcanode)]
            }
            if (edgelen > branch.length.threshold)
                return(edgelen)
            else return(NA)
        }
        else {
            return(NA)
        }
    }
    else {
        return(NA)
    }
}
}
\keyword{ Node-density-effect }