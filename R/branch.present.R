# Missing exclusion of branches by branch support

branch.present <- function(locus.tree, sp.tree, sptredge = 1, branch.support.threshold = 0, branch.length.threshold = 1e-6){
	
	subrootchildren <- Descendants(sp.tree, sp.tree$edge[sptredge, 1], type = "children")
	
	all.sisters <- sp.tree$tip.label[Descendants(sp.tree, sp.tree$edge[sptredge, 1], type = "tips")[[1]]]
        other.sister <- all.sisters[-which(all.sisters %in% sp.tree$tip.label[Descendants(sp.tree, sp.tree$edge[sptredge, 2], type = "tips")[[1]]])]
	
	if(!any(other.sister %in% locus.tree$tip.label)) return(NA)
        if(!is.monophyletic(locus.tree, locus.tree$tip.label[which(locus.tree$tip.label %in% all.sisters)])) return(NA)
        if(!is.monophyletic(locus.tree, locus.tree$tip.label[which(locus.tree$tip.label %in% other.sister)])) return(NA)
	
	if(length(subrootchildren) == 3){
		otherchildren <- subrootchildren[-which(subrootchildren == sp.tree$edge[sptredge, 2])]
		other.sister1 <- sp.tree$tip.label[Descendants(sp.tree, otherchildren[1], type = "tips")[[1]]]
		if(!any(other.sister1 %in% locus.tree$tip.label)) return(NA)
		other.sister2 <- sp.tree$tip.label[Descendants(sp.tree, otherchildren[2], type = "tips")[[1]]]
		if(!any(other.sister2 %in% locus.tree$tip.label)) return(NA)
		if(!is.monophyletic(locus.tree, locus.tree$tip.label[which(locus.tree$tip.label %in% other.sister1)])) return(NA)
		if(!is.monophyletic(locus.tree, locus.tree$tip.label[which(locus.tree$tip.label %in% other.sister2)])) return(NA)
	}
	
	# If branch is a tip
	if(sp.tree$edge[sptredge, 2] <= Ntip(sp.tree)){
		if(!sp.tree$tip.label[sp.tree$edge[sptredge, 2]] %in% locus.tree$tip.label) return(NA)
		tipname <- sp.tree$tip.label[sp.tree$edge[sptredge, 2]]
		edgelen <- locus.tree$edge.length[which(locus.tree$edge[,2] == which(locus.tree$tip.label == tipname))]
		return(edgelen)
	}
	
	sisteredges <- Descendants(sp.tree, sp.tree$edge[sptredge, 2], type = "children")
	spedgetaxa1 <- sp.tree$tip.label[Descendants(sp.tree, sisteredges[1], type = "tips")[[1]]]
	spedgetaxa2 <- sp.tree$tip.label[Descendants(sp.tree, sisteredges[2], type = "tips")[[1]]]
	if(any(spedgetaxa1 %in% locus.tree$tip.label) && any(spedgetaxa2 %in% locus.tree$tip.label)){
		loctax1 <- locus.tree$tip.label[which(locus.tree$tip.label %in% spedgetaxa1)]
		loctax2 <- locus.tree$tip.label[which(locus.tree$tip.label %in% spedgetaxa2)]
		if(is.monophyletic(locus.tree, locus.tree$tip.label[which(locus.tree$tip.label %in% c(loctax1, loctax2))])){
			if(!is.monophyletic(locus.tree, loctax1) || !is.monophyletic(locus.tree, loctax2)) return(NA)
			
			locus.tree <- root(locus.tree, outgroup = locus.tree$tip.label[which(!locus.tree$tip.label %in% c(loctax1, loctax2))], resolve.root = T)
			
			mrcanode <- getMRCA(locus.tree, c(loctax1, loctax2))
			
			# Remove branch with low support - only takes nodelabels that are single values!
        		if(!is.null(locus.tree$node.label)) if(!is.na(suppressWarnings(as.numeric(locus.tree$node.label[mrcanode - Ntip(locus.tree)])))) if(as.numeric(locus.tree$node.label[mrcanode - Ntip(locus.tree)]) < branch.support.threshold) return(NA)
			
			edgelen <- locus.tree$edge.length[which(locus.tree$edge[, 2] == mrcanode)]
			if(length(edgelen) == 0){
				alldaughters <- Descendants(locus.tree, mrcanode, "children")
				mrcanode <- alldaughters[which(alldaughters != getMRCA(locus.tree, loctax1) && alldaughters != getMRCA(locus.tree, loctax2))]
				edgelen <- locus.tree$edge.length[which(locus.tree$edge[, 2] == mrcanode)]
			}
			if(length(edgelen) < branch.length.threshold) return(NA)
			return(edgelen)
		} else {
		        return(NA)
		}
	} else {
	        return(NA)
	}
}