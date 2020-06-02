# Plot basic data in clockspace, e.g., number of taxa, branch support, tree length, topological distance from species tree (plus other user-selected variables). other.data must be a matrix or data.frame where each column is a numeric variable.

write.clocks.plots <- function(groupclocks, loctrs, sptr, other.data = NULL, default.variables = c("ntaxa", "nmissingbr", "brsup", "trlen", "rttdist", "topdist"), pdf.file = "clock.diagnosis"){
	
	# Extract default variables and create their colours
	
	if("bsd.clock.space" %in% names(groupclocks)) mds <- T else mds <- F
        if("pca.clock.space" %in% names(groupclocks)) pca <- T else pca <- F
	
	if(mds & pca){
	       varstoplot <- data.frame(as.numeric(groupclocks$bsd.clustering[[1]][, groupclocks$bsd.best.k[1]]), as.numeric(groupclocks$weighted.bsd.clustering[[1]][, groupclocks$weighted.bsd.best.k[1]]), as.numeric(groupclocks$pca.clustering[[1]][, groupclocks$pca.best.k[1]]), as.numeric(groupclocks$weighted.pca.clustering[[1]][, groupclocks$weighted.pca.best.k[1]]))
	       colnames(varstoplot) <- c("Relative rates clusters (MDS)", "Weighted rates clusters (MDS)", "Relative rates clusters (PCA)", "Weighted rates clusters (PCA)")
	} else if(mds){
	       varstoplot <- data.frame(as.numeric(groupclocks$bsd.clustering[[1]][, groupclocks$bsd.best.k[1]]), as.numeric(groupclocks$weighted.bsd.clustering[[1]][, groupclocks$weighted.bsd.best.k[1]]))
               colnames(varstoplot) <- c("Relative rates clusters (MDS)", "Weighted rates clusters (MDS)")
	} else if(pca){
	       varstoplot <- data.frame(as.numeric(groupclocks$pca.clustering[[1]][, groupclocks$pca.best.k[1]]), as.numeric(groupclocks$weighted.pca.clustering[[1]][, groupclocks$weighted.pca.best.k[1]]))
	       colnames(varstoplot) <- c("Relative rates clusters (PCA)", "Weighted rates clusters (PCA)")
	}
	

	colnames(varstoplot) <- c("Relative rates clusters", "Weighted rates clusters")
		
	ntaxa <- sapply(loctrs, Ntip)
	if("ntaxa" %in% default.variables) varstoplot["N taxa"] <- topo.colors(10)[as.numeric(cut(ntaxa, breaks = 10))]
	
	nmissingbr <- sapply(1:length(loctrs), function(x) sum(is.na(groupclocks[[1]][x,])) / Nedge(loctrs[[x]]))
	if("nmissingbr" %in% default.variables) varstoplot["N missing species tree branches"] <- topo.colors(10)[as.numeric(cut(nmissingbr, breaks = 10))]

	brsup <- sapply(loctrs, function(x) mean(as.numeric(x$node.label), na.rm = T))
	if("brsup" %in% default.variables) varstoplot["Branch support"] <- topo.colors(10)[as.numeric(cut(brsup, breaks = 10))]

	trlen <- sapply(loctrs, function(x) sum(x$edge.length, na.rm = T) / Ntip(x))
	if("trlen" %in% default.variables) varstoplot["Tree length"] <- topo.colors(10)[as.numeric(cut(trlen, breaks = 10))]
	
	rtts <- lapply(loctrs, function(x) pathnode(midpoint(x))[[1]])
	rttdist <- sapply(rtts, function(x) sd(x, na.rm = T) / mean(x, na.rm = T))
	if("rttdist" %in% default.variables) varstoplot["Midpoint-root-to-tip CV"] <- topo.colors(10)[as.numeric(cut(rttdist, breaks = 10))]

	topdist <- sapply(loctrs, function(x) RF.dist(x, drop.tip(sptr, sptr$tip.label[which(!sptr$tip.label %in% x$tip.label)])))
	if("topdist" %in% default.variables) varstoplot["Distance to species tree"] <- topo.colors(10)[as.numeric(cut(topdist, breaks = 10))]

	# Create colours for other variables
	
	uniquevals <- list()
	
	if(!is.null(other.data)){
		for(i in 1:ncol(other.data)){
		      if(class(other.data[,i]) == "numeric"){
				varstoplot[colnames(other.data)[i]] <- topo.colors(10)[as.numeric(cut(other.data[,i], breaks = 10))]
		      } else {
				other.data[,i] <- as.character(other.data[,i])
				uniquevals[[ncol(varstoplot) + 1]] <- unique(other.data[,i])
				varcols <- other.data[,i]
				for(j in 1:length(uniquevals[[ncol(varstoplot) + 1]])) varcols[which(other.data[,i] == uniquevals[[ncol(varstoplot) + 1]][j])] <- as.numeric(j + 1)
				varstoplot <- cbind(varstoplot, data.frame(as.numeric(varcols)))
				colnames(varstoplot)[ncol(varstoplot)] <- colnames(other.data)[i]
		      }
		}
	}
	
	# Create pdf for MDS:
	
	if(mds){

	pdf(paste0(pdf.file, ".MDS.pdf"), useDingbats = F, width = 7, height = if(pca) (length(varstoplot) - 4) * 4 else (length(varstoplot) - 2) * 4)
        if(pca) par(mfrow = c(length(varstoplot) - 3, 2)) else par(mfrow = c(length(varstoplot) - 1, 2))
	
	# First plot the clustering results
	plot(groupclocks$bsd.clock.space[[1]]$points, pch = 19, col = as.numeric(groupclocks$bsd.clustering[[1]][, groupclocks$bsd.best.k[1]]) + 1, xlab = "Dim 1", ylab = "Dim 2", main = "Relative rates clusters (MDS)")
	legend("topright", legend = paste0("k = ", groupclocks$bsd.best.k[1]))
	plot(groupclocks$weighted.bsd.clock.space[[1]]$points, pch = 19, col = as.numeric(groupclocks$weighted.bsd.clustering[[1]][, groupclocks$weighted.bsd.best.k[1]]) + 1, xlab = "Dim 1", ylab = "Dim 2", main = "Residual rates clusters (MDS)")
	legend("topright", legend = paste0("k = ", groupclocks$weighted.bsd.best.k[1]))
	
	# Plot the new variable results
	if(pca) colstart <- 5 else colstart <- 3
	for(i in colstart:length(varstoplot)){
              plot(groupclocks$bsd.clock.space[[1]]$points, pch = 19, col = varstoplot[,i], xlab = "Dim 1", ylab = "Dim 2", main = paste(colnames(varstoplot)[i], "\n(relative rates)"))
	      if(length(uniquevals) > 0) if(length(uniquevals) >= i) if(!is.null(uniquevals[[i]])) legend("topleft", legend = uniquevals[[i]], col = 2:(length(uniquevals[[i]]) + 1), pch = 19)
              plot(groupclocks$weighted.bsd.clock.space[[1]]$points, pch = 19, col = varstoplot[,i], xlab = "Dim 1", ylab = "Dim 2", main = paste(colnames(varstoplot)[i], "\n(residual rates)"))
	}

	dev.off()
	
	}

	if(pca){
	
	pdf(paste0(pdf.file, ".PCA.pdf"), useDingbats = F, width = 7, height = if(mds) (length(varstoplot) - 4) * 4 else (length(varstoplot) - 2) * 4)
        if(mds) par(mfrow = c(length(varstoplot) - 3, 2)) else par(mfrow = c(length(varstoplot) - 1, 2))
	
	# First plot the clustering results
	topload <- head(order(groupclocks$pca.clock.space[[1]][[2]][,1]), 5)
        ordpc2 <- order(groupclocks$pca.clock.space[[1]][[2]][,2])
        topload <- c(topload, head(ordpc2[-which(ordpc2 %in% topload)], 5))
	lx <- cbind(groupclocks$pca.clock.space[[1]][[2]][topload,1])
	ly <- cbind(groupclocks$pca.clock.space[[1]][[2]][topload,2])
	rownames(lx) <- rownames(ly) <- gsub("V", "br ", rownames(lx))
	if(any(abs(c(lx, ly)) < 1e-4)){
		nonzeroloads <- which(abs(lx) > 1e-4 & abs(ly) > 1e-4)
		lx <- cbind(lx[nonzeroloads,])
		ly <- cbind(ly[nonzeroloads,])
	}
	plot(groupclocks$pca.clock.space[[1]][[1]][,1], groupclocks$pca.clock.space[[1]][[1]][,2], pch = 19, col = as.numeric(groupclocks$pca.clustering[[1]][, groupclocks$pca.best.k[1]]) + 1, xlab=paste0("PC 1 (", round(groupclocks$pca.clock.space[[1]][[3]]$importance[2]*100, 1), "%)"), ylab=paste0("PC 2 (", round(groupclocks$pca.clock.space[[1]][[3]]$importance[5]*100, 1), "%)"), main = "Relative rates clusters\nand top branch loadings (PCA)", xlim = c(min(c(range(lx), range(groupclocks$pca.clock.space[[1]][[1]][,1]))), max(c(range(lx), range(groupclocks$pca.clock.space[[1]][[1]][,1])))), ylim = c(min(c(range(ly), range(groupclocks$pca.clock.space[[1]][[1]][,2]))), max(c(range(ly), range(groupclocks$pca.clock.space[[1]][[1]][,2])))))
	legend("topright", legend = paste0("k = ", groupclocks$pca.best.k[1]))
	abline(v=0, lty=2, col="grey50"); abline(h=0, lty=2, col="grey50")
	arrows(x0=0, x1=lx, y0=0, y1=ly, length=0.1, lwd=1)
	text(lx, ly, labels=rownames(lx), pos = c(3,1,3,1), offset=0.3, cex=1)
	
	topload <- head(order(groupclocks$weighted.pca.clock.space[[1]][[2]][,1]), 5)
	ordpc2 <- order(groupclocks$weighted.pca.clock.space[[1]][[2]][,2])
	topload <- c(topload, head(ordpc2[-which(ordpc2 %in% topload)], 5))
	lxW <- cbind(groupclocks$weighted.pca.clock.space[[1]][[2]][topload,1])
    lyW <- cbind(groupclocks$weighted.pca.clock.space[[1]][[2]][topload,2])
    rownames(lxW) <- rownames(lyW) <- gsub("V", "br ", rownames(lxW))
    if(any(abs(c(lxW, lyW)) < 1e-4)){
                nonzeroloads <- which(abs(lxW) > 1e-4 & abs(lyW) > 1e-4)
                lxW <- cbind(lxW[nonzeroloads,])
                lyW <- cbind(lyW[nonzeroloads,])
    }
    plot(groupclocks$weighted.pca.clock.space[[1]][[1]][,1], groupclocks$weighted.pca.clock.space[[1]][[1]][,2], pch = 19, col = as.numeric(groupclocks$weighted.pca.clustering[[1]][, groupclocks$weighted.pca.best.k[1]]) + 1, xlab=paste0("PC 1 (", round(groupclocks$weighted.pca.clock.space[[1]][[3]]$importance[2]*100, 1), "%)"), ylab=paste0("PC 2 (", round(groupclocks$weighted.pca.clock.space[[1]][[3]]$importance[5]*100, 1), "%)"), main = "Residual rates clusters\nand top branch loadings (PCA)", xlim = c(min(c(range(lxW), range(groupclocks$weighted.pca.clock.space[[1]][[1]][,1]))), max(c(range(lxW), range(groupclocks$weighted.pca.clock.space[[1]][[1]][,1])))), ylim = c(min(c(range(lyW), range(groupclocks$weighted.pca.clock.space[[1]][[1]][,2]))), max(c(range(lyW), range(groupclocks$weighted.pca.clock.space[[1]][[1]][,2])))))
	legend("topright", legend = paste0("k = ", groupclocks$weighted.pca.best.k[1]))
	abline(v=0, lty=2, col="grey50"); abline(h=0, lty=2, col="grey50")
    arrows(x0=0, x1=lxW, y0=0, y1=lyW, length=0.1, lwd=1)
    text(lxW, lyW, labels=rownames(lxW), pos = c(3,1,3,1), offset=0.1, cex=1)

	# Plot the new variable results
	if(mds) colstart <- 5 else colstart <- 3
	for(i in colstart:length(varstoplot)){
              plot(groupclocks$pca.clock.space[[1]][[1]], pch = 19, col = varstoplot[,i], main = paste(colnames(varstoplot)[i], "\n(relative rates)"), xlab = "PC 1", ylab = "PC 2", xlim = c(min(c(range(lx), range(groupclocks$pca.clock.space[[1]][[1]][,1]))), max(c(range(lx), range(groupclocks$pca.clock.space[[1]][[1]][,1])))), ylim = c(min(c(range(ly), range(groupclocks$pca.clock.space[[1]][[1]][,2]))), max(c(range(ly), range(groupclocks$pca.clock.space[[1]][[1]][,2])))))
	      abline(v=0, lty=2, col="grey50"); abline(h=0, lty=2, col="grey50")
              if(length(uniquevals) > 0) if(length(uniquevals) >= i) if(!is.null(uniquevals[[i]])) legend("topleft", legend = uniquevals[[i]], col = 2:(length(uniquevals[[i]]) + 1), pch = 19)
	      arrows(x0=0, x1=lx, y0=0, y1=ly, length=0.1, lwd=1)
              text(lx, ly, labels=rownames(lx), pos = c(3,1,3,1), offset=0.3, cex=1)
              plot(groupclocks$weighted.pca.clock.space[[1]][[1]], pch = 19, col = varstoplot[,i], main = paste(colnames(varstoplot)[i], "\n(residual rates)"), xlab = "PC 1", ylab = "PC 2", xlim = c(min(c(range(lxW), range(groupclocks$weighted.pca.clock.space[[1]][[1]][,1]))), max(c(range(lxW), range(groupclocks$weighted.pca.clock.space[[1]][[1]][,1])))), ylim = c(min(c(range(lyW), range(groupclocks$weighted.pca.clock.space[[1]][[1]][,2]))), max(c(range(lyW), range(groupclocks$weighted.pca.clock.space[[1]][[1]][,2])))))
	      abline(v=0, lty=2, col="grey50"); abline(h=0, lty=2, col="grey50")
	      arrows(x0=0, x1=lxW, y0=0, y1=lyW, length=0.1, lwd=1)
    	      text(lxW, lyW, labels=rownames(lxW), pos = c(3,1,3,1), offset=0.1, cex=1)
        }
	
	dev.off()
	
	# Plot species trees with branches showing loadings in PC1 and PC2
	if(is.rooted(sptr)) sptr <- unroot(sptr)
	pdf(paste0(pdf.file, ".branchLoadings.pdf"), useDingbats = F, width = 7, height = Ntip(sptr) * 0.5)
	plot(sptr, main = paste0("Colours by nomalized PC 1 loadings (", round(groupclocks$pca.clock.space[[1]][[3]]$importance[2]*100, 1), "%)\nrelative rates (unrooted tree)"), edge.width = 5, edge.color = topo.colors(10)[as.numeric(cut(scale(groupclocks$pca.clock.space[[1]][[2]][,1]), breaks = 10))])
        edgelabels(frame = "circle", bg = "white")
        plot(sptr, main = paste0("Colours by normalized PC 2 loadings (", round(groupclocks$pca.clock.space[[1]][[3]]$importance[5]*100, 1), "%)\nrelative rates (unrooted tree)"), edge.width = 5, edge.color = topo.colors(10)[as.numeric(cut(scale(groupclocks$pca.clock.space[[1]][[2]][,2]), breaks = 10))])
        edgelabels(frame = "circle", bg = "white")
        sptrPC1 <- sptrPC2 <- sptr
        sptrPC1$edge.length <- groupclocks$pca.clock.space[[1]][[2]][,1]
        sptrPC2$edge.length <- groupclocks$pca.clock.space[[1]][[2]][,2]

	plot(sptr, main = paste0("Colours by nomalized PC 1 loadings (", round(groupclocks$weighted.pca.clock.space[[1]][[3]]$importance[2]*100, 1), "%)\nresidual rates (unrooted tree)"), edge.width = 5, edge.color = topo.colors(10)[as.numeric(cut(scale(groupclocks$weighted.pca.clock.space[[1]][[2]][,1]), breaks = 10))])
	edgelabels(frame = "circle", bg = "white")
	plot(sptr, main = paste0("Colours by normalized PC 2 loadings (", round(groupclocks$weighted.pca.clock.space[[1]][[3]]$importance[5]*100, 1), "%)\nresidual rates (unrooted tree)"), edge.width = 5, edge.color = topo.colors(10)[as.numeric(cut(scale(groupclocks$weighted.pca.clock.space[[1]][[2]][,2]), breaks = 10))])
	edgelabels(frame = "circle", bg = "white")
	sptrPC1resid <- sptrPC2resid <- sptr
	sptrPC1resid$edge.length <- groupclocks$weighted.pca.clock.space[[1]][[2]][,1]
	sptrPC2resid$edge.length <- groupclocks$weighted.pca.clock.space[[1]][[2]][,2]
	dev.off()
	}
	
	# Return table with colours
	if(pca) reslist <- c(groupclocks, list(species.tree.branch.loadings.pc1 = sptrPC1, species.tree.branch.loadings.pc2 = sptrPC2, species.tree.weighted.branch.loadings.pc1 = sptrPC1resid, species.tree.weighted.branch.loadings.pc2 = sptrPC2resid))
	reslist <- c(reslist, list(variable.colours = varstoplot))
	
	cat("Output includes:", fill = T)
        for(i in 1:length(reslist)) cat(paste0(i, ". ", names(reslist)[i]), fill = T)
	return(reslist)
}