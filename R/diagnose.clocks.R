# Plot basic data in clockspace, e.g., number of taxa, branch support, tree length, topological distance from species tree (plus other user-selected variables). other.data must be a matrix or data.frame where each column is a numeric variable.

diagnose.clocks <- function(groupclocks, loctrs, sptr, other.data = NULL, default.variables = c("ntaxa", "nmissingbr", "brsup", "trlen", "rttdist", "topdist"), pdf.file = "clock.diagnosis.pdf"){
	
	# Extract default variables and create their colours
	
	varstoplot <- data.frame("Raw rates clusters" = as.numeric(groupclocks$raw.clustering[[1]][, groupclocks$raw.best.k[1]]), "Relative rates clusters" <- as.numeric(groupclocks$bsd.clustering[[1]][, groupclocks$bsd.best.k[1]]), "Weighted rates clusters" <- as.numeric(groupclocks$weighted.bsd.clustering[[1]][, groupclocks$weighted.bsd.best.k[1]]))
	
	grPal <- colorRampPalette(c('lightblue', 'blue','darkblue'))
	
	ntaxa <- sapply(loctrs, Ntip)
	if("ntaxa" %in% default.variables) varstoplot["N taxa"] <- grPal(10)[as.numeric(cut(ntaxa, breaks = 10))]
	
	nmissingbr <- sapply(1:length(loctrs), function(x) sum(is.na(groupclocks[[1]][x,])) / Nedge(loctrs[[x]]))
	if("nmissingbr" %in% default.variables) varstoplot["N missing species tree branches"] <- grPal(10)[as.numeric(cut(nmissingbr, breaks = 10))]

	brsup <- sapply(loctrs, function(x) mean(as.numeric(x$node.label), na.rm = T))
	if("brsup" %in% default.variables) varstoplot["Branch support"] <- grPal(10)[as.numeric(cut(brsup, breaks = 10))]

	trlen <- sapply(loctrs, function(x) sum(x$edge.length, na.rm = T) / Ntip(x))
	if("trlen" %in% default.variables) varstoplot["Tree length"] <- grPal(10)[as.numeric(cut(trlen, breaks = 10))]
	
	rtts <- lapply(loctrs, function(x) pathnode(midpoint(x))[[1]])
	rttdist <- sapply(rtts, function(x) sd(x, na.rm = T) / mean(x, na.rm = T))
	if("rttdist" %in% default.variables) varstoplot["Midpoint-root-to-tip CV"] <- grPal(10)[as.numeric(cut(rttdist, breaks = 10))]

	topdist <- sapply(loctrs, function(x) RF.dist(x, drop.tip(sptr, sptr$tip.label[which(!sptr$tip.label %in% x$tip.label)])))
	if("topdist" %in% default.variables) varstoplot["Distance to species tree"] <- grPal(10)[as.numeric(cut(topdist, breaks = 10))]

	# Create colours for other variables
	
	uniquevals <- list()
	
	if(!is.null(other.data)){
		for(i in 1:ncol(other.data)){
		      if(class(other.data[,i]) == "numeric"){
				varstoplot[colnames(other.data)[i]] <- grPal(10)[as.numeric(cut(other.data[,i], breaks = 10))]
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
	
	# Create pdf where each plot is 4x4 in size.
	pdf(pdf.file, useDingbats = F, width = 12, height = (length(varstoplot) - 2) * 4)
	par(mfrow = c(length(varstoplot) - 2, 3))
	
	# First plot the clustering results
	plot(groupclocks$raw.clock.space[[1]]$points, pch = 19, col = as.numeric(groupclocks$raw.clustering[[1]][, groupclocks$raw.best.k[1]]) + 1, xlab = "Dim 1", ylab = "Dim 2", main = "Raw rates clusters")
	legend("topright", legend = paste0("k = ", groupclocks$raw.best.k[1]))
	plot(groupclocks$bsd.clock.space[[1]]$points, pch = 19, col = as.numeric(groupclocks$bsd.clustering[[1]][, groupclocks$bsd.best.k[1]]) + 1, xlab = "Dim 1", ylab = "Dim 2", main = "Relative rates clusters")
	legend("topright", legend = paste0("k = ", groupclocks$bsd.best.k[1]))
	plot(groupclocks$weighted.bsd.clock.space[[1]]$points, pch = 19, col = as.numeric(groupclocks$weighted.bsd.clustering[[1]][, groupclocks$weighted.bsd.best.k[1]]) + 1, xlab = "Dim 1", ylab = "Dim 2", main = "Residual rates clusters")
	legend("topright", legend = paste0("k = ", groupclocks$weighted.bsd.best.k[1]))
	
	# Plot the new variable results
	for(i in 4:length(varstoplot)){
	      plot(groupclocks$raw.clock.space[[1]]$points, pch = 19, col = varstoplot[,i], xlab = "Dim 1", ylab = "Dim 2", main = paste(colnames(varstoplot)[i], "\n(raw rates)"))
	      if(length(uniquevals) > 0) if(length(uniquevals) >= i) if(!is.null(uniquevals[[i]])) legend("topleft", legend = uniquevals[[i]], col = 2:(length(uniquevals[[i]]) + 1), pch = 19)
              plot(groupclocks$bsd.clock.space[[1]]$points, pch = 19, col = varstoplot[,i], xlab = "Dim 1", ylab = "Dim 2", main = paste(colnames(varstoplot)[i], "\n(relative rates)"))
              plot(groupclocks$weighted.bsd.clock.space[[1]]$points, pch = 19, col = varstoplot[,i], xlab = "Dim 1", ylab = "Dim 2", main = paste(colnames(varstoplot)[i], "\n(residual rates)"))
	}

	dev.off()
	
	# Return table with colours
	reslist <- c(groupclocks, list(variable.colours = varstoplot))
	return(reslist)
}