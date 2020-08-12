clock.space <- function(ratesmat, sptr, pca = T, mds = F, sp.time.tree = T, log.branches = F, pca.permutations = 100, mean.scaling.brlen = 0.05, ncore = 1, make.plots = F, sammon.correction = F, verbose = T){

	# If the package mice is to be used for imputation (possibly to add unwanted signal) add the following two arguments: N.imputations = 1, prop.sample.for.imputation = 0.1,
	
	if(sp.time.tree){
		rootlineages <- lapply(Descendants(sptr, getMRCA(sptr, sptr$tip.label), type = "children"), function(x) sptr$tip.label[Descendants(sptr, x, type = "tips")[[1]]])
		if(is.rooted(sptr)) sptr <- unroot(sptr)
		rootbr <- which(apply(sptr$edge, 1, function(x) all(c(mrca.phylo(sptr, rootlineages[[1]]), mrca.phylo(sptr, rootlineages[[2]])) %in% x)))
	} else if(is.rooted(sptr)){
	        sptr <- unroot(sptr)
	}
	
	# Impute data (uses mean for avoiding distortions). Branches with all NAs are turned into the global mean (no contribution to space); branches with a single value also have all NAs replaced with global mean, while branches with more than a single value have NAs replaced by the branch mean.
	raw.rates.matrix <- apply(ratesmat$raw.rates.matrix, 2, function(x){
		 if(log.branches){
			if(any(x == 0, na.rm = T)) x[x == 0] <- 1e-6
			x <- log(x)
		 }
		 if(sum(is.na(x)) < (length(x) - 2) & any(x < 0, na.rm = T)) x <- x + abs(min(x, na.rm = T))
		 return(x)
	})
	
	meanlen <- mean(as.numeric(raw.rates.matrix), na.rm = T)
	impudats <- apply(raw.rates.matrix, 2, function(x){
		 nas <- is.na(x)
		 if(any(nas)){
			if(all(nas)) return(rep(meanlen, length(x)))
			if(sum(nas) == (length(x) - 1)) x[nas] <- rep(meanlen, (length(x) - 1)) else x[nas] <- mean(x, na.rm = T)
			return(x)
		 } else {
		   	return(x)
		 }
	})
	
	if(sp.time.tree) impudats[,rootbr] <- meanlen
	
	# Create weighted matrices (excludes lineage effects)
	impudats.weighted <- weigh.clocks(impudats)
	
	# Might need to change to make only 1 or 2 matrices, as opposed to the whopping 3 matrices (1 raw rates matrix and 2 bsd matrices).
	
	if(mds){
		cat("Making pairwise comparisons of clocks across loci", fill = T)
		rawmatfun <- function(y){
	       	       trs <- apply(y, 1, function(x){
	       	       trtemp <- sptr
		       trtemp$edge.length <- x
	       	       return(trtemp)
	               })
	               class(trs) <- "multiPhylo"
	       	       trtempKF <- KF.dist(trs)
	       	       return(trtempKF)
		}
		rawmat <- rawmatfun(impudats)
		bsdmat <- bsd.matrix(impudats, sptr, mean.brlen = mean.scaling.brlen, ncore = ncore)
		bsdmats.weighted <- bsd.matrix(impudats.weighted, sptr, mean.brlen = mean.scaling.brlen, ncore = ncore)
		cat("Finished making pairwise comparisons", fill = T)
	
		# Dimensionality reduction using MDS (use the function sammon for then sammon's non-linear mapping)
		require(MASS)
		mdsdata.bsd <- if(sammon.correction) sammon(as.dist(bsdmat)) else cmdscale(as.dist(bsdmat), eig = T)
		mdsdata.weighted <- if(sammon.correction) sammon(as.dist(bsdmats.weighted)) else cmdscale(as.dist(bsdmats.weighted), eig = T)
	}
	if(pca){
		pcadata.raw <- pca.perm.test(impudats, perm = pca.permutations)
        	pcadata.weighted <- pca.perm.test(impudats.weighted, perm = pca.permutations)
	}
	
	if(make.plots){
		if(mds){
		    dev.new(width=8, height=4, unit="in")
                    par(mfrow = c(1,2))
		    plot(mdsdata.bsd$points, pch = 19, xlab = "Dim 1", ylab = "Dim 2", main = "Relative rates (MDS)")
		    plot(mdsdata.weighted$points, pch = 19, xlab = "Dim 1", ylab = "Dim 2", main = "Residual rates (MDS)")
		}
		if(pca){
		    dev.new(width=16, height=4, unit="in")
                    par(mfrow = c(1,4))
	topload <- head(order(apply(pcadata.raw[[1]]$rotation[,1:2]^2, 1, function(x) sqrt(sum(x))), decreasing = T), 10)
        lx <- cbind(pcadata.raw[[1]]$rotation[topload,1])
        ly <- cbind(pcadata.raw[[1]]$rotation[topload,2])
        rownames(lx) <- rownames(ly) <- gsub("V", "br ", rownames(lx))
        if(any(abs(c(lx, ly)) < 1e-4)){
                nonzeroloads <- which(abs(lx) > 1e-4 & abs(ly) > 1e-4)
                lx <- cbind(lx[nonzeroloads,])
                ly <- cbind(ly[nonzeroloads,])
        }
        plot(pcadata.raw[[1]]$x[,1:2], pch = 19, xlab=paste0("PC 1 (", round(summary(pcadata.raw[[1]])$importance[2]*100, 1), "%)"), ylab=paste0("PC 2 (", round(summary(pcadata.raw[[1]])$importance[5]*100, 1), "%)"), main = "Raw rates (PCA)", xlim = c(min(c(range(lx), range(pcadata.raw[[1]]$x[,1]))), max(c(range(lx), range(pcadata.raw[[1]]$x[,1])))), ylim = c(min(c(range(ly), range(pcadata.raw[[1]]$x[,2]))), max(c(range(ly), range(pcadata.raw[[1]]$x[,2])))))
        abline(v=0, lty=2, col="grey50"); abline(h=0, lty=2, col="grey50")
        arrows(x0=0, x1=lx, y0=0, y1=ly, length=0.1, lwd=1)
        text(lx, ly, labels=rownames(lx), pos = c(3,1,3,1), offset=0.3, cex=1)

	nulldat <- rbind(cbind(as.vector(pcadata.raw[[2]]), rep(1:length(pcadata.raw[[1]]$sdev), each = pca.permutations)), cbind(summary(pcadata.raw[[1]])$importance[2,], 1:length(pcadata.raw[[1]]$sdev)))
                plot(nulldat[,c(2, 1)], pch = c(rep(20, nrow(nulldat) - length(pcadata.raw[[1]]$sdev)), rep(17, length(pcadata.raw[[1]]$sdev))), col = c(rep("darkgrey", nrow(nulldat) - length(pcadata.raw[[1]]$sdev)), rep("red", length(pcadata.raw[[1]]$sdev))), main = "Residual rates PCA scree", ylab = "Proportion variance explained", xlab = "Principal component")
                lines(tail(nulldat[,c(2, 1)], length(pcadata.raw[[1]]$sdev)), col = 2)
                legend("topright", pch = c(17, 20), col = c("red", "darkgrey"), legend = c("Empirical", "Permuted"))
		    
	topload <- head(order(apply(pcadata.weighted[[1]]$rotation[,1:2]^2, 1, function(x) sqrt(sum(x))), decreasing = T), 10)
        lxW <- cbind(pcadata.weighted[[1]]$rotation[topload,1])
    lyW <- cbind(pcadata.weighted[[1]]$rotation[topload,2])
    rownames(lxW) <- rownames(lyW) <- gsub("V", "br ", rownames(lxW))
    if(any(abs(c(lxW, lyW)) < 1e-4)){
                nonzeroloads <- which(abs(lxW) > 1e-4 & abs(lyW) > 1e-4)
                lxW <- cbind(lxW[nonzeroloads,])
                lyW <- cbind(lyW[nonzeroloads,])
    }
        	plot(pcadata.weighted[[1]]$x[,1:2], pch = 19, xlab=paste0("PC 1 (", round(summary(pcadata.weighted[[1]])$importance[2]*100, 1), "%)"), ylab=paste0("PC 2 (", round(summary(pcadata.weighted[[1]])$importance[5]*100, 1), "%)"), main = "Residual rates (PCA)", xlim = c(min(c(range(lxW), range(pcadata.weighted[[1]]$x[,1]))), max(c(range(lxW), range(pcadata.weighted[[1]]$x[,1])))), ylim = c(min(c(range(lyW), range(pcadata.weighted[[1]]$x[,2]))), max(c(range(lyW), range(pcadata.weighted[[1]]$x[,2])))))
        	abline(v=0, lty=2, col="grey50"); abline(h=0, lty=2, col="grey50")
    		arrows(x0=0, x1=lxW, y0=0, y1=lyW, length=0.1, lwd=1)
    		text(lxW, lyW, labels=rownames(lxW), pos = c(3,1,3,1), offset=0.1, cex=1)
		
		nulldat <- rbind(cbind(as.vector(pcadata.weighted[[2]]), rep(1:length(pcadata.weighted[[1]]$sdev), each = pca.permutations)), cbind(summary(pcadata.weighted[[1]])$importance[2,], 1:length(pcadata.weighted[[1]]$sdev)))
		plot(nulldat[,c(2, 1)], pch = c(rep(20, nrow(nulldat) - length(pcadata.weighted[[1]]$sdev)), rep(17, length(pcadata.weighted[[1]]$sdev))), col = c(rep("darkgrey", nrow(nulldat) - length(pcadata.weighted[[1]]$sdev)), rep("red", length(pcadata.weighted[[1]]$sdev))), main = "Residual rates PCA scree", ylab = "Proportion variance explained", xlab = "Principal component")
		lines(tail(nulldat[,c(2, 1)], length(pcadata.weighted[[1]]$sdev)), col = 2)
		legend("topright", pch = c(17, 20), col = c("red", "darkgrey"), legend = c("Empirical", "Permuted"))
		
		}
	}
	
	# Return all imputed rates matrices and rates spaces.
	reslist <- c(ratesmat, list(imputed.clocks = impudats, weighted.imputed.clocks = impudats.weighted))
	if(mds) reslist <- c(reslist, list(bsd.pairwise.matrix = bsdmat, weighted.bsd.pairwise.matrix = bsdmats.weighted, bsd.clock.space = mdsdata.bsd, weighted.bsd.clock.space = mdsdata.weighted))
	if(pca) reslist <- c(reslist, list(pca.clock.space = pcadata.raw, weighted.pca.clock.space = pcadata.weighted))
	if(verbose){ cat("Output includes:", fill = T)
        for(i in 1:length(reslist)) cat(paste0(i, ". ", names(reslist)[i]), fill = T) }
	return(reslist)
}