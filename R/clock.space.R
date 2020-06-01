clock.space <- function(ratesmat, sptr, pca = T, mds = F, log.branches = T, mean.scaling.brlen = 0.05, ncore = 1, make.plots = F, sammon.correction = F, verbose = T){

	# If the package mice is to be used for imputation (possibly to add unwanted signal) add the following two arguments: N.imputations = 1, prop.sample.for.imputation = 0.1,
	
	if(is.rooted(sptr)) sptr <- unroot(sptr)
	# Impute data (uses mean for avoiding distortions). Branches with all NAs are turned into the global mean (no contribution to space); branches with a single value also have all NAs replaced with global mean, while branches with more than a single value have NAs replaced by the branch mean.
	raw.rates.matrix <- apply(ratesmat$raw.rates.matrix, 2, function(x){
		 if(log.branches){
			if(any(x == 0, na.rm = T)) x[x == 0] <- 1e-6
			x <- log(x)
		 }
		 if(sum(is.na(x)) < (length(x) - 2) & any(x < 0, na.rm = T)) x <- x + abs(min(x, na.rm = T))
		 return(x)
	})
	
	medianlen <- mean(as.numeric(raw.rates.matrix), na.rm = T)
	impudats <- list(apply(raw.rates.matrix, 2, function(x){
		 nas <- is.na(x)
		 if(any(nas)){
			if(all(nas)) return(rep(medianlen, length(x)))
			if(sum(nas) == (length(x) - 1)) x[nas] <- rep(medianlen, (length(x) - 1)) else x[nas] <- mean(x, na.rm = T)
			return(x)
		 } else {
		   	return(x)
		 }
	}))
	
	# Create weighted matrices (excludes lineage effects)
	impudats.weighted <- lapply(impudats, weigh.clocks)
	
	# Might need to change to make only 1 or 2 matrices, as opposed to the whopping 3 matrices (1 raw rates matrix and 2 bsd matrices).
	
	if(mds){
		cat("Making pairwise comparisons of clocks across loci", fill = T)
		rawmat <- lapply(impudats, function(y){
	       	       trs <- apply(y, 1, function(x){
	       	       trtemp <- sptr
		       trtemp$edge.length <- x
	       	       return(trtemp)
	               })
	               class(trs) <- "multiPhylo"
	       	       trtempKF <- KF.dist(trs)
	       	       return(trtempKF)
		})
		bsdmat <- lapply(impudats, function(x) bsd.matrix(x, sptr, mean.brlen = mean.scaling.brlen, ncore = ncore))
		bsdmats.weighted <- lapply(impudats.weighted, function(x) bsd.matrix(x, sptr, mean.brlen = mean.scaling.brlen, ncore = ncore))
		cat("Finished making pairwise comparisons", fill = T)
	
		# Dimensionality reduction using MDS (use the function sammon for then sammon's non-linear mapping)
		require(MASS)
		mdsdata.bsd <- lapply(bsdmat, function(x) if(sammon.correction) sammon(as.dist(x)) else cmdscale(as.dist(x), eig = T))
		mdsdata.weighted <- lapply(bsdmats.weighted, function(x) if(sammon.correction) sammon(as.dist(x)) else cmdscale(as.dist(x), eig = T))
	}
	if(pca){
		pcadata.raw <- lapply(impudats, function(x){
                    pcs1n2 <- prcomp(x, center = T)
                    res <- list(points = pcs1n2$x[,1:2], rotation = pcs1n2$rotation, summary = summary(pcs1n2))
                    return(res)
        	})
        	pcadata.weighted <- lapply(impudats.weighted, function(x){
                    pcs1n2 <- prcomp(x, center = T)
                    res    <- list(points = pcs1n2$x[,1:2], rotation = pcs1n2$rotation, summary = summary(pcs1n2))
                    return(res)
        	})
	}
	
	if(make.plots){
		if(mds & pca){
		       dev.new(width=8, height=8, unit="in")
		       par(mfrow = c(2,2))
		} else {
		       dev.new(width=8, height=4, unit="in")
		       par(mfrow = c(1,2))
		}
		if(mds){
		    plot(mdsdata.bsd[[1]]$points, pch = 19, xlab = "Dim 1", ylab = "Dim 2", main = "Relative rates (MDS)")
		    plot(mdsdata.weighted[[1]]$points, pch = 19, xlab = "Dim 1", ylab = "Dim 2", main = "Residual rates (MDS)")
		}
		if(pca){
		    
		    topload <- head(order(pcadata.raw[[1]][[2]][,1]), 5)
        ordpc2 <- order(pcadata.raw[[1]][[2]][,2])
        topload <- c(topload, head(ordpc2[-which(ordpc2 %in% topload)], 5))
        lx <- cbind(pcadata.raw[[1]][[2]][topload,1])
        ly <- cbind(pcadata.raw[[1]][[2]][topload,2])
        rownames(lx) <- rownames(ly) <- gsub("V", "br ", rownames(lx))
        if(any(abs(c(lx, ly)) < 1e-4)){
                nonzeroloads <- which(abs(lx) > 1e-4 & abs(ly) > 1e-4)
                lx <- cbind(lx[nonzeroloads,])
                ly <- cbind(ly[nonzeroloads,])
        }
        plot(pcadata.raw[[1]]$points, pch = 19, xlab=paste0("PC 1 (", round(pcadata.raw[[1]][[3]]$importance[2]*100, 1), "%)"), ylab=paste0("PC 2 (", round(pcadata.raw[[1]][[3]]$importance[5]*100, 1), "%)"), main = "Raw rates (PCA)", xlim = c(min(c(range(lx), range(pcadata.raw[[1]][[1]][,1]))), max(c(range(lx), range(pcadata.raw[[1]][[1]][,1])))), ylim = c(min(c(range(ly), range(pcadata.raw[[1]][[1]][,2]))), max(c(range(ly), range(pcadata.raw[[1]][[1]][,2])))))
        abline(v=0, lty=2, col="grey50"); abline(h=0, lty=2, col="grey50")
        arrows(x0=0, x1=lx, y0=0, y1=ly, length=0.1, lwd=1)
        text(lx, ly, labels=rownames(lx), pos = c(3,1,3,1), offset=0.3, cex=1)
		    
		    topload <- head(order(pcadata.weighted[[1]][[2]][,1]), 5)
        ordpc2 <- order(pcadata.weighted[[1]][[2]][,2])
        topload <- c(topload, head(ordpc2[-which(ordpc2 %in% topload)], 5))
        lxW <- cbind(pcadata.weighted[[1]][[2]][topload,1])
    lyW <- cbind(pcadata.weighted[[1]][[2]][topload,2])
    rownames(lxW) <- rownames(lyW) <- gsub("V", "br ", rownames(lxW))
    if(any(abs(c(lxW, lyW)) < 1e-4)){
                nonzeroloads <- which(abs(lxW) > 1e-4 & abs(lyW) > 1e-4)
                lxW <- cbind(lxW[nonzeroloads,])
                lyW <- cbind(lyW[nonzeroloads,])
    }
        	plot(pcadata.weighted[[1]]$points, pch = 19, xlab=paste0("PC 1 (", round(pcadata.weighted[[1]][[3]]$importance[2]*100, 1), "%)"), ylab=paste0("PC 2 (", round(pcadata.weighted[[1]][[3]]$importance[5]*100, 1), "%)"), main = "Residual rates (PCA)", xlim = c(min(c(range(lxW), range(pcadata.weighted[[1]][[1]][,1]))), max(c(range(lxW), range(pcadata.weighted[[1]][[1]][,1])))), ylim = c(min(c(range(lyW), range(pcadata.weighted[[1]][[1]][,2]))), max(c(range(lyW), range(pcadata.weighted[[1]][[1]][,2])))))
        	abline(v=0, lty=2, col="grey50"); abline(h=0, lty=2, col="grey50")
    		arrows(x0=0, x1=lxW, y0=0, y1=lyW, length=0.1, lwd=1)
    		text(lxW, lyW, labels=rownames(lxW), pos = c(3,1,3,1), offset=0.1, cex=1)
		}
	}
	
	# Return all imputed rates matrices and rates spaces.
	reslist <- c(ratesmat, list(imputed.clocks = impudats))
	if(mds) reslist <- c(reslist, list(bsd.pairwise.matrix = bsdmat, weighted.bsd.pairwise.matrix = bsdmats.weighted, bsd.clock.space = mdsdata.bsd, weighted.bsd.clock.space = mdsdata.weighted))
	if(pca) reslist <- c(reslist, list(pca.clock.space = pcadata.raw, weighted.pca.clock.space = pcadata.weighted))
	if(verbose){ cat("Output includes:", fill = T)
        for(i in 1:length(reslist)) cat(paste0(i, ". ", names(reslist)[i]), fill = T) }
	return(reslist)
}