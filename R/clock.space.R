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
		    plot(pcadata.raw[[1]]$points, pch = 19, xlab = "Dim 1", ylab = "Dim 2", main = "Raw rates (PCA)")
                    plot(pcadata.weighted[[1]]$points, pch = 19, xlab = "Dim 1", ylab = "Dim 2", main = "Residual rates (PCA)")
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