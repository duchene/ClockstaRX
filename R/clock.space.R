# Impute missing rates (or fill using median branch lengths), and perform dimensionality reduction to plot rate space. Requires packages mice and MASS.

clock.space <- function(ratesmat, sptr, pca = T, mds = F, mean.scaling.brlen = 0.05, ncore = 1, make.plot = F){

	# If the package mice is to be used for imputation (possibly to add unwanted signal) add the following two arguments: N.imputations = 1, prop.sample.for.imputation = 0.1,
	
	if(is.rooted(sptr)) sptr <- unroot(sptr)

	predmat <- matrix(0, nrow = nrow(ratesmat[[1]]), ncol = nrow(ratesmat[[1]]))
	for(i in 1:nrow(predmat)) predmat[sample(c(1:nrow(predmat))[-i], round(nrow(ratesmat[[1]]) * prop.sample.for.imputation)),i] <- 1
	
	# Proper rates imputation (possibly increases distortion of rates space)
	#print("Data imputation is starting")
	#imputation <- suppressWarnings(mice(t(ratesmat[[1]]), m = N.imputations, maxit = 30, method = 'norm', predictorMatrix = predmat, seed = 123))
	#print("Data imputation completed.")
	#impudats <- lapply(1:N.imputations, function(x) t(complete(imputation, x)))
	
	# Alternative to imputation by filling missing data with the median. This might minimize distortion of clock space
	medianlen <- mean(as.numeric(ratesmat[[1]]), na.rm = T)
	impudats <- list(apply(ratesmat[[1]], 2, function(x){
		 nas <- is.na(x)
		 if(any(nas)){
			if(all(nas)) return(rep(medianlen, length(x)))
			x[nas] <- mean(x, na.rm = T)
		 	return(x)
		 } else {
		   	return(x)
		 }
	}))
	
	# Create weighted matrices (excludes lineage effects)
	impudats.weighted <- lapply(impudats, weigh.clocks)
	
	# Might need to change to make only 1 or 2 matrices, as opposed to the whopping 3 matrices (1 raw rates matrix and 2 bsd matrices).
	
	if(mds){
		print("Making pairwise comparisons...")
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
		print("Finished making pairwise comparisons.")
	
		# Dimensionality reduction using MDS and then sammon's non-linear mapping
		require(MASS)
		mdsdata.bsd <- lapply(bsdmat, function(x) sammon(as.matrix(x)))
		mdsdata.weighted <- lapply(bsdmats.weighted, function(x) sammon(as.matrix(x)))
	}
	if(pca){
		pcadata.raw <- lapply(impudats, function(x){
                    pcs1n2 <- prcomp(x, center = T)
                    res <- list(points = pcs1n2$x[,1:2], summary = summary(pcs1n2))
                    return(res)
        	})
        	pcadata.weighted <- lapply(impudats.weighted, function(x){
                    pcs1n2 <- prcomp(x, center = T)
                    res    <- list(points = pcs1n2$x[,1:2], summary = summary(pcs1n2))
                    return(res)
        	})
	}
	
	if(make.plot){
		dev.new()
		if(mds && pca) par(mfrow = c(2,2)) else par(mfrow = c(1,2))
		if(mds){
		    plot(mdsdata.bsd[[1]]$points, pch = 19, xlab = "Dim 1", ylab = "Dim 2", main = "Relative rates (Sammon)")
		    plot(mdsdata.weighted[[1]]$points, pch = 19, xlab = "Dim 1", ylab = "Dim 2", main = "Residual rates (Sammon)")
		}
		if(pca){
		    plot(pcadata.raw[[1]]$points, pch = 19, xlab = "Dim 1", ylab = "Dim 2", main = "Raw rates (PCA)")
                    plot(pcadata.weighted[[1]]$points, pch = 19, xlab = "Dim 1", ylab = "Dim 2", main = "Residual rates (PCA)")
		
		}
	}
	
	# Return all imputed rates matrices and rates spaces.
	reslist <- c(ratesmat, list(imputed.clocks = impudats))
	if(mds) reslist <- c(reslst, list(bsd.pairwise.matrix = bsdmat, weighted.bsd.pairwise.matrix = bsdmats.weighted, bsd.clock.space = mdsdata.bsd, weighted.bsd.clock.space = mdsdata.weighted))
	if(pca) reslist <- c(reslist, list(pca.clock.space = pca.raw, weighted.pca.clock.space = pcadata.weighted))
	return(reslist)
}
