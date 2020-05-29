group.clocks <- function(clockspace, boot.samps = 50, kmax = 10, make.plots = F, verbose = T){

	#Clustering. The criterion for clustering is the following: if there is a significant decline in gap from one cluster cluster to two, then there is a single cluster. Otherwhise, the largest increase in gap leads to the optimum number of clusters (see Duchene et al. 2018, Syst. Biol.).
	
	if("bsd.clock.space" %in% names(clockspace)) mds <- T else mds <- F
	if("pca.clock.space" %in% names(clockspace)) pca <- T else pca <- F
	
	nimput <- 1:length(clockspace[[grep("clock[.]space", names(clockspace))[1]]])
	
	# Identify best number of clusters for each rates space.
	require(cluster)
	
	if(mds){
	
	clusdatbsd <- lapply(nimput, function(x) clusGap(clockspace$bsd.clock.space[[x]]$points, B = boot.samps, FUNcluster = pam, K.max = kmax))
	gapstablebsd <- lapply(nimput, function(x){
                     tab <- clusdatbsd[[x]]$Tab
             	     tab[is.infinite(tab)] <- NA
                     return(tab)
    	})
	npartbsd <- sapply(nimput, function(x) gapCR(gapstablebsd[[x]][,3], gapstablebsd[[x]][,4]))
	
	clusdatbsdW <- lapply(nimput, function(x) clusGap(clockspace$weighted.bsd.clock.space[[x]]$points, B = boot.samps, FUNcluster = pam, K.max = kmax))
        gapstablebsdW <- lapply(nimput, function(x){
                     tab <- clusdatbsdW[[x]]$Tab
                     tab[is.infinite(tab)] <- NA
                     return(tab)
        })
        npartbsdW <- sapply(nimput, function(x) gapCR(gapstablebsdW[[x]][,3], gapstablebsdW[[x]][,4]))
    
	# Collect cluster data at all numbers of clusters, for each rates space.
    	clust.bsd <- lapply(nimput, function(y) do.call("cbind", lapply(1:kmax, function(x) pam(clockspace$bsd.clock.space[[y]]$points, k = x)$clustering)))
    	clust.bsdW <- lapply(nimput, function(y) do.call("cbind", lapply(1:kmax, function(x) pam(clockspace$weighted.bsd.clock.space[[y]]$points, k = x)$clustering)))
    	
	}
    
	if(pca){
	
	clusdatpca <- lapply(nimput, function(x) clusGap(clockspace$pca.clock.space[[x]]$points, B = boot.samps, FUNcluster = pam, K.max = kmax))
        gapstablepca <- lapply(nimput, function(x){
                     tab <- clusdatpca[[x]]$Tab
                     tab[is.infinite(tab)] <- NA
                     return(tab)
        })
        npartpca <- sapply(nimput, function(x) gapCR(gapstablepca[[x]][,3], gapstablepca[[x]][,4]))

        clusdatpcaW <- lapply(nimput, function(x) clusGap(clockspace$weighted.pca.clock.space[[x]]$points, B = boot.samps, FUNcluster = pam, K.max = kmax))
        gapstablepcaW <- lapply(nimput, function(x){
                     tab <- clusdatpcaW[[x]]$Tab
                     tab[is.infinite(tab)] <- NA
                     return(tab)
        })
        npartpcaW <- sapply(nimput, function(x) gapCR(gapstablepcaW[[x]][,3], gapstablepcaW[[x]][,4]))
	
	clust.pca <- lapply(nimput, function(y) do.call("cbind", lapply(1:kmax, function(x) pam(clockspace$pca.clock.space[[y]]$points, k = x)$clustering)))
    	clust.pcaW <- lapply(nimput, function(y) do.call("cbind", lapply(1:kmax, function(x) pam(clockspace$weighted.pca.clock.space[[y]]$points, k = x)$clustering)))

    	}
    
	#rownames(clust.raw) <- rownames(clust.bsd) <- rownames(clust.bsdW) <- rownames(clockspace$imputed.clocks[[1]])
    	#colnames(clust.raw) <- colnames(clust.bsd) <- colnames(clust.bsdW) <- paste0("k=", 1:kmax)
	
	if(make.plots){
		if(mds & pca){
		        dev.new(width=16, height=8, unit="in")
			par(mfrow = c(2, 4))
		} else {
		        dev.new(width=16, height=4, unit="in")
			par(mfrow = c(1, 4))
		}
		if(mds){
			plot(clockspace$bsd.clock.space[[1]]$points, pch = 19, col = as.numeric(clust.bsd[[1]][,npartbsd[1]]) + 1, xlab = "Dim 1", ylab = "Dim 2", main = "Relatve rates (MDS)")
			legend("topright", legend = paste0("k = ", npartbsd[1]))
			plot(clusdatbsd[[1]], main = c("Support for k clusters\nof relative rates (MDS)"))
			plot(clockspace$weighted.bsd.clock.space[[1]]$points, pch = 19, col = as.numeric(clust.bsdW[[1]][,npartbsdW[1]]) + 1, xlab = "Dim 1", ylab = "Dim 2", main = "Residual rates (MDS)")
			legend("topright", legend = paste0("k = ", npartbsdW[1]))
			plot(clusdatbsdW[[1]], main = c("Support for k clusters\nof residual rates (MDS)"))
		}
		if(pca){
			plot(clockspace$pca.clock.space[[1]]$points, pch = 19, col = as.numeric(clust.pca[[1]][,npartpca[1]]) + 1, xlab = "Dim 1", ylab = "Dim 2", main = "Relatve rates (PCA)")
                     	legend("topright", legend = paste0("k = ", npartpca[1]))
		     	plot(clusdatpca[[1]], main = c("Support for k clusters\nof relative rates (PCA)"))
			plot(clockspace$weighted.pca.clock.space[[1]]$points, pch = 19, col = as.numeric(clust.pcaW[[1]][,npartpcaW[1]]) + 1, xlab = "Dim 1", ylab = "Dim 2", main = "Residual rates (PCA)")
                 	legend("topright", legend = paste0("k = ", npartpcaW[1]))
			plot(clusdatpcaW[[1]], main = c("Support for k clusters\nof residual rates (PCA)"))
		}
		
		# If using the package mice for imputation, then add plot summarising the clustering findings of different imputations.

	}
	
	# Return clustering data and best numbers of clusters for each rates matrix.
	reslist <- clockspace
	if(mds) reslist <- c(reslist, list(bsd.cluster.support = clusdatbsd, weighted.bsd.cluster.support = clusdatbsdW, bsd.clustering = clust.bsd, weighted.bsd.clustering = clust.bsdW, bsd.best.k = npartbsd, weighted.bsd.best.k = npartbsdW))
	if(pca) reslist <- c(reslist, list(pca.cluster.support = clusdatpca, weighted.pca.cluster.support = clusdatpcaW, pca.clustering = clust.pca, weighted.pca.clustering = clust.pcaW, pca.best.k = npartpca, weighted.pca.best.k = npartpcaW))
	if(verbose){ cat("Output includes:", fill = T)
        for(i in 1:length(reslist)) cat(paste0(i, ". ", names(reslist)[i]), fill = T) }
	return(reslist)
}