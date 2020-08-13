group.clocks <- function(clockspace, boot.samps = 50, kmax = 10, make.plots = F, verbose = T){

	#Clustering. The criterion for clustering is the following: if there is a significant decline in gap from one cluster cluster to two, then there is a single cluster. Otherwhise, the largest increase in gap leads to the optimum number of clusters (see Duchene et al. 2018, Syst. Biol.).
	
	if("bsd.clock.space" %in% names(clockspace)) mds <- T else mds <- F
	if("pca.clock.space" %in% names(clockspace)) pca <- T else pca <- F
	
	# Identify best number of clusters for each rates space.
	require(cluster)
	
	if(mds){
	
	clusdatbsd <- clusGap(clockspace$bsd.clock.space$points, B = boot.samps, FUNcluster = pam, K.max = kmax)
	gapstablebsd <- clusdatbsd$Tab
        gapstablebsd[is.infinite(gapstablebsd)] <- NA
        
	npartbsd <- gapCR(gapstablebsd[,3], gapstablebsd[,4])
	
	clusdatbsdW <- clusGap(clockspace$weighted.bsd.clock.space$points, B = boot.samps, FUNcluster = pam, K.max = kmax)
        gapstablebsdW <- clusdatbsdW$Tab
        gapstablebsdW[is.infinite(gapstablebsdW)] <- NA
        
        npartbsdW <- gapCR(gapstablebsdW[,3], gapstablebsdW[,4])
    
	# Collect cluster data at all numbers of clusters, for each rates space.
    	clust.bsd <- do.call("cbind", lapply(1:kmax, function(x) pam(clockspace$bsd.clock.space$points, k = x)$clustering))
    	clust.bsdW <- do.call("cbind", lapply(1:kmax, function(x) pam(clockspace$weighted.bsd.clock.space$points, k = x)$clustering))
    	
	}
    
	if(pca){
	
	clusdatpca <- clusGap(clockspace$pca.clock.space[[1]]$x[,1:2], B = boot.samps, FUNcluster = pam, K.max = kmax)
        gapstablepca <- clusdatpca$Tab
        gapstablepca[is.infinite(gapstablepca)] <- NA
        
        npartpca <- gapCR(gapstablepca[,3], gapstablepca[,4])

        clusdatpcaW <- clusGap(clockspace$weighted.pca.clock.space[[1]]$x[,1:2], B = boot.samps, FUNcluster = pam, K.max = kmax)
        gapstablepcaW <- clusdatpcaW$Tab
        gapstablepcaW[is.infinite(gapstablepcaW)] <- NA
        
        npartpcaW <- gapCR(gapstablepcaW[,3], gapstablepcaW[,4])
	
	clust.pca <- do.call("cbind", lapply(1:kmax, function(x) pam(clockspace$pca.clock.space[[1]]$x[,1:2], k = x)$clustering))
    	clust.pcaW <- do.call("cbind", lapply(1:kmax, function(x) pam(clockspace$weighted.pca.clock.space[[1]]$x[,1:2], k = x)$clustering))

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
			plot(clockspace$bsd.clock.space$points, pch = 19, col = as.numeric(clust.bsd[,npartbsd[1]]) + 1, xlab = "Dim 1", ylab = "Dim 2", main = "Relatve rates (MDS)")
			legend("topright", legend = paste0("k = ", npartbsd))
			plot(clusdatbsd, main = c("Support for k clusters\nof relative rates (MDS)"))
			plot(clockspace$weighted.bsd.clock.space$points, pch = 19, col = as.numeric(clust.bsdW[,npartbsdW[1]]) + 1, xlab = "Dim 1", ylab = "Dim 2", main = "Residual rates (MDS)")
			legend("topright", legend = paste0("k = ", npartbsdW))
			plot(clusdatbsdW, main = c("Support for k clusters\nof residual rates (MDS)"))
		}
		if(pca){
			
	lam <- clockspace$pca.clock.space[[1]]$sdev * sqrt(nrow(clockspace$pca.clock.space[[1]]$x))
	topload <- head(order(apply(clockspace$pca.clock.space[[1]]$rotation[,1:2]^2, 1, function(x) sqrt(sum(x))), decreasing = T), 10)
        lx <- cbind(t(t(clockspace$pca.clock.space[[1]]$rotation[topload,1]) * lam[1]))
        ly <- cbind(t(t(clockspace$pca.clock.space[[1]]$rotation[topload,2]) * lam[2]))
        rownames(lx) <- rownames(ly) <- gsub("V", "br ", rownames(lx))
        if(any(abs(c(lx, ly)) < 1e-4)){
                nonzeroloads <- which(abs(lx) > 1e-4 & abs(ly) > 1e-4)
                lx <- cbind(lx[nonzeroloads,])
                ly <- cbind(ly[nonzeroloads,])
        }
			plot(clockspace$pca.clock.space[[1]]$x[,1:2], pch = 19, col = as.numeric(clust.pca[,npartpca[1]]) + 1, xlab=paste0("PC 1 (", round(summary(clockspace$pca.clock.space[[1]])$importance[2]*100, 1), "%)"), ylab=paste0("PC 2 (", round(summary(clockspace$pca.clock.space[[1]])$importance[5]*100, 1), "%)"), main = "Relatve rates (PCA)", xlim = c(min(c(range(lx), range(clockspace$pca.clock.space[[1]]$x[,1]))), max(c(range(lx), range(clockspace$pca.clock.space[[1]]$x[,1])))), ylim = c(min(c(range(ly), range(clockspace$pca.clock.space[[1]]$x[,2]))), max(c(range(ly), range(clockspace$pca.clock.space[[1]]$x[,2])))))
            legend("topright", legend = paste0("k = ", npartpca))
            abline(v=0, lty=2, col="grey50"); abline(h=0, lty=2, col="grey50")
        arrows(x0=0, x1=lx, y0=0, y1=ly, length=0.1, lwd=1)
        text(lx, ly, labels=rownames(lx), pos = c(3,1,3,1), offset=0.3, cex=1)
		    
		    plot(clusdatpca, main = c("Support for k clusters\nof relative rates (PCA)"))
			
	lamW <- clockspace$weighted.pca.clock.space[[1]]$sdev * sqrt(nrow(clockspace$weighted.pca.clock.space[[1]]$x))
	topload <- head(order(apply(clockspace$weighted.pca.clock.space[[1]]$rotation[,1:2]^2, 1, function(x) sqrt(sum(x))), decreasing = T), 10)
        lxW <- cbind(t(t(clockspace$weighted.pca.clock.space[[1]]$rotation[topload,1]) * lamW[1]))
    lyW <- cbind(t(t(clockspace$weighted.pca.clock.space[[1]]$rotation[topload,2]) * lamW[2]))
    rownames(lxW) <- rownames(lyW) <- gsub("V", "br ", rownames(lxW))
    if(any(abs(c(lxW, lyW)) < 1e-4)){
                nonzeroloads <- which(abs(lxW) > 1e-4 & abs(lyW) > 1e-4)
                lxW <- cbind(lxW[nonzeroloads,])
                lyW <- cbind(lyW[nonzeroloads,])
    }
			plot(clockspace$weighted.pca.clock.space[[1]]$x[,1:2], pch = 19, col = as.numeric(clust.pcaW[,npartpcaW[1]]) + 1, xlab=paste0("PC 1 (", round(summary(clockspace$weighted.pca.clock.space[[1]])$importance[2]*100, 1), "%)"), ylab=paste0("PC 2 (", round(summary(clockspace$weighted.pca.clock.space[[1]])$importance[5]*100, 1), "%)"), main = "Residual rates (PCA)", xlim = c(min(c(range(lxW), range(clockspace$weighted.pca.clock.space[[1]]$x[,1]))), max(c(range(lxW), range(clockspace$weighted.pca.clock.space[[1]]$x[,1])))), ylim = c(min(c(range(lyW), range(clockspace$weighted.pca.clock.space[[1]]$x[,2]))), max(c(range(lyW), range(clockspace$weighted.pca.clock.space[[1]]$x[,2])))))
            legend("topright", legend = paste0("k = ", npartpcaW))
            abline(v=0, lty=2, col="grey50"); abline(h=0, lty=2, col="grey50")
    arrows(x0=0, x1=lxW, y0=0, y1=lyW, length=0.1, lwd=1)
    text(lxW, lyW, labels=rownames(lxW), pos = c(3,1,3,1), offset=0.1, cex=1)
			
			plot(clusdatpcaW, main = c("Support for k clusters\nof residual rates (PCA)"))
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