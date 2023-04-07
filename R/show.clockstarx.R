show.clockstarx <- function(x) {
	cat("This ClockstaRX object includes:\n")
	for(i in 1:length(x)) cat(paste0(i, ". ", names(x)[i]), fill = T)
	if(any(grepl("weighted", names(x)))) cat("Elements with the name 'weighted' involve 'residual' rates.\n")
	if(any(grepl("pca", names(x)))){
		if(x$pca.clock.space$pPCs[1] < 0.01) cat("Significant PCs of relative (raw) rates include:\n") else cat("No PCs significantly describe variation in relative (raw) rates.\n")
		for(i in 1:length(x$pca.clock.space$pPCs)){
			if(x$pca.clock.space$pPCs[i] > 0.01) break
			cat(paste0("PC", i, " with ", length(which(x$pca.clock.space$pPCs[,i] < 0.01)), " significant branches.\n"))
		}
		if(x$weighted.pca.clock.space$pPCs[1] < 0.01) cat("Significant PCs of relative (raw) rates include:\n") else cat("No PCs significantly describe variation in residual (weighted) rates.\n")
		for(i in 1:length(x$weighted.pca.clock.space$pPCs)){
                        if(x$weighted.pca.clock.space$pPCs[i] > 0.01) break
                        cat(paste0("PC", i, " with ", length(which(x$weighted.pca.clock.space$pPCs[,i] < 0.01)), " significant branches.\n"))
                }
		cat("The full results of tests of clocks are found in the elements pPhi, pPsi, pPCs, and pIL in the pca.clock.space objects.\n")
	}
}
