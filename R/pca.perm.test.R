pca.perm.test <- function(tab, perm = 100){
	      
	      tabrows <- nrow(tab)
	      empPCA <- prcomp(tab, center = T)
              perPCA <- lapply(1:perm, function(x) prcomp(apply(tab, 2, function(y) y[sample(tabrows)]), center = T))

	      #First, test the whole PCA analysis
	      pcaPhi <- function(x) sqrt(abs((log(sum(((x$sdev)^2)^2, na.rm = T)) - length(x$sdev)) / (length(x$sdev)*(length(x$sdev)-1))))
	      pcaPsi <- function(x) sum((((x$sdev)^2) - mean((x$sdev)^2, na.rm = T))^2, na.rm = T)
	      
	      empPhi <- pcaPhi(empPCA)
	      empPsi <- pcaPsi(empPCA)
	      perPhi <- sapply(perPCA, pcaPhi)
	      perPsi <- sapply(perPCA, pcaPsi)
	      
	      pPhi <- pnorm((empPhi - mean(perPhi, na.rm = T)) / sd(perPhi, na.rm = T), lower.tail = T)
	      pPsi <- pnorm((empPsi - mean(perPsi, na.rm = T)) / sd(perPsi, na.rm = T), lower.tail = F)
	      
	      #Second, test each PC
	      
	      pcnull <- do.call("rbind", lapply(perPCA, function(y) y$sdev^2))

	      pPC <- sapply(1:length(empPCA$sdev), function(x){
		  pcPval <- pnorm((empPCA$sdev[x]^2 - mean(pcnull[,x], na.rm = T)) / sd(pcnull[,x], na.rm = T), lower.tail = F)
		  return(pcPval)
	      })

	      permimp <- do.call("rbind", lapply(perPCA, function(y) summary(y)$importance[2,]))

	      #Third test each variable in each PC

	      pcaIL <- function(x) sapply(1:length(x$sdev), function(y) x$rotation[,y]^2 * ((x$sdev[y]^2)^2))
	      
	      empIL <- pcaIL(empPCA)
	      perIL <- lapply(perPCA, pcaIL)
	      
	      pIL <- matrix(NA, nrow(empPCA$rotation), ncol(empPCA$rotation))
	      for(i in 1:ncol(empPCA$rotation)){
	      	    for(j in 1:nrow(empPCA$rotation)){
		    	  nullIL <- sapply(perIL, function(x) x[j, i])
		    	  pIL[j, i] <- pnorm((empIL[j, i] - mean(nullIL, na.rm = T)) / sd(nullIL, na.rm = T), lower.tail = F)
		    }
	      }
	      
	      colnames(pcnull) <- names(pPC) <- colnames(pIL) <- paste0("PC", 1:length(empPCA$sdev))
	      
	      res <- list(empPCA = empPCA, permImportance = permimp, pPhi = pPhi, pPsi = pPsi, pPCs = pPC, pIL = pIL)
	      return(res)
}