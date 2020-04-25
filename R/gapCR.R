# This function selects the best number of clusters according to the following criterion: If there is a significant decline in gap from one cluster to two, then there is a single cluster; otherwise, the largest increase in gap leads to the optimum number of clusters. 

gapCR <- function(gaps, SEs){

      if(gaps[1] > (gaps[2] + SEs[2])){
	optimgap <- 1
      } else {
      	if(any(gaps == Inf)){ 
		    gaps <- gaps[1:(which(gaps == Inf)[1] - 1)]
		    SEs <- SEs[1:length(gaps)]
	}
      	gapdiffs <- sapply(1:(length(gaps)-1), function(i) gaps[i+1] - gaps[i])
	optimgap <- which(gapdiffs == max(gapdiffs)) + 1
	if((gaps[(optimgap - 1)] + SEs[(optimgap - 1)]) > gaps[optimgap]){
		print("The biggest difference in gaps is not significant! Returning 1")
		return(1)
	}
      }
      return(optimgap)
}