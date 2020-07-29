diagnose.clocks <- function(loctrs, sptr, sp.time.tree = T, branch.support.threshold = 0, log.branches = F, pca = T, mds = F, pca.permutations = 10, mean.scaling.brlen = 0.05, ncore = 1, make.plots = F, sammon.correction = F, clustering.boot.samps = 50, kmax = 10, other.data = NULL, pdf.file = "clock.diagnosis"){
	cat("Collecting clocks", fill = T)
	crx <- collect.clocks(loctrs, sptr, sp.time.tree = T, branch.support.threshold, verbose = F)
	cat("Creating representation of clock space", fill = T)
	crx <- clock.space(crx, sptr, pca, mds, log.branches, pca.permutations, mean.scaling.brlen, ncore, make.plots, sammon.correction, verbose = F)
	cat("Grouping clocks", fill = T)
	crx <- group.clocks(crx, boot.samps = clustering.boot.samps, kmax, make.plots, verbose = !make.plots)
	if(make.plots){ cat("Saving results PDFs", fill = T)
		crx <- write.clocks.plots(crx, loctrs, sptr, other.data, default.variables = c("ntaxa", "nmissingbr", "brsup", "trlen", "rttdist", "topdist"), pdf.file)
	}
    	return(crx)
}