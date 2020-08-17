diagnose.clocks <- function(loctrs, sptr, sp.time.tree = T, branch.support.threshold = 0, branch.length.threshold = 1e-6, log.branches = F, pca = T, mds = F, pca.permutations = 100, mean.scaling.brlen = 0.05, ncore = 1, sammon.correction = F, clustering.boot.samps = 50, kmax = 10, other.data = NULL, pdf.file = "clock.diagnosis"){
	cat("Collecting clocks", fill = T)
	crx <- collect.clocks(loctrs, sptr, sp.time.tree = T, branch.support.threshold, branch.length.threshold, verbose = F)
	cat("Creating representation of clock space", fill = T)
	crx <- clock.space(crx, sptr, pca, mds, sp.time.tree, log.branches, pca.permutations, mean.scaling.brlen, ncore, sammon.correction, verbose = F, pdf.file)
	cat("Grouping clocks", fill = T)
	crx <- group.clocks(crx, boot.samps = clustering.boot.samps, kmax, verbose = is.null(pdf.file), pdf.file)
	if(!is.null(pdf.file)){ cat("Saving results PDFs", fill = T)
		crx <- write.clocks.plots(crx, loctrs, sptr, other.data, default.variables = c("ntaxa", "nmissingbr", "brsup", "trlen", "rttdist", "topdist"), pdf.file)
	}
    	return(crx)
}