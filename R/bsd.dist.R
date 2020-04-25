bsd.dist <- function(tree1 , tree2, mean.brlen = 0.05){

        # In every case the tree that is rescaled is tree2, which is the shortest
	if(sum(tree1$edge.length) < sum(tree2$edge.length)){
		tree3 <- tree1
		tree1 <- tree2
		tree2 <- tree3
	}

        # This is the objective function for the optimisation implementation with optim()
        tree.dist.opt <- function(x){
            tree3 <- tree2
            tree3$edge.length <- tree2$edge.length*x
            return(KF.dist(tree1, tree3))
        }

        # Optimisation of the tree distance with a starting value of 0
        opt.dist <- optim(0, fn=tree.dist.opt, method = "Brent", lower = 0, upper = 50)

        min.bdi <- opt.dist$value
        scaling <- opt.dist$par

        # Scaling for tree2
        tree2.scaled <- tree2
        tree2.scaled$edge.length <- tree2$edge.length * scaling

        # The trees are scaled so that the mean branch length is 0.05
        # This is an arbitrary value, but is useful so that the total tree lengths are
        # in similar scales
        root.scaling <- mean.brlen / mean(c(tree1$edge.length[tree1$edge.length > 0.00001] , tree2.scaled$edge.length[tree2.scaled$edge.length > 0.00001]))

    	tree1.root.scaled <- tree1
    	tree2.root.scaled <- tree2.scaled

    	tree1.root.scaled$edge.length <- tree1$edge.length * root.scaling
    	tree2.root.scaled$edge.length <- tree2.scaled$edge.length * root.scaling

    	min.bdi.root.scaled <- KF.dist(tree1.root.scaled, tree2.root.scaled)
        
        return(min.bdi.root.scaled)
}