
# this function takes posterior of trees and data
# it returns a list with pruned data and random sample of possible
# tip values if multiple points are available.
getData <- function(trees, dat){
  tree.genera <- trees[[1]]$tip.label
  good.genera <- unique(dat$X[which(dat$X %in% tree.genera)])
  # get total possible matched data
  missing <- trees[[1]]$tip.label[!trees[[1]]$tip.label %in% dat$X]
  trees.pruned <- list()
  tree.depths <- c()
  for(i in 1:100){
    cur.tree <- drop.tip(trees[[i]], tip = missing)
    tree.depths[i] <- max(branching.times(cur.tree))
    cur.tree$edge.length <-  cur.tree$edge.length / max(branching.times(cur.tree))
    trees.pruned[[i]] <- cur.tree
  }
}