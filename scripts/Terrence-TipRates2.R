GetTipRates2 <- function(tree, Q, tip.states){
  recon <- asr_mk_model( tree = tree,
                         tip_states = tip.states,
                         transition_matrix = Q,
                         Nstates = ncol(Qfilled),
                         include_ancestral_likelihoods = TRUE)
  est <- c()
  for(i in 1:nrow(recon$ancestral_likelihoods)){
    est[i] <- which.max(as.vector(recon$ancestral_likelihoods[i,]))
  }
  ## first get the node numbers of the tips
  tips <- tree$tip.label
  nodes<-sapply(tips,function(x,y) which(y==x),y=tree$tip.label)
  ## then get the edge lengths for those nodes
  edge.lengths<-setNames(tree$edge.length[sapply(nodes,
                                                 function(x,y) which(y==x),y=tree$edge[,2])],names(nodes))
  
  
  tiprates <- nodepulls <- c()
  for(i in 1:length(tree$tip.label)){
    nodepulls[i] <- getParent(tree, nodes[i])-length(tree$tip.label)
  }
  tip.changes <- abs(est[nodepulls]-tip.states)
  names(tip.changes) <- tree$tip.label
  tiprate <- tip.changes/edge.lengths
  return(tiprate)
}