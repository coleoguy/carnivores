# This function will generate a Q matrix which has all the parameters that 
# describes chromosoe number evolution
getQ <- function (data = NULL, 
                  hyper = NULL,
                  polyploidy = NULL) {
  parMat <- matrix(0, ncol(data), ncol(data))
  colnames(parMat) <- colnames(data)
  rownames(parMat) <- colnames(parMat)
  split <- ncol(parMat)/2
  if (hyper == T) 
    chroms <- as.numeric(colnames(data)[1:split])
  if (hyper == F) 
    chroms <- as.numeric(colnames(data))
  if (hyper == F) {
    print("Constraining model to simple chromevol version")
    for (i in 1:(nrow(parMat) - 1)) {
      if ((chroms[i] * 2) <= max(chroms)) 
        parMat[i, which(chroms == (chroms[i] * 2))] <- 5
      if (constrain$drop.demi == F) {
        if ((ceiling(chroms[i] * 1.5)) <= max(chroms)) {
          x <- chroms[i] * 1.5
          if (x%%1 == 0) 
            parMat[i, which(chroms == x)] <- 10
          if (x%%1 != 0) 
            parMat[i, which(chroms %in% c(floor(x), ceiling(x)))] <- 11
        }
      }
      parMat[i, (i + 1)] <- 1
      parMat[(i + 1), i] <- 2
    }
  }
  if (hyper == T & polyploidy == T) {
    print("Creating rate matrix for chosen chromosome model")
    for (i in 1:(split - 1)) {
      if ((chroms[i] * 2) <= max(chroms)) 
        parMat[i, (which(chroms[i] * 2 == chroms) + split)] <- 5
      if ((ceiling(chroms[i] * 1.5)) <= max(chroms)) {
        x <- chroms[i] * 1.5
        if (x%%1 == 0) 
          parMat[i, (which(chroms == x) + split)] <- 10
        if (x%%1 != 0) 
          parMat[i, (which(chroms %in% c(floor(x), ceiling(x))) + 
                       split)] <- 11
      }
      parMat[i, (i + 1)] <- 1
      parMat[(i + 1), i] <- 2
    }
    for (i in (split + 1):(nrow(parMat) - 1)) {
      if ((chroms[i - split] * 2) <= max(chroms)) 
        parMat[i, (which(chroms[i - split] * 2 == chroms) + 
                     split)] <- 6
      if ((ceiling(chroms[i - split] * 1.5)) <= max(chroms)) {
        x <- chroms[i - split] * 1.5
        if (x%%1 == 0) 
          parMat[i, (which(chroms == x) + split)] <- 12
        if (x%%1 != 0) 
          parMat[i, (which(chroms %in% c(floor(x), ceiling(x))) + 
                       split)] <- 13
      }
      parMat[i, (i - split)] <- 7
      if (i == (nrow(parMat) - 1)) 
        parMat[(i + 1), (i + 1 - split)] <- 7
      parMat[i, (i + 1)] <- 3
      parMat[(i + 1), i] <- 4
    }
  }
  if (hyper == T & polyploidy == F) {
    print("Creating rate matrix for chosen chromosome model")
    for (i in 1:(split - 1)) {
      if ((chroms[i] * 2) <= max(chroms)) 
        parMat[i, which(chroms == (chroms[i] * 2))] <- 5
      if ((ceiling(chroms[i] * 1.5)) <= max(chroms)) {
        x <- chroms[i] * 1.5
        if (x%%1 == 0) 
          parMat[i, which(chroms == x)] <- 10
        if (x%%1 != 0) 
          parMat[i, which(chroms %in% c(floor(x), ceiling(x)))] <- 11
      }
      parMat[i, (i + split)] <- 8
      if (i == (split - 1)) 
        parMat[(i + 1), (i + 1 + split)] <- 8
      parMat[i, (i + 1)] <- 1
      parMat[(i + 1), i] <- 2
    }
    for (i in (split + 1):(nrow(parMat) - 1)) {
      if ((chroms[i - split] * 2) <= max(chroms)) 
        parMat[i, (which(chroms[i - split] * 2 == chroms) + 
                     split)] <- 6
      if ((ceiling(chroms[i - split] * 1.5)) <= max(chroms)) {
        x <- chroms[i - split] * 1.5
        if (x%%1 == 0) 
          parMat[i, (which(chroms == x) + split)] <- 12
        if (x%%1 != 0) 
          parMat[i, (which(chroms %in% c(floor(x), ceiling(x))) + 
                       split)] <- 13
      }
      parMat[i, (i - split)] <- 9
      if (i == (nrow(parMat) - 1)) 
        parMat[(i + 1), (i + 1 - split)] <- 9
      parMat[i, (i + 1)] <- 3
      parMat[(i + 1), i] <- 4
    }
  }
  return(parMat)
}


# This function will fill the Q matrix produced by getQ function
fillQ <- function(data = NULL, 
                  Q = NULL,
                  hyper = NULL,
                  polyploidy = NULL){
  # here there are rates that we are not interested at
  # discard them from the table
  Q[!(Q %in% c(1,2,3,4,8,9))] <- 0
  
  # fill the qmat
  Q[Q == 1] <- mean(chromplus$asc1)
  Q[Q == 2] <- mean(chromplus$desc1)
  Q[Q == 3] <- mean(chromplus$asc2)
  Q[Q == 4] <- mean(chromplus$desc2)
  Q[Q == 8] <- mean(chromplus$tran12)
  Q[Q == 9] <- mean(chromplus$tran21)
  
  # fill the diagonal so that row sums are zero
  diag(Q) <- -rowSums(Q)
  
  return(Q)
}

# this function will assign a unique value for each state in the chromosome 
# matrix
getChromStates <- function(chrom.mat){
  # get the number of states
  states <- 1:ncol(chrom.mat)
  # assign names to states
  names(states) <- colnames(chrom.mat)
  return(states)
}

# this function will get the tip states in the chromplus format
getTipStates <- function(tree = NULL, 
                         chrom.mat = NULL){
  tipstates <- c()
  for(i in 1:Ntip(tree)){
    tipstates[i] <- as.numeric(names(which(chrom.mat[tree$tip.label[i],] == 1)))
  }
  # assign names for tip states
  names(tipstates) <- tree$tip.label
  return(tipstates)
}

# This function will compute the tip rates of a tree given a Q matrix and 
# tip states
GetTipRates <- function(tree = NULL, 
                        Q = NULL,
                        tip.states = NULL,
                        tip.probability = NULL,
                        hyper=NULL){
  # time the code
  tm.start <- Sys.time()
  
  #### ---- define inputs ---- ####
  
  # tree > a rooted tree of class phylo
  # Q > a transition matrix
  # tip.states > a vector of size Ntips specifying the state of each tip
  # tip.probability > matrix output from datatoMatrix function in ChromPlus
  # hyper > logical. whether data have a binary hyper state
  
  #### ---- define inputs ---- ####
  
  # set defaults for missing arguments
  if(is.null(hyper)){
    hyper <- F
  }
  
  # set tip conversion as false.
  # this is used if the tip states are needed to be converted to integers
  tipConversion <- F
  
  #### perform checks ####
  
  # following checks are done if a hyper state is present in data.
  # if not following checks are ignored
  
  # additionally following checks are performed correct data in several ways
  # tip states are not provided but tip probabilities are provided
  # here tip states will be calculated based on the tip probability matrix
  # assumes a single state for a given tip
  # tip states are provided but not numeric / integers
  # tip states will be converted to a numeric vector based on the 
  # tip probability matrix
  # Column names and row names of the Q matrix are not integers
  # column and row names are converted to integers
  
  if(hyper == T){
    
    # if a hyper trait is present then tip probabilities are required.
    # this is the output from datatomatrix function in ChromPlus
    if(is.null(tip.probability)){
      stop("tip probabilities are not provided")
    }
    
    # next check is to see if the tip states are integers
    # if tip states are not numeric then it is converted to integers
    if(is.numeric(tip.states) == F){
      if(is.null(tip.probability) == F){
        cat("\n Tip states empty or are not integers. Converting of tip states to integers\n")
        # assign states a value from 1 to the number of states
        states <- 1:ncol(tip.probability)
        # assign names to states
        names(states) <- colnames(tip.probability)
        # change tip probability column names accordingly
        colnames(tip.probability) <- 1:ncol(tip.probability)
        # make a new tip state vector
        tip.state.int <- vector(mode = "numeric", length = length(tip.states))
        # convert tip states to integers
        # if tip states are not provided then get the tip states based on the
        # tip probability matrix
        if(is.null(tip.states)){
          cat("\n Tip states are empty. Assigning tip states based on the provided tip probability matrix\n")
          for(i in 1:Ntip(tree)){
            tip.state.int[i] <- names(which(tip.probability[tree$tip.label[i],] == 1))
          }
          # assign names for tip states
          names(tip.state.int) <- tree$tip.label
          tipConversion <- T
          # if tip states are provided then assign the proper integer value to
          # tip states
        }else{
          cat("\n Tip states are provided but are not integers. Converting tip states to integers\n")
          for(i in 1:length(tip.state.int)){
            hit <- which(names(states) == tip.states[i])
            tip.state.int[i] <- as.numeric(states[hit])
            tipConversion <- T
          }
        }
      }
    }
    # come back to this code chunk and make corrections
    # idea here is that to make sure all the states in the tip states are
    # present in the Q matrix
    # # check to see if the Q matrix column and raw names need to be converted
    # if(sum(tip.state.int %in% colnames(Q)) != length(tip.state.int)){
    #   cat("\n States in Q matrix are not integers. Converting of Qmatrix states to integers\n")
    #   colnames(Q) <- rownames(Q) <- 1:ncol(Q)
    # }
    # # check if the state names in the Q matrix and the state names assigned 
    # # using the chromosome probability matrix matches or not
    # if(sum(states %in% colnames(Q)) != length(states)){
    #   stop("All the states are not present in the Q matrix")
    # }
  }
  # if the tip conversion was done use the new tip states for calculation of
  # tip rates
  if(tipConversion == T){
    tip.states <- tip.state.int
    tip.states <- as.numeric(tip.states)
  }
  # reconstruct the ancestral states
  recon <- asr_mk_model(tree = tree,
                        tip_states = tip.states,
                        transition_matrix = Q,
                        Nstates = ncol(Q),
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
  
  if(hyper == T && tipConversion == T){
    # assign correct chromosome number for tips and nodes
    est.chroms <- tip.chroms <- c()
    # get the correct value for chromosome number at the nodes
    for(i in 1:length(est)){
      hit.est <- names(which(states == est[i]))
      est.chroms[i] <- as.numeric(gsub("h", "", hit.est))
    }
    # get the correct value for chromosome number at the tips
    for(i in 1:length(tip.states)){
      hit.tip <- names(which(states == tip.states[i]))
      tip.chroms[i] <- as.numeric(gsub("h", "", hit.tip))
    }
    tip.changes <- abs(est.chroms[nodepulls]-tip.chroms)
  }else{
    tip.changes <- abs(est[nodepulls]-tip.states)
  }  
  # assign names for tip changes
  names(tip.changes) <- tree$tip.label
  # get the tip rates
  tiprate <- tip.changes/edge.lengths
  tm.end <- Sys.time()
  cat("\n Process ended in", difftime(time2 = tm.start, time1 = tm.end, units = "secs")[[1]], "seconds.\n")
  return(tiprate)
}

GetTipRates2 <- function(tree, Q, tip.states){
  recon <- asr_mk_model( tree = tree,
                         tip_states = tip.states,
                         transition_matrix = Q,
                         Nstates = ncol(Q),
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
