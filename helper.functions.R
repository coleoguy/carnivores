# Terrence Sylvester
# pradakshanas@gmail.com
# April 29, 2020
# helper functions

# whenwver there are multiple chromosome values or a range of chromosome values
# are given this function will sample a single chromosome number from the
# given range. Infput for this function is a data table of chromosome numbers

chromSampler <- function(dat){
  # this function will sample a single chromosome number from a given
  # range of chromosome numbers.
  for(i in 1:nrow(dat)){
    if(dat$notesMale2N[i] != ""){
      foo.a <- as.list(unlist(strsplit(dat$notesMale2N[i], split = ",")))
      foo.x <- vector(mode = "list", length = length(foo.a))
      for(j in 1:length(foo.a)){
        if((0 < gregexpr(pattern = "-", foo.a[[j]])[[1]])){
          x <- unlist(strsplit(x = foo.a[[j]], split = "-"))
          foo.x[[j]] <- seq(from = as.numeric(x[1]),
                            to = as.numeric(x[2]),
                            by = 2)
        }else{
          foo.x[[j]] <- as.numeric(foo.a[[j]])
        }
        dat$male2N[i] <-  sample(unlist(foo.x), size = 1)
      }
    }
  }
  return(dat)
}

sp.matches <- function(dat, trees, hyper = F){
  cat("This function assumes that in the dataset provided second column is genera names and third column is species names\n")
  if(hyper == T){
    cat("hyper parameter is set to true. Expecting fourth and fifth columns to be chromosome number and the binary trait respectively\n")
  }
  dat.new <- dat
  
  if(class(trees) == "multiPhylo"){
    cat("randomly sampling a single tree out of the provided list of trees\n")
    tree <- trees[[sample(length(trees),1)]]
  }else{
    tree <- trees
  }
  # get the species level matches
  dat.new.s.m <- dat.new[dat.new[,3] %in%  tree$tip.label,]
  
  # from these species level matches if there are multiple hits for a single
  # species get only 1 hit
  unique.sp.matches <- unique(dat.new.s.m[,3])
  # make an empty data frame to store values unique species matches
  dat.new.s.m.unique <- as.data.frame(matrix(data = NA,
                                             nrow = length(unique(dat.new.s.m[,3])),
                                             ncol = ncol(dat.new)))
  # give column names
  colnames(dat.new.s.m.unique) <- colnames(dat.new)
  # fill in the data frame
  for(i in 1:length(unique.sp.matches)){
    dat.new.temp <- dat.new.s.m[dat.new.s.m[,3] %in% unique.sp.matches[i],]
    dat.new.s.m.unique[i,1] <- sample(dat.new.temp[,1],1) 
    dat.new.s.m.unique[i,2] <- sample(dat.new.temp[,2],1)
    dat.new.s.m.unique[i,3] <- sample(dat.new.temp[,3],1)
    # chromosome number
    if(length(dat.new.temp[,4]) > 1){
      dat.new.s.m.unique[i,4] <- sample(dat.new.temp[,4],1)
    }else{
      dat.new.s.m.unique[i,4] <- dat.new.temp[,4]  
    }
    # hyper parameter
    if(hyper == T){
      if(length(dat.new.temp[,5]) > 1){
        dat.new.s.m.unique[i,5] <- sample(dat.new.temp[,5],1)
      }else{
        dat.new.s.m.unique[i,5] <- dat.new.temp[,5]  
      }
    }
  }
  # subset the original dataset by the species which had no hit to the tree 
  dat.new.unhit <- dat.new[!(dat.new[,3] %in% tree$tip.label),]
  # get the genera names of the species had matches to the tree dataset
  hit.genera <- c(unique(dat.new.s.m[,2]))
  # get the genera names that have no species level matches to the tree dataset
  unhit.gen <- unique(dat.new[,2][!(dat.new[,2] %in% hit.genera)])
  # from the dataset that have species which have no species level matches remove
  # species that share the genus with the species that have a match with the tree
  dat.new.unhit <- dat.new.unhit[dat.new.unhit[,2] %in% unhit.gen,]
  # get the names of the species from the tree that have no species level match
  tree.unhit <- tree$tip.label[!(tree$tip.label %in% unique.sp.matches)]
  # get the genus names of the species hat have no species level match with the
  # trait dataset
  tree.unhit.gen <- c()
  for(i in 1:length(tree.unhit)){
    tree.unhit.gen[i] <- unlist(strsplit(tree.unhit[i], split = "_"))[[1]]
  }
  # get the index of the genera which have species level matches
  tree.hit.gen.num <- which(tree.unhit.gen %in% hit.genera)
  # get the names of the genera that appear on the phylogeny that have no species level match
  tree.true.unhit.gen <- tree.unhit.gen[!(tree.unhit.gen %in% hit.genera)]
  # get genera level matches
  dat.new.g.m <- dat.new.unhit[dat.new.unhit[,2] %in%  tree.true.unhit.gen,]
  # get the names of unique genera level matches
  unique.gn.matches <- unique(dat.new.g.m[,2])
  # from these genera level matches if there are multiple hits for a single
  # genera get only 1 hit
  dat.new.g.m.unique <- as.data.frame(matrix(data = NA,
                                             nrow = length(unique(dat.new.g.m[,2])),
                                             ncol = ncol(dat.new)))
  # give column names
  colnames(dat.new.g.m.unique) <- colnames(dat.new)
  # fill in the data frame
  dat.new.temp <- c()
  for(i in 1:length(unique.gn.matches)){
    dat.new.temp <- dat.new.g.m[dat.new.g.m[,2] %in% unique.gn.matches[i],]
    
    dat.new.g.m.unique[i,1] <- sample(dat.new.temp[,1],1) 
    dat.new.g.m.unique[i,2] <- sample(dat.new.temp[,2],1)
    dat.new.g.m.unique[i,3] <- paste(dat.new.g.m.unique[,2][i], "_sp.", sep = "")
    # chromosome number
    if(length(dat.new.temp[,4]) > 1){
      dat.new.g.m.unique[i,4] <- sample(dat.new.temp[,4],1)
    }else{
      dat.new.g.m.unique[i,4] <- dat.new.temp[,4]  
    }
    # hyper parameter
    if(hyper == T){
      if(length(dat.new.temp[,5]) > 1){
        dat.new.g.m.unique[i,5] <- sample(dat.new.temp[,5],1)
      }else{
        dat.new.g.m.unique[i,5] <- dat.new.temp[,5]  
      }
    }
  }
  # get the index of genera level matches
  tree..unhit.hit.gen.num <- which(tree.unhit.gen %in% dat.new.g.m.unique[,2])
  # get the names of the species as they appear on the phylogeny which have a
  # genera level match
  tree..unhit.hit.sp.names <- tree.unhit[tree..unhit.hit.gen.num]
  # this will be used to rename the names in the tree
  # make a data frame to store the species names of those species that have a
  # genera level match
  unhit.temp.dat.new <- as.data.frame(matrix(data = NA,
                                             nrow = length(tree..unhit.hit.sp.names),
                                             ncol = 2))
  # give column names
  colnames(unhit.temp.dat.new) <- c("gen", "sp")
  # fill in the species names 
  unhit.temp.dat.new$sp <- tree..unhit.hit.sp.names
  # give them the proper naming for genera level match
  for(i in 1:length(tree..unhit.hit.sp.names)){
    unhit.temp.dat.new$gen[i] <- unlist(strsplit(tree..unhit.hit.sp.names[i], split = "_"))[[1]]
  }
  # some genera are duplicated. get the unique genera name and randomly sample
  # a single species to represent these genera
  tree.gen.matches <- as.data.frame(matrix(data = NA,
                                           nrow = length(unique(unhit.temp.dat.new$gen)),
                                           ncol = 2))
  # give column names
  colnames(tree.gen.matches) <- c("gen", "sp")
  # fill in the data table  
  tree.gen.matches$gen <- unique(unhit.temp.dat.new$gen)
  for(i in 1:nrow(tree.gen.matches)){
    tree.gen.matches$sp[i] <- sample(unhit.temp.dat.new$sp[unhit.temp.dat.new$gen == tree.gen.matches$gen[i]],1)
  }
  # rename the genera
  tree.gen.matches$gen <- paste(tree.gen.matches$gen, "_sp.", sep = "")
  # for tree name correction
  tree.names <- as.data.frame(matrix(data = NA,
                                     nrow = (nrow(dat.new.s.m.unique) + nrow(dat.new.g.m.unique)),
                                     ncol = 2))
  colnames(tree.names) <- c("species_name", "name_on_tree")
  # fill in
  tree.names$species_name[1:nrow(dat.new.s.m.unique)] <- tree.names$name_on_tree[1:nrow(dat.new.s.m.unique)] <- dat.new.s.m.unique[,3]
  tree.names$species_name[(nrow(dat.new.s.m.unique)+1):nrow(tree.names)] <- tree.gen.matches$gen
  tree.names$name_on_tree[(nrow(dat.new.s.m.unique)+1):nrow(tree.names)] <- tree.gen.matches$sp
  # get the finalized dataset
  finaldat.new <- rbind(dat.new.s.m.unique, dat.new.g.m.unique)
  # get the results
  results <- list(finaldat.new, tree.names)
  names(results) <- c("chroms", "name_corrections")
  cat("done")
  return(results)
  
}

buildQ <- function (pars, limits){ 
  print("building q-matrix")
  if (length(pars) != 10) 
    stop("pars should have length of 10")
  if (limits[2] < 100) 
    pad <- 2
  if (limits[2] >= 100) 
    pad <- 3
  if (limits[2] < 10) 
    pad <- 1
  parMat <- matrix(0, 2 * length(limits[1]:limits[2]), 
                   2 * length(limits[1]:limits[2]))
  colnames(parMat) <- sprintf(paste("%0", pad, "d", 
                                    sep = ""), 1:ncol(parMat))
  rownames(parMat) <- colnames(parMat)
  split <- ncol(parMat)/2
  chroms <- limits[1]:limits[2]
  for (i in 1:(split - 1)) {
    parMat[i, (i + 1)] <- pars[1]
    parMat[(i + 1), i] <- pars[3]
    if ((chroms[i] * 2) <= max(chroms)) 
      parMat[i, which(chroms == (chroms[i] * 2))] <- parMat[i, 
                                                            which(chroms == (chroms[i] * 2))] + pars[7]
    if ((ceiling(chroms[i] * 1.5)) <= max(chroms)) {
      x <- chroms[i] * 1.5
      if (x%%1 == 0) 
        parMat[i, which(chroms == x)] <- parMat[i, 
                                                which(chroms == x)] + pars[5]
      if (x%%1 != 0) 
        parMat[i, which(chroms %in% c(floor(x), ceiling(x)))] <- parMat[i, 
                                                                        which(chroms %in% c(floor(x), ceiling(x)))] + 
          (pars[5]/2)
    }
    parMat[i, (i + split)] <- pars[9]
    if (i == (split - 1)) 
      parMat[(i + 1), (i + 1 + split)] <- pars[9]
  }
  for (i in (split + 1):(nrow(parMat) - 1)) {
    parMat[i, (i + 1)] <- pars[2]
    parMat[(i + 1), i] <- pars[4]
    if ((chroms[i - split] * 2) <= max(chroms)) 
      parMat[i, (which(chroms[i - split] * 2 == chroms) + 
                   split)] <- parMat[i, (which(chroms[i - split] * 
                                                 2 == chroms) + split)] + pars[8]
    if ((ceiling(chroms[i - split] * 1.5)) <= max(chroms)) {
      x <- chroms[i - split] * 1.5
      if (x%%1 == 0) 
        parMat[i, (which(chroms == x) + split)] <- parMat[i, 
                                                          (which(chroms == x) + split)] + pars[6]
      if (x%%1 != 0) 
        parMat[i, (which(chroms %in% c(floor(x), ceiling(x))) + 
                     split)] <- parMat[i, (which(chroms %in% c(floor(x), 
                                                               ceiling(x))) + split)] + (pars[6]/2)
    }
    parMat[i, (i - split)] <- pars[10]
    if (i == (nrow(parMat) - 1)) 
      parMat[(i + 1), (i + 1 - split)] <- pars[10]
  }
  return(parMat)
}

# edit the constrainMUSSE function in package ChromePlus so that it can constrain 
# the direction of the binary trait transision
constrainMuSSE.custom <- function (data,
                                   lik,
                                   hyper = T,
                                   polyploidy = F,
                                   s.lambda = T, 
                                   s.mu = T,
                                   verbose = F,
                                   oneway = F,
                                   direction = NULL,
                                   constrain = list(drop.poly = F,
                                                    drop.demi = F,
                                                    symmetric = F,
                                                    nometa = F,
                                                    meta = "ARD")) 
{
  if (length(constrain) < 5) {
    if (is.null(constrain$drop.pol)) 
      constrain$drop.poly = F
    if (is.null(constrain$drop.demi)) 
      constrain$drop.demi = F
    if (is.null(constrain$symmetric)) 
      constrain$symmetric = F
    if (is.null(constrain$nometa)) 
      constrain$nometa = F
    if (is.null(constrain$meta)) 
      constrain$meta = "ARD"
  }
  if(oneway == T){
    if(is.null(direction)){
      stop("Direction of the binary trait is not provided. It should be either tran12 or tran21")
    }
    if(direction %in% c("tran12", "tran21") == F){
      stop("Please provide the direction in the correct format. It should be either tran12 or tran21")
    }
  }
  if (ncol(data) < 100) 
    pad <- 2
  if (ncol(data) >= 100) 
    pad <- 3
  if (ncol(data) < 10) 
    pad <- 1
  parMat <- matrix(0, ncol(data), ncol(data))
  colnames(parMat) <- sprintf(paste("%0", pad, "d", 
                                    sep = ""), 1:ncol(parMat))
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
      if ((ceiling(chroms[i] * 1.5)) <= max(chroms)) {
        x <- chroms[i] * 1.5
        if (x%%1 == 0) 
          parMat[i, which(chroms == x)] <- 10
        if (x%%1 != 0) 
          parMat[i, which(chroms %in% c(floor(x), ceiling(x)))] <- 11
      }
      parMat[i, (i + 1)] <- 1
      parMat[(i + 1), i] <- 2
    }
  }
  if (hyper == T & polyploidy == T) {
    print("Constraining model where ploidy is a meta state and different rates of chromosome evolution are possible based on being polyploid or diploid")
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
    print("Constraining model with a hyper state that may have different rates of chromsome number evolution")
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
  rate.table <- as.data.frame(matrix(, nrow(parMat) * ncol(parMat), 
                                     3))
  rate.table[, 1] <- rep(as.character(row.names(parMat)), each = ncol(parMat))
  rate.table[, 2] <- rep(as.character(colnames(parMat)), nrow(parMat))
  rate.table[, 3] <- as.character(c(t(parMat)))
  rate.table <- rate.table[rate.table[, 1] != rate.table[2], 
  ]
  rate.table[rate.table[, 3] == 1, 3] <- "asc1"
  rate.table[rate.table[, 3] == 2, 3] <- "desc1"
  rate.table[rate.table[, 3] == 3, 3] <- "asc2"
  rate.table[rate.table[, 3] == 4, 3] <- "desc2"
  rate.table[rate.table[, 3] == 5, 3] <- "pol1"
  rate.table[rate.table[, 3] == 6, 3] <- "pol2"
  rate.table[rate.table[, 3] == 7, 3] <- "redip"
  rate.table[rate.table[, 3] == 8, 3] <- "tran12"
  rate.table[rate.table[, 3] == 9, 3] <- "tran21"
  rate.table[rate.table[, 3] == 10, 3] <- "dem1"
  rate.table[rate.table[, 3] == 11, 3] <- ".5*dem1"
  rate.table[rate.table[, 3] == 12, 3] <- "dem2"
  rate.table[rate.table[, 3] == 13, 3] <- ".5*dem2"
  if (constrain$nometa == T) {
    rate.table[rate.table[, 3] == "asc2", 3] <- "asc1"
    rate.table[rate.table[, 3] == "desc2", 3] <- "desc1"
    rate.table[rate.table[, 3] == "pol2", 3] <- "pol1"
    rate.table[rate.table[, 3] == "dem2", 3] <- "dem1"
    rate.table[rate.table[, 3] == ".5*dem2", 3] <- ".5*dem1"
  }
  if (constrain$drop.poly == T) {
    rate.table[rate.table[, 3] == "pol1", 3] <- "0"
    rate.table[rate.table[, 3] == "pol2", 3] <- "0"
  }
  if (constrain$drop.demi == T) {
    rate.table[rate.table[, 3] == "dem1", 3] <- "0"
    rate.table[rate.table[, 3] == ".5*dem1", 3] <- "0"
    rate.table[rate.table[, 3] == "dem2", 3] <- "0"
    rate.table[rate.table[, 3] == ".5*dem2", 3] <- "0"
  }
  if (constrain$symmetric == T) {
    rate.table[rate.table[, 3] == "desc1", 3] <- "asc1"
    rate.table[rate.table[, 3] == "desc2", 3] <- "asc2"
    rate.table[rate.table[, 3] == "pol2", 3] <- "pol1"
    rate.table[rate.table[, 3] == 7, 3] <- "redip"
    rate.table[rate.table[, 3] == 8, 3] <- "tran12"
    rate.table[rate.table[, 3] == 9, 3] <- "tran21"
    rate.table[rate.table[, 3] == "dem2", 3] <- "dem1"
    rate.table[rate.table[, 3] == ".5*dem2", 3] <- ".5*dem1"
  }
  if (constrain$meta == "SYM") {
    rate.table[rate.table[, 3] == "tran21", 3] <- "tran12"
  }
  if(oneway == T){
    if(direction == "tran12"){
      rate.table[rate.table[, 3] == "tran21", 3] <- "0"  
    }
    if(direction == "tran21"){
      rate.table[rate.table[, 3] == "tran12", 3] <- "0"  
    }
  }
  formulae <- vector(mode = "character", length = nrow(rate.table))
  for (i in 1:nrow(rate.table)) {
    formulae[i] <- paste("q", rate.table[i, 1], rate.table[i, 
                                                           2], " ~ ", rate.table[i, 3], collapse = "", 
                         sep = "")
  }
  lambda <- mu <- vector()
  for (i in 1:nrow(parMat)) {
    if (hyper == F) {
      lambda <- c(lambda, paste("lambda", colnames(parMat)[i], 
                                " ~ lambda1", sep = ""))
      mu <- c(mu, paste("mu", colnames(parMat)[i], 
                        " ~ mu1", sep = ""))
    }
    if (hyper == T & s.lambda == T) {
      lambda <- c(lambda, paste("lambda", colnames(parMat)[i], 
                                " ~ lambda1", sep = ""))
    }
    if (hyper == T & s.mu == T) {
      mu <- c(mu, paste("mu", colnames(parMat)[i], 
                        " ~ mu1", sep = ""))
    }
    if (hyper == T & s.lambda == F) {
      if (i <= split) {
        lambda <- c(lambda, paste("lambda", colnames(parMat)[i], 
                                  " ~ lambda1", sep = ""))
      }
      if (i > split) {
        lambda <- c(lambda, paste("lambda", colnames(parMat)[i], 
                                  " ~ lambda2", sep = ""))
      }
    }
    if (hyper == T & s.mu == F) {
      if (i <= split) {
        mu <- c(mu, paste("mu", colnames(parMat)[i], 
                          " ~ mu1", sep = ""))
      }
      if (i > split) {
        mu <- c(mu, paste("mu", colnames(parMat)[i], 
                          " ~ mu2", sep = ""))
      }
    }
  }
  extras <- c("asc1", "desc1", "asc2", "desc2", 
              "pol1", "pol2", "dem1", "dem2", 
              "redip", "tran12", "tran21", "lambda1", 
              "mu1", "lambda2", "mu2")
  lik.con <- constrain(lik, formulae = c(formulae, lambda, 
                                         mu), extra = extras)
  colnames(parMat) <- rownames(parMat) <- colnames(data)
  if (verbose == T) 
    return(list(lik.con, parMat))
  if (verbose == F) 
    return(lik.con)
}

# this function will get the radis of each concentric circle (still under construction)
getRadius <- function(scale = NULL,
                      width = NULL,
                      tree = NULL,
                      tip.labels = FALSE,
                      trait.values = NULL,
                      classes = NULL){
  # get the data from last plot phylogeny
  obj <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  # get max hight of the tree
  h <- max(nodeHeights(tree))
  # get the point where the bar starts
  sw <- strwidth("l")
  # get the bar hights
  x <- trait.values * scale
  # get the width of the bars
  w <- width
  
  # make breaks
  # get the  circle limit
  theta <- atan(obj$yy[which.max(trait.values)]/obj$xx[which.max(trait.values)])
  if(obj$xx[which.max(trait.values)] > 0){
    s <- 1
  }else{
    s <- -1
  }
  # get starting X and Y values
  dx <- s * h * cos(theta) + s * cos(theta) * sw
  dy <- s * h * sin(theta) + s * sin(theta) * sw
  x1 <- dx + (w/2) * cos(pi/2 - theta) - s * min(0,min(x)) * cos(theta)
  y1 <- dy - (w/2) * sin(pi/2 - theta) - s * min(0,min(x)) * sin(theta)
  x2 <- dx - (w/2) * cos(pi/2 - theta) - s * min(0,min(x)) * cos(theta)
  y2 <- dy + (w/2) * sin(pi/2 - theta) - s * min(0,min(x)) * sin(theta)
  # get where each circle is
  divisions <- seq(from = 0, to = 1, length.out = (classes + 1))[-1]
  # get ending X and Y values
  x3 <- y3 <- x4 <- y4 <- vector(mode = "numeric", length = length(classes))
  for(i in 1:classes){
    x3[i] <- s * x[which.max(trait.values)] * divisions[i] * cos(theta) + x2
    y3[i] <- s * x[which.max(trait.values)] * divisions[i] * sin(theta) + y2
    x4[i] <- s * x[which.max(trait.values)] * divisions[i] * cos(theta) + x1
    y4[i] <- s * x[which.max(trait.values)] * divisions[i] * sin(theta) + y1
    
  }
  # get the radius of each circle
  radi <- vector(mode = "numeric", length = classes)
  for(i in 1:classes){
    xMax <- max(abs(c(x1, x2, x3[i], x4[i]))) 
    yMax <- max(abs(c(y1, y2, y3[i], y4[i])))
    # pythogores
    radi[i] <- sqrt(xMax^2 + yMax^2) 
  }
  # get lables of each circle
  labs <- round(max(trait.values),0)
  names(radi) <- labs * divisions
  return(radi)
}
