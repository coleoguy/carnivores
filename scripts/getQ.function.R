# This function will generate a Q matrix which has all the parameters that 
# describes chromosome number evolution
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