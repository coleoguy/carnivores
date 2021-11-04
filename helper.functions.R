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
