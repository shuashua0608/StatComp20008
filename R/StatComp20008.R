#' @title Bart Simpson density function
#' @description Bart Simpson density function f(x)=1/2 phi(x;0,1)+1/10 sum_{j=0}^4 phi(x;j/2-1,1/10)
#' @param  x the chosen interval
#' @return Bart Simpson function\code{n}
#' @examples
#' \dontrun{
#' x<-seq(-1,1,0.01)
#' plot(x,BS(x))
#' }
#' @export
BS <- function(x){
  # probability of BS density
  f <- 1/2 * dnorm(x, 0, 1)
  for (j in 0:4) {
    f <- f + 1/10 * dnorm(x, j/2-1, 1/10)
  }
  return(f)
}

#' @title Generate samples from Bart Simpson density and draw histogram
#' @description Get random points generated from Bart Simpson density function and draw histogram accordingly.
#' @param n number of desired random sample points
#' @param h value chosen for bandwidth
#' @return  Bart Simpson histogram\code{n}
#' @export
HBS<- function(n,h){
  # random points from BS density
  u <- runif(n)
  y <- u
  i <- which(u > 0.5)
  #index for those generated from N(0,1)
  y[i] <- rnorm(length(i), 0, 1)
  for (j in 0:4) {
    i <- which(u > j * 0.1 & u <= (j+1) * 0.1)
    #index for those generated from N(j/2-1,1/10^2)
    y[i] <- rnorm(length(i), j/2 -1, 1/10)
  }
  x<-seq(-4,4,h)
  return(hist(y, breaks = x, probability = TRUE))
}


#' @title Use naive density estimator to estimate BS density,
#' @description Use naive density estimator to estimate density,
#' @param  x estimate points
#' @param  y random samples
#' @param  h bandwidth
#' @return estimated density at x \code{n}
#' @export
NBS <- function(x, y, h){
  m <- length(x)
  n <- length(y)
  ye<- rep(0, m)
  for (i in 1:n){
    ye <- ye + as.numeric((x >= y[i] - h) & (x < y[i] + h))
  }
  ye <- ye / (2*h*n)
  return(ye)
}


