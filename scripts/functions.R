
### Estimate cluster centers with Lloyds algorithm initialized by spectral clustering
#
# Input:
# y numeric matrix with observations, where nrow(y) equal to the number of observations and ncol(y) equal to the dimension of an individual observations
# K integer which gives the number of cluster centers
#
# Output:
# Matrix with rows equal to the K estimated centers each of dimension M
clustLloySpec <- function(y, K){
  
  y_svd <- svd(t(y), nu = K, nv = 0)
  hy <- t(y_svd$u %*% t(y_svd$u) %*% t(y))
  
  ind <- sample(1:nrow(hy), K)
  centers <- hy[ind,]
  clusters_init <- kmeans(hy, centers)
  error <- sum((fitted(clusters_init) - hy)^2)
  
  error_best <- error
  centers_best <- centers
  
  for(i in 1:(K*200) ){
    ind1 <- sample(1:nrow(hy),1)
    ind2 <- sample(1:K,1)
    
    if(min(colSums( (t(centers) - hy[ind1, ])^2 )) != 0){
      centers[ind2,] <- hy[ind1,]
      clusters_init <- kmeans(hy, centers)
      error <- sum((fitted(clusters_init) - hy)^2)
      
      if(error < error_best){
        centers_best <- centers
        error_best <- error
      }
    }
    
  }
  
  clusters <- kmeans(y, centers_best, iter.max = max(10, ceiling(4* log(n)) ))$centers
  attr(clusters, "init") <- centers_best

  return(clusters)
}


### Estimates the mixing weights from a linear model with unknown finite alphabet design
#
# Input:
# y numeric matrix with observations, where nrow(y) equal to the number of observations and ncol(y) equal to the dimension of an individual observation
# m integer with the number of sources
# al numeric vector with the alphabet
#
# Output:
# omega numeric matrix with m rows and entries in (0,1) and each of the M columns sums up to one.
estOmega <- function(y, m, al = c(0,1)){
  n <- nrow(y)
  M <- ncol(y)
  K <- length(al)^m
  al <- sort(al)
  

  if (K > n){
    stop("Number of clusters is larger than number of observations.")
  }
  
  y <- (y - al[1])/(al[2] - al[1]) #renormalize so that alphabet is of the form al = (0,1,...)
  al <- (al - al[1])/(al[2] - al[1])

  val <- clustLloySpec(y, K) #initalize cluster centers with spectral clustering
  val <- val[order(rowSums(abs(val)^2)),] #reorder centers by their norm
  
  val <- val[-1,]
  
  omega <- matrix(nrow = m, ncol = M)
  omega[1,] <- val[1,]
  m_run <- 1

  for (i in 2:m) {
    
    alM <-  as.matrix(expand.grid(c(rep(list(al), m_run - 1), list(al[-1]))))
    val_run <- alM %*% omega[1:m_run, ]
    
    for(j in 1:nrow(val_run)){
      ind_best <- which.min( colSums( (val_run[j, ] - t(val))^2) )
      val <- val[-ind_best, ]
    }

    omega[m_run + 1, ] <- val[1, ]
    m_run <- m_run + 1
  }
    
  return(omega)
}


### Generates a random weight matrix 
#
# Input:
# m integer with number of sources
# M integer with dimension of weights
#
# Output:
# omega numeric matrix with m rows and entries in (0,1) and each of the M columns sums up to one.
rOmega <- function(m,M){
  ans <- matrix(nrow=m, ncol=M)
  for (i in 1:M)
    ans[,i] <- (diff(sort(c(0,stats::runif(m-1),1))))
  ord <- order(rowSums(ans))
  ans <- ans[ord, ]
  as.matrix(ans, nrow = m, ncol = M)
}
