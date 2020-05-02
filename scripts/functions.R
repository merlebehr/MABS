
### Estimate cluster centers with Lloyds algorithm initialized by spectral clustering
#
# Input:
# y numeric matrix with observations, where nrow(y) equal to the number of observations and ncol(y) equal to the dimension of an individual observations
# K integer which gives the number of cluster centers
#
# Output:
# Matrix with rows equal to the K estimated centers each of dimension M
clustLloySpec <- function(y, K){
  
  library(mgcv)
  
  y_svd <- svd(t(y), nu = K, nv = 0)
  hy <- t(y_svd$u %*% t(y_svd$u) %*% t(y))
  hy <- uniquecombs(hy)
  
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
  
  clusters <- kmeans(y, centers_best, iter.max = max(10, ceiling(4* log(nrow(y))) ))$centers
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
  
  baseline <- val[1,]
  
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
    
  ind <- order(rowSums(omega), decreasing = T)
  omega <- omega[ind,]
  
  attr(omega, "baseline") <- baseline
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


### Lloyd's algorithms for iterative updates
lloyd <- function (y, m, al, omega) {
  library(lsei)
  
  tsh   <-  1e-3
  M <- ncol(y)
  n <- nrow(y)
  
  omdiff <- tsh
  biasV <- rep(0, M)
  

  
  maxIt   <- 100 # 10e2
  countIt <- 0
  Al <- as.matrix(expand.grid(rep(list(al), m)))
  
  while (omdiff >= tsh && countIt <= maxIt) {
    
    val <- Al %*% omega[1:m,]
    pi <- sapply(1:nrow(y), function(x) which.min(colSums((y[x,] - t(val))^2)))
    

    omegaOld <- omega
    
    A = Al[pi,]
    E = rbind(diag(m),-diag(m), rep(-1,ncol(A)))
    f = c(rep(0,m),rep(-1,m), -1)

    aux <- sapply(1:M, function(x) lsei(a = A, b = y[,x], e = E, f = f))
    omega <- aux[1:m, ]
   
    omdiff  <- max(abs(omega- omegaOld))
    countIt <- countIt + 1
  }
  
  if (omdiff >= tsh)
    warning(paste("Lloyds algorithm not converged with omega diff =",
                  round(omdiff, digits = 6), "and threshold =", tsh))
  
  ord   <- order(rowSums(omega),decreasing=TRUE)
  omega <- omega[ord,]
  pi    <- sapply(pi, function(x) which(colSums((Al[x,ord] - t(Al))^2) == 0))
  
  estOm <- omega
  estA <- Al[pi,]
  colnames(estOm) <- colnames(y)
  rownames(estA) <- rownames(y)
  return(list(estOm = estOm, estA = estA))
}




