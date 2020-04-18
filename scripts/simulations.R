#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
print('start running job')

library(combinat)
source('functions.R')

args <- as.numeric(args)
print("args = ")
print(args)

m <- args[1]
sigma <- args[2]/100
al <- 0:args[3]
N <- args[4]

print(paste('m = ', m, 'sigma = ', sigma, 'al = 0:', max(al), 'N = ', N))


K <- length(al)^m
Am <- as.matrix(expand.grid(rep(list(al), m)))

MN <- c(10, 20, 100)
nN <- c(50, 60, 70, 80, 90, seq(100, 500, 50))

mseN <- vector('list', length(MN) * length(nN))

ptm <- proc.time()
for(l in 1:length(MN)){
  M <- MN[l]
  print(paste('M = ', M, 'iteration', l, 'out of', length(MN)))
  
  for(j in 1:length(nN)){
    
    n <- nN[j]
    print(paste('n = ', n, 'iteration', j, 'out of', length(nN)))
    
    set.seed(1)
    mse <- numeric(N) + 0
    for(k in 1:N){
      trOm <- rOmega(m , M)
      trA <- matrix(sample(al, n * m, replace = T), ncol = m)
      y <- trA %*% trOm + rnorm(n * M, sd = sigma)
      
      estOm <- estOmega(y, m, al)
      
      perm <- permn(1:m)
      mse_perm <- numeric(length(perm))
      for(i in 1:length(perm)){
        mse_perm[i] <- max(rowSums((estOm[perm[[i]],] - trOm)^2))/M
      }
      mse[k] <- min(mse_perm)
    }
    mseN[[length(nN)*(l - 1) + j]] <- mse
  }
}
time <- round(proc.time() - ptm,2)
print(paste("total run time was", time[1], time[2], time[3]))

save(file = paste0("../results/mse_m_",m, "_sigma_", 100 * sigma, "_al_", max(al), "_N_", N, ".Rdata"), mseN, MN, nN, m, sigma, al, N)

