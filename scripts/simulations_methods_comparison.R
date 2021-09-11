#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE) 

#arguments passed to the script should be: 
# 1. m (number of souces), 
# 2. 100 * sigma (standard deviation), 
# 3. maximum alphabet value k (alphabet = 0,...,k),
# 4. number of Monte Carlo runs
# 5. random seed

print('start running job')

library(combinat)
library(NMF)
library(mixtools)
source('functions.R')

args <- as.numeric(args)
print("args = ")
print(args)

m <- args[1]
sigma <- args[2]/100
al <- 0:args[3]
N <- args[4]
seed <- args[5]

print(paste('m = ', m, 'sigma = ', sigma, 'al = 0:', max(al), 'N = ', N, 'seed = ', seed))


K <- length(al)^m
Am <- as.matrix(expand.grid(rep(list(al), m)))

MN <- c(20, 50)
nN <- c(50, 60, 70, 80, 90, seq(100, 500, 50), seq(600, 1000, 200), 1500, 2000)


mseN <- vector('list', length(MN) * length(nN))
accAN <- vector('list', length(MN) * length(nN))
mseImN <- vector('list', length(MN) * length(nN))
runTmN <- vector('list', length(MN) * length(nN))


mseN_nmf <- vector('list', length(MN) * length(nN))
accAN_nmf <- vector('list', length(MN) * length(nN))
mseImN_nmf <- vector('list', length(MN) * length(nN))
runTmN_nmf <- vector('list', length(MN) * length(nN))


mseImN_em <- vector('list', length(MN) * length(nN))
runTmN_em <- vector('list', length(MN) * length(nN))


ptm <- proc.time()
for(l in 1:length(MN)){
  M <- MN[l]
  print(paste('M = ', M, 'iteration', l, 'out of', length(MN)))
  
  for(j in 1:length(nN)){
    
    n <- nN[j]
    print(paste('n = ', n, 'iteration', j, 'out of', length(nN)))
    
    set.seed(seed)
    mse <- numeric(N) + 0
    accA <- numeric(N) + 0
    mseIm <- numeric(N) + 0
    runT <- numeric(N) + 0
    
    
    mse_nmf <- numeric(N) + 0
    accA_nmf <- numeric(N) + 0
    mseIm_nmf <- numeric(N) + 0
    runT_nmf <- numeric(N) + 0
    
    mseIm_em <- numeric(N) + 0
    runT_em <- numeric(N) + 0
    
    
    for(k in 1:N){
      trOm <- rOmega(m , M)
      trA <- matrix(sample(al, n * m, replace = T), ncol = m)
      y <- trA %*% trOm + rnorm(n * M, sd = sigma)
      
      y_pos = apply(y, c(1,2), function(x) pmin(1, pmax(10e-5, x)))
      
      timeNow = proc.time()
      ans_nmf = nmf(y_pos, m)
      estOm_nmf = ans_nmf@fit@H
      estA_nmf = ans_nmf@fit@W
      runT_nmf[k] = max(proc.time() - timeNow)
      
      timeNow = proc.time()
      ans_em = try(mvnormalmixEM(y, k = length(al)^m, arbvar = FALSE), silent = T)
      if(class(ans_em) == "mixEM"){
        estIm_em = ans_em$mu[apply(ans_em$posterior, 1, function(x) which.max(x))]
        estIm_em = as.data.frame(do.call(rbind, estIm_em))
      }else{
        estIm_em = NA #t(matrix(rep(colMeans(y), nrow(y)), ncol(y)))
      }
      runT_em[k] = max(proc.time() - timeNow)
      
      
      timeNow = proc.time()
      estOm <- estOmega(y, m, al)
      valEst <- Am %*% estOm
      indEstA <- sapply(1:nrow(y), function(i) which.min(colSums( (y[i,] - t(valEst) )^2 ) ))
      estA <- Am[indEstA, ]
      runT[k] = max(proc.time() - timeNow)
      
      
      perm <- permn(1:m)
      
      mse_perm <- numeric(length(perm))
      accA_perm <- numeric(length(perm))

      mse_nmf_perm <- numeric(length(perm))
      accA_nmf_perm <- numeric(length(perm))

      for(i in 1:length(perm)){
        mse_perm[i] <- max(rowSums((estOm[perm[[i]],] - trOm)^2))/M
        accA_perm[i] <- max(colSums((estA[,perm[[i]]] - trA)^2))/n
        
        mse_nmf_perm[i] <- max(rowSums((estOm_nmf[perm[[i]],] - trOm)^2))/M
        accA_nmf_perm[i] <- max(colSums((estA_nmf[,perm[[i]]] - trA)^2))/n
        
      }
      mse[k] <- min(mse_perm)
      accA[k] <- min(accA_perm)
      mseIm[k] <- sum((estA %*% estOm - trA %*% trOm)^2)/(M * n)
      
      mse_nmf[k] <- min(mse_nmf_perm)
      accA_nmf[k] <- min(accA_nmf_perm)
      mseIm_nmf[k] <- sum((estA_nmf %*% estOm_nmf - trA %*% trOm)^2)/(M * n)
      
      mseIm_em[k] <- sum((estIm_em - trA %*% trOm)^2)/(M * n)
      
    }
    mseN[[length(nN)*(l - 1) + j]] <- mse
    accAN[[length(nN)*(l - 1) + j]] <- accA
    mseImN[[length(nN)*(l - 1) + j]] <- mseIm
    runTmN[[length(nN)*(l - 1) + j]] <- runT
    
    
    mseN_nmf[[length(nN)*(l - 1) + j]] <- mse_nmf
    accAN_nmf[[length(nN)*(l - 1) + j]] <- accA_nmf
    mseImN_nmf[[length(nN)*(l - 1) + j]] <- mseIm_nmf
    runTmN_nmf[[length(nN)*(l - 1) + j]] <- runT_nmf
    
    
    mseImN_em[[length(nN)*(l - 1) + j]] <- mseIm_em
    runTmN_em[[length(nN)*(l - 1) + j]] <- runT_em
    
    
  }
}
time <- round(proc.time() - ptm,2)
print(paste("total run time was", time[1], time[2], time[3]))

save(file = paste0("../results/comparison_m_",m, "_sigma_", 100 * sigma, "_al_", max(al), "_N_", N, "_seed_", seed, ".Rdata"), 
     mseN, mseN_nmf,
     accAN, accAN_nmf,
     mseImN, mseImN_nmf,mseImN_em,
     runTmN, runTmN_nmf, runTmN_em,
     MN, nN, m, sigma, al, N)

