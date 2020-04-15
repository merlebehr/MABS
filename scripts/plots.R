
m <- 2
al <- 0:1
sigma <- 0.01
N <- 2
file_name <- paste0("../results/mse_m_",m, "_sigma_", 100 * sigma, "_al_", max(al), "_N_", N, ".Rdata")
if(file.exists(file_name)){
  load(file = file_name)
}else{
  stop("No such file.")
}

mseNm <- (matrix(sapply(mseN, mean), nrow = length(nN)))


par(mfrow = c(1,1))
ind <- 1:length(nN)
plot(nN[ind], mseNm[ind,1], type = "b", ylim = range(mseNm[ind,]), xlab = "n", ylab = "MSE") 
for(i in 2:ncol(mseNm)){
  lines(nN[ind], (mseNm[ind,i]), type = "b") 
  y <- mseNm[ind,i]
  x <- 1/nN[ind]
  fit <- lm(y ~ x)
  lines(nN[ind], fit$fitted.values, col = "red", lwd = 2)
  
}

