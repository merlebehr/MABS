
m <- 2
al <- 0:2
sigma <- 1
plot_fit <- T

file_name <- file_name <- paste0("mse_m_",m, "_sigma_", 100 * sigma, "_al_", max(al))
file_path <- paste0("../results/")
file_name_all <- list.files(path =  file_path, pattern = file_name)

if(length(file_name_all) > 0){
  load(file = paste0(file_path, file_name_all[1]) )
  mseN_c <- mseN
  accAN_c <- accAN
  for(i in 2:(length(file_name_all) - 0) ){
    load(file = paste0(file_path, file_name_all[i]) )
    if(length(mseN_c) == length(mseN)){
      mseN_c <- mapply(c, mseN_c, mseN, SIMPLIFY=FALSE)
      accAN_c <- mapply(c, accAN_c, accAN, SIMPLIFY=FALSE)
    }
  }
}else{
  stop("No such file.")
}

mc_runs <- length(mseN_c[[1]])
print(paste("In total there are ", mc_runs, "MC runs."))

mseNm <- (matrix(sapply(mseN_c, mean), nrow = length(nN)))
accANm <- (matrix(sapply(accAN_c, mean), nrow = length(nN)))

pdf(file = paste0("../results/plots/plot_",file_name,".pdf"), width = 12, height = 7)
par(mfrow = c(1,2), mar = c(5,5,1,1), oma = c(0,0,3,0), cex = 1.5)
ind <- 1:length(nN)
limY <- range(mseNm[,])
colM <- c("black", "blue", "orange", "darkgreen", "red", "brown", "purple")
plot(nN, mseNm[,1], type = "l", ylim = limY, xlab = "n", ylab = expression(paste("||",hat(Omega) - Omega,"||"^2,""["Inf,2"]/M)), col = colM[1], lwd = 2) 

if(max(al) == 1){
  al_name = "{0,1}"
}else if (max(al) == 2){
  al_name = "{0,1,2}"
}else{
  al_name = paste0("{", min(al), "...", max(al), "}")
}

title(main = paste0("m = ",m, ", alphabet = ", al_name, ", sigma = ",sigma), outer = T) #", number of MC runs = ", mc_runs

if(plot_fit){
  y <- mseNm[ind,1]
  x <- 1/nN[ind]
  fit <- lm(y ~ x - 1)
  lines(nN[ind], fit$fitted.values, col = colM[1], lwd = 1, lty = 2)
  lines(nN, mseNm[,1], col = colM[1], lwd = 2) 
}

for(i in 2:ncol(mseNm)){
  if(plot_fit){
    y <- mseNm[ind,i]
    fit <- lm(y ~ x -1)
    lines(nN[ind], fit$fitted.values, col = colM[i], lwd = 1, lty = 2)
  }
  lines(nN, (mseNm[,i]), col = colM[i], lwd = 2) 
}

name_ylab <- expression(paste('mean( ', F ==  hat(F), ' )'))

expression(paste("||",hat(Omega) - Omega,"||"^2,""["Inf,2"]/M))

ylab_name <- expression("( " * Sigma * ""["i,j"] * "1" * ""[ "{ F" * ""["i,j"] * " = " * hat(F) * ""["i,j"] * "}" ] * " ) / (nm)")


plot(nN[ind], accANm[ind,1], type = "l", ylim = c(0,1), xlab = "n",ylab = ylab_name, col = colM[1], lwd = 2) 
for(i in 2:ncol(mseNm)){
  lines(nN[ind], (accANm[ind,i]), type = "l", col = colM[i], lwd = 2) 
  legend(x = min(nN[ind]) +  0.45 *( max(nN[ind]) - min(nN[ind])), y = 0.5, 
         legend = paste("M = ", MN), col = colM, lty = rep(1, length(MN)),
         x.intersp = 0.8, y.intersp = 0.8)
  abline(h = 1, col = "gray")
}
dev.off()

