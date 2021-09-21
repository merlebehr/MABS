
m <- 2
al <- 1
sigma <- 1

file_name <- file_name <- paste0("comparison_m_",m, "_sigma_", 100 * sigma, "_al_", max(al))
file_path <- paste0("../results/")
file_name_all <- list.files(path =  file_path, pattern = file_name)

if(length(file_name_all) > 0){
  load(file = paste0(file_path, file_name_all[1]) )
  mseN_c <- mseN
  accAN_c <- accAN
  mseImN_c <- mseImN
  runTmN_c <- runTmN
  
  mseN_nmf_c <- mseN_nmf
  accAN_nmf_c <- accAN_nmf
  mseImN_nmf_c <- mseImN_nmf
  runTmN_nmf_c <- runTmN_nmf
  
  mseImN_em_c <- mseImN_em
  runTmN_em_c <- runTmN_em
  
  
  if(length(file_name_all) > 1){
    for(i in 2:(length(file_name_all) - 0) ){
      load(file = paste0(file_path, file_name_all[i]) )
      if(length(mseN_c) == length(mseN)){
        mseN_c <- mapply(c, mseN_c, mseN, SIMPLIFY=FALSE)
        accAN_c <- mapply(c, accAN_c, accAN, SIMPLIFY=FALSE)
        mseImN_c <- mapply(c, mseImN_c, mseImN, SIMPLIFY=FALSE)
        runTN_c <- mapply(c, runTmN_c, runTmN, SIMPLIFY=FALSE)
        
        
        mseN_nmf_c <- mapply(c, mseN_nmf_c, mseN_nmf, SIMPLIFY=FALSE)
        accAN_nmf_c <- mapply(c, accAN_nmf_c, accAN_nmf, SIMPLIFY=FALSE)
        mseImN_nmf_c <- mapply(c, mseImN_nmf_c, mseImN_nmf, SIMPLIFY=FALSE)
        runTN_nmf_c <- mapply(c, runTmN_nmf_c, runTmN_nmf, SIMPLIFY=FALSE)
        
        
        mseImN_em_c <- mapply(c, mseImN_em_c, mseImN_em, SIMPLIFY=FALSE)
        runTN_em_c <- mapply(c, runTmN_em_c, runTmN_em, SIMPLIFY=FALSE)
        
        
      }
    }
  }
}else{
  stop("No such file.")
}

mc_runs <- length(mseN_c[[1]])
print(paste("In total there are ", mc_runs, "MC runs."))

indN_range = 8:length(nN)

mseNm <- (matrix(sapply(mseN_c, mean), nrow = length(nN)))
accANm <- (matrix(sapply(accAN_c, mean), nrow = length(nN)))
mseImNm <- (matrix(sapply(mseImN_c, mean), nrow = length(nN)))
runTNm <- (matrix(sapply(runTmN_c, mean), nrow = length(nN)))


mseNm_nmf <- (matrix(sapply(mseN_nmf_c, mean), nrow = length(nN)))
accANm_nmf <- (matrix(sapply(accAN_nmf_c, mean), nrow = length(nN)))
mseImNm_nmf <- (matrix(sapply(mseImN_nmf_c, mean), nrow = length(nN)))
runTNm_nmf <- (matrix(sapply(runTmN_nmf_c, mean), nrow = length(nN)))


mseImNm_em <- (matrix(sapply(mseImN_em_c, mean, na.rm = TRUE), nrow = length(nN)))
runTNm_em <- (matrix(sapply(runTmN_em_c, mean), nrow = length(nN)))



pdf(file = paste0("../results/plots/plot_comparison_",file_name,".pdf"), width = 12, height = 12)
par(mfrow = c(2,2), mar = c(5,5,1,1), oma = c(0,0,3,0), cex = 1.5)

limY <- range(c(mseNm[indN_range,], mseNm_nmf[indN_range,]))
colM <- c("black", "blue", "orange", "darkgreen", "red", "brown", "purple")
plot(nN[indN_range], mseNm[indN_range,1], type = "l", ylim = limY, xlab = "n", 
     ylab = expression(paste("||",hat(Omega) - Omega,"||"^2,""["Inf,2"]/M)), col = colM[1], lwd = 2) 
lines(nN, mseNm_nmf[,1], col = colM[1], lwd = 2, lty = 2)

if(max(al) == 1){
  al_name = "{0,1}"
}else if (max(al) == 2){
  al_name = "{0,1,2}"
}else{
  al_name = paste0("{", min(al), "...", max(al), "}")
}

title(main = paste0("m = ",m, ", alphabet = ", al_name, ", sigma = ",sigma), outer = T) #", number of MC runs = ", mc_runs

for(i in 2:ncol(mseNm)){
  lines(nN, (mseNm[,i]), col = colM[i], lwd = 2) 
  lines(nN, (mseNm_nmf[,i]), col = colM[i], lwd = 2, lty = 2) 
  
}



#ylab_name <- expression("( " * Sigma * ""["i,j"] * "1" * ""[ "{ F" * ""["i,j"] * " = " * hat(F) * ""["i,j"] * "}" ] * " ) / (nM)")
ylab_name <- expression(paste("||",hat(F) - F,"||"^2,""["2,Inf"]/n))
limY <- range(c(accANm[indN_range,], accANm_nmf[indN_range,]))

plot(nN[indN_range], accANm[indN_range,1], type = "l", ylim = limY, xlab = "n",ylab = ylab_name, col = colM[1], lwd = 2) 
lines(nN, accANm_nmf[,1], col = colM[1], lwd = 2, lty = 2) 

for(i in 2:ncol(mseNm)){
  lines(nN, (accANm[,i]), type = "l", col = colM[i], lwd = 2) 
  lines(nN, (accANm_nmf[,i]), type = "l", col = colM[i], lwd = 2, lty = 2) 
  abline(h = 1, col = "gray")
}




limY <- range(c(mseImNm[indN_range,], mseImNm_nmf[indN_range,], mseImNm_em[indN_range,]), na.rm = T)
colM <- c("black", "blue", "orange", "darkgreen", "red", "brown", "purple")
plot(nN[indN_range], mseImNm[indN_range,1], type = "l", ylim = limY, xlab = "n", 
     ylab = expression(paste("||",hat(F) * hat(Omega) - F * Omega,"||"^2,""[""]/(Mn))), col = colM[1], lwd = 2) 
lines(nN, mseImNm_nmf[,1], col = colM[1], lwd = 2, lty = 2)
lines(nN, mseImNm_em[,1], col = colM[1], lwd = 2, lty = 3)


for(i in 2:ncol(mseNm)){
  lines(nN, (mseImNm[,i]), col = colM[i], lwd = 2) 
  lines(nN, (mseImNm_nmf[,i]), col = colM[i], lwd = 2, lty = 2) 
  lines(nN, (mseImNm_em[,i]), col = colM[i], lwd = 2, lty = 3) 
  
}


limY <- range(c(runTNm[indN_range,], runTNm_nmf[indN_range,], runTNm_em[indN_range[1:2], 1:2]), na.rm = T) #runTNm_em
plot(nN[indN_range], runTNm[indN_range,1], type = "l", ylim = limY, xlab = "n", 
     ylab = "runtime in seconds", col = colM[1], lwd = 2) 
lines(nN, runTNm_nmf[,1], col = colM[1], lwd = 2, lty = 2)
lines(nN, runTNm_em[,1], col = colM[1], lwd = 2, lty = 3)


for(i in 2:ncol(runTNm)){
  lines(nN, (runTNm[,i]), col = colM[i], lwd = 2) 
  lines(nN, (runTNm_nmf[,i]), col = colM[i], lwd = 2, lty = 2) 
  lines(nN, (runTNm_em[,i]), col = colM[i], lwd = 2, lty = 3) 
  
}

legend(x = min(nN[indN_range]) +  0.45 *( max(nN[indN_range]) - min(nN[indN_range])), y = 0.8 * max(limY), 
       legend = c(paste("M = ", MN), c("MABS", "NMF", "EM")), 
       col = c(colM[1:length(MN)], rep("darkgray", 3)), 
       lty = c(rep(1, length(MN)), 1:3),
       x.intersp = 0.8, y.intersp = 0.8)

dev.off()

