data <- readRDS(file = '../data/haplotypeData.rds')
source('functions.R')

y <- data[[1]] 
omega <- data[[2]]
f <- data[[3]]
al <- 0:1


set.seed(1)
ans <- estOmega(y, m = 5)
baseline <- attr(ans, "baseline")

y <- t(t(y) - baseline)
estOm <- estOmega(y, m = 5)

Am <- as.matrix(expand.grid(rep(list(al), nrow(estOm))))
valEst <- Am %*% estOm
indEstA <- sapply(1:nrow(y), function(i) which.min(colSums( (y[i,] - t(valEst) )^2 ) ))
estA <- Am[indEstA, ]


colM <- c("black", "blue", "orange", "darkgreen", "red", "brown", "purple")
ind <- order(rowSums(omega), decreasing = T)[1:5]



pdf(file = "../results/plots/haplotypes_weights.pdf", width = 7, height = 7)
plot(estOm[1,], type = "l", ylim = range(c(estOm, omega)), col = colM[1], lwd = 2, xlab = "j", 
     ylab = expression(paste(omega, ""["i,j"])))
lines(omega[ind[1], ], col = colM[1], lty = 2, lwd = 2)
for(i in 2:nrow(estOm)){
  lines(estOm[i,], col = colM[i], lwd = 2)
  lines(omega[ind[i], ], col = colM[i], lty = 2, lwd = 2)
}
legend(x = 1.5, y = 0.55, 
       legend = paste("i = ", 1:5), col = colM[1:5], lty = rep(1, 5),
       x.intersp = 0.8, y.intersp = 0.8, lwd =2)
dev.off()


print(colMeans(estA == f[,ind]))

pdf(file = "../results/plots/haplotypes_structure.pdf", width = 7, height = 7)
par(mfrow = c(2,1), mar = c(4,4,3,1))
image(y = 1:ncol(estA), x = 1:nrow(estA), f[,ind], ylab = "true haplotypes", xlab = "", 
      col = c("white", "black"))

image(y = 1:ncol(estA), x = 1:nrow(estA),  abs(estA - f[,ind]), ylab = "residuals ", xlab = "SNPs", 
      col = c("white", "black"))
dev.off()
