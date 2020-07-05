library(readxl)
library(tidyverse)
library(combinat)
library(Clomial)
source('functions.R')

# Example from Clomial R package as in 
# "Inferring clonal composition from multiple sections of a breast cancer", Zare et al., 2014
set.seed(1)
data(breastCancer)
Dc <- breastCancer$Dc
Dt <- breastCancer$Dt

# m = 2
ClomialResult <-Clomial(Dc=Dc,Dt=Dt,maxIt=20,C=3,binomTryNum=2)
chosen <- choose.best(models=ClomialResult$models)
M1 <- chosen$bestModel
print("Genotypes:")
round(M1$Mu)
estA_clomial <- M1$Mu
print("Clone frequencies:")
round(M1$P,2)
estOm_clomial <- M1$P

# m = 3
ClomialResult <-Clomial(Dc=Dc,Dt=Dt,maxIt=20,C=4,binomTryNum=2)
chosen <- choose.best(models=ClomialResult$models)
M1 <- chosen$bestModel
print("Genotypes:")
round(M1$Mu)
estA_clomial_m3 <- M1$Mu
print("Clone frequencies:")
round(M1$P,2)
estOm_clomial_m3 <- M1$P

print("MSE Zare:")
sqrt(mean((Dt/Dc - estA_clomial %*% estOm_clomial)^2))
sqrt(mean((Dt/Dc - estA_clomial_m3 %*% estOm_clomial_m3)^2))




# load full data set
# data can be downloaded from "https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003703#s5"
data <- vector('list', 14)
for(i in 1:length(data)){
  data[[i]] <- read_excel("../data/pcbi_data2.xls", sheet = i)
}

#extract somatic variants targeted by Zare et al
ind_somVar <- which(!is.na(data[[1]][,4]))
names_somVar <- unlist(data[[1]][ind_somVar,4])

#extract major variant from normal tissue
ind_major <- data[[14]] %>% select(A,T,C,G) %>% apply(1, which.max) 

#extract total counts per site
total_data <- sapply(c(14, 2:13), function(i) sapply(1:length(ind_major), function(j)  unlist(data[[i]]$total[j])))

#extract major count per site
major_data <- sapply(c(14,2:13), function(i) sapply(1:length(ind_major), function(j)  unlist(data[[i]][j, ind_major[j] + 1])))

#compute alternative allele frequency as 1 - major allele frequency
af_data <- (total_data - major_data)/total_data

rownames(af_data) <- unlist(data[[1]][,1])
colnames(af_data) <- c("N", "P1-1", "P1-2", "P1-3", "P1-4", "P2-1", "P2-2", "P2-3", "P2-4", "P3-1", "P3-2", "M1-1", "M1-2")

y <- af_data[ind_somVar,]

# compare with breastCancer data set form Clomial R Package
order_sites <- c(1:5,8,9,10,11,6,7,13,12)
af_clomial <- (Dt/Dc)[,order_sites]
estOm_clomial <- estOm_clomial[, order_sites]
colnames(estOm_clomial) <-colnames(y)
estOm_clomial_m3 <- estOm_clomial_m3[, order_sites]
colnames(estOm_clomial_m3) <-colnames(y)

order_variants <- match(rownames(y), rownames(af_clomial))
af_clomial <- af_clomial[order_variants,]
estA_clomial <- estA_clomial[order_variants, ]
estA_clomial_m3 <- estA_clomial_m3[order_variants, ]

# data from Clomial R Package coincides with y
round(y - af_clomial, 2)
round(y,2)
rownames(y) <- names_somVar

# extract row indeces which are not from somatic variants, have no NAs, and are homozygous for normal tissue
ind_noNA <- setdiff(which(apply(af_data, 1, function(x) !any(is.na(x)) & x[13] == 0)), ind_somVar)
n_add <- 20 #number of randomly selected non-somatic variants
ind_add <- sample(ind_noNA, n_add)

# add randomly selected non-somatic variants
y <- rbind(y, af_data[ind_add, ])



# two sources and binary alphabet
m <- 2
al <- 0:1
Am <- as.matrix(expand.grid(rep(list(al), m)))

# estimate weights
estOm <- estOmega(y = y, m = m, al = al)

# estimate sources
valEst <- Am %*% estOm
indEstA <- sapply(1:nrow(y), function(i) which.min(colSums( (y[i,] - t(valEst) )^2 ) ))
estA <- Am[indEstA, ]

round(estOm,2) 

# compare mse 
print("Total mse:")
sqrt(mean((y - estA %*% estOm)^2))

print("Mse on somatic variants:")
sqrt(mean((y[1:length(ind_somVar)] - (estA %*% estOm)[1:length(ind_somVar),] )^2))

print("Mse from clomial")
sqrt(mean((y[1:length(ind_somVar)] - estA_clomial %*% estOm_clomial)^2))


ans <- lloyd(y, m, al = al, omega = estOm)
round(ans$estOm, 2) 
ans$estA

print("Total mse with Lloyds:")
sqrt(mean((y - (ans$estA %*% ans$estOm))^2))

print("Mse on somatic variants with Lloyds:")
sqrt(mean((y[1:length(ind_somVar),] - (ans$estA %*% ans$estOm)[1:length(ind_somVar),])^2))

# Results from Zara et al manuscript
estOm_zara <- rbind(c(1, 47, 43, 57, 75, 77, 51,61, 55, 44, 63, 77, 84),
                    c(0, 51, 55, 43, 25, 21, 48, 38, 44, 54, 36, 0, 0),
                    c(0,2,3,0,0,1,1, 0,0,1,2, 22, 16))/100

estA_zara <- t(matrix(c(1,0,  1,0,  1,1,  1,0,  1,0,  1,0,  1,0,  1,0,  1,0,  0,1,  1,1,  0,1,  1,0,  0,1,  0,1, 0,1,  0,1), nrow = 2))
estA_zara <- cbind(rep(0, nrow(estA_zara)), estA_zara)

print("Mse from Zara et al manuscript")
sqrt(mean((y[1:length(ind_somVar)] - estA_zara %*% estOm_zara)^2))


cbind(estA_zara[,-1], estA_clomial[, c(3,2)], ans$estA[1:17, c(2,1)])
round(y[1:17, ],2)



###### make plots
m <- 2
al <- c(0,0.5,1)
estOm <- estOmega(y = y, m = m, al = al)
ans <- lloyd(y, m, al = al, omega = estOm)
mse <- sqrt(mean((y[1:length(ind_somVar),] - (ans$estA %*% ans$estOm)[1:length(ind_somVar),])^2))
res <- y[1:length(ind_somVar),] - (ans$estA %*% ans$estOm)[1:length(ind_somVar),]
res_all <- y - (ans$estA %*% ans$estOm)


pdf(file = paste0("../results/plots/cancer_example_m_", m, "al_length_", length(al),"_mse_",round(mse,3)*1000, ".pdf"), width = 15, height = 9) 
layout(rbind(1,2), heights = c(2,1))
col_clone <- c("black", "red", "blue")
pch_clones <- c(19, 17, 15)

plot(ans$estOm[1,], axes = F, xlab = "location", ylab = "frequency", col = col_clone[1], pch = pch_clones[1], ylim = c(0,0.65))
for(i in 2:nrow(ans$estOm)){
  lines(ans$estOm[i,], col = col_clone[i], pch = pch_clones[i], type = "p")
}
axis(1, at = 1:ncol(ans$estOm), labels = colnames(ans$estOm))
axis(2, at = seq(0,1,0.1))
legend(1,0.64, legend = paste("clone (C)", 1:nrow(ans$estOm)), col = col_clone, 
       lty = rep(1, nrow(ans$estOm)), lwd = rep(2, nrow(ans$estOm)))
box()

if(length(al) == 2){
  col_fill <- c("white", "gray")
}else{
  col_fill <- c("white", "gray", "black")
}

if(m == 2){
  at_clones <- c(0,1)
}
if(m == 3){
  at_clones <- c(0,0.5,1)
}
at_location <- seq(0,16/17,1/17)

image(at_location, at_clones, ans$estA[1:17,], col = col_fill, axes = F, xlab = "SNP", ylab = "clone")
axis(1, at =at_location, labels = names_somVar)
#axis(1, at = 16/17, labels = names_somVar[17])
axis(2, at = at_clones, labels = paste("C", 1:ncol(ans$estA)))
if(m == 3){
  axis(2, at = at_clones[2], labels = "C2")
}

box()
dev.off()

data <- data.frame(res = as.numeric(res), gene = rep(rownames(res), ncol(res) ), 
                   location = rep(colnames(res), each = nrow(res)))
ggplot(data, aes(location, gene)) +
  geom_raster(aes(fill = res)) +
  scale_fill_gradientn(colours=c("blue","white", "red"), limits = c(-0.25,0.25))

ggsave( paste0("../results/plots/residuals_cancer_example_", length(al), "_m", m, ".png"))


#------------ simulate with estimate as ground truth

trOm <- ans$estOm
trA <- ans$estA
n <- nrow(ans$estA)
M <- ncol(trOm)
sigma <- sd(res_all[1:length(ind_somVar),]) 
sigma0 <- sd(res_all[(length(ind_somVar) + 1):nrow(res_all),])

set.seed(5)
K <- 100
estOm_K <- vector('list', K)
estA_K <- vector('list', K)
for(k in 1:K){
  noise <-  rbind( matrix(rnorm(length(ind_somVar) * M, sd = sigma), ncol = ncol(trOm)), 
                   matrix(rnorm((nrow(res_all) - length(ind_somVar)) * M, sd = sigma0), nrow = n - length(ind_somVar), ncol = ncol(trOm))) # sample(res_all, n * M, replace = T)
  y <- trA %*% trOm +  noise
  round(trOm,2)
  estOm <- estOmega(y = y, m = m, al = al)
  round(estOm,2)
  ansSim <- lloyd(y, m, al = al, omega = estOm)
  round(ansSim$estOm,2)
  
  estOm <- ansSim$estOm
  estA <- ansSim$estA
  
  perm <- permn(1:m)
  mse_perm <- numeric(length(perm))
  accA_perm <- numeric(length(perm))
  for(i in 1:length(perm)){
    mse_perm[i] <- max(rowSums((estOm[perm[[i]],] - trOm)^2))/M
    accA_perm[i] <- mean(estA[, perm[[i]]] == trA)
  }
  i_best <- which.max(accA_perm)
  estOm_K[[k]] <- estOm[perm[[i_best]],]
  estA_K[[k]] <- estA[, perm[[i_best]]]
}

estOm_q <- apply(simplify2array(estOm_K), 1:2, quantile, prob = c(0.05, 0.95))
estOm_q2 <- apply(simplify2array(estOm_K), 1:2, quantile, prob = c(0.25, 0.75))

estA_acc <- apply(simplify2array(estA_K), 1:2, function(x) mean(x == median(x)))




pdf(file = paste0("../results/plots/cancer_example_error_sim_m_", m, "al_length_", length(al),"_mse_",round(mse,3)*1000, ".pdf"), width = 15, height = 9) 
layout(rbind(1,2), heights = c(2,1))
eps <- 0.15
col_clone <- c("black", "red", "blue")
pch_clones <- c(19, 17, 15)

plot(trOm[1,], axes = F, xlab = "location", ylab = "frequency", col = col_clone[1], pch = pch_clones[1], ylim = c(0,0.85))
arrows(x0=1:M, y0 = estOm_q[1,1,], x1=1:M, y1=estOm_q[2,1,], code=3, angle=90, length=0, col = col_clone[1])
arrows(x0=1:M, y0 = estOm_q2[1,1,], x1=1:M, y1=estOm_q2[2,1,], code=3, angle=90, length=0.1, col = col_clone[1])

for(i in 2:nrow(ans$estOm)){
  x <- 1:M + (i - 1)*eps
  lines(x, trOm[i,], col = col_clone[i], pch = pch_clones[i], type = "p")
  arrows(x0=x, y0 = estOm_q[1,i,], x1=x, y1=estOm_q[2,i,], code=3, angle=90, length=0, col = col_clone[i])
  arrows(x0=x, y0 = estOm_q2[1,i,], x1=x, y1=estOm_q2[2,i,], code=3, angle=90, length=0.1, col = col_clone[i])
  
}
axis(1, at = 1:ncol(ans$estOm), labels = colnames(ans$estOm))
axis(2, at = seq(0,1,0.1))
box()



if(m == 2){
  at_clones <- c(0,1)
}
if(m == 3){
  at_clones <- c(0,0.5,1)
}

at_location <- seq(0,16/17,1/17)

image(at_location, at_clones, ans$estA[1:17,], col = col_fill, axes = F, xlab = "SNP", ylab = "clone")
text(x = rep(at_location,m), y = rep(at_clones, each = 17), round(estA_acc[1:17,],2), col = "red", cex = 1)
axis(1, at = at_location, labels = names_somVar)
axis(2, at = at_clones, labels = paste("C", 1:ncol(ans$estA)))
if(m == 3){
  axis(2, at = at_clones[2], labels = "C2")
}

box()
dev.off()


