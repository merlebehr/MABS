library(readxl)
library(tidyverse)
library(Clomial)
source('functions.R')

data <- vector('list', 14)
for(i in 1:length(data)){
  data[[i]] <- read_excel("../data/pcbi_data2.xls", sheet = i)
}

ind_target <- which(!is.na(data[[1]][,4]))
names_target <- unlist(data[[1]][ind_target,4])


ind_major <- data[[14]] %>% select(A,T,C,G) %>% apply(1, which.max) 

af_data <- sapply(c(14, 2:13), function(i) sapply(1:length(ind_major), function(j) 
  max(0, 1 - as.numeric(unlist(data[[i]][j, ind_major[j] + 1] / data[[i]]$total[j])))))

total_data <- sapply(2:14, function(i) sapply(1:length(ind_major), function(j)  unlist(data[[i]]$total[j])))
minor_data <- sapply(2:14, function(i) sapply(1:length(ind_major), function(j)  unlist(data[[i]]$total[j] - data[[i]][j, ind_major[j] + 1])))



rownames(af_data) <- unlist(data[[1]][,1])
colnames(af_data) <- c("P11", "P12", "P13", "P14", "P21", "P22", "P23", "P24", "P31", "P32", "M11", "M12", "N11")


ind_noNA <- which(apply(af_data, 1, function(x) !any(is.na(x)) & x[13] == 0))
ind_noNA <- setdiff(ind_noNA, ind_target)
n_add <- 20
ind_add <- sample(ind_noNA, n_add)
#ind_add <- c(ind_add, sample(ind_noNA, 30, prob = rowSums(af_data[ind_noNA,])))



y <- af_data[ind_target,]
rownames(y) <- names_target

y <- rbind(y, af_data[ind_add, ])

y <- rbind(y_test, af_data[ind_add, ])



m <- 2
al <- c(0,0.5,1)
#al <- 0:1
Am <- as.matrix(expand.grid(rep(list(al), m)))

estOm <- estOmega(y = y, m = m, al = al)
valEst <- Am %*% estOm
indEstA <- sapply(1:nrow(y), function(i) which.min(colSums( (y[i,] - t(valEst) )^2 ) ))
estA <- Am[indEstA, ]
rownames(estA) <- rownames(y)
colnames(estOm) <- colnames(y)

round(y,2)
round(estOm,2) 
estA

sqrt(mean((y - estA %*% estOm)^2))

ans <- lloyd(y, m, al = al, omega = estOm)
round(ans$estOm, 2) 
ans$estA

sqrt(mean((y - (ans$estA %*% ans$estOm))^2))
sqrt(mean((y[1:17,] - (ans$estA %*% ans$estOm)[1:17,])^2))


# plot(y[3,] - y[11,], ylim = c(0,0.5))
# lines(colSums(ans$estOm)/2, type = "p", col = "red")



set.seed(1)
data(breastCancer)
Dc <- breastCancer$Dc
Dt <- breastCancer$Dt
y_test <- Dt/Dc
ClomialResult <-Clomial(Dc=Dc,Dt=Dt,maxIt=20,C=3,binomTryNum=2)
chosen <- choose.best(models=ClomialResult$models)
M1 <- chosen$bestModel
print("Genotypes:")
round(M1$Mu)
estA <- M1$Mu
estA
print("Clone frequencies:")
round(M1$P,2)
estOm <- M1$P

sqrt(mean((y_test - estA %*% estOm)^2))
sqrt(mean((y[1:17,] - estA %*% estOm)^2))

