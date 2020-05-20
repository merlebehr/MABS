library(tidyverse)
library(ggmap)
source("functions.R")

al <- 0:2
M <- 3
seed <- 3


m <- 2
omRange <- seq(0,1,0.005)
AsbOm <- matrix(nrow = length(omRange), ncol = length(omRange))

if(M > 2){
  set.seed(seed)
  omRest <- rOmega(m, M-2)
}

for(i in 1:length(omRange)){
  for(j in 1:length(omRange)){
    omega <- cbind(c(omRange[i], 1 - omRange[i]), c(omRange[j], 1 - omRange[j]))
    if(M > 2){
      omega <- cbind(omega, omRest)
    }
    AsbOm[i,j] <- asbOm(omega, al)
  }
}

data <- data.frame(asb = as.numeric(AsbOm), omega11 = rep(omRange, length(omRange)), omega21 = rep(omRange, each = length(omRange)))
ggplot(data, aes(omega11, omega21)) +
  geom_raster(aes(fill = asb)) +
  scale_fill_gradientn(colours=c("white","black"), limits = c(0,1))

ggsave( paste0("../results/plots/asb_plot_al_", max(al),"_M_",M,"_seed_",seed,".png"))


