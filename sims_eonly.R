

rm(list = ls(all=TRUE))

set.seed(4)


library(MASS)
library(ppcor)
source('MRest.R')

make_geno <- function(nid, nsnp, af)
{
  return(matrix(rbinom(nid * nsnp, 2, af), nid, nsnp))
}


reps = 2000
n = 150000 #number of individuals
l = 150    #number of SNPs (for exposures - total)
lo = 10 #number of SNPs for outcome

results = NULL

#define effects outside the repetitions so they are consistent across the simulations
effs_g <- rnorm(l,0,sqrt(0.15/l))  
effs_g2 <- 0.3*effs_g + rnorm((l),0,sqrt(0.15/l))
effs_out <- rnorm(lo,0,sqrt(0.3/l))
effs_c1 <- 0.5
effs_c2 <- 0.5

for(i in 1:reps){
  
  
  #variables
  g <- make_geno(n,(l+lo),0.5)
  Sigma = matrix(c(1, 0.8,0.8,0.8,1,0.8,0.8,0.8,1), nrow = 3)
  u <- mvrnorm(n, c(1,1,1), Sigma)
  ua = u[,1]
  u2a = u[,2]
  u3a = u[,3]
  
  gb <- make_geno(n,(l+lo),0.5)
  u <- mvrnorm(n, c(1,1,1), Sigma)
  ub = u[,1]
  u2b= u[,2]
  u3b = u[,3]
  
  effs_g3 <- rnorm(l,0,sqrt(0.15/l))
  
  L1 <- g[,1:l]%*%effs_g 
  L2 <- g[,1:l]%*%effs_g2
  L3 <- g[,1:l]%*%effs_g3
  
  x1 <- L1 + effs_c1*ua + rnorm(n,0,1)
  x2 <- L2 + effs_c2*u2a + 0.1*x1 + 0.9*rnorm(n,0,1) 
  x3 <- L3 + effs_c2*u3a + 0.1*x2 + 0.9*rnorm(n,0,1) 
  
  L1b <- gb[,1:l]%*%effs_g 
  L2b <- gb[,1:l]%*%effs_g2 
  L3b <- gb[,1:l]%*%effs_g3
  
  x1b <- L1b + effs_c1*ub + rnorm(n,0,1)
  x2b <- L2b + effs_c2*u2b + 0.1*x1b + 0.9*rnorm(n,0,1) 
  x3b <- L3b + effs_c2*u3b + 0.1*x2b + 0.9*rnorm(n,0,1) 
  
  #outcome
  y <- 0.2*x1b + 0.3*x2b + 0.2*x3b + gb[,(l+1):(l+lo)]%*%effs_out + 0.3*effs_c1*ub + 0.3*effs_c2*u2b + 0.3*effs_c2*u3b
  
  res <- MRest()
  rese <- data.frame("e", res)
  colnames(rese)[1] <- ("sim")
  
  resd$beta1_u <- 0.2 + 0.1*(0.3+0.1*0.2) + cor(effs_g, effs_g2)*0.3 + cor(effs_g, effs_g3)*0.2
  resd$beta2_u <- 0.3 + 0.1*0.2 + cor(effs_g, effs_g2)*0.2 + cor(effs_g2, effs_g3)*0.2
  
  resd$beta1_m <- 0.2 + (pcor.test(effs_g,effs_g3,effs_g2)$estimate)*0.2
  resd$beta2_m <- 0.3 + (pcor.test(effs_g2,effs_g3,effs_g)$estimate)*0.2 +  0.1*0.2 
  
  results <- rbind(results,rese)
  
}


save(results, file="sim_output_e.Rda")

