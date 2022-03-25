

rm(list = ls(all=TRUE))

set.seed(4)


library(MASS)
library(ppcor)
source('MRest_x1x3.R')

make_geno <- function(nid, nsnp, af)
{
  return(matrix(rbinom(nid * nsnp, 2, af), nid, nsnp))
}


reps = 1000
n = 150000 #number of individuals (should be 150000)
l = 150    #number of SNPs (for exposures - total)
lo = 10 #number of SNPs for outcome

beta1 = 0
beta2 = 0.3
beta3 = 0

results = NULL

#define effects outside the repetitions so they are consistent across the simulations
effs_g <- rnorm(l,0,sqrt(0.15/l))  
effs_g2 <- 0.05*effs_g + rnorm((l),0,sqrt(0.15/l))
effs_g3 <- 0*effs_g + 0.5*effs_g2 + 0.5*rnorm(l,0,sqrt(0.15/l))

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
  y <- beta1*x1b + beta2*x2b + beta3*x3b + gb[,(l+1):(l+lo)]%*%effs_out + 0.3*effs_c1*ub + 0.3*effs_c2*u2b + 0.3*effs_c2*u3b+ rnorm(n,0,1)
  
  res <- MRest_x1x3()
  res_x1x3b <- data.frame("x1x3b", res)
  colnames(res_x1x3b)[1] <- ("sim")

  res_x1x3b$beta1_u <- beta1 + 0.1*(beta2+0.1*beta3) + cor(effs_g, effs_g2)*beta2 + cor(effs_g, effs_g3)*beta3
  res_x1x3b$beta3_u <- beta3 + cor(effs_g, effs_g3)*(beta1+0.1*beta2) + cor(effs_g2, effs_g3)*beta2
  
  res_x1x3b$beta1_m <- beta1 + (pcor.test(effs_g,effs_g2,effs_g3)$estimate)*beta2
  res_x1x3b$beta3_m <- beta3 + (pcor.test(effs_g2,effs_g3,effs_g)$estimate)*beta2
  
  results <- rbind(results,res_x1x3b)
  
}


save(results, file="sim_output_x1x3b.Rda")

