library(simulateGP)

source('MRest_weighting.R')
reps = 100

n = 300000 #number of individuals (300,000)
l = 150    #number of SNPs (for exposures - total)
lo=10

resultsA = NULL

params <- expand.grid( proportion = c(0.20, 0.50,1), 
                       type = c("limited", "weight"))

params<-params[-c(6),]
params<-params[-c(5),]

for(i in 1:reps){
  
  for(j in 1:nrow(params)){
  
  
  #effects
  effs_g <- rnorm(l,0,sqrt(0.1/l))  
  effs_g2 <- 0.25*effs_g + 0.75*rnorm(l,0,sqrt(0.1/l))
  effs_out <- rnorm(lo,0,sqrt(0.25/l))
  effs_c1 <- 0.5
  effs_c2 <- 0.5
  
  #variables
  g <- make_geno(n,(l+lo),0.5)
  ua = rnorm(n,0,1)
  u2a = rnorm(n,0,1)
  
  gb <- make_geno(n,(l+lo),0.5)
  ub = rnorm(n,0,1)
  u2b = rnorm(n,0,1)
  
  
  ##a - SNPs have different effects on X1 and X2;
  # This is the standard reference setup where everything should work.
  
  x1 <- g[,1:l]%*%effs_g + effs_c1*ua + rnorm(n,0,1)
  x2 <- g[,1:l]%*%effs_g2 + effs_c2*u2a + 0.5*x1 + 0.5*rnorm(n,0,1) 
  
  x1b <- gb[,1:l]%*%effs_g + effs_c1*ub + rnorm(n,0,1)
  x2b <- gb[,1:l]%*%effs_g2 + effs_c2*u2b + 0.5*x1b + 0.5*rnorm(n,0,1) 
  
  #outcome
  y <- 0.2*x1b + 0.3*x2b + gb[,(l+1):(l+lo)]%*%effs_out + effs_c1*ub + effs_c2*u2b
  
  res <-  MRweight(params$type[j],params$proportion[j])
  resa <- data.frame(params$type[j],params$proportion[j], "diff_effects", res)
  colnames(resa)[1] <- ("type")
  colnames(resa)[2] <- ("proportion")
  colnames(resa)[3] <- ("datagen")
  
  
  ##b - SNPs have same effect on X1 and X2
  x1 <- g[,1:l]%*%effs_g + effs_c1*ua + rnorm(n,0,1)
  x2 <- g[,1:l]%*%effs_g + effs_c2*u2a + 0.5*x1 + 0.5*rnorm(n,0,1) 
  #exposures
  x1b <- gb[,1:l]%*%effs_g + effs_c1*ub + rnorm(n,0,1)
  x2b <- gb[,1:l]%*%effs_g + effs_c2*u2b + 0.5*x1b + 0.5*rnorm(n,0,1) 
  #outcome
  y <- 0.2*x1b + 0.3*x2b + gb[,(l+1):(l+lo)]%*%effs_out + effs_c1*ub + effs_c2*u2b + rnorm(n,0,1)
  
  res <- MRweight(params$type[j],params$proportion[j])
  resb <- data.frame(params$type[j],params$proportion[j], "same_effects", res)
  colnames(resb)[1] <- ("type")
  colnames(resb)[2] <- ("proportion")
  colnames(resb)[3] <- ("datagen")
  
  resultsA <- rbind(resultsA,resa, resb)
  }
}