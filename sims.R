

results = NULL

for(i in 1:reps){
  #effects
  effs_g <- rnorm(l,0,sqrt(0.15/l))  
  effs_g2 <- 0.25*effs_g + 0.75*rnorm(l,0,sqrt(0.15/l))
  effs_out <- rnorm(lo,0,sqrt(0.3/l))
  effs_c1 <- 0.5
  effs_c2 <- 0.5
  
  #variables
  g <- make_geno(n,(l+lo),0.5)
  Sigma = matrix(c(1, 0.8,0.8,0.8,1,0.8,0.8,0.8,1), nrow = 3)
  u <- mvrnorm(n, c(0,0,0), Sigma)
  ua = u[,1]
  u2a = u[,2]
  u3a = u[,3]
  
  gb <- make_geno(n,(l+lo),0.5)
  u <- mvrnorm(n, c(0,0,0), Sigma)
  ub = u[,1]
  u2b= u[,2]
  u3b = u[,3]
  
  ##a - x1 and x2 are associated with different latent periods
  # This is the standard reference setup where everything should work.
  
  L1 <- g[,1:l]%*%effs_g 
  L2 <- g[,1:l]%*%effs_g2
  
  x1 <- L1 + effs_c1*ua + rnorm(n,0,1)
  x2 <- L2 + effs_c2*u2a + 0.1*x1 + 0.9*rnorm(n,0,1) 
  
  L1b <- gb[,1:l]%*%effs_g 
  x1b <- L1b + effs_c1*ub + rnorm(n,0,1)
  L2b <- gb[,1:l]%*%effs_g2 
  x2b <- L2b + effs_c2*u2b + 0.1*x1b + 0.9*rnorm(n,0,1) 
  
  #outcome
  y <- 0.2*x1b + 0.3*x2b + gb[,(l+1):(l+lo)]%*%effs_out + 0.5*effs_c1*ub + 0.5*effs_c2*u2b
  
  res <- MRest()
  resa <- data.frame("a", res)
  colnames(resa)[1] <- ("sim")
  
  
  ##b - x1 and x2 are associated with the same latent period
  L1 <- g[,1:l]%*%effs_g 
  L2 <- g[,1:l]%*%effs_g2
  
  x1 <- L1 + effs_c1*ua + rnorm(n,0,1)
  x2 <- L1 + effs_c2*u2a + 0.1*x1 + 0.9*rnorm(n,0,1) 
  
  L1b <- gb[,1:l]%*%effs_g 
  x1b <- L1b + effs_c1*ub + rnorm(n,0,1)
  L2b <- gb[,1:l]%*%effs_g2 
  x2b <- L1b + effs_c2*u2b + 0.1*x1b + 0.9*rnorm(n,0,1) 
  
  y <- 0.2*x1b + 0.3*x2b + gb[,(l+1):(l+lo)]%*%effs_out + 0.5*effs_c1*ub + 0.5*effs_c2*u2b
  
  res <- MRest()
  resb <- data.frame("b", res)
  colnames(resb)[1] <- ("sim")
  
  
  #c. model where y -> x2
  
  #exposures
  L1 <- g[,1:l]%*%effs_g 
  L2 <- g[,1:l]%*%effs_g2
  
  x1 <- L1 + effs_c1*ua + rnorm(n,0,1)
  ya <- 0.2*x1 +g[,(l+1):(l+lo)]%*%effs_out + effs_c1*ua + effs_c2*u2a + rnorm(n,0,1)
  x2 <- L2 + effs_c2*u2a + 0.1*x1 + 0.5*ya + 0.5*rnorm(n,0,1) 

  L1b <- gb[,1:l]%*%effs_g 
  x1b <- L1b + effs_c1*ub + rnorm(n,0,1)
  L2b <- gb[,1:l]%*%effs_g2 
 
  y <- 0.2*x1b +gb[,(l+1):(l+lo)]%*%effs_out + effs_c1*ub + effs_c2*u2b + rnorm(n,0,1)
  x2b <- L2b + effs_c2*u2b + 0.1*x1b + 0.5*y + 0.5*rnorm(n,0,1) 
  
  res <- MRest()
  resc <- data.frame("c", res)
  colnames(resc)[1] <- ("sim")
  
  #d. model where there are three periods but only two included.
  
  effs_g3 <- 0.1*effs_g + 0.25*effs_g2 + 0.65*rnorm(l,0,sqrt(0.15/l))
  
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
  resd <- data.frame("d", res)
  colnames(resd)[1] <- ("sim")
  
  #e. same as d but no correlation between the genetic effects on latent period 3 and the others
  
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
  
  
  results <- rbind(results,resa, resb, resc, resd, rese)
  
}
