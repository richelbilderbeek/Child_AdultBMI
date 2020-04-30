
resultsC = NULL

for(i in 1:reps){
  
  ##SNPs have different effects on X1 and X2; Confounder has the same effect
  #Y -> X2 
  
  #variables
  g <- make_geno(n,l,0.5)
  u = rnorm(n,0,1)
  
  
  #effects
  effs_g <- (choose_effects(l, 0.05))
  effs_g2 <- (choose_effects(l, 0.05))
  effs_c1 <- (choose_effects(1, 0.25))
  #effs_12 <- choose_effects(1, 0.5)
  
  #i. model where y -> x2
  
  #exposures
  x1 <- make_phen(c(effs_g, effs_c1), cbind(g, u))
  
  #outcome 
  y <- make_phen(c(0.3, effs_c1), cbind(x1, u))
  x2 <- make_phen(c(effs_g2, effs_c1, 0.5, 0.2), cbind(g, u, x1, y))
  
  res <- MRest()
  res1 <- data.frame("C.i", res)
  colnames(res1)[1] <- ("sim")
  
  #ii. same model but x2 -> y and x1 !-> y
  
  #exposures
  x1 <- make_phen(c(effs_g, effs_c1), cbind(g, u))
  x2 <- make_phen(c(effs_g2, effs_c1, 0.5), cbind(g, u, x1))
  
  #outcome 
  y <- make_phen(c(0, 0.2, effs_c1), cbind(x1, x2, u))
  
  
  res <- MRest()
  res2 <- data.frame("C.ii", res)
  colnames(res2)[1] <- ("sim")
  
  #iii. same model but change between x1 and x2 -> y
  
  x1 <- make_phen(c(effs_g, effs_c1), cbind(g, u))
  x2 <- make_phen(c(effs_g2, effs_c1, 0.5), cbind(g, u, x1))
  
  ch <- x2 - x1
  
  #outcome 
  y <- make_phen(c(0.3, 0.2, effs_c1), cbind(x1, ch, u))
  
  
  res <- MRest()
  res3 <- data.frame("C.iii", res)
  colnames(res3)[1] <- ("sim")

  
  #iv. Y is observed binary but actual level affects x2
  
  x1 <- make_phen(c(effs_g, effs_c1), cbind(g, u))
  
  #outcome 
  y_cont <- make_phen(c(0.3, effs_c1), cbind(x1, u))
  y <- ifelse(y_cont >= mean(y_cont), 1, 0)
  x2 <- make_phen(c(effs_g2, effs_c1, 0.5, 0.2), cbind(g, u, x1, y_cont))
  
  res <- MRest()
  res4 <- data.frame("C.iv", res)
  colnames(res4)[1] <- ("sim")
  
  resultsC <- rbind(resultsC,res1,res2,res3,res4)
  
}
