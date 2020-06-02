
resultsC = NULL

for(i in 1:reps){
  
  ##SNPs have different effects on X1 and X2; Confounder has the same effect

  #effects
  effs_g <- (choose_effects(l, 0.05))
  effs_g2 <- 0.25*effs_g + 0.75*(choose_effects(l, 0.05))
  effs_out <- (choose_effects(lo, 0.05))
  effs_c1 <- (choose_effects(1, 0.15))
  effs_c2 <- (choose_effects(1, 0.15))
  
  
  #variables
  g <- make_geno(n,(l+lo),0.5)
  u = rnorm(n,0,1)
  u2 <- rnorm(n,0,1)
  
  gb <- make_geno(n,(l+lo),0.5)
  ub = rnorm(n,0,1)
  u2b <- rnorm(n,0,1)
  
  #i. model where y -> x2
  
  #exposures
  x1 <- make_phen(c(effs_g, effs_c1), cbind(g[,1:l], u))
  
  #outcome 
  ya <- make_phen(c(0.3, effs_out, effs_c1, effs_c2), cbind(x1, g[,(l+1):(l+lo)], u, u2))
  x2 <- make_phen(c(effs_g2, effs_c2, 0.5, 0.25), cbind(g[,1:l], u2, x1, ya))
  
  #exposures
  x1b <- make_phen(c(effs_g, effs_c1), cbind(gb[,1:l], ub))
  
  #outcome 
  y <- make_phen(c(0.3, effs_out, effs_c1, effs_c2), cbind(x1b, gb[,(l+1):(l+lo)], ub, u2b))
  x2b <- make_phen(c(effs_g2, effs_c2, 0.5, 0.25), cbind(gb[,1:l], u2b, x1b, y))
  
  res <- MRest()
  res1 <- data.frame("C.i", res)
  colnames(res1)[1] <- ("sim")
  
  #ii. same model but x2 -> y and x1 !-> y
  
  #exposures
  x1 <- make_phen(c(effs_g, effs_c1), cbind(g[,1:l], u))
  x2 <- make_phen(c(effs_g2, effs_c2, 0.5), cbind(g[,1:l], u2, x1))
  
  #outcome
  x1b <- make_phen(c(effs_g, effs_c1), cbind(gb[,1:l], ub))
  x2b <- make_phen(c(effs_g2, effs_c2, 0.5), cbind(gb[,1:l], u2b, x1b))
  y <- make_phen(c(0, 0.3, effs_out, effs_c1, effs_c2), cbind(x1b, x2b, gb[,(l+1):(l+lo)], ub, u2b))
  
  
  res <- MRest()
  res2 <- data.frame("C.ii", res)
  colnames(res2)[1] <- ("sim")
  
  #iii. same model but change between x1 and x2 -> y
  
  x1 <- make_phen(c(effs_g, effs_c1), cbind(g[,1:l], u))
  x2 <- make_phen(c(effs_g2, effs_c2, 0.5), cbind(g[,1:l], u2, x1))
  
  ch <- x2 - x1
  
  #outcome 
  x1b <- make_phen(c(effs_g, effs_c1), cbind(gb[,1:l], ub))
  x2b <- make_phen(c(effs_g2, effs_c2, 0.5), cbind(gb[,1:l], u2b, x1b))
  
  chb <- x2b - x1b
  
  y <- make_phen(c(0.2, -0.4, effs_out, effs_c1), cbind(x1b, chb, gb[,(l+1):(l+lo)], ub))
  
  
  res <- MRest()
  res3 <- data.frame("C.iii", res)
  colnames(res3)[1] <- ("sim")

  ##same but only an effect of the change and no effect of level. 
  
  x1 <- make_phen(c(effs_g, effs_c1), cbind(g[,1:l], u))
  x2 <- make_phen(c(effs_g2, effs_c2, 0.5), cbind(g[,1:l], u2, x1))
  
  ch <- x2 - x1
  
  #outcome 
  x1b <- make_phen(c(effs_g, effs_c1), cbind(gb[,1:l], ub))
  x2b <- make_phen(c(effs_g2, effs_c2, 0.5), cbind(gb[,1:l], u2b, x1b))
  
  chb <- x2b - x1b
  
  y <- make_phen(c(0.0, -0.4, effs_out, effs_c1), cbind(x1b, chb, gb[,(l+1):(l+lo)], ub))
  
  
  res <- MRest()
  res3b <- data.frame("C.iiib", res)
  colnames(res3b)[1] <- ("sim")
  
  #iv. Y is observed binary but actual level affects x2 Y -> x2
  
  x1 <- make_phen(c(effs_g, effs_c1), cbind(g[,1:l], u))
  y_conta <- make_phen(c(0.3, effs_out, effs_c1, effs_c2), cbind(x1, g[,(l+1):(l+lo)], u, u2))
  x2 <- make_phen(c(effs_g2, effs_c2, 0.5, 0.5), cbind(g[,1:l], u2, x1, y_conta))
  
 
  #outcome 
  x1b <- make_phen(c(effs_g, effs_c1), cbind(g[,1:l], u))
  y_contb <- make_phen(c(0.3, effs_out, effs_c1, effs_c2), cbind(x1b, gb[,(l+1):(l+lo)], ub, u2b))
  y <- as.numeric(y_contb>quantile(y_contb, c(0.75)))
  x2b <- make_phen(c(effs_g2, effs_c2, 0.5, 0.5), cbind(gb[,1:l], u2b, x1b, y_contb))
  
  res <- MRest()
  res4 <- data.frame("C.iv", res)
  colnames(res4)[1] <- ("sim")
  
  resultsC <- rbind(resultsC,res1,res2,res3,res3b,res4)
  
}
