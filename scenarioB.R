
resultsB = NULL
 a = 0.2
 b = 0.2
 c = 0.2
 d = 0.5
 e = 0.5

for(i in 1:reps){
  
  #effects
  effs_g1 <- choose_effects(l1, 0.025)
  effs_ga <- (choose_effects(l, 0.025))
  effs_gb <- (choose_effects(l, 0.05))
  effs_gc <- (choose_effects(l, 0.05))
  effs_c1 <- (choose_effects(1, 0.25))
  
  
   ##i - have snps for x1
  
  #g1 from a separate sample 
  g_1 <- make_geno(n,l,0.5)
  g1_1 <- make_geno(n,l1,0.5)
  u_1 = rnorm(n,0,1)
  x1 <- make_phen(c(effs_g1, effs_ga, effs_c1), cbind(g1_1, g_1, u_1))
  
  #variables
  g <- make_geno(n,l,0.5)
  g1 <- make_geno(n,l1,0.5)
  u = rnorm(n,0,1)
 
  #exposures
  x1a <- make_phen(c(effs_g1, effs_ga, effs_c1), cbind(g1, g, u))
  x2 <- make_phen(c(effs_gb, effs_c1, d), cbind(g, u, x1a))
  x3 <- make_phen(c(effs_gc, effs_c1, d), cbind(g, u, x2))
  
  #outcome
  
  
  #variables
  gb <- make_geno(n,l,0.5)
  g1b <- make_geno(n,l1,0.5)
  ub = rnorm(n,0,1)
  
  
  #exposures
  x1b <- make_phen(c(effs_g1, effs_ga, effs_c1), cbind(g1b, gb, ub))
  x2b <- make_phen(c(effs_gb, effs_c1, d), cbind(gb, ub, x1b))
  x3b <- make_phen(c(effs_gc, effs_c1, d), cbind(gb, ub, x2b))
  
  y <- make_phen(c(a, b, c, effs_c1), cbind(x1b, x2b, x3b, ub))
 
  res1 <- MRest_3()
  
  resultsB <- rbind(resultsB,res1)
   
}
  