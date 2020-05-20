
resultsA = NULL

for(i in 1:reps){

##i - SNPs have same effect on X1 and X2; But different confounders

#variables
g <- make_geno(n,(l+lo),0.5)
u = rnorm(n,0,1)
u2 = rnorm(n,0,1)

#effects
effs_g <- (choose_effects(l, 0.05))
effs_g2 <- 0.25*effs_g + 0.75*(choose_effects(l, 0.05))
effs_out <- (choose_effects(lo, 0.05))
effs_c1 <- (choose_effects(1, 0.25))
effs_c2 <- (choose_effects(1, 0.25))

#exposures
x1 <- make_phen(c(effs_g, effs_c1), cbind(g[,1:l], u))
x2 <- make_phen(c(effs_g, effs_c2, 0.5), cbind(g[,1:l], u2, x1))

#outcome
y <- make_phen(c(0.1, 0.3, effs_out, effs_c1, effs_c2), cbind(x1, x2, g[,(l+1):(l+lo)], u, u2))

res <- MRest()
res1 <- data.frame("A.i", res)
colnames(res1)[1] <- ("sim")

#ii. SNPs only affect x2 through x1, different confounders, both have an effect on the outcome. 

#variables
g <- make_geno(n,(l+lo),0.5)
u = rnorm(n,0,1)
u2 = rnorm(n,0,1)

#effects
effs_g <- (choose_effects(l, 0.05))
effs_g2 <- 0.25*effs_g + 0.75*(choose_effects(l, 0.05))
effs_out <- (choose_effects(lo, 0.05))
effs_c1 <- (choose_effects(1, 0.25))
effs_c2 <- (choose_effects(1, 0.25))


#exposures
x1 <- make_phen(c(effs_g, effs_c1), cbind(g[,1:l], u))
x2 <- make_phen(c(effs_c2, 0.5), cbind(u2, x1))

#outcome
y <- make_phen(c(0.1, 0.3, effs_out, effs_c1, effs_c2), cbind(x1, x2, g[,(l+1):(l+lo)], u, u2))

res <- MRest()
res2 <- data.frame("A.ii", res)
colnames(res2)[1] <- ("sim")

##iii - SNPs have different effects on X1 and X2; Confounder has the same effect

#variables
g <- make_geno(n,(l+lo),0.5)
u = rnorm(n,0,1)

#effects
effs_g <- (choose_effects(l, 0.05))
effs_g2 <- 0.25*effs_g + 0.75*(choose_effects(l, 0.05))
effs_out <- (choose_effects(lo, 0.05))
effs_c1 <- (choose_effects(1, 0.25))

#exposures
x1 <- make_phen(c(effs_g, effs_c1), cbind(g[,1:l], u))
x2 <- make_phen(c(effs_g2, effs_c1, 0.5), cbind(g[,1:l], u, x1))

#outcome
y <- make_phen(c(0.1, 0.3, effs_out, effs_c1), cbind(x1, x2, g[,(l+1):(l+lo)], u))


res <- MRest()
res3 <- data.frame("A.iii", res)
colnames(res3)[1] <- ("sim")

##iv - As iii but all binary (underlying cont latent variable)

#variables
g <- make_geno(n,(l+lo),0.5)
u = rnorm(n,0,1)

#effects
effs_g <- (choose_effects(l, 0.05))
effs_g2 <- 0.25*effs_g + 0.75*(choose_effects(l, 0.05))
effs_out <- (choose_effects(lo, 0.05))
effs_c1 <- (choose_effects(1, 0.25))

#exposures
x1_c <- make_phen(c(effs_g, effs_c1), cbind(g[,1:l], u))
x2_c <- make_phen(c(effs_g2, effs_c1, 0.5), cbind(g[,1:l], u, x1))
#outcome
y_c <- make_phen(c(0.1, 0.3, effs_out, effs_c1), cbind(x1, x2, g[,(l+1):(l+lo)], u))

x1 <- as.numeric(x1_c>quantile(x1_c, c(0.75)))
x2 <- as.numeric(x2_c>quantile(x2_c, c(0.75)))
y <- as.numeric(y_c>quantile(y_c, c(0.75)))

#estimation using LPM
res <- MRest()
res4 <- data.frame("A.iv_a", res)
colnames(res4)[1] <- ("sim")

#estimation using logit
res <- MRest_OR()
res5 <- data.frame("A.iv_b", res)
colnames(res5)[1] <- ("sim")

resultsA <- rbind(resultsA,res1, res2, res3, res4, res5)

}
