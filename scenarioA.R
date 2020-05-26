
resultsA = NULL

for(i in 1:reps){
#effects
effs_g <- (choose_effects(l, 0.05))
effs_g2 <- 0.25*effs_g + 0.75*(choose_effects(l, 0.05))
effs_out <- (choose_effects(lo, 0.05))
effs_c1 <- (choose_effects(1, 0.25))
effs_c2 <- (choose_effects(1, 0.25))

#variables
g <- make_geno(n,(l+lo),0.5)
ua = rnorm(n,0,1)
u2a = rnorm(n,0,1)

gb <- make_geno(n,(l+lo),0.5)
ub = rnorm(n,0,1)
u2b = rnorm(n,0,1)


##i - SNPs have same effect on X1 and X2; But different confounders

x1 <- make_phen(c(effs_g, effs_c1), cbind(g[,1:l], ua))
x2 <- make_phen(c(effs_g, effs_c2, 0.5), cbind(g[,1:l], u2a, x1))

#exposures
x1b <- make_phen(c(effs_g, effs_c1), cbind(gb[,1:l], ub))
x2b <- make_phen(c(effs_g, effs_c2, 0.5), cbind(gb[,1:l], u2b, x1b))

#outcome
y <- make_phen(c(0.1, 0.3, effs_out, effs_c1, effs_c2), cbind(x1b, x2b, gb[,(l+1):(l+lo)], ub, u2b))

res <- MRest()
res1 <- data.frame("A.i", res)
colnames(res1)[1] <- ("sim")

#ii. SNPs only affect x2 through x1, different confounders, both have an effect on the outcome. 

#exposures
x1 <- make_phen(c(effs_g, effs_c1), cbind(g[,1:l], ua))
x2 <- make_phen(c(effs_c2, 0.5), cbind(u2a, x1))

x1b <- make_phen(c(effs_g, effs_c1), cbind(gb[,1:l], ub))
x2b <- make_phen(c(effs_c2, 0.5), cbind(u2b, x1b))


#outcome
y <- make_phen(c(0.1, 0.3, effs_out, effs_c1, effs_c2), cbind(x1b, x2b, gb[,(l+1):(l+lo)], ub, u2b))

res <- MRest()
res2 <- data.frame("A.ii", res)
colnames(res2)[1] <- ("sim")

##iii - SNPs have different effects on X1 and X2; Confounder has the same effect

#exposures
x1 <- make_phen(c(effs_g, effs_c1), cbind(g[,1:l], ua))
x2 <- make_phen(c(effs_g2, effs_c1, 0.5), cbind(g[,1:l], ua, x1))

x1b <- make_phen(c(effs_g, effs_c1), cbind(gb[,1:l], ub))
x2b <- make_phen(c(effs_g2, effs_c1, 0.5), cbind(gb[,1:l], ub, x1))

#outcome
y <- make_phen(c(0.1, 0.3, effs_out, effs_c1), cbind(x1b, x2b, gb[,(l+1):(l+lo)], ub))


res <- MRest()
res3 <- data.frame("A.iii", res)
colnames(res3)[1] <- ("sim")

##iv - As iii but all binary (underlying cont latent variable)



#exposures
x1_c <- make_phen(c(effs_g, effs_c1), cbind(g[,1:l], ua))
x2_c <- make_phen(c(effs_g2, effs_c1, 0.5), cbind(g[,1:l], ua, x1_c))

x1b_c <- make_phen(c(effs_g, effs_c1), cbind(gb[,1:l], ub))
x2b_c <- make_phen(c(effs_g2, effs_c1, 0.5), cbind(gb[,1:l], ub, x1b_c))

#outcome
y_c <- make_phen(c(0.1, 0.3, effs_out, effs_c1), cbind(x1b_c, x2b_c, gb[,(l+1):(l+lo)], ub))

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
