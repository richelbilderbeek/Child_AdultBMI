
resultsA = NULL

for(i in 1:reps){

##i - SNPs have same effect on X1 and X2; But different confounders

#variables
g <- make_geno(n,l,0.5)
u = rnorm(n,0,1)
u2 = rnorm(n,0,1)

#effects
effs_g <- (choose_effects(l, 0.05))
effs_c1 <- (choose_effects(1, 0.25))
effs_c2 <- (choose_effects(1, 0.25))
#effs_12 <- choose_effects(1, 0.5)

#exposures
x1 <- make_phen(c(effs_g, effs_c1), cbind(g, u))
x2 <- make_phen(c(effs_g, effs_c2, 0.5), cbind(g, u2, x1))

#outcome
y <- make_phen(c(0.1, 0.3, effs_c1, effs_c2), cbind(x1, x2, u, u2))

res <- MRest()
res1 <- data.frame("A.i", res)
colnames(res1)[1] <- ("sim")

#ii. SNPs only affect x2 through x1, different confounders, both have an effect on the outcome. 

#variables
g <- make_geno(n,l,0.5)
u = rnorm(n,0,1)
u2 = rnorm(n,0,1)

#effects
effs_g <- (choose_effects(l, 0.05))
effs_c1 <- (choose_effects(1, 0.25))
effs_c2 <- (choose_effects(1, 0.25))
#effs_12 <- choose_effects(1, 0.5)

#exposures
x1 <- make_phen(c(effs_g, effs_c1), cbind(g, u))
x2 <- make_phen(c(effs_c2, 0.5), cbind(u2, x1))

#outcome
y <- make_phen(c(0.1, 0.3, effs_c1, effs_c2), cbind(x1, x2, u, u2))

res <- MRest()
res2 <- data.frame("A.ii", res)
colnames(res2)[1] <- ("sim")

##iii - SNPs have different effects on X1 and X2; Confounder has the same effect

#variables
g <- make_geno(n,l,0.5)
u = rnorm(n,0,1)


#effects
effs_g <- (choose_effects(l, 0.05))
effs_g2 <- (choose_effects(l, 0.05))
effs_c1 <- (choose_effects(1, 0.25))
#effs_12 <- choose_effects(1, 0.5)

#exposures
x1 <- make_phen(c(effs_g, effs_c1), cbind(g, u))
x2 <- make_phen(c(effs_g2, effs_c1, 0.5), cbind(g, u, x1))

#outcome
y <- make_phen(c(0.1, 0.3, effs_c1), cbind(x1, x2, u))


res <- MRest()
res3 <- data.frame("A.iii", res)
colnames(res3)[1] <- ("sim")

resultsA <- rbind(resultsA,res1, res2, res3)

}
