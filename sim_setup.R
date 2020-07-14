
resultsA = NULL

for(i in 1:reps){
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
x2 <- g[,1:l]%*%effs_g2 + effs_c2*u2a + 0.5*x1 + rnorm(n,0,1) 

x1b <- gb[,1:l]%*%effs_g + effs_c1*ub + rnorm(n,0,1)
x2b <- gb[,1:l]%*%effs_g2 + effs_c2*u2b + 0.5*x1b + rnorm(n,0,1) 

#outcome
y <- 0.2*x1b + 0.3*x2b + gb[,(l+1):(l+lo)]%*%effs_out + effs_c1*ub + effs_c2*u2b

res <- MRest()
resa <- data.frame("a", res)
colnames(resa)[1] <- ("sim")


##b - SNPs have same effect on X1 and X2
x1 <- g[,1:l]%*%effs_g + effs_c1*ua + rnorm(n,0,1)
x2 <- g[,1:l]%*%effs_g + effs_c2*u2a + 0.5*x1 + rnorm(n,0,1) 


#exposures
x1b <- gb[,1:l]%*%effs_g + effs_c1*ub + rnorm(n,0,1)
x2b <- gb[,1:l]%*%effs_g + effs_c2*u2b + 0.5*x1b + rnorm(n,0,1) 

#outcome
y <- 0.2*x1b + 0.3*x2b + gb[,(l+1):(l+lo)]%*%effs_out + effs_c1*ub + effs_c2*u2b + rnorm(n,0,1)

res <- MRest()
resb <- data.frame("b", res)
colnames(resb)[1] <- ("sim")

#c. SNPs only affect x2 through x1 

#exposures
x1 <- g[,1:l]%*%effs_g + effs_c1*ua + rnorm(n,0,1)
x2 <- effs_c2*u2a + 0.5*x1 + rnorm(n,0,1) 

x1b <- gb[,1:l]%*%effs_g + effs_c1*ub + rnorm(n,0,1)
x2b <- effs_c2*u2b + 0.5*x1b + rnorm(n,0,1) 


#outcome
y <- 0.2*x1b + 0.3*x2b + gb[,(l+1):(l+lo)]%*%effs_out + effs_c1*ub + effs_c2*u2b + rnorm(n,0,1)

res <- MRest()
resc <- data.frame("c", res)
colnames(resc)[1] <- ("sim")



#d. model where y -> x2

#exposures
x1 <- g[,1:l]%*%effs_g + effs_c1*ua + rnorm(n,0,1)
#outcome 
ya <- 0.2*x1 +g[,(l+1):(l+lo)]%*%effs_out + effs_c1*ua + effs_c2*u2a + rnorm(n,0,1)
x2 <- g[,1:l]%*%effs_g2 + effs_c2*u2a + 0.5*x1 + 0.5*ya + rnorm(n,0,1) 


#exposures
x1b <- gb[,1:l]%*%effs_g + effs_c1*ub + rnorm(n,0,1)
#outcome 
y <- 0.2*x1b +gb[,(l+1):(l+lo)]%*%effs_out + effs_c1*ub + effs_c2*u2b + rnorm(n,0,1)
x2b <- gb[,1:l]%*%effs_g2 + effs_c2*u2b + 0.5*x1b + 0.5*y + rnorm(n,0,1) 


#res <- MRest(TRUE)      # nb this all of the SNPs for the outcome as well and will select them for estimation if they are significantly associated with x2
res <- MRest()
resd <- data.frame("d", res)
colnames(resd)[1] <- ("sim")


#e. change between x1 and x2 -> y

x1 <- g[,1:l]%*%effs_g + effs_c1*ua + rnorm(n,0,1)
x2 <- g[,1:l]%*%effs_g2 + effs_c2*u2a + 0.5*x1 + rnorm(n,0,1) 

ch <- x2 - x1

#outcome 
x1b <- gb[,1:l]%*%effs_g + effs_c1*ub + rnorm(n,0,1)
x2b <- gb[,1:l]%*%effs_g2 + effs_c2*u2b + 0.5*x1b + rnorm(n,0,1) 

chb <- x2b - x1b

y <- 0.2*x1b + -0.4*chb + gb[,(l+1):(l+lo)]%*%effs_out + effs_c1*ub + effs_c2*u2b + rnorm(n,0,1)


res <- MRest()
rese <- data.frame("e", res)
colnames(rese)[1] <- ("sim")

##f. same as e but only an effect of the change and no effect of level. 

x1 <- g[,1:l]%*%effs_g + effs_c1*ua + rnorm(n,0,1)
x2 <- g[,1:l]%*%effs_g2 + effs_c2*u2a + 0.5*x1 + rnorm(n,0,1) 

ch <- x2 - x1

#outcome 
x1b <- gb[,1:l]%*%effs_g + effs_c1*ub + rnorm(n,0,1)
x2b <- gb[,1:l]%*%effs_g2 + effs_c2*u2b + 0.5*x1b + rnorm(n,0,1) 

chb <- x2b - x1b

y <- -0.4*chb + gb[,(l+1):(l+lo)]%*%effs_out + effs_c1*ub + effs_c2*u2b + rnorm(n,0,1)


res <- MRest()
resf <- data.frame("f", res)
colnames(resf)[1] <- ("sim")


##g. setup where most of the snps have the same effect but a few have a different effect
l_q = round(3*(l/4))
effs_g <- rnorm(l,0,sqrt(0.1/l))  
effs_g2 <- c(effs_g[1:l_q],
             rnorm(l-l_q,0,sqrt(0.1/l)))

x1 <- g[,1:l]%*%effs_g + effs_c1*ua + rnorm(n,0,1)
x2 <- g[,1:l]%*%effs_g2 + effs_c2*u2a + 0.5*x1 + rnorm(n,0,1) 

x1b <- gb[,1:l]%*%effs_g + effs_c1*ub + rnorm(n,0,1)
x2b <- gb[,1:l]%*%effs_g2 + effs_c2*u2b + 0.5*x1b + rnorm(n,0,1) 

#outcome
y <- 0.2*x1b + 0.3*x2b + gb[,(l+1):(l+lo)]%*%effs_out + effs_c1*ub + effs_c2*u2b

res <- MRest()
resg <- data.frame("g", res)
colnames(resg)[1] <- ("sim")

resultsA <- rbind(resultsA,resa, resb, resc, resd, rese, resf, resg)

}
