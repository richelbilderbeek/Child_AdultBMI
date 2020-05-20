
MRest <- function()
{

 mvdat <- make_mvdat(list(x1, x2), y, g)
 mvdat$minp <- apply(mvdat$exposure_pval, 1, min)
 mvdat <- as.data.frame(mvdat)
 dat1 <- subset(mvdat, mvdat$exposure_pval.x1 < 5e-08)
 unix1 <- summary(lm(dat1$outcome_beta ~ -1 + dat1$exposure_beta.x1, weights = 1/(dat1$outcome_se^2)))
 
 F_1 <- (dat1$exposure_beta.x1/dat1$exposure_se.x1)^2
 
 dat2 <- subset(mvdat, mvdat$exposure_pval.x2 < 5e-08)
 unix2 <- summary(lm(dat2$outcome_beta ~ -1 + dat2$exposure_beta.x2, weights = 1/(dat2$outcome_se^2)))
 F_2 <- (dat2$exposure_beta.x2/dat2$exposure_se.x2)^2

  mv <- subset(mvdat, mvdat$minp < 5e-08)
  mvmr <- summary(lm(mv$outcome_beta ~ -1 + mv$exposure_beta.x1 + mv$exposure_beta.x2, weights = 1/(mv$outcome_se^2)))
  
  rho = cor(x1,x2)
  sig12 = rho*mv$exposure_se.x1*mv$exposure_se.x2
  
  delta1 <- lm(mv$exposure_beta.x1~ -1 + mv$exposure_beta.x2)$coefficients[1]
  delta2 <- lm(mv$exposure_beta.x2~ -1 + mv$exposure_beta.x1)$coefficients[1]
  
  Qind_1 <- ((mv$exposure_beta.x1 - delta1*mv$exposure_beta.x2)^2)/(mv$exposure_se.x1^2 + (delta1^2)*mv$exposure_se.x2^2 - 2*delta1*sig12)
  Qind_2 <- ((mv$exposure_beta.x2 - delta2*mv$exposure_beta.x1)^2)/(mv$exposure_se.x2^2 + (delta2^2)*mv$exposure_se.x1^2 - 2*delta2*sig12)
  
  snps_outx1 <- as.numeric(sum((mvdat$exposure_pval.x1[(l+1):(l+lo)]<5e-08)))
  snps_outx2 <- as.numeric(sum((mvdat$exposure_pval.x2[(l+1):(l+lo)]<5e-08)))
  
out <- data.frame(rho, unix1$coefficients[1,1], nrow(dat1), mean(F_1), unix2$coefficients[1,1], nrow(dat2), mean(F_2), 
                      mvmr$coefficients[1,1], sum(Qind_1)/(nrow(mv)-1), mvmr$coefficients[2,1], nrow(mv), sum(Qind_2)/(nrow(mv)-1), 
                      snps_outx1, snps_outx2)

colnames(out) <- c("rho", "uni_x1_b", "uni_x1_nsnp","F_x1" , "uni_x2_b", "uni_x2_nsnp", "F_x2", 
                       "mv_x1_b", "CF_x1", "mv_x2_b", "mv_nsnp", "CF_x2", "snps_outx1", "snps_outx2")
  
return(out)

}

MRest_OR <- function()
{
  mvdat <- data.frame()
  
  
  # mvdat <- make_mvdat(list(x1, x2), y, g)
  for(s in 1:(l+lo)){
    x1_res <- summary(glm(x1~g[,s], family = binomial))
    x2_res <- summary(glm(x2~g[,s], family = binomial))
    y_res <- summary(glm(y~g[,s], family = binomial))
    mvdat[s,1] <- exp(x1_res$coefficient[2,1])
    mvdat[s,2] <- exp(x2_res$coefficient[2,1])
    mvdat[s,3] <- exp(y_res$coefficient[2,1])
    
    mvdat[s,4] <- (x1_res$coefficient[2,2])
    mvdat[s,5] <- (x2_res$coefficient[2,2])
    mvdat[s,6] <- (y_res$coefficient[2,2])
    
    mvdat[s,7] <- (x1_res$coefficient[2,4])
    mvdat[s,8] <- (x2_res$coefficient[2,4])
    mvdat[s,9] <- (y_res$coefficient[2,4])
    
  }
  
  colnames(mvdat) <- c("exposure_beta.x1", "exposure_beta.x2", "outcome_beta","exposure_se.x1", "exposure_se.x2", 
                       "outcome_se","exposure_pval.x1", "exposure_pval.x2", "outcome_pval" )
  
  
  
  mvdat$minp <- apply(mvdat[,c("exposure_pval.x1","exposure_pval.x2")], 1, min)
  mvdat <- as.data.frame(mvdat)
  dat1 <- subset(mvdat, mvdat$exposure_pval.x1 < 5e-08)
  unix1 <- summary(lm(dat1$outcome_beta ~ -1 + dat1$exposure_beta.x1, weights = 1/(dat1$outcome_se^2)))
  
  F_1 <- (log(dat1$exposure_beta.x1)/(dat1$exposure_se.x1))^2
  
  dat2 <- subset(mvdat, mvdat$exposure_pval.x2 < 5e-08)
  unix2 <- summary(lm(dat2$outcome_beta ~ -1 + dat2$exposure_beta.x2, weights = 1/(dat2$outcome_se^2)))
  F_2 <- (log(dat2$exposure_beta.x2)/(dat2$exposure_se.x2))^2
  
  mv <- subset(mvdat, mvdat$minp < 5e-08)
  mvmr <- summary(lm(mv$outcome_beta ~ -1 + mv$exposure_beta.x1 + mv$exposure_beta.x2, weights = 1/(mv$outcome_se^2)))
  
  rho = cor(x1,x2)
  sig12 = rho*mv$exposure_se.x1*mv$exposure_se.x2
  
  delta1 <- lm(mv$exposure_beta.x1~ -1 + mv$exposure_beta.x2)$coefficients[1]
  delta2 <- lm(mv$exposure_beta.x2~ -1 + mv$exposure_beta.x1)$coefficients[1]
  
  Qind_1 <- ((mv$exposure_beta.x1 - delta1*mv$exposure_beta.x2)^2)/(mv$exposure_se.x1^2 + (delta1^2)*mv$exposure_se.x2^2 - 2*delta1*sig12)
  Qind_2 <- ((mv$exposure_beta.x2 - delta2*mv$exposure_beta.x1)^2)/(mv$exposure_se.x2^2 + (delta2^2)*mv$exposure_se.x1^2 - 2*delta2*sig12)
  
  snps_outx1 <- as.numeric(sum((mvdat$exposure_pval.x1[(l+1):(l+lo)]<5e-08)))
  snps_outx2 <- as.numeric(sum((mvdat$exposure_pval.x2[(l+1):(l+lo)]<5e-08)))
  
  out <- data.frame(rho, unix1$coefficients[1,1], nrow(dat1), mean(F_1), unix2$coefficients[1,1], nrow(dat2), mean(F_2), 
                    mvmr$coefficients[1,1], sum(Qind_1)/(nrow(mv)-1), mvmr$coefficients[2,1], nrow(mv), sum(Qind_2)/(nrow(mv)-1),
                    snps_outx1, snps_outx2)
  
  colnames(out) <- c("rho", "uni_x1_b", "uni_x1_nsnp","F_x1" , "uni_x2_b", "uni_x2_nsnp", "F_x2", 
                     "mv_x1_b", "CF_x1", "mv_x2_b", "mv_nsnp", "CF_x2",  "snps_outx1", "snps_outx2")
  
  return(out)
  
}


MRest_3 <- function()
{
  
  #x1 only SNPs
  dat <- get_effs(x1, y, g1)
  dat <- subset(dat, dat$pval.exposure < 5e-08)
  n_g1 <- nrow(dat)
  unix1_y <- summary(lm(dat$beta.outcome ~ -1 + dat$beta.exposure, weights = 1/(dat$se.outcome^2)))
  
  dat <- get_effs(x1, x2, g1)
  dat <- subset(dat, dat$pval.exposure < 5e-08)
  unix1_x2 <- summary(lm(dat$beta.outcome ~ -1 + dat$beta.exposure, weights = 1/(dat$se.outcome^2)))
  
  dat <- get_effs(x1, x3, g1)
  dat <- subset(dat, dat$pval.exposure < 5e-08)
  unix1_x3 <- summary(lm(dat$beta.outcome ~ -1 + dat$beta.exposure, weights = 1/(dat$se.outcome^2)))
  
  F_1 <- (dat$beta.exposure/dat$se.exposure)^2
  
  
  #SNPs for all exposures
  dat <- get_effs(x1, y, g)
  dat <- subset(dat, dat$pval.exposure < 5e-08)
  n_g <- nrow(dat)
  unix1_yp <- summary(lm(dat$beta.outcome ~ -1 + dat$beta.exposure, weights = 1/(dat$se.outcome^2)))
  
  dat <- get_effs(x1, x2, g)
  dat <- subset(dat, dat$pval.exposure < 5e-08)
  unix1_x2p <- summary(lm(dat$beta.outcome ~ -1 + dat$beta.exposure, weights = 1/(dat$se.outcome^2)))
  
  dat <- get_effs(x1, x3, g)
  dat <- subset(dat, dat$pval.exposure < 5e-08)
  unix1_x3p <- summary(lm(dat$beta.outcome ~ -1 + dat$beta.exposure, weights = 1/(dat$se.outcome^2)))
  
  F_1p <- (dat$beta.exposure/dat$se.exposure)^2
  
  #g_all <- cbind(g1, g)
  
  mvdat <- make_mvdat(list(x1, x2, x3), y, g)
  mvdat$minp <- apply(mvdat$exposure_pval, 1, min)
  mvdat <- as.data.frame(mvdat)
  mv <- subset(mvdat, mvdat$minp < 5e-08)
  mvmr <- summary(lm(mv$outcome_beta ~ -1 + mv$exposure_beta.x1 + mv$exposure_beta.x2 + mv$exposure_beta.x3, weights = 1/(mv$outcome_se^2)))
  
  rho12 = cor(x1,x2)
  rho13 = cor(x1,x3)
  rho23 = cor(x2,x3)
  sig12 = rho12*mv$exposure_se.x1*mv$exposure_se.x2
  sig13 = rho13*mv$exposure_se.x1*mv$exposure_se.x3
  sig23 = rho23*mv$exposure_se.x2*mv$exposure_se.x3
  
  delta12 <- lm(mv$exposure_beta.x1~ -1 + mv$exposure_beta.x2 + mv$exposure_beta.x3)$coefficients[1]
  delta13 <- lm(mv$exposure_beta.x1~ -1 + mv$exposure_beta.x2 + mv$exposure_beta.x3)$coefficients[2]
  
  delta21 <- lm(mv$exposure_beta.x2~ -1 + mv$exposure_beta.x1 + mv$exposure_beta.x3)$coefficients[1]
  delta23 <- lm(mv$exposure_beta.x2~ -1 + mv$exposure_beta.x1 + mv$exposure_beta.x3)$coefficients[2]
  
  delta31 <- lm(mv$exposure_beta.x3~ -1 + mv$exposure_beta.x1 + mv$exposure_beta.x2)$coefficients[1]
  delta32 <- lm(mv$exposure_beta.x3~ -1 + mv$exposure_beta.x1 + mv$exposure_beta.x2)$coefficients[2]
  
  Qind_1 <- ((mv$exposure_beta.x1 - delta12*mv$exposure_beta.x2 - delta13*mv$exposure_beta.x3)^2)/
    (mv$exposure_se.x1^2 + (delta12^2)*mv$exposure_se.x2^2 + (delta13^2)*mv$exposure_se.x3^2 
                        - 2*delta12*sig12 - 2*delta13*sig13 + delta12*delta13*sig23)
 
  Qind_2 <- ((mv$exposure_beta.x2 - delta21*mv$exposure_beta.x1 - delta23*mv$exposure_beta.x3)^2)/
    (mv$exposure_se.x2^2 + (delta21^2)*mv$exposure_se.x1^2 + (delta23^2)*mv$exposure_se.x3^2 
     - 2*delta21*sig12 - 2*delta23*sig23 + delta21*delta23*sig13)
  
  Qind_3 <- ((mv$exposure_beta.x3 - delta32*mv$exposure_beta.x2 - delta31*mv$exposure_beta.x1)^2)/
    (mv$exposure_se.x3^2 + (delta32^2)*mv$exposure_se.x2^2 + (delta31^2)*mv$exposure_se.x1^2 
     - 2*delta32*sig23 - 2*delta31*sig13 + delta32*delta31*sig12)
  
  
  out <- data.frame(n_g1, unix1_y$coefficients[1,1], unix1_x2$coefficients[1,1], unix1_x3$coefficients[1,1], mean(F_1),
                    n_g, unix1_yp$coefficients[1,1], unix1_x2p$coefficients[1,1], unix1_x3p$coefficients[1,1], mean(F_1p),
                    nrow(mv), mvmr$coefficients[1,1], sum(Qind_1)/(nrow(mv)-2), mvmr$coefficients[2,1], sum(Qind_2)/(nrow(mv)-2),
                    mvmr$coefficients[3,1], sum(Qind_3)/(nrow(mv)-2))
  
  colnames(out) <- c("uni_x1_nsnp", "uni_x1y_b", "uni_x1x2_b", "uni_x1x3_b", "Fx1",
                     "pl_x1_nsnp", "pl_x1y_b", "pl_x1x2_b",  "pl_x1x3_b", "pl_Fx1",
                     "mv_nsnp", "mv_x1_b", "CF_x1",  "mv_x2_b", "CF_x2", "mv_x3_b", "CF_x3")
  
  return(out)
  
}

