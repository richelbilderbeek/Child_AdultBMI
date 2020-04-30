
MRest <- function()
{

 mvdat <- make_mvdat(list(x1, x2), y, g)
 mvdat$minp <- apply(mvdat$exposure_pval, 1, min)
 mvdat <- as.data.frame(mvdat)
 dat1 <- subset(mvdat, mvdat$exposure_pval.x1 < 5e-08)
 unix1 <- summary(lm(dat1$outcome_beta ~ -1 + dat1$exposure_beta.x1, weights = 1/(dat1$outcome_se^2)))
 
 dat2 <- subset(mvdat, mvdat$exposure_pval.x2 < 5e-08)
 unix2 <- summary(lm(dat2$outcome_beta ~ -1 + dat2$exposure_beta.x2, weights = 1/(dat2$outcome_se^2)))


  mv <- subset(mvdat, mvdat$minp < 5e-08)
  mvmr <- summary(lm(mv$outcome_beta ~ -1 + mv$exposure_beta.x1 + mv$exposure_beta.x2, weights = 1/(mv$outcome_se^2)))

out <- data.frame(unix1$coefficients[1,1], unix1$coefficients[1,2], nrow(dat1),unix2$coefficients[1,1], unix2$coefficients[1,2], nrow(dat2), 
                      mvmr$coefficients[1,1], mvmr$coefficients[1,2],
                      mvmr$coefficients[2,1], mvmr$coefficients[2,2], nrow(mv))

colnames(out) <- c("uni_x1_b", "uni_x1_se", "uni_x1_nsnp", "uni_x2_b", "uni_x2_se", "uni_x2_nsnp", 
                       "mv_x1_b", "mv_x1_se", "mv_x2_b", "mv_x2_se", "mv_nsnp")
  
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
  
  
  #g_all <- cbind(g1, g)
  
  mvdat <- make_mvdat(list(x1, x2, x3), y, g)
  mvdat$minp <- apply(mvdat$exposure_pval, 1, min)
  mvdat <- as.data.frame(mvdat)
  mv <- subset(mvdat, mvdat$minp < 5e-08)
  mvmr <- summary(lm(mv$outcome_beta ~ -1 + mv$exposure_beta.x1 + mv$exposure_beta.x2 + mv$exposure_beta.x3, weights = 1/(mv$outcome_se^2)))
  
  out <- data.frame(n_g1, unix1_y$coefficients[1,1], unix1_y$coefficients[1,2], 
                    unix1_x2$coefficients[1,1], unix1_x2$coefficients[1,2], 
                    unix1_x3$coefficients[1,1], unix1_x3$coefficients[1,2],
                    n_g, unix1_yp$coefficients[1,1], unix1_yp$coefficients[1,2], 
                    unix1_x2p$coefficients[1,1], unix1_x2p$coefficients[1,2], 
                    unix1_x3p$coefficients[1,1], unix1_x3p$coefficients[1,2],
                    nrow(mv), mvmr$coefficients[1,1], mvmr$coefficients[1,2], 
                    mvmr$coefficients[2,1], mvmr$coefficients[2,2],
                    mvmr$coefficients[3,1], mvmr$coefficients[3,2])
  
  colnames(out) <- c("uni_x1_nsnp", "uni_x1y_b", "uni_x1y_se", "uni_x1x2_b", "uni_x1x2_se", "uni_x1x3_b", "uni_x1x3_se",
                     "pl_x1_nsnp", "pl_x1y_b", "pl_x1y_se", "pl_x1x2_b", "pl_x1x2_se", "pl_x1x3_b", "pl_x1x3_se",
                     "mv_nsnp", "mv_x1_b", "mv_x1_se",  "mv_x2_b", "mv_x2_se", "mv_x3_b", "mv_x3_se")
  
  return(out)
  
}

