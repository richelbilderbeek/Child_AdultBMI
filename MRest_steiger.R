
MRest <- function()
{
  
  mvdat <- data.frame()
  for(s in 1:(l+lo)){
    x1_res <- summary(lm(x1~g[,s]))
    x2_res <- summary(lm(x2~g[,s]))
    y_res <- summary(lm(y~gb[,s]))
    mvdat[s,1] <- (x1_res$coefficient[2,1])
    mvdat[s,2] <- (x2_res$coefficient[2,1])
    mvdat[s,3] <- (y_res$coefficient[2,1])
    
    mvdat[s,4] <- (x1_res$coefficient[2,2])
    mvdat[s,5] <- (x2_res$coefficient[2,2])
    mvdat[s,6] <- (y_res$coefficient[2,2])
    
    mvdat[s,7] <- (x1_res$coefficient[2,4])
    mvdat[s,8] <- (x2_res$coefficient[2,4])
    mvdat[s,9] <- (y_res$coefficient[2,4])
    
  }
  
  colnames(mvdat) <- c("exposure_beta.x1", "exposure_beta.x2", "outcome_beta","exposure_se.x1", "exposure_se.x2", 
                       "outcome_se","exposure_pval.x1", "exposure_pval.x2", "outcome_pval" )
  
 #steiger filtering
  mvdat<-filter(mvdat, !(exposure_pval.x2<5e-8 & outcome_pval < exposure_pval.x2))

 dat1 <- subset(mvdat, mvdat$exposure_pval.x1 < 5e-08)

 unix1 <- summary(lm(dat1$outcome_beta ~ -1 + dat1$exposure_beta.x1, weights = 1/(dat1$outcome_se^2)))
 
 F_1 <- (dat1$exposure_beta.x1/dat1$exposure_se.x1)^2
 
 dat2 <- subset(mvdat, mvdat$exposure_pval.x2 < 5e-08)

 unix2 <- summary(lm(dat2$outcome_beta ~ -1 + dat2$exposure_beta.x2, weights = 1/(dat2$outcome_se^2)))
 F_2 <- (dat2$exposure_beta.x2/dat2$exposure_se.x2)^2
 
 mvdat$minp <- apply(mvdat[,c("exposure_pval.x1","exposure_pval.x2")], 1, min)
  mv <- subset(mvdat, mvdat$minp < 5e-08)
  mvmr <- summary(lm(mv$outcome_beta ~ -1 + mv$exposure_beta.x1 + mv$exposure_beta.x2, weights = 1/(mv$outcome_se^2)))
  
  rho = cor(x1,x2)
  sig12 = as.vector(rho)*mv$exposure_se.x1*mv$exposure_se.x2
  
  delta1 <- lm(mv$exposure_beta.x1~ -1 + mv$exposure_beta.x2)$coefficients[1]
  delta2 <- lm(mv$exposure_beta.x2~ -1 + mv$exposure_beta.x1)$coefficients[1]
  
  Qind_1 <- ((mv$exposure_beta.x1 - delta1*mv$exposure_beta.x2)^2)/(mv$exposure_se.x1^2 + (delta1^2)*mv$exposure_se.x2^2 - 2*delta1*sig12)
  Qind_2 <- ((mv$exposure_beta.x2 - delta2*mv$exposure_beta.x1)^2)/(mv$exposure_se.x2^2 + (delta2^2)*mv$exposure_se.x1^2 - 2*delta2*sig12)
  
  snps_outx1 <- as.numeric(sum((mvdat$exposure_pval.x1[(l+1):(l+lo)]<5e-08)))
  snps_outx2 <- as.numeric(sum((mvdat$exposure_pval.x2[(l+1):(l+lo)]<5e-08)))
  
 
colnames(out) <- c("rho", "L1_L2_cor", "uni_x1_b","uni_x1_se", "uni_x1_nsnp","F_x1" , 
                   "uni_x2_b","uni_x2_se", "uni_x2_nsnp", "F_x2", 
                       "mv_x1_b",     "mv_x1_se", "CF_x1", 
                   "mv_x2_b","mv_x2_se",  "CF_x2","mv_nsnp", 
                   "mv_nsnp_res", "mv_x1_b_res","mv_x1_se_res", "CF_x1_res", 
                   "mv_x2_b_res",   "mv_x2_se_res", "CF_x2_res", 
                   "snps_outx1", "snps_outx2")
  
return(out)

}

