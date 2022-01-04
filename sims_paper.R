
rm(list = ls(all=TRUE))

#args <- commandArgs(T)
#set.seed(as.numeric(args[1])*100)
#job_id <- (as.numeric(args[1]))

library(MASS)
library(simulateGP)
source('MRest.R')

reps = 1
n = 150000 #number of individuals
l = 150    #number of SNPs (for exposures - total)
lo = 10 #number of SNPs for outcome

#source('sims.R')
source('sims_aonly.R')
#message("filesave", sprintf("sim_output%s.Rda", job_id))

#save(resultsA, file=sprintf("sim_output%s.Rda", job_id))
