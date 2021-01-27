
rm(list = ls(all=TRUE))

args <- commandArgs(T)
set.seed(as.numeric(args[1])*100)
job_id <- (as.numeric(args[1]))

library(MASS)
library(simulateGP)
source('MRest.R')

reps = 50
n = 150000 #number of individuals
l = 150    #number of SNPs (for exposures - total)
lo = 50 #number of SNPs for outcome

source('sims.R')

message("filesave", sprintf("sim_output%s.Rda", job_id))

save(resultsA, file=sprintf("sim_output%s.Rda", job_id))
