
rm(list = ls(all=TRUE))

args <- commandArgs(T)
set.seed(as.numeric(args[1]))
job_id <- (as.numeric(args[1]))

library(simulateGP)
library(TwoSampleMR)
source('MRest.R')

#clock<-proc.time(); clock;

reps = 5
n = 200000 #number of individuals
l = 150    #number of SNPs (for exposures - total)
lo = 50 #number of SNPs for outcome

source('scenarioA.R')

l1 <- 100 #number of additional SNPs for X1
source('scenarioB.R')

source('scenarioC.R')


#clock_dur<-proc.time()-clock; clock_dur;
#(time <- paste0("Computation time = ", round(clock_dur[3]/60, 2)," minutes"))


message("filesave", sprintf("sim_output_A%s.Rda", job_id))

save(resultsA, file=sprintf("sim_output_A%s.Rda", job_id))
save(resultsB, file=sprintf("sim_output_B%s.Rda", job_id))
save(resultsC, file=sprintf("sim_output_C%s.Rda", job_id))
