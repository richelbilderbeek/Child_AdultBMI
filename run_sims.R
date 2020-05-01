
rm(list = ls(all=TRUE))

#args <- commandArgs(T)
#set.seed(as.numeric(args[1]))
#job_id <- (as.numeric(args[1]))
job_id <- 2


install.packages("devtools")
devtools::install_github("MRCIEU/TwoSampleMR")
library(simulateGP)
library(TwoSampleMR)
source('MRest.R')


reps = 1
n=250000 #number of individuals
l=150    #number of SNPs (total)

source('scenarioA.R')

l1 <- 50 #number of additional SNPs for X1
source('scenarioB.R')

source('scenarioC.R')


message("filesave", sprintf("sim_output_A%s.Rda", job_id))

save(resultsA, file=sprintf("sim_output_A%s.Rda", job_id))
save(resultsB, file=sprintf("sim_output_B%s.Rda", job_id))
save(resultsC, file=sprintf("sim_output_C%s.Rda", job_id))
