
rm(list = ls(all=TRUE))


load("sim_output_A1.Rda")
resultsA_all <- resultsA
load("sim_output_B1.Rda")
resultsB_all <- resultsB
load("sim_output_C1.Rda")
resultsC_all <- resultsC

for (i in 1:40){
  load("sim_output_A1.Rda")
  load("sim_output_B1.Rda")
  load("sim_output_C1.Rda")
  
  resultsA_all <- rbind(resultsA_all, resultsA)
  resultsB_all <- rbind(resultsB_all, resultsB)
  resultsC_all <- rbind(resultsC_all, resultsC)
  
}

save(resultsA_all, file=("sim_output_A.Rda"))
save(resultsB_all, file=("sim_output_B.Rda"))
save(resultsC_all, file=("sim_output_C.Rda"))
