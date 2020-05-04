
rm(list = ls(all=TRUE))

load("sim_output_A1.Rda")
resultsA_all <- resultsA
load("sim_output_B1.Rda")
resultsB_all <- resultsB
load("sim_output_C1.Rda")
resultsC_all <- resultsC


file.namesA <- dir(pattern ="sim_output_A")
file.namesB <- dir(pattern ="sim_output_B")
file.namesC <- dir(pattern ="sim_output_C")

for (i in 2:40){
  
  load(file.namesA[i])
  load(file.namesB[i])
  load(file.namesC[i])
  
  resultsA_all <- rbind(resultsA_all, resultsA)
  resultsB_all <- rbind(resultsB_all, resultsB)
  resultsC_all <- rbind(resultsC_all, resultsC)
  
}

save(resultsA_all, file=("Results_A.Rda"))
save(resultsB_all, file=("Tesults_B.Rda"))
save(resultsC_all, file=("Results_C.Rda"))
