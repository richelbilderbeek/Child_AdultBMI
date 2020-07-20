
rm(list = ls(all=TRUE))

resultsA_all <- NULL
file.namesA <- dir(pattern ="sim_output_A")


for (i in 1:length(file.namesA)){
  
  load(file.namesA[i])
  
  resultsA_all <- rbind(resultsA_all, resultsA)
  
}

save(resultsA_all, file=("Results_A.Rda"))

