#load in libraries to concatenate all results into a single csv file
library(dplyr)
library(readr)

#read in the  list of files
results <- list.files(path = "/Users/michelle/Desktop/Carninvore/results_25small75large", 
                      full.names = TRUE)
#use regex to sort by human number order
fileNum <- as.numeric(gsub(
  '^/Users/michelle/Desktop/Carninvore/results_25small75large/tree.carn([0123456789]*)\\.csv$',
  '\\1', results))

#apply sorting to the list of files
files.ordered <- results[order(fileNum)]

#concatenate all results in human order
carn.data.all <- bind_rows(lapply(files.ordered, read_csv))

#remove first column from reading in
carn.data.all <- carn.data.all[,-c(1)]

#write out a csv with all model results
write.csv(carn.data.all, "results_25small75large/carn_data_all.csv", row.names = F)
