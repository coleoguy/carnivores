#load in libraries to concatenate all results into a single csv file
library(dplyr)
library(readr)

#use to dplyr and readr to concatenate csv
carn.data.all <-
  #give a list of files
  list.files(path = "/Users/michelle/Desktop/Carninvore/results/", 
             full.names = TRUE) %>%
  #read all csv files in
  lapply(read_csv) %>%
  #bind tigether all rows of all csv files
  bind_rows 

#remove first column from reading in
carn.data.all <- carn.data.all[,-c(1)]

#write out a csv with all model results
write.csv(carn.data.all, "results/carn_data_all.csv", row.names = F)