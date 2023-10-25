#read in the csv file with data
pd <- read.csv("../data/pop_density.csv")
#look at a histagram of the data
hist(pd$Pop.Density.Correction, breaks =200)
#adds a line to the histogram where the median is
abline(v = median(pd$Pop.Density.Correction, na.rm = T), col = "red")
#calculates the median of the range data
median <- median(pd$Pop.Density.Correction, na.rm = T)
#assigns the range values into two groups for the model.
pd[,5] <- pd[,4] < median
#changes T/F to 1/0
pd$V5[pd$V5 == "TRUE"] <- 1 
write.csv(pd, "carn.pd.discrete.csv", row.names = F)
