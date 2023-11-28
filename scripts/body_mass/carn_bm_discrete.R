#read in the csv file with data
bm <- read.csv("../data/body_mass.csv")
#look at a histagram of the data
hist(bm$AdultBodyMass_g, breaks =200)
#adds a line to the histogram where the median is
abline(v = median(bm$AdultBodyMass_g, na.rm = T), col = "red")
#calculates the median of the range data
median <- median(bm$AdultBodyMass_g, na.rm = T)
#assigns the range values into two groups for the model.
bm[,3] <- bm[,2] < median
#changes T/F to 1/0
bm$V3[bm$V3 == "TRUE"] <- 1 
write.csv(bm, "../data/body_mass/carn_bm_discrete.csv", row.names = F)
