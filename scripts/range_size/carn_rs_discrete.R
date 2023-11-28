#read in the csv file with data
range <- read.csv("range_size.csv")
#look at a histagram of the data
hist(range$x, breaks =200)
#adds a line to the histogram where the median is
abline(v = median(range$x), col = "red")
#calculates the median of the range data
median <- median(range$x)
#assigns the range values into two groups for the model.
range[,3] <- range[,2] < median
#changes T/F to 1/0
range$V3[range$V3 == "TRUE"] <- 1 
write.csv(range, "carn_rs_discrete.csv", row.names = F)
