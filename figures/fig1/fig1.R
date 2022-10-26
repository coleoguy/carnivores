#Michelle Jonika
#september 7, 2021

#creates a phylogenetic plot visualizing chromosome number and range size as a
#discrete trait

###LOAD IN DATA NEEDED###-------------------------------------------------------

#load in chromosome number and binary trait data
data <- read.csv("../data/data_sets/dataset1.csv")

#load in tree data
trees <- read.nexus("../data/carnivorapruned.nex")
#reduce to one tree from dataset that is needed
tree <- trees[[1]]


#empty vector to store correct order in 
neworder <- c()
#loop through to store the correct order of the data ito match the tree tips
for(i in 1:length(tree$tip.label)){
  neworder[i] <- which(tree$tip.label[i] == data$species)
}
#reorder the data into a new data frame
data_order <- data[neworder,]

#remove unecessary variables from environment
rm(trees, i, neworder, data)

###PLOTTING###------------------------------------------------------

#make a plot of the chromosome nuber data
plot(y = 1:110, x = data_order$hap.chrom, 
     xlab = "Haploid Chromosome Number",
     xlim = c(14,41), 
     pch = 21,
     cex = 0.9,
     col = "black",
     bg = c("#FDE725FF", "#39568CFF")[data_order$range.size + 1])

#export as PDF 8.5"x11"

#add a legend to the plot
legend(x = "top", 
       legend = c("Small Range Size", "Large Range Size"), 
       pch = 22, 
       pt.cex = 2, 
       box.col = "transparent", 
       pt.bg = c("#FDE725FF", "#39568CFF"))

#export as PDF 8.5"x11"






