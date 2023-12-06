#Michelle Jonika
#September 7, 2021

#creates a phylogenetic plot visualizing chromosome number and range size as a
#discrete trait

###LOAD IN DATA NEEDED###-------------------------------------------------------

#load in chromosome number and binary trait data
load("../../data/range_size/datalists_rangesize.RData")
#0 = small; 1 = large pop size

#load in tree data
trees <- read.nexus("../../data/range_size/carnivora_rs_pruned.nex")
#reduce to one tree from dataset that is needed
tree <- trees[[1]]
#save the order of the tree tip labels in order from 1:length of data rather
#than in phylogenetic order
test <- untangle(tree, "read.tree")
rm(trees, tree)


#loop through to make sure all of the datalists are in the same order as each
#of the trees
data_all <- c()
for(i in 1:length(datalist[[1]])){
  #empty vector to store correct order in 
  neworder <- c()
  #loop through to store the correct order of the data to match the tree tips
  for(j in 1:length(test$tip.label)){
    neworder[j] <- which(test$tip.label[j] == datalist[[1]]$species)
  }
  #reorder the data into a new data frame
  data_all <- datalist[[1]][neworder,]
}
rm(datalist, i, j, neworder)

#create a named vector with the chromosome data 
data_chrom <- data_all$hap.chrom 
names(data_chrom) <- data_all$species

###PLOTTING###------------------------------------------------------

#make a plot of the chromosome nuber data
plot(y = 1:110,
     x = data_chrom, 
     xlab = "Haploid Chromosome Number",
     xlim = c(15, 40), 
     pch = 21,
     cex = 0.9,
     col = "black",
     bg = c("#FDE725FF", "#39568CFF")[data_all$range.size + 1])

#export as PDF 8.5"x11"



