library(ape)
species <- read.nexus("carnivora.nex")[[1]]$tip.label
chroms <- read.csv("mammal_chroms_update.csv", as.is=T)
# we pared down this dataset to just location and
# species name
ranges <- read.delim("carnivorerange.csv",
                    header = T,
                    sep="\t",
                    encoding="UTF-8",
                    as.is = T)
ranges <- ranges[,c(10,22,23)]
write.csv(ranges, file="carnivorerange-reducedcolumns.csv")

# find the species in both chromosome and tree datasets
good.species <- c()
for(i in 1:nrow(chroms)){
  if(chroms$tree.name[i] %in% species){
    good.species <- c(good.species, chroms$tree.name[i])
  }
}

# example of subsetting the data
# goal here is to keep only the chromosome data for
# species that are in the phylogeny


good.table <- chroms[which(chroms$tree.name %in% good.species ), ]
rm(chroms, species, i)

# read the table that has coordinates
ranges <- read.csv("carnivorerange-reducedcolumns.csv", as.is=T)

# keep all of the data that is for one of the species in good.species
good.species2 <- c()
for(i in 1:length(good.species)){
  good.species2[i] <- paste(strsplit(good.species[i], "_")[[1]], collapse=" ")
}

new.data <- ranges[ranges$species %in% good.species2, ]

write.csv(new.data, "pruned.ranges.csv")

#Need to find the species present in good.species but not in new.data
#all.equal(good.species, new.data$species, ... = ,
#          check.attributes = TRUE, use.names = TRUE)
#isTRUE(all.equal(good.table$species, new.data$species))
#colnames(new.data)

good.species2[!good.species2 %in% new.data$species] 




