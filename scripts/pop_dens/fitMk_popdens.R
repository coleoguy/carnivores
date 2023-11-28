trees <- read.nexus("../data/pop_dens/carnivora_pd_pruned.nex")
pd <- read.csv("../data/pop_dens/pop_dens.csv")

tree <- trees[[1]]

x <- pd$Pop.Density.Correction
names(x) <- pd$species

mk <- fitMk(tree = tree,
            x = x, 
            model="SYM")
