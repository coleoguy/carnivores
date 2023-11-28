bm <- read.csv("../data/body_mass.csv")

rs <- read.csv("../data/range_size.csv")

for(i in 1:nrow(bm)){
  hit <- which(rs$X == bm$species[i])
  if(length(hit) > 1){
    hit <- sample(hit, 1)
  }
  bm$range_size[i] <- rs[hit, 2]
}

plot(bm$AdultBodyMass_g ~ bm$range_size)
#both variables fail a normality test (p-value significant) so take that 
#into account with the results of the correlation test
shapiro.test(bm$AdultBodyMass_g)
shapiro.test(bm$range_size)

#correlation test
#r = -0.01
cor <- cor.test(bm$AdultBodyMass_g, 
                bm$range_size, 
                method="pearson")
