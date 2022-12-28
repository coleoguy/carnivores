# Michelle Jonika
# December 27, 2022

# compares published range size estimates to those that we estimated for our 
# analyses

comps <- read.csv("../data/carn.range.comp.csv")
comps <- comps[-c(2,4,13,19,23),-c(3:4,6)]

cor.test(x = comps$km.2,
         y = comps$Final.Published)
#cor estimate: 0.6085696