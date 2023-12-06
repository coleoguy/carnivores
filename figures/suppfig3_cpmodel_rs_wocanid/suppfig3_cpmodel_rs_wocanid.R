# Michelle Jonika
# 6 December
# This code makes the figure for the analysis of chromosome number in
# carnivores with range size as the binary trait

#load in libraries needed
library(chromePlus)

#read in chromplus data for range size
data <- read.csv("../../results/range_size/rs_wocanid.csv")

#read in chromplus data for range size
fission <- temp$asc2 - temp$asc1
fusion <- temp$desc2 - temp$desc1
data_munge <- data.frame(fission, fusion)

plotChromeplus(data = data_munge,
               colors = c("#FDE725FF", "#39568CFF"),
               main_title = "",
               x_title = "rate difference (per MY)\n small - large range size",
               alpha_geom = 0.75,
               alpha_line = 0.45)

#save as PDF 6"x6"