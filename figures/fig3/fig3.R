#### PACKAGES ####
library(viridis)
library(ggplot2)

#### LOAD DATA ####
load("../results/scalar.outputs.RData")

#### PLOT ####

#set basic theme
theme_basic <- theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.background= element_blank(),
                     axis.line = element_line(colour = "black"),
                     axis.text.x = element_text(size = 15),
                     axis.text.y = element_text(size = 15),
                     axis.title = element_text(face = "bold",
                                               size = 15),
                     legend.title = element_blank(),
                     plot.title = element_text(face = "bold",
                                               size = 17,
                                               hjust=0.5))

# barplot
scalar.state <- ggplot(data=scalars.summarized ,
                       aes(x=scalar.state,y=proportion,fill=family))+
  geom_bar(stat = "identity",position=position_dodge())+
  geom_bar(stat = "identity",position=position_dodge())+
  #ggtitle(paste0("Overlap in mean scalars across posterior"))+
  scale_y_continuous("Proportion of clade")+
  scale_fill_viridis_d()+
  xlab("Rates of evolution")+
  theme_basic

#save
ggsave(scalar.state,
       filename = paste0("./figures/fig3.pdf"),
       width = 14,
       height = 7,
       units = "in")
ggsave(scalar.state,
       filename = paste0("./figures/fig3.jpg"),
       width = 14,
       height = 7,
       units = "in")
ggsave(scalar.state,
       filename = paste0("./figures/fig3.png"),
       width = 14,
       height = 7,
       units = "in")





