#### PACKAGES ####
library(viridis)
library(ggplot2)

#### LOAD DATA ####
load("../../results/summarized_scalar_outputs.RData")

#### PLOT ####
# back-to-back barplot
scalar.state <- ggplot(data = scalars.summarized,
                       aes(x = family,
                           y = proportion, 
                           fill = scalar.state))+
  geom_bar(data = subset(scalars.summarized,
                         scalar.state == "Elevated"),
           stat = "identity", 
           fill = "#FDE725FF")+
  geom_bar(data = subset(scalars.summarized,
                         scalar.state == "Depressed"),
           aes(y = -proportion),
           stat = "identity",
           fill = "#39568CFF")+
  xlab("")+
  ylab("Proportion of Branches in Each State")+
  labs(fill = "State of Chromosome \n\ Number Evolution")+
  geom_hline(yintercept=0, 
             linetype="dashed", 
             color = "black") +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", angle = 45, hjust = 1, size = 15),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(color = "black", size = 15))

#save
ggsave(scalar.state,
       filename = paste0("fig3_scalar_rates.pdf"),
       width = 14,
       height = 7,
       units = "in")
ggsave(scalar.state,
       filename = paste0("fig3_scalar_rates.jpg"),
       width = 14,
       height = 7,
       units = "in")
ggsave(scalar.state,
       filename = paste0("fig3_scalar_rates.png"),
       width = 14,
       height = 7,
       units = "in")





