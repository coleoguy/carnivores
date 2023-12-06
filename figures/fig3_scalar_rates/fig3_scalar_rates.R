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
           stat = "identity")+
  geom_bar(data = subset(scalars.summarized,
                         scalar.state == "Depressed"),
           aes(y = -proportion),
           stat = "identity")+
  xlab("")+
  ylab("Proportion of Family in Each State")+
  labs(fill = "State of Chromosome \n\ Number Evolution")+
  geom_hline(yintercept=0, 
             linetype="dashed", 
             color = "black")

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





