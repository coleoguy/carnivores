#script to compare rate values between the 25/75 and 75/25 models


## LOAD IN 25/75 DATA---------------------------------------------------------##
#load in data
carn.25 <- read.csv("results_25small75large/carn_data_all.csv")

#empty list to store
dat_25 <- list()
#vector that stores the beginning row value for each tree
x  <- seq(1, 245001, by=2500)

counter <- 1
#loop that stores each tree's rate data into a separate list
for(i in x){
  dat_25[[counter]] <- carn.25[i:(i+2499), ]
  counter <- counter + 1
}

rm(carn.25, counter, i, x)

## LOAD IN 75/25 DATA---------------------------------------------------------##
#load in data
carn.75 <- read.csv("results_75small25large/carn_data_all.csv")

#empty list to store
dat_75 <- list()
#vector that stores the beginning row value for each tree
x  <- seq(1, 245001, by=2500)

counter <- 1
#loop that stores each tree's rate data into a separate list
for(i in x){
  dat_75[[counter]] <- carn.75[i:(i+2499), ]
  counter <- counter + 1
}

rm(carn.75, counter, i, x)

## FORMAT DATA TO DROP BURNIN-------------------------------------------------##

#drop observations 1:2000 for burnin
for(i in 1:99){
  dat_25[[i]]<- dat_25[[i]][-c(1:2000),]
}

#drop observations 1:2000 for burnin
for(i in 1:99){
  dat_75[[i]]<- dat_75[[i]][-c(1:2000),]
}

rm(i)

model25 <- dat_25[[1]]
model75 <- dat_75[[1]]

long_data_25 <- pivot_longer(model25, cols = c(2:7), names_to = "type", 
                             values_to = "rate")
long_data_25 <- long_data_25[ ,-c(1:2)]

long_data_75 <- pivot_longer(model75, cols = c(2:7), names_to = "type", 
                             values_to = "rate")
long_data_75 <- long_data_75[ ,-c(1:2)]

ggplot(long_data_25 + aes(type, rate)) +
  geom_jitter(cex = .5, alpha = .1,
              position = position_jitterdodge(.3)) +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Parameter") + ylab("Rate")


stat_summary(aes(x = type, y = rate),
             fun.data = "mean_hdci", fun.args = list(mult=1),
             size = 0.4, position = position_jitterdodge(0),
             inherit.aes = FALSE)

ggplot(model75, aes(x=type, y=rate, color=order)) +
  geom_jitter(cex=.5, alpha=.1,position=position_jitterdodge(.3)) +
  theme_bw()+
  theme(legend.position = "none")+
  xlab("Parameter") + ylab("Rate") +
  stat_summary(aes(x=type, y=rate, fill=order), fun.data="mean_hdci", fun.args = list(mult=1),
               size = 0.4, position = position_jitterdodge(0),
               inherit.aes = FALSE)