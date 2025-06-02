library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggridges)
library(ggpubr)

params <- read.csv('Files/Params.csv') #read in parameters compiled from literature
df <- data.frame(matrix(ncol = 6, nrow = 1000)) #create empty dataframe for parameter distributions
colnames(df) <- c('Fecundity', 'Egg', 'parr', '0_1', '1_2', '2_plus') #assign column names to dataframe

# Loop through each parameter and generate lognormal values for each parameter
for (i in 1:nrow(params)) {
  temp <- rlnorm(1000, meanlog = log(params$Mean[i]^2 / sqrt(params$SD[i]^2 + params$Mean[i]^2)),
                 sdlog = sqrt(log(1 + params$SD[i]^2 / params$Mean[i]^2)))
  df[, i] <- temp
}

df[df > .99 & df < 10 | df <= 0] <- NA #assign NAs to any 'probabilities' >= 1 or < 0
df <- df %>%
  select(-1) %>%
  na.omit() #get rid of NAs

#user defined numbers
sims <- 200 #number of simulations
redd_num <- 1000 #number of redds for adult approach

#Create empty data.frames for JPE simulations
JPE_adult <- data.frame(Adult = matrix(ncol = 1, nrow = sims))
JPE_juvenile <- data.frame(Juvenile = matrix(ncol = 1, nrow = sims))
juv <- data.frame(value = matrix(ncol = 1, nrow = sims))
#for loop to calculate juveniles
for(i in 1:sims){
  temp <- redd_num* #calculates the number of juveniles for each simulation
    sample(df[,1], 1)* 
    sample(df[,2], 1)*
    sample(df[,3], 1)*
    sample(df[,4], 1)
  juv[i,] <- temp #compiles juvenile calculations into dataframe
  for(i in 1:nrow(juv)){ #calculates adult approach from the juvenile calculations
    JPE <- juv$value[i]*
      sample(df[,5], 1)* 
      sample(df[,6], 1)
    JPE_adult[i,] <- JPE #compiles into a dataframe
  }
  for(i in 1:sims){ #calculates juvenile approach based on median of juvenile estimates in juv dataframe
    JPE_juv <- median(juv$value)*
      sample(df[,5], 1)* 
      sample(df[,6], 1)
    JPE_juvenile[i,] <- JPE_juv #compiles into a dataframe
  }
}

#combining and cleaning the two dataframes for graphing
JPE_all <- JPE_adult %>% bind_cols(JPE_juvenile) %>% rename('Adult' = 1, 'Juvenile' = 2) %>%
  gather(key = 'Approach', value = 'JPE', 1:2)

#used to automate scale on graphs
x_max <- plyr::round_any(max(JPE_all$JPE), 5000, f = ceiling)

#boxplot to illustrate differences in two approaches
jpe_graph <- ggplot(JPE_all, aes(x = Approach, y = JPE, fill = Approach)) + 
  geom_boxplot(alpha = 0.5) +
  theme(legend.position = 'none',
        plot.margin = margin(0.75,0.25,0.25,0.25, unit = 'cm'),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin=margin(r=12)),
        axis.text.x = element_blank()) +
  scale_y_continuous(breaks = seq(0,x_max, if_else(x_max > 25000, 10000, 5000)), limits = c(0,x_max),
                     labels = scales::comma) +
  coord_flip()
jpe_graph

#density histogram to illustrate differences in two approaches
hist_lognorm <- ggplot(JPE_all, aes(x = JPE, y = Approach, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE, scale = 1) +
  scale_fill_gradient(low = "white", high = "green4",
                      name = "Tail Prob.") +
  theme(plot.margin = margin(0.25,0.25,0,0.25, unit = 'cm'),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.title.x = element_text(margin=margin(t=12)),
        axis.title.y = element_text(margin=margin(r=12)),
        legend.position=c(.5,.85),
        legend.direction = 'horizontal',
        legend.box.background = element_rect(colour = "black")) +
  scale_x_continuous(breaks = seq(0,x_max, if_else(x_max > 25000, 10000, 5000)), limits = c(0,x_max),
                     labels = scales::comma)
hist_lognorm

#combine both plots
final <- ggarrange(jpe_graph, hist_lognorm, ncol = 1, nrow = 2, align = 'v')
final

#graph of the parameter distributions
dist <- df %>% rename('Egg to Fry Survival' = 1, 'Fry to Parr Survival' = 2, 
                      'Parr to 1+ Survival' = 3, '1+ to 2+ Survival' = 4, 
                      '2+ to Smolt and Survive' = 5) %>% 
  gather(key = 'Parameter', value = 'value', 1:5) %>% 
  mutate(Parameter = factor(Parameter, levels = c('Egg to Fry Survival', 
                                                  'Fry to Parr Survival', 'Parr to 1+ Survival', 
                                                  '1+ to 2+ Survival', '2+ to Smolt and Survive'))) %>%
  ggplot(aes(x = value)) +
  geom_histogram(fill = 'steelblue3') + 
  facet_wrap(~ Parameter, ncol = 2) +
  labs(x = 'Probability', y = 'Frequency') +
  theme_bw() +
  theme(legend.position = 'none',
        plot.margin = margin(0.25,0.25,0,0.25, unit = 'cm'),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 8),
        axis.title.x = element_text(margin=margin(t=12)),
        axis.title.y = element_text(margin=margin(r=12)),
        strip.text = element_text(size = 12),
        panel.spacing.x = unit(6, "mm"))
dist

#save final plots
ggsave(final, file = 'Graphs/final.png', unit = 'px', width = 2500, height = 2000)
ggsave(dist, file = 'Graphs/dist.png', unit = 'in', width = 7, height = 6)

#save parameter distributions and JPE simulation estimates
write.csv(JPE_all, file = 'Files/JPE_sims.csv', row.names = FALSE)
write.csv(df, file = 'Files/Param_dist.csv', row.names = FALSE)
