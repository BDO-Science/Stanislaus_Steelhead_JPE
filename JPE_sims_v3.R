library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggridges)
library(patchwork)
options(scipen = 999)
params <- read.csv('Files/Params.csv') #read in parameters compiled from literature
age <- read.csv('Files/age_structure.csv') #read in age structure data
df <- data.frame(matrix(ncol = 6, nrow = 1000)) #create empty dataframe for parameter distributions
colnames(df) <- c('Fecundity', 'Egg', 'parr', '0_1', '1_2', '2_plus') #assign column names to dataframe

# Loop through each parameter and generate lognormal values for each parameter
for (i in 1:nrow(params)) {
  temp <- rlnorm(1000, meanlog = log(params$Mean[i]^2 / sqrt(params$SD[i]^2 + params$Mean[i]^2)),
                 sdlog = sqrt(log(1 + params$SD[i]^2 / params$Mean[i]^2)))
  df[, i] <- temp
}

df[df > .99 & df < 10 | df <= 0] <- NA #assign NAs to any 'probabilities' >= 1 or < 0
df <- df %>% na.omit()

#user defined numbers
sims <- 200 #number of simulations
sizeProp <- data.frame(Age = c(1,2,3,4,5,6), #setting up size-age structure
                       minFL = c(104,181,228,324,720,720), 
                       maxFL = c(284,310,400,690,720,720)) %>%
  left_join(age, by = 'Age') %>%
  group_by(Year) %>%
  mutate(prop = prop.table(Count)) 

#create empty dataframes for sim results
JPE_adult <- data.frame(Adult = matrix(ncol = 1, nrow = sims))
JPE_juvenile <- data.frame(Juvenile = matrix(ncol = 1, nrow = sims))
juv <- data.frame(value = matrix(ncol = 1, nrow = sims))

for(i in 1:sims){
  year <- sample(min(sizeProp$Year):max(sizeProp$Year), 1)
  females <- sample(4121:12656, 1)
  eggs <- sizeProp %>% 
    filter(Year == year) %>%
    mutate(females = round(prop*females,0)) %>%
    uncount(females) %>%
    rowwise() %>%
    mutate(FL = sample(minFL:maxFL, 1, replace = TRUE)) %>%
    ungroup() %>%
    mutate(FL = if_else(Age %in% c(5,6), 720, FL)) %>%
    mutate(Fecundity = round((0.213 *FL)^2.43,0)) %>%
    summarize(sum = sum(Fecundity)) %>%
    pull()
  temp <- eggs* #calculates the number of juveniles for each simulation
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

JPE_all <- JPE_adult %>% 
  bind_cols(JPE_juvenile) %>% 
  rename('Adult' = 1, 'Juvenile' = 2) %>%
  gather(key = 'Approach', value = 'JPE', 1:2)

jpe_boxplot <- ggplot(JPE_all, aes(x = Approach, y = JPE)) + 
  geom_boxplot(width = 0.3) +
  theme_minimal() +
  theme(legend.position = 'none',
        plot.margin = margin(0.75,0.25,0.25,0.25, unit = 'cm'),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank()) +
  scale_y_continuous(breaks = seq(0, max(JPE_all$JPE), 500000), labels = scales::comma) +
  coord_flip()
jpe_boxplot

jpe_hist <- ggplot(JPE_all, aes(x = JPE, y = Approach)) +
  geom_density_ridges(stat = "binline", scale = 0.95) +
  scale_x_continuous(breaks = seq(0, max(JPE_all$JPE), 500000), labels = scales::comma) +
  labs(x = 'Juvenile Production Estimate') +
  theme_minimal() +
  theme(legend.position = 'none',
        plot.margin = margin(0.75,0.25,0.25,0.25, unit = 'cm'),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.title.x = element_text(margin=margin(t=12)))
jpe_hist

graph <- jpe_boxplot/jpe_hist
graph

summary <- JPE_all %>%
  group_by(Approach) %>%
  summarize(
    min = min(JPE), 
    q1 = quantile(JPE, 0.25), 
    median = median(JPE), 
    mean = mean(JPE), 
    q3 = quantile(JPE, 0.75), 
    max = max(JPE),
    ci_lower = mean - qt(0.975, df = n() - 1) * (sd(JPE) / sqrt(n())),
    ci_upper = mean + qt(0.975, df = n() - 1) * (sd(JPE) / sqrt(n()))
  ) %>%
  mutate(across(min:ci_upper, ~ scales::comma(round(., 0)))) %>%
  mutate(Range = paste0(min,' - ',max),
         'Median (1st and 3rd Quantiles)' = paste0(median,' (',q1, ' - ',q3,')'),
         'Mean (95% CI)' = paste0(mean, ' (',ci_lower,' - ', ci_upper, ')')) %>%
  select(1,10,11,12) %>%
  t()

write.csv(summary, file = 'summary.csv')
ggsave(graph, file = 'graph.png', units = 'in', width = 8, height = 6)
