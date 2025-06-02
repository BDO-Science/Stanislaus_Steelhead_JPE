library(tidyverse)


#######import and establish necessary data
age <- read.csv('Files/age_structure.csv') #data frame of count of fish by age class in the Stanislaus for 3 years
n_sim <- 500
alpha <- 0.213 #for fecundity function from hodge et al. 2016
beta <- 2.43 #for fecundity function from hodge et al. 2016


#######calculate egg production
sizeProp <- data.frame(Age = c(1,2,3,4,5,6), #setting up size-age structure using FL data from Cramer study
                       minFL = c(104,181,228,324,720,720), 
                       maxFL = c(284,310,400,690,720,720)) %>%
  left_join(age, by = 'Age') %>% #join age data to this
  group_by(Year) %>%
  mutate(prop = prop.table(Count)) #calculate proportion of each age by year

egglist <- list() #empty list for storing egg production for each simulation
for(i in 1:n_sim){
  year <- sample(min(sizeProp$Year):max(sizeProp$Year), 1) #randomly sample a year from the sizeProp data frame
  females <- sample(2075:2450, 1) *.56 #randomly sample from spawner distribution for 2022 and multiple by .56 for females
  eggs <- sizeProp %>% 
    filter(Year == year) %>%
    mutate(females = round(prop*females,0)) %>% #calculate the proportion of females for each age
    uncount(females) %>% #uncount function to give each individual female their own row of data
    rowwise() %>% #indicating to mutate by individual rows
    mutate(FL = sample(minFL:maxFL, 1, replace = TRUE)) %>% #sample from min and max FL for a given age
    ungroup() %>% #ungrouping to speed it up for some reason
    mutate(FL = if_else(Age %in% c(5,6), 720, FL)) %>% #cleaning up some NAs
    mutate(Fecundity = round((alpha *FL)^beta,0)) %>% #calculating with fecundity function
    summarize(sum = sum(Fecundity)) #summing total eggs for all females
  egglist[[i]] <- eggs
}

egg_prod <- bind_rows(egglist) #put total egg production into dataframe