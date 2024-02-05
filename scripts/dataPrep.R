# Load libraries ----

library(tidyverse)

# Load data ----

chns <- read.csv("./data/CHNS.csv", header=T)

decomp <- read.csv("./data/decomp_12.3.20.csv", header = TRUE) %>%
  mutate(
    Site = as.factor(Site),
    Batch = as.factor(Batch),
    Source_Cond = as.factor(Source_Cond),
    Site_Cond = as.factor(Site_Cond),
    LdryMass = log(dryMass),
    LinitialMass = log(initialMass)
  ) %>%
  filter(numDaysField < 400)

decomp$batch <- as.factor(substr(decomp$Batch, start=2, stop=2))
decomp <- droplevels(decomp[-which(is.na(decomp$dryMass)),])
decomp$months <- cut(decomp$numDaysField, breaks=c(-1,25,50,100,200,400), labels=c(0,1,2,4,8))