# Load libraries

library(tidyverse)

# Load data ----

source("./scripts/dataPrep.R")

# Prop. mass remaining ----

initialMass0 <- decomp %>%
  filter(months == 0) %>%
  summarise(mean = mean(dryMass)) %>%
  as.numeric()

decomp %>%
  filter(months != 0) %>%
  mutate(propMassRemain = dryMass/initialMass0) %>%
  group_by(Source_Cond, Site_Cond, months) %>%
  summarise(mean = mean(propMassRemain),
            sd = sd(propMassRemain)) %>%
  sort_df("months") %>%
  pivot_wider(names_from = 1:2,
              values_from = 4:5, values_fn = mean) %>%
  write.csv(file = "./out/propMassRemain.csv", quote = FALSE)

chns %>%
  group_by(Source) %>%
  summarise(C = mean(C),
            N = mean(N),
            CN = mean(CN))

plots <- read.csv("../NSSF2_main/Data/NSSF2env_190729.csv") %>%
  left_join(read.csv("../NSSF2_main/Data/NSSF2texture_190729.csv") %>%
              group_by(plot) %>%
              summarise_at(-1, mean),
            by = "plot")

decomp %>%
  select(Site) %>%
  unique() %>%
  mutate(plot = paste0("Q", Site)) %>%
  left_join(plots, by = "plot") %>%
  group_by(type1) %>%
  summarise_at(c(-1,-2), mean) %>%
  mutate(NP = N/P)
