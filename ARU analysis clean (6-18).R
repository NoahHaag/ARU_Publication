setwd("C:/Users/Noah/OneDrive/Desktop/Northeastern/Kimbro lab/My research/Artificial reef meta")
# Library ####
library(ggplot2)
library(RColorBrewer)
library(ggthemes)
library(forcats)
library(dplyr)
library(tidyverse)
library(reshape2)
library(ggpmisc)
library(ggpubr)
library(lme4)
library(car)
library(bbmle)
library(multcompView)
library(ggiraphExtra)
library(ggiraph)
library(stringr)
library(cowplot)
library(patchwork)
library(scales)
library(renv)
library(gg.gap)

# Functions ####
## Standard error
SE <- function(vec){
  sd(vec, na.rm= T) / sqrt(sum(!is.na(vec)))
}

# Load data ####
reefData <- read.csv("Artificial reefs (4-30-26).csv",  fileEncoding="latin1")
papers <- read.csv("List of papers.csv")

# Begin data cleaning ####
reefData[reefData == ''] <- NA

# make sure metrics are numeric 
reefData$Fish.Biomass <- suppressWarnings(as.numeric(reefData$Fish.Biomass))
reefData$Fish.Abundance <- suppressWarnings(as.numeric(reefData$Fish.Abundance))

# make IDs a factor
reefData$ID <- as.factor(reefData$ID)

# create log response ratio data frame
effect.df <- reefData[,c("ID",
                         "Authors",
                         "Latitude",
                         "Longitude",
                         "Ocean.Basin",
                         "Rugosity",
                         "Vertical.relief",
                         "Openings",
                         "General.Reef.Type",
                         "Specific.Reef.Type",
                         "Substrate.type",
                         "Experimental.group",
                         "N",
                         "Averageable",
                         "Fish.Abundance",
                         "Abundance.upper.error",
                         "Abundance.lower.error",
                         "Abundance.SD",
                         "Abundance.SE",
                         "Fish.Biomass",
                         "Biomass.upper.error",
                         "Biomass.lower.error",
                         "Biomass.SD",
                         "Biomass.SE",
                         "Fish.Species.Richness",
                         "Species.richness.upper.error",
                         "Species.richness.lower.error",
                         "Species.richness.SD",
                         "Species.richness.SE")]

# begin cleaning log response ratio data frame

## Standardize Na/None
effect.df[effect.df == 'None'] <- NA

## remove rows which include no experimental group
effect.df <- subset(effect.df, !is.na(Experimental.group))


### experimental3 (from effect.df)
experimental3 <- subset(effect.df, Experimental.group == "Experimental")

## Abundance
abund_mean <- with(experimental3,
                   ave(Fish.Abundance, ID, Rugosity, Vertical.relief, Openings,
                       FUN = function(x) if (all(is.na(x))) NA else mean(x, na.rm = TRUE)))
experimental3$Fish.Abundance <- ifelse(is.na(abund_mean),
                                       experimental3$Fish.Abundance,
                                       abund_mean)

abund_sd <- experimental3 %>%
  filter(is.numeric(N) & is.numeric(Abundance.SD) & !is.na(N) & !is.na(Abundance.SD)) %>%
  group_by(ID, Rugosity, Vertical.relief, Openings) %>%
  summarise(upper = (sum((N - 1) * (Abundance.SD^2))),
            lower = (sum(N - 1, na.rm = TRUE)),
            pooledSD = sqrt(upper/lower))

data_withpooled_SD <- experimental3 %>%
  left_join(abund_sd, by = c("ID","Rugosity", "Vertical.relief", "Openings"))

experimental3$Abundance.SD <- data_withpooled_SD$pooledSD


## Biomass
biom_mean <- with(experimental3,
                  ave(Fish.Biomass, ID, Rugosity, Vertical.relief, Openings,
                      FUN = function(x) if (all(is.na(x))) NA else mean(x, na.rm = TRUE)))
experimental3$Fish.Biomass <- ifelse(is.na(biom_mean),
                                     experimental3$Fish.Biomass,
                                     biom_mean)

bio_sd <- experimental3 %>%
  group_by(ID, Rugosity, Vertical.relief, Openings) %>%
  summarise(upper = (sum((N - 1) * (Biomass.SD^2))),
            lower = (sum(N - 1)),
            pooledSD = sqrt(upper/lower))

data_withpooled_SD <- experimental3 %>%
  left_join(bio_sd, by = c("ID","Rugosity", "Vertical.relief", "Openings"))

experimental3$Biomass.SD <- data_withpooled_SD$pooledSD


## Species richness
rich_mean <- with(experimental3,
                  ave(Fish.Species.Richness, ID, Rugosity, Vertical.relief, Openings,
                      FUN = function(x) if (all(is.na(x))) NA else mean(x, na.rm = TRUE)))
experimental3$Fish.Species.Richness <- ifelse(is.na(rich_mean),
                                              experimental3$Fish.Species.Richness,
                                              rich_mean)

rich_sd <- experimental3 %>%
  group_by(ID, Rugosity, Vertical.relief, Openings) %>%
  summarise(upper = (sum((N - 1) * (Species.richness.SD^2))),
            lower = (sum(N - 1)),
            pooledSD = sqrt(upper/lower))

data_withpooled_SD <- experimental3 %>%
  left_join(rich_sd, by = c("ID","Rugosity", "Vertical.relief", "Openings"))

experimental3$Species.species.SD <- data_withpooled_SD$pooledSD


## Remove duplicate rows (same as distinct)

# experimental3 <- experimental3[!duplicated(experimental3), ]

experimental3 <- distinct(experimental3,
                          across(-c(Authors, Latitude, Longitude, Specific.Reef.Type)),
                          .keep_all = TRUE)

## --- Subset control (soft sediment) group ---
control3 <- subset(effect.df, Experimental.group == "Control")

## --- Replace each value with its mean per paper (ID) ---
control3$Fish.Abundance <- with(control3,
                                ave(Fish.Abundance, ID, FUN = function(x) mean(x, na.rm = TRUE))
)

abund_sd <- control3 %>%
  group_by(ID) %>%
  summarise(upper = (sum((N - 1) * (Abundance.SD^2))),
            lower = (sum(N - 1)),
            pooledSD = sqrt(upper/lower))

data_withpooled_SD <- control3 %>%
  left_join(abund_sd, by = c("ID"))

control3$Abundance.SD <- data_withpooled_SD$pooledSD

control3$Fish.Biomass <- with(control3,
                              ave(Fish.Biomass, ID, FUN = function(x) mean(x, na.rm = TRUE))
)

bio_sd <- control3 %>%
  group_by(ID) %>%
  summarise(upper = (sum((N - 1) * (Biomass.SD^2))),
            lower = (sum(N - 1)),
            pooledSD = sqrt(upper/lower))

data_withpooled_SD <- control3 %>%
  left_join(bio_sd, by = c("ID"))

control3$Biomass.SD <- data_withpooled_SD$pooledSD

control3$Fish.Species.Richness <- with(control3,
                                       ave(Fish.Species.Richness, ID, FUN = function(x) mean(x, na.rm = TRUE))
)

rich_sd <- control3 %>%
  group_by(ID) %>%
  summarise(upper = (sum((N - 1) * (Species.richness.SD^2))),
            lower = (sum(N - 1)),
            pooledSD = sqrt(upper/lower))

data_withpooled_SD <- control3 %>%
  left_join(rich_sd, by = c("ID"))

control3$Species.richness.SD <- data_withpooled_SD$pooledSD

control3 <- distinct(control3,
                     across(-c(Authors, Latitude, Longitude, Specific.Reef.Type)),
                     .keep_all = TRUE)

## each paper has one control value, and multiple experimental values
effect.df.better3 <- merge(experimental3, control3, by.x = c("ID"), by.y = c("ID"))

## NEW CALCULATIONS ####

effect.df.better3$abundance.mean.difference <- effect.df.better3$Fish.Abundance.x - effect.df.better3$Fish.Abundance.y
effect.df.better3$biomass.mean.difference <- effect.df.better3$Fish.Biomass.x - effect.df.better3$Fish.Biomass.y
effect.df.better3$richness.mean.difference <- effect.df.better3$Fish.Species.Richness.x - effect.df.better3$Fish.Species.Richness.y

effect.df.better3$abundance.s2Pooled <- ((effect.df.better3$N.x-1)*(effect.df.better3$Abundance.SD.x^2) + (effect.df.better3$N.y-1)*(effect.df.better3$Abundance.SD.y^2)) / ((effect.df.better3$N.x-1 ) + (effect.df.better3$N.y-1))
effect.df.better3$biomass.s2Pooled <- ((effect.df.better3$N.x-1)*(effect.df.better3$Biomass.SD.x^2) + (effect.df.better3$N.y-1)*(effect.df.better3$Biomass.SD.y^2)) / ((effect.df.better3$N.x-1 ) + (effect.df.better3$N.y-1))
effect.df.better3$richness.s2Pooled <- ((effect.df.better3$N.x-1)*(effect.df.better3$Species.richness.SD.x^2) + (effect.df.better3$N.y-1)*(effect.df.better3$Species.richness.SD.y^2)) / ((effect.df.better3$N.x-1 ) + (effect.df.better3$N.y-1))

effect.df.better3$abundance.SyPooled <- sqrt(effect.df.better3$abundance.s2Pooled / ((effect.df.better3$N.x * effect.df.better3$N.y) / (effect.df.better3$N.x + effect.df.better3$N.y)))
effect.df.better3$biomass.SyPooled <- sqrt(effect.df.better3$biomass.s2Pooled / ((effect.df.better3$N.x * effect.df.better3$N.y) / (effect.df.better3$N.x + effect.df.better3$N.y)))
effect.df.better3$richness.SyPooled <- sqrt(effect.df.better3$richness.s2Pooled / ((effect.df.better3$N.x * effect.df.better3$N.y) / (effect.df.better3$N.x + effect.df.better3$N.y)))

effect.df.better3$abundance.CohenD <- effect.df.better3$abundance.mean.difference / effect.df.better3$abundance.SyPooled
effect.df.better3$biomass.CohenD <- effect.df.better3$biomass.mean.difference / effect.df.better3$biomass.SyPooled
effect.df.better3$richness.CohenD <- effect.df.better3$richness.mean.difference / effect.df.better3$richness.SyPooled

effect.df.better3$abundance.CohenD[is.na(effect.df.better3$abundance.CohenD) | effect.df.better3$abundance.CohenD=="Inf" | effect.df.better3$abundance.CohenD=="-Inf"] = NA
effect.df.better3$richness.CohenD[is.na(effect.df.better3$richness.CohenD) | effect.df.better3$richness.CohenD=="Inf" | effect.df.better3$richness.CohenD=="-Inf"] = NA
effect.df.better3$biomass.CohenD[is.na(effect.df.better3$biomass.CohenD) | effect.df.better3$biomass.CohenD=="Inf" | effect.df.better3$biomass.CohenD=="-Inf"] = NA

### Save for back transformation if needed ####
effect.df.better3$abundance.CohenD.original <- effect.df.better3$abundance.CohenD
effect.df.better3$biomass.CohenD.original <- effect.df.better3$biomass.CohenD
effect.df.better3$richness.CohenD.original <- effect.df.better3$richness.CohenD

### Log transformations #####
effect.df.better3$abundance.CohenD <- log(effect.df.better3$abundance.CohenD + abs(min(effect.df.better3$abundance.CohenD, na.rm = TRUE)) + 1)
effect.df.better3$biomass.CohenD <- log(effect.df.better3$biomass.CohenD + abs(min(effect.df.better3$biomass.CohenD, na.rm = TRUE)) + 1)
effect.df.better3$richness.CohenD <- log(effect.df.better3$richness.CohenD + abs(min(effect.df.better3$richness.CohenD, na.rm = TRUE)) + 1)

effect.df.better3 <- effect.df.better3 %>% filter(!is.na(abundance.CohenD) | !is.na(biomass.CohenD) | !is.na(richness.CohenD))

write.csv(effect.df.better3, "LRR_data_frame.csv")
# Create data frames for summary statistics and model selection ####

# find means for each papers LRR
IDs <- data.frame(unique(effect.df.better3$ID))
colnames(IDs) <- "ID"

# mean metric by paper 
abundanceData <- subset(effect.df.better3, !is.na(effect.df.better3$abundance.CohenD))
biomassData <- subset(effect.df.better3, !is.na(effect.df.better3$biomass.CohenD))
richnessData <- subset(effect.df.better3, !is.na(effect.df.better3$richness.CohenD))

# Supplemental Information ####

## Supplemental Table 1 #####

# abundance metrics for Supplemental Table 1A
abundanceData %>%
  summarise(
    n = n_distinct(ID),
    t_crit = qt(0.975, df = n - 1),
    mean = mean(abundance.CohenD.original, na.rm = TRUE),
    sd = sd(abundance.CohenD.original, na.rm = TRUE),
    se = sd / sqrt(n),
    margin = t_crit * se,
    ci_lower = mean - margin,
    ci_upper = mean + margin
  )

# biomass metrics for Supplemental Table 1A
biomassData %>%
  summarise(
    n = n_distinct(ID),
    t_crit = qt(0.975, df = n - 1),
    mean = mean(biomass.CohenD.original, na.rm = TRUE),
    sd = sd(biomass.CohenD.original, na.rm = TRUE),
    se = sd / sqrt(n),
    margin = t_crit * se,
    ci_lower = mean - margin,
    ci_upper = mean + margin
  )

# richness metrics for Supplemental Table 1A
richnessData %>%
  summarise(
    n = n_distinct(ID),
    t_crit = qt(0.975, df = n - 1),
    mean = mean(richness.CohenD.original, na.rm = TRUE),
    sd = sd(richness.CohenD.original, na.rm = TRUE),
    se = sd / sqrt(n),
    margin = t_crit * se,
    ci_lower = mean - margin,
    ci_upper = mean + margin
  )

## create a data frame with latitudinal zones based on latitude of experiment ####
septData5 <- effect.df.better3
septData5$Ocean.Basin.x[septData5$Ocean.Basin.x == "Mediterranean "]<- "North Atlantic"
septData5$LatZone <- ifelse(septData5$Latitude.x <= 23 & septData5$Latitude.x >= -23, "Tropical", "Subtropical")
septData5$LatZone <- ifelse(septData5$LatZone != "Subtropical",
                            septData5$LatZone,
                            ifelse(septData5$Latitude.x <= 35 & septData5$Latitude.x >= -35, "Subtropical", "Temperate"))
septData5$LatZone <- ifelse(septData5$LatZone != "Temperate",
                            septData5$LatZone,
                            ifelse(septData5$Latitude.x <= 50 & septData5$Latitude.x >= -50, "Temperate", "Subarctic"))
septData5$Ocean.Basin.x <- ifelse(septData5$Ocean.Basin.x != "North Atlantic",
                                  septData5$Ocean.Basin.x,
                                  ifelse((septData5$Longitude.x <= -50),
                                         "Atlantic", "Atlantic"))
septData5$Ocean.Basin.x <- ifelse(septData5$Ocean.Basin.x != "South Atlantic",
                                  septData5$Ocean.Basin.x,
                                  ifelse((septData5$Longitude.x <= -50),
                                         "Atlantic", "Atlantic"))

septData5$Ocean.Basin.x <- ifelse(septData5$Ocean.Basin.x != "North Pacific",
                                  septData5$Ocean.Basin.x,
                                  ifelse((septData5$Longitude.x <= -50),
                                         "Pacific", "Pacific"))
septData5$Ocean.Basin.x <- ifelse(septData5$Ocean.Basin.x != "South Pacific",
                                  septData5$Ocean.Basin.x,
                                  ifelse((septData5$Longitude.x <= -50),
                                         "Pacific", "Pacific"))

# septData5$LatZone <- ifelse(septData5$LatZone == "Tropical" | septData5$LatZone == "Subtropical",
#                yes = "Low Latitude",
#                no = "High Latitude")

## Supplemental Table 2 ####

septData5 %>%
  group_by(Ocean.Basin.x, General.Reef.Type.x) %>%
  filter(!is.na(richness.CohenD)) %>%
  filter(Ocean.Basin.x == "Pacific") %>%
  summarize(
    numPapers = n_distinct(ID)
  )
 
septData5 %>%
  group_by(LatZone, General.Reef.Type.x) %>%
  filter(!is.na(richness.CohenD)) %>%
  filter(LatZone == "Tropical") %>%
  summarize(
    numPapers = n_distinct(ID)
  )

# Summary Statistics ####

abundance_summary <- abundanceData %>%
  tidyr::pivot_longer(
    cols = c(Rugosity.x, Vertical.relief.x, Openings.x),
    names_to = "Factor",
    values_to = "Level"
  ) %>%
  group_by(Factor, Level) %>%
  summarise(
    n = n_distinct(ID),
    t_crit = qt(0.975, df = n - 1), 
    mean = mean(abundance.CohenD, na.rm = TRUE),
    sd = sd(abundance.CohenD, na.rm = TRUE),
    se = sd / sqrt(n),
    margin = t_crit * se,
    ci_lower = mean - margin,
    ci_upper = mean + margin,
    .groups = "drop"
  )

abundance_summary

biomass_summary <- biomassData %>%
  tidyr::pivot_longer(
    cols = c(Rugosity.x, Vertical.relief.x, Openings.x),
    names_to = "Factor",
    values_to = "Level"
  ) %>%
  group_by(Factor, Level) %>%
  summarise(
    n = n_distinct(ID),
    t_crit = qt(0.975, df = n - 1),
    mean = mean(biomass.CohenD, na.rm = TRUE),
    sd = sd(biomass.CohenD, na.rm = TRUE),
    se = sd / sqrt(n),
    margin = t_crit * se,
    ci_lower = mean - margin,
    ci_upper = mean + margin,
    .groups = "drop"
  )

biomass_summary

richness_summary <- richnessData %>%
  tidyr::pivot_longer(
    cols = c(Rugosity.x, Vertical.relief.x, Openings.x),
    names_to = "Factor",
    values_to = "Level"
  ) %>%
  group_by(Factor, Level) %>%
  summarise(
    n = n_distinct(ID),
    t_crit = qt(0.975, df = n - 1),
    mean = mean(richness.CohenD, na.rm = TRUE),
    sd = sd(richness.CohenD, na.rm = TRUE),
    se = sd / sqrt(n),
    margin = t_crit * se,
    ci_lower = mean - margin,
    ci_upper = mean + margin,
    .groups = "drop"
  )

richness_summary

# Summary statistics for latitude zones 
bioExplain <- data.frame((effect.df.better3$Rugosity.x), effect.df.better3$Vertical.relief.x, effect.df.better3$Openings.x, (effect.df.better3$ID))
colnames(bioExplain) <- c("Rugosity", "Vertical Relief", "Openings","ID")

es.LatZone_Lat <- septData5 |>
  select(ID, LatZone, General.Reef.Type.x, 
         abundance.CohenD, biomass.CohenD, richness.CohenD) |>
  rename(LatZone = LatZone, 
         Reef_Type = General.Reef.Type.x,
         meanAbundance = abundance.CohenD,
         meanBiomass = biomass.CohenD,
         meanRichness = richness.CohenD) |>
  distinct(ID, LatZone, Reef_Type, .keep_all = TRUE)

LatZone_Stats <- es.LatZone_Lat %>%
  tidyr::pivot_longer(
    cols = c(meanAbundance, meanBiomass, meanRichness),
    names_to = "Metric",
    values_to = "EffectSize"
  ) %>%
  filter(!is.na(EffectSize)) %>%
  group_by(LatZone, Metric) %>%
  summarise(
    n = n_distinct(ID),
    t_crit = qt(0.975, df = n - 1),
    mean = mean(EffectSize, na.rm = TRUE),
    sd = sd(EffectSize, na.rm = TRUE),
    se = SE(EffectSize),
    margin = t_crit * se,
    ci_lower = mean - margin,
    ci_upper = mean + margin,
    .groups = "drop"
  )

LatZone_Stats

# summary statistics for ocean basin 

es.LatZone_Basin <- septData5 |>
  select(ID, Ocean.Basin.x, General.Reef.Type.x, 
         abundance.CohenD, biomass.CohenD, richness.CohenD) |>
  rename(OceanBasin = Ocean.Basin.x, 
         Reef_Type = General.Reef.Type.x,
         meanAbundance = abundance.CohenD,
         meanBiomass = biomass.CohenD,
         meanRichness = richness.CohenD) |>
  distinct(ID, OceanBasin, Reef_Type, .keep_all = TRUE)



OceanBasin_Stats <- es.LatZone_Basin %>%
  tidyr::pivot_longer(
    cols = c(meanAbundance, meanBiomass, meanRichness),
    names_to = "Metric",
    values_to = "EffectSize"
  ) %>%
  filter(!is.na(EffectSize)) %>%
  group_by(OceanBasin, Metric) %>%
  summarise(
    n = n_distinct(ID),
    t_crit = qt(0.975, df = n - 1),
    mean = mean(EffectSize, na.rm = TRUE),
    sd = sd(EffectSize, na.rm = TRUE),
    se = SE(EffectSize),
    margin = t_crit * se,
    ci_lower = mean - margin,
    ci_upper = mean + margin,
    .groups = "drop"
  )

OceanBasin_Stats

# Model Selection ####

## Latitudinal Zone vs Ocean Basin models #####
abund_septData5 <- septData5 %>% 
  filter(General.Reef.Type.x == "Designed reefs")


latModelAbund <- list()
latModelAbund[[1]] <- lmer(abundance.CohenD ~ 1 + (1|ID), data = abund_septData5)
latModelAbund[[2]] <- lmer(abundance.CohenD ~ Ocean.Basin.x + (1|ID), data = abund_septData5)
# latModelAbund[[3]] <- lmer(abundance.CohenD ~ General.Reef.Type.x + (1|ID), data = abund_septData5)
# latModelAbund[[4]] <- lmer((abundance.CohenD) ~ Ocean.Basin.x + General.Reef.Type.x + (1|ID), data = abund_septData5)

latModelAbund_AIC <- AICctab(latModelAbund,
        weights=TRUE,
        delta=TRUE,
        base = T,
        sort = F) 

latModelAbund_AIC
shapiro.test(resid(latModelAbund[[which.max(latModelAbund_AIC$weight)]]))
qqnorm(resid(latModelAbund[[which.max(latModelAbund_AIC$weight)]]))
qqline(resid(latModelAbund[[which.max(latModelAbund_AIC$weight)]]))

# Latitude zone biomass with full data set

# Latitude zone biomass without subtropical and subarctic, removed for limited data
biomassModelData <- subset(septData5, c(septData5$LatZone != "subtropical" & septData5$LatZone != "subarctic"))
latModelBio2 <- list()
latModelBio2[[1]] <- lmer((biomass.CohenD) ~ 1 + (1|ID), data = biomassModelData)
# latModelBio2[[2]] <- lmer((biomass.CohenD) ~ Ocean.Basin.x + (1|ID), data = biomassModelData)
latModelBio2[[2]] <- lmer((biomass.CohenD) ~ General.Reef.Type.x + (1|ID), data = biomassModelData)
# latModelBio2[[4]] <- lmer((biomass.CohenD) ~ Ocean.Basin.x + General.Reef.Type.x + (1|ID), data = biomassModelData)

latModelBio2_AIC<-AICctab(latModelBio2,
        weights=TRUE,
        delta=TRUE,
        base = T,
        sort = F)
latModelBio2_AIC
shapiro.test(resid(latModelBio2[[which.max(latModelBio2_AIC$weight)]]))
qqnorm(resid(latModelBio2[[which.max(latModelBio2_AIC$weight)]]))
qqline(resid(latModelBio2[[which.max(latModelBio2_AIC$weight)]]))

# latitude zone and ocean basin with subtropical and Indian ocean removed for lack of data
richnessModelData <- subset(septData5, c(septData5$LatZone != "subarctic" & septData5$Ocean.Basin.x != "Indian Ocean"))
latModelRich2 <- list()
latModelRich2[[1]] <- lmer((richness.CohenD) ~ 1 + (1|ID), data = richnessModelData)
# latModelRich2[[2]] <- lmer((richness.CohenD) ~ Ocean.Basin.x + (1|ID), data = richnessModelData)
latModelRich2[[2]] <- lmer((richness.CohenD) ~ General.Reef.Type.x + (1|ID), data = richnessModelData)
# latModelRich2[[4]] <- lmer((richness.CohenD) ~ Ocean.Basin.x + General.Reef.Type.x + (1|ID), data = richnessModelData)


latModelRich2_AIC<-AICctab(latModelRich2,
        weights=TRUE,
        delta=TRUE,
        base = T,
        sort = F)
latModelRich2_AIC

shapiro.test(resid(latModelRich2[[which.max(latModelRich2_AIC$weight)]]))
qqnorm(resid(latModelRich2[[which.max(latModelRich2_AIC$weight)]]))
qqline(resid(latModelRich2[[which.max(latModelRich2_AIC$weight)]]))

## Ecological Predictor Models ####

### Abundance ####
## only additive models due to data imbalances
abundanceData_model <- abundanceData %>%
  filter(General.Reef.Type.x == "Designed reefs")

ecoModelAbund <- list()
ecoModelAbund[[1]] <- lmer(abundance.CohenD~1 + (1|ID), data = abundanceData_model)
##### Single factor
ecoModelAbund[[2]] <- lmer(abundance.CohenD ~ Rugosity.x + (1|ID), data = abundanceData_model)
ecoModelAbund[[3]] <- lmer(abundance.CohenD ~ Vertical.relief.x + (1|ID), data = abundanceData_model)
ecoModelAbund[[4]] <- lmer(abundance.CohenD ~ Openings.x + (1|ID), data = abundanceData_model)
##### Double factor
ecoModelAbund[[5]] <- lmer(abundance.CohenD ~ Rugosity.x + Vertical.relief.x + (1|ID), data = abundanceData_model)
ecoModelAbund[[6]] <- lmer(abundance.CohenD ~ Rugosity.x + Openings.x + (1|ID), data = abundanceData_model)
ecoModelAbund[[7]] <- lmer(abundance.CohenD ~ Vertical.relief.x + Openings.x + (1|ID), data = abundanceData_model)
##### Triple factor
ecoModelAbund[[8]] <- lmer(abundance.CohenD ~ Rugosity.x + Vertical.relief.x + Openings.x + (1|ID), data = abundanceData_model)

ecoModelAbund_AIC<- AICctab(ecoModelAbund,
        weights=TRUE,
        delta=TRUE,
        base = T, sort = F)
ecoModelAbund_AIC

shapiro.test(resid(ecoModelAbund[[which.max(ecoModelAbund_AIC$weight)]]))
qqnorm(resid(ecoModelAbund[[which.max(ecoModelAbund_AIC$weight)]]))
qqline(resid(ecoModelAbund[[which.max(ecoModelAbund_AIC$weight)]]))

### Biomass ####
## only additive models due to data imbalances 
biomassData_model <- biomassData %>%
  filter(General.Reef.Type.x == "Designed reefs")

ecoModelBio <- list()
ecoModelBio[[1]] <- lmer(biomass.CohenD~1 + (1|ID), data = biomassData_model)
##### Single factor
ecoModelBio[[2]] <- lmer(biomass.CohenD ~ Rugosity.x + (1|ID), data = biomassData_model)
ecoModelBio[[3]] <- lmer(biomass.CohenD ~ Vertical.relief.x + (1|ID), data = biomassData_model)
# ecoModelBio[[4]] <- lmer(biomass.CohenD ~ Openings.x + (1|ID), data = biomassData)
##### Double factor
ecoModelBio[[4]] <- lmer(biomass.CohenD ~ Rugosity.x + Vertical.relief.x + (1|ID), data = biomassData_model)
# ecoModelBio[[6]] <- lmer(biomass.CohenD ~ Rugosity.x + Openings.x + (1|ID), data = biomassData)
# ecoModelBio[[5]] <- lmer(biomass.CohenD ~ Vertical.relief.x + Openings.x + (1|ID), data = biomassData_model)
##### Triple factor
# ecoModelBio[[8]] <- lmer(biomass.CohenD ~ Rugosity.x + Vertical.relief.x + Openings.x + (1|ID), data = biomassData)

ecoModelBio_AIC<- AICtab(ecoModelBio,
       weights=TRUE,
       delta=TRUE,
       base = T, 
       sort = F)
ecoModelBio_AIC

shapiro.test(resid(ecoModelBio[[which.max(ecoModelBio_AIC$weight)]]))
qqnorm(resid(ecoModelBio[[which.max(ecoModelBio_AIC$weight)]]))
qqline(resid(ecoModelBio[[which.max(ecoModelBio_AIC$weight)]]))

### Richness ####
## only additive models due to data imbalances
richnessData_model <- richnessData %>%
  filter(General.Reef.Type.x == "Designed reefs")

ecoModelRich <- list()
ecoModelRich[[1]] <- lmer(richness.CohenD~1 + (1|ID), data = richnessData_model)
##### Single factor
ecoModelRich[[2]] <- lmer(richness.CohenD ~ Rugosity.x + (1|ID), data = richnessData_model)
ecoModelRich[[3]] <- lmer(richness.CohenD ~ Vertical.relief.x + (1|ID), data = richnessData_model)
ecoModelRich[[4]] <- lmer(richness.CohenD ~ Openings.x + (1|ID), data = richnessData_model)
##### Double factor
ecoModelRich[[5]] <- lmer(richness.CohenD ~ Rugosity.x + Vertical.relief.x + (1|ID), data = richnessData_model)
ecoModelRich[[6]] <- lmer(richness.CohenD ~ Rugosity.x + Openings.x + (1|ID), data = richnessData_model)
ecoModelRich[[7]] <- lmer(richness.CohenD ~ Vertical.relief.x + Openings.x + (1|ID), data = richnessData_model)
##### Triple factor
ecoModelRich[[8]] <- lmer(richness.CohenD ~ Rugosity.x + Vertical.relief.x + Openings.x + (1|ID), data = richnessData_model)

ecoModelRich_AIC<- AICtab(ecoModelRich,
       weights=TRUE,
       delta=TRUE,
       base = T, 
       sort = F)
ecoModelRich_AIC
shapiro.test(resid(ecoModelRich[[which.max(ecoModelRich_AIC$weight)]]))
qqnorm(resid(ecoModelRich[[which.max(ecoModelRich_AIC$weight)]]))
qqline(resid(ecoModelRich[[which.max(ecoModelRich_AIC$weight)]]))


abundanceData <- subset(effect.df.better3, !is.na(effect.df.better3$abundance.CohenD))
biomassData <- subset(effect.df.better3, !is.na(effect.df.better3$biomass.CohenD))
richnessData <- subset(effect.df.better3, !is.na(effect.df.better3$richness.CohenD))

# Figures ####
## Figure 2 ####

papers$Kept..Yes.or.no.[papers$Kept..Yes.or.no. == "No "]<- "No"
papers$Kept..Yes.or.no.[papers$Kept..Yes.or.no. == "Yes "]<- "Yes"
papers$Kept..Yes.or.no.[papers$Kept..Yes.or.no. == ""]<- NA
papers2 <- subset(papers, !is.na(Kept..Yes.or.no.))

## Correct naming convention for plotting
papers2 <- subset(papers2, Kept..Yes.or.no. == "No" | Kept..Yes.or.no. == "Yes" |Kept..Yes.or.no. == "Maybe" |Kept..Yes.or.no. == "Hardbottom")
papers2$Kept..Yes.or.no.[papers2$Kept..Yes.or.no. == "Maybe"]<- "Natural control"
papers2$Kept..Yes.or.no.[papers2$Kept..Yes.or.no. == "No"]<- "Removed for other reason"
papers2$Kept..Yes.or.no.[papers2$Kept..Yes.or.no. == "Yes"]<- "Bare control"

papers2$Kept..Yes.or.no.[str_detect(papers2$Reason.for.removal, "No control|No controls")]<- "No control"

length(unique(papers2$Kept..Yes.or.no.))

sum((papers2$Kept..Yes.or.no.=="Bare control"), na.rm = T)
sum((papers2$Kept..Yes.or.no.=="Natural control"), na.rm = T)
sum((papers2$Kept..Yes.or.no.=="No control"), na.rm = T)
sum((papers2$Kept..Yes.or.no.=="Removed for other reason"), na.rm = T)
sum((papers2$Kept..Yes.or.no.=="Hardbottom"), na.rm = T)


papers2$count <- 1
papers2 <- papers2[order(papers2$Publication.Year),]
papers2$cum_sum <- ave(papers2$count, papers2$Kept..Yes.or.no., FUN=cumsum)

papers2$Kept..Yes.or.no. <- as.factor(papers2$Kept..Yes.or.no.)

papers3 <- subset(papers2, Kept..Yes.or.no. != "Removed for other reason")

sum((papers3$Kept..Yes.or.no.=="Bare control"), na.rm = T)
sum((papers3$Kept..Yes.or.no.=="Natural control"), na.rm = T)
sum((papers3$Kept..Yes.or.no.=="No control"), na.rm = T)
sum((papers2$Kept..Yes.or.no.=="Removed for other reason"), na.rm = T)
sum((papers3$Kept..Yes.or.no.=="Hardbottom"), na.rm = T)

Figure2<- ggplot(data = papers3)+
  geom_line(aes(y = cum_sum,
                x = Publication.Year,
                color = Kept..Yes.or.no.),
            size = 2.5
  )+
  theme(legend.position = c(0.18,0.871),
        legend.background = element_rect(fill= NA),
        legend.key = element_rect(fill = NA),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "grey95"),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        axis.title = element_text(size = 23),
        axis.text = element_text(size = 23))+
  guides(color=guide_legend(ncol= 1,
                            keyheight = 1,
                            keywidth = 1,
                            reverse = TRUE,
                            title = ""))+
  scale_color_grey(labels = c( "Soft-sediment control","Hardbottom control", "Natural control", "No control"))+
  xlab("Publication year")+
  ylab("Cumulative number of papers on artificial reefs")
Figure2

## Figure 3 ####


# abundanceData %>%
#   group_by(ID)%>%
# summarise(
#     abundance.CohenD = mean(abundance.CohenD, na.rm = TRUE),
#     abundance.CohenD.original = mean(abundance.CohenD.original, na.rm = TRUE),
#   ) -> abundanceData
# 
# biomassData %>%
#   group_by(ID)%>%
#   summarise(
#     biomass.CohenD = mean(biomass.CohenD, na.rm = TRUE),
#     biomass.CohenD.original = mean(biomass.CohenD.original, na.rm = TRUE),
#   ) -> biomassData
# 
# richnessData %>%
#   group_by(ID)%>%
#   summarise(
#     richness.CohenD = mean(richness.CohenD, na.rm = TRUE),
#     richness.CohenD.original = mean(richness.CohenD.original, na.rm = TRUE),
#   ) -> richnessData


Figure3 <- ggplot()+
  geom_point(data = abundanceData,
             aes(x = paste("Abundance"),
                 y = mean(abundanceData$abundance.CohenD.original, na.rm = TRUE),
                 color = 'Abundance'),
             size = 5,
             shape = ifelse(mean(abundanceData$abundance.CohenD.original, na.rm = TRUE) - (qt(0.975, df=(length(unique(abundanceData$ID))-1))*(sd(abundanceData$abundance.CohenD.original, na.rm = TRUE)/sqrt(n_distinct(abundanceData$ID)))) > 0,
                            15, 0),
             stroke = 1
             )+
  geom_errorbar(data = abundanceData,
                aes(x = paste("Abundance"),
                    ymin = mean(abundanceData$abundance.CohenD.original, na.rm = TRUE) - (qt(0.975, df=(length(unique(abundanceData$ID))-1))*(sd(abundanceData$abundance.CohenD.original, na.rm = TRUE)/sqrt(n_distinct(abundanceData$ID)))),
                    ymax = mean(abundanceData$abundance.CohenD.original, na.rm = TRUE) + (qt(0.975, df=(length(unique(abundanceData$ID))-1))*(sd(abundanceData$abundance.CohenD.original, na.rm = TRUE)/sqrt(n_distinct(abundanceData$ID)))),
                    color = 'Abundance'),
                size = 1,
                width = 0
  )+
  geom_point(data = biomassData,
             aes(x = paste("Biomass"),
                 y = mean(biomassData$biomass.CohenD.original, na.rm = TRUE),
                 color = 'Biomass'),
             size = 5, 
             shape = ifelse(mean(biomassData$biomass.CohenD.original, na.rm = TRUE) - (qt(0.975, df=(length(unique(biomassData$ID))-1))*(sd(biomassData$biomass.CohenD.original, na.rm = TRUE)/sqrt(n_distinct(biomassData$ID)))) > 0,
                            16, 1), 
             stroke = 1)+
  geom_errorbar(data = biomassData,
                aes(x = paste("Biomass"),
                    ymin = mean(biomassData$biomass.CohenD.original, na.rm = TRUE) - (qt(0.975, df=(length(unique(biomassData$ID))-1))*(sd(biomassData$biomass.CohenD.original, na.rm = TRUE)/sqrt(n_distinct(biomassData$ID)))),
                    ymax = mean(biomassData$biomass.CohenD.original, na.rm = TRUE) + (qt(0.975, df=(length(unique(biomassData$ID))-1))*(sd(biomassData$biomass.CohenD.original, na.rm = TRUE)/sqrt(n_distinct(biomassData$ID)))),
                    color = 'Biomass'),
                size = 1,
                width = 0
  )+
  geom_point(data = richnessData,
             aes(x = paste("Richness"),
                 y = mean(richnessData$richness.CohenD.original, na.rm = TRUE),
                 color = 'Richness'),
             size = 5, 
             shape = ifelse(mean(richnessData$richness.CohenD.original, na.rm = TRUE) - (qt(0.975, df=(length(unique(richnessData$ID))-1))*(sd(richnessData$richness.CohenD.original, na.rm = TRUE)/sqrt(n_distinct(richnessData$ID)))) > 0,
                            17, 2), 
             stroke = 1)+
  geom_errorbar(data = richnessData,
                aes(x = paste("Richness"),
                    ymin = mean(richnessData$richness.CohenD.original, na.rm = TRUE) - (qt(0.975, df=(length(unique(richnessData$ID))-1))*(sd(richnessData$richness.CohenD.original, na.rm = TRUE)/sqrt(n_distinct(richnessData$ID)))),
                    ymax = mean(richnessData$richness.CohenD.original, na.rm = TRUE) + (qt(0.975, df=(length(unique(richnessData$ID))-1))*(sd(richnessData$richness.CohenD.original, na.rm = TRUE)/sqrt(n_distinct(richnessData$ID)))),
                    color = 'Richness'),
                size = 1,
                width = 0 
  )+
  scale_color_manual(values = c('#e0e002','#73DFFF', '#38A800'))+
  # scale_x_discrete(limits = c(paste("Abundance (", length(unique(abundanceData$ID)),")", sep = ""), paste("Biomass (", length(unique(biomassData$ID)), ")", sep = ""),paste("Richness (", length(unique(richnessData$ID)),")", sep = "")))+
  theme(legend.position = "none")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  xlab(NULL)+
  ylab("Standard Effect Size")+
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid = element_line(color = "grey95"),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    axis.title = element_text(size = 23),
    axis.text = element_text(size = 23))+
  geom_text(data = abundanceData,
            aes(x = paste("Abundance"),
                y = mean(abundanceData$abundance.CohenD.original, na.rm = TRUE) + (qt(0.975, df=(length(unique(abundanceData$ID))-1))*(sd(abundanceData$abundance.CohenD.original, na.rm = TRUE)/sqrt(n_distinct(abundanceData$ID)))) + 5),
            label = paste("n = ", length(unique(abundanceData$ID))," (", length(abundanceData$abundance.CohenD), ")", sep = ""),
            size = 5)+
  geom_text(data = biomassData,
            aes(x = paste("Biomass"),
                y = mean(biomassData$biomass.CohenD.original, na.rm = TRUE) + (qt(0.975, df=(length(unique(biomassData$ID))-1))*(sd(biomassData$biomass.CohenD.original, na.rm = TRUE)/sqrt(n_distinct(biomassData$ID)))) + 5),
            label = paste("n = ", length(unique(biomassData$ID)), " (", length(biomassData$biomass.CohenD), ")" , sep = ""),
            size = 5)+
  geom_text(data = richnessData,
            aes(x = paste("Richness"),
                y = mean(richnessData$richness.CohenD.original, na.rm = TRUE) + (qt(0.975, df=(length(unique(richnessData$ID))-1))*(sd(richnessData$richness.CohenD.original, na.rm = TRUE)/sqrt(n_distinct(richnessData$ID)))) + 5),
            label = paste("n =", length(unique(richnessData$ID)), " (", length(richnessData$richness.CohenD), ")", sep = ""),
            size = 5)

Figure3

ggsave(filename = "Figure3.png",
       path = "Figures/",
       plot = Figure3,
       width = 27,
       height = 22,
       scale = 1,
       units = "cm",
       dpi = 300)

## Figure 4 ####
Figure4Table<- septData5 %>%
  filter(!is.na(abundance.CohenD) & General.Reef.Type.x == "Designed reefs") %>%
  group_by(Ocean.Basin.x) %>% 
  summarise(
    n = n_distinct(ID), 
    N_effect_size = length(abundance.CohenD),
    t_crit = qt(0.975, df = n - 1),
    abundance_offset = abs(min(effect.df.better3$abundance.CohenD.original, na.rm = TRUE)) + 1,
    transformed_mean_abund = mean(abundance.CohenD, na.rm = TRUE),
    transformed_sd_abund = sd(abundance.CohenD, na.rm = TRUE),
    transformed_se_abund = transformed_sd_abund / sqrt(n),
    transfomed_margin_abund = t_crit * transformed_se_abund,
    transformed_ci_lower_abund = transformed_mean_abund - transfomed_margin_abund,
    transformed_ci_upper_abund = transformed_mean_abund + transfomed_margin_abund,
    
    back_transformed_mean_abund = exp(transformed_mean_abund) - abundance_offset,
    ci_lower_abund = exp(transformed_ci_lower_abund) - abundance_offset,
    ci_upper_abund = exp(transformed_ci_upper_abund) - abundance_offset,
    ci_crosses_zero = replace_na(ci_lower_abund <= 1, FALSE),
    shape_val = ifelse(ci_crosses_zero, 0, 15)
    
  )

NewFigure4<- ggplot(data = Figure4Table)+
  geom_point(aes(x = Ocean.Basin.x, 
                 y = back_transformed_mean_abund,
                 color = Ocean.Basin.x,
                 shape = shape_val),
             size = 5,
             stroke = 1)+
  geom_errorbar(aes(x = Ocean.Basin.x, 
                    ymin = ci_lower_abund,
                    ymax = ci_upper_abund,
                    color = Ocean.Basin.x),
                size = 1,
                width = 0)+
  geom_text(aes(x = Ocean.Basin.x, 
                y = ci_upper_abund + 7.5,
                label = paste0("n = ", n, " (", N_effect_size, ")", sep = "")),
            size = 5,
            color = "black")+
  scale_color_manual(values = c('#e0e002','#b3b302', '#868601'),
                     breaks=c("Atlantic", "Indian Ocean", "Pacific"))+
  scale_x_discrete(labels = c("Atlantic Ocean", "Indian Ocean", "Pacific Ocean"),
                   guide = guide_axis(n.dodge = 2))+
  scale_shape_identity() +
  ylab("Standard Effect Size: \n Abundance")+
  xlab(NULL)+
  geom_hline(yintercept = 0, linetype = "dashed")+ 
  guides(color = guide_legend(
    ncol = 1,
    keyheight = 1,
    keywidth = 1,
    reverse = F,
    title = "Metric",
    override.aes = list(shape = c(15), linetype = "blank")
  ))+
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid = element_line(color = "grey95"),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    axis.title = element_text(size = 23),
    axis.text = element_text(size = 23))


ggsave(filename = "NewFigure4.png",
       path = "Figures/",
       plot = NewFigure4,
       width = 27,
       height = 27,
       units = "cm",
       bg = "white")


## Figure 5 ########

abundance_offset = abs(min(effect.df.better3$abundance.CohenD.original, na.rm = TRUE)) + 1
richness_offset = abs(min(effect.df.better3$richness.CohenD.original, na.rm = TRUE)) + 1
biomass_offset = abs(min(effect.df.better3$biomass.CohenD.original, na.rm = TRUE)) + 1

### Panel 1 #######
Newpanel1 <- abundanceData %>%
  filter(General.Reef.Type.x == "Designed reefs") %>%
  tidyr::pivot_longer(
    cols = c(Rugosity.x),
    names_to = "Factor",
    values_to = "Level"
  ) %>%
  group_by(Factor, Level, General.Reef.Type.x) %>%
  summarise(
    abundance_offset = abundance_offset,
    n = n_distinct(ID),
    N_effect_size = length(abundance.CohenD),
    Reef_Type = General.Reef.Type.x,
    t_crit = qt(0.975, df = n - 1), 
    transformed_mean = mean(abundance.CohenD, na.rm = TRUE),
    transformed_sd = sd(abundance.CohenD, na.rm = TRUE),
    transformed_se = transformed_sd / sqrt(n),
    transformed_margin = t_crit * transformed_se,
    transformed_ci_lower = transformed_mean - transformed_margin,
    transformed_ci_upper = transformed_mean + transformed_margin,
    
    mean = exp(transformed_mean) - abundance_offset,
    ci_lower = exp(transformed_ci_lower) - abundance_offset,
    ci_upper = exp(transformed_ci_upper) - abundance_offset,
    
    shape_type = ifelse(ci_lower > 0, 15, 0),
    .groups = "drop"
  ) %>%
  distinct(Factor, Level, Reef_Type, .keep_all = TRUE) %>%
  
  ggplot(aes(x = Level, y = mean, color = Level)) +
  geom_point(aes(shape = as.factor(shape_type)), size = 5, stroke = 1) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Low rugosity" = "#e0e002",
                                "High rugosity" = "#e0e002")) +
  scale_shape_manual(values = c("15" = 15, "0" = 0),
                     labels = c("Low rugosity" = "Low Rugosity", "High rugosity" = "High Rugosity")) +  # map numeric to shape
  scale_x_discrete(limits = c("Low rugosity", "High rugosity"),
                   labels = scales::label_wrap(5)(c("Low Rugosity", "High Rugosity"))) +
  ylim(-5,50) +
  labs(
    x = NULL,
    y = "Standard Effect Size:\n Abundance",
    color = "Rugosity",
    shape = "Significance"
  ) +
  theme(
    legend.position = c("none"),
    legend.background = element_rect(fill = NA),
    legend.key = element_rect(fill = NA),
    panel.background = element_rect(fill = "white"),
    panel.grid = element_line(color = "grey95"),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    axis.title = element_text(size = 23),
    axis.text = element_text(size = 23),
    plot.title = element_text(hjust = 0.5)
  )+
  
  geom_text(aes(y = ci_upper + 2.5,
                label = paste0("n = ", n, " (", N_effect_size, ")", sep = "")),
            size = 5,
            color = "black"
  )

Newpanel1

### Panel 2 #######

Newpanel2 <- abundanceData %>%
  filter(General.Reef.Type.x == "Designed reefs") %>%
  tidyr::pivot_longer(
    cols = c(Vertical.relief.x),
    names_to = "Factor",
    values_to = "Level"
  ) %>%
  group_by(Factor, Level, General.Reef.Type.x) %>%
  summarise(
    abundance_offset = abundance_offset,
    n = n_distinct(ID),
    N_effect_size = length(abundance.CohenD),
    Reef_Type = General.Reef.Type.x,
    t_crit = qt(0.975, df = n - 1), 
    transformed_mean = mean(abundance.CohenD, na.rm = TRUE),
    transformed_sd = sd(abundance.CohenD, na.rm = TRUE),
    transformed_se = transformed_sd / sqrt(n),
    transformed_margin = t_crit * transformed_se,
    transformed_ci_lower = transformed_mean - transformed_margin,
    transformed_ci_upper = transformed_mean + transformed_margin,
    
    mean = exp(transformed_mean) - abundance_offset,
    ci_lower = exp(transformed_ci_lower) - abundance_offset,
    ci_upper = exp(transformed_ci_upper) - abundance_offset,
    
    shape_type = ifelse(ci_lower > 0, 15, 0),
    .groups = "drop"
  ) %>%
  distinct(Factor, Level, Reef_Type, .keep_all = TRUE) %>%
  
  ggplot(aes(x = Level, y = mean, color = Level)) +
  geom_point(aes(shape = as.factor(shape_type)), size = 5, stroke = 1) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Low vertical relief" = "#e0e002",
                                "High vertical relief" = "#e0e002")) +
  scale_shape_manual(values = c("15" = 15, "0" = 0),
                     labels = c("Low vertical relief" = "Low Vertical Relief", "High vertical relief" = "High Vertical Relief")) +  # map numeric to shape
  scale_x_discrete(limits = c("Low vertical relief", "High vertical relief"),
                   labels = scales::label_wrap(5)(c("Low Vertical Relief", "High Vertical Relief"))) +
  ylim(-5,50) +
  labs(
    x = NULL,
    y = NULL,
    color = "Vertical Relief",
    shape = "Significance"
  ) +
  theme(
    legend.position = c("none"),
    legend.background = element_rect(fill = NA),
    legend.key = element_rect(fill = NA),
    panel.background = element_rect(fill = "white"),
    panel.grid = element_line(color = "grey95"),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    axis.title = element_text(size = 23),
    axis.text = element_text(size = 23),
    plot.title = element_text(hjust = 0.5)
  )+
  
  geom_text(aes(y = ci_upper + 2.5,
                label = paste0("n = ", n, " (", N_effect_size, ")", sep = "")),
            size = 5,
            color = "black"
  )

Newpanel2

row1 <- plot_grid(Newpanel1 + Newpanel2, align = "h")

### Panel 3 #######
Newpanel3 <- richnessData %>%
  filter(General.Reef.Type.x == "Designed reefs") %>%
  tidyr::pivot_longer(
    cols = c(Rugosity.x),
    names_to = "Factor",
    values_to = "Level"
  ) %>%
  group_by(Factor, Level, General.Reef.Type.x) %>%
  summarise(
    richness_offset = richness_offset,
    n = n_distinct(ID),
    N_effect_size = length(richness.CohenD),
    Reef_Type = General.Reef.Type.x,
    t_crit = qt(0.975, df = n - 1), 
    transformed_mean = mean(richness.CohenD, na.rm = TRUE),
    transformed_sd = sd(richness.CohenD, na.rm = TRUE),
    transformed_se = transformed_sd / sqrt(n),
    transformed_margin = t_crit * transformed_se,
    transformed_ci_lower = transformed_mean - transformed_margin,
    transformed_ci_upper = transformed_mean + transformed_margin,
    
    mean = exp(transformed_mean) - richness_offset,
    ci_lower = exp(transformed_ci_lower) - richness_offset,
    ci_upper = exp(transformed_ci_upper) - richness_offset,
    
    shape_type = ifelse(ci_lower > 0, 17, 2),
    .groups = "drop"
  ) %>%
  distinct(Factor, Level, Reef_Type, .keep_all = TRUE) %>%
  
  ggplot(aes(x = Level, y = mean, color = Level)) +
  geom_point(aes(shape = as.factor(shape_type)), size = 5, stroke = 1) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Low rugosity" = "#38A800",
                                "High rugosity" = "#38A800")) +
  scale_shape_manual(values = c("17" = 17, "2" = 2),
                     labels = c("Low rugosity" = "Low Rugosity", "High rugosity" = "High Rugosity")) +  # map numeric to shape
  scale_x_discrete(limits = c("Low rugosity", "High rugosity"),
                   labels = scales::label_wrap(5)(c("Low Rugosity", "High Rugosity"))) +
  ylim(-5,100) +
  labs(
    x = NULL,
    y = "Standard Effect Size:\n Richness",
    color = "Rugosity",
    shape = "Significance"
  ) +
  theme(
    legend.position = c("none"),
    legend.background = element_rect(fill = NA),
    legend.key = element_rect(fill = NA),
    panel.background = element_rect(fill = "white"),
    panel.grid = element_line(color = "grey95"),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    axis.title = element_text(size = 23),
    axis.text = element_text(size = 23),
    plot.title = element_text(hjust = 0.5)
  )+
  geom_text(aes(y = ci_upper + 5,
                label = paste0("n = ", n, " (", N_effect_size, ")", sep = "")),
            size = 5,
            color = "black"
  )

Newpanel3

### Panel 4 #############
Newpanel4 <- richnessData %>%
  filter(General.Reef.Type.x == "Designed reefs") %>%
  tidyr::pivot_longer(
    cols = c(Vertical.relief.x),
    names_to = "Factor",
    values_to = "Level"
  ) %>%
  group_by(Factor, Level, General.Reef.Type.x) %>%
  summarise(
    richness_offset = richness_offset,
    n = n_distinct(ID),
    N_effect_size = length(richness.CohenD),
    Reef_Type = General.Reef.Type.x,
    t_crit = qt(0.975, df = n - 1), 
    transformed_mean = mean(richness.CohenD, na.rm = TRUE),
    transformed_sd = sd(richness.CohenD, na.rm = TRUE),
    transformed_se = transformed_sd / sqrt(n),
    transformed_margin = t_crit * transformed_se,
    transformed_ci_lower = transformed_mean - transformed_margin,
    transformed_ci_upper = transformed_mean + transformed_margin,
    
    mean = exp(transformed_mean) - richness_offset,
    ci_lower = exp(transformed_ci_lower) - richness_offset,
    ci_upper = exp(transformed_ci_upper) - richness_offset,
    
    shape_type = ifelse(ci_lower > 0, 17, 2),
    .groups = "drop"
  ) %>%
  distinct(Factor, Level, Reef_Type, .keep_all = TRUE) %>%
  
  ggplot(aes(x = Level, y = mean, color = Level)) +
  geom_point(aes(shape = as.factor(shape_type)), size = 5, stroke = 1) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Low vertical relief" = "#38A800",
                                "High vertical relief" = "#38A800")) +
  scale_shape_manual(values = c("17" = 17, "2" = 2),
                     labels = c("Low vertical relief" = "Low Vertical Relief", "High vertical relief" = "High Vertical Relief")) +  # map numeric to shape
  scale_x_discrete(limits = c("Low vertical relief", "High vertical relief"),
                   labels = scales::label_wrap(15)(c("Low Vertical Relief", "High Vertical Relief"))) +
  ylim(-5,60) +
  labs(
    x = NULL,
    y = "Standard Effect Size:\n Species Richness",
    color = "Vertical Relief",
    shape = "Significance"
  ) +
  theme(
    legend.position = c("none"),
    legend.background = element_rect(fill = NA),
    legend.key = element_rect(fill = NA),
    panel.background = element_rect(fill = "white"),
    panel.grid = element_line(color = "grey95"),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    axis.title = element_text(size = 23),
    axis.text = element_text(size = 23),
    plot.title = element_text(hjust = 0.5)
  )+
  geom_text(aes(y = ci_upper + 2.5,
                label = paste0("n = ", n, " (", N_effect_size, ")", sep = "")),
            size = 5,
            color = "black"
  )

Newpanel4

### Panel 5 #############
Newpanel5 <- richnessData %>%
  filter(General.Reef.Type.x == "Designed reefs") %>%
  tidyr::pivot_longer(
    cols = c(Openings.x),
    names_to = "Factor",
    values_to = "Level"
  ) %>%
  group_by(Factor, Level, General.Reef.Type.x) %>%
  summarise(
    richness_offset = richness_offset,
    n = n_distinct(ID),
    N_effect_size = length(richness.CohenD),
    Reef_Type = General.Reef.Type.x,
    t_crit = qt(0.975, df = n - 1), 
    transformed_mean = mean(richness.CohenD, na.rm = TRUE),
    transformed_sd = sd(richness.CohenD, na.rm = TRUE),
    transformed_se = transformed_sd / sqrt(n),
    transformed_margin = t_crit * transformed_se,
    transformed_ci_lower = transformed_mean - transformed_margin,
    transformed_ci_upper = transformed_mean + transformed_margin,
    
    mean = exp(transformed_mean) - richness_offset,
    ci_lower = exp(transformed_ci_lower) - richness_offset,
    ci_upper = exp(transformed_ci_upper) - richness_offset,
    
    shape_type = ifelse(ci_lower > 0, 17, 2),
    .groups = "drop"
  ) %>%
  distinct(Factor, Level, Reef_Type, .keep_all = TRUE) %>%
  
  ggplot(aes(x = Level, y = mean, color = Level)) +
  geom_point(aes(shape = as.factor(shape_type)), size = 5, stroke = 1) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("No openings" = "#38A800",
                                "With openings" = "#38A800")) +
  scale_shape_manual(values = c("17" = 17, "2" = 2),
                     labels = c("No openings" = "No Openings", "With openings" = "With Openings")) +  # map numeric to shape
  scale_x_discrete(limits = c("No openings", "With openings"),
                   labels = scales::label_wrap(5)(c("No Openings", "With Openings"))) +
  # ylim(-5, 65) +
  labs(
    x = NULL,
    y = NULL,
    color = "Porosity",
    shape = "Significance"
  ) +
  theme(
    legend.position = c("none"),
    legend.background = element_rect(fill = NA),
    legend.key = element_rect(fill = NA),
    panel.background = element_rect(fill = "white"),
    panel.grid = element_line(color = "grey95"),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    axis.title = element_text(size = 23),
    axis.text = element_text(size = 23),
    plot.title = element_text(hjust = 0.5)
  )+
  geom_text(aes(y = ci_upper + 5,
                label = paste0("n = ", n, " (", N_effect_size, ")", sep = "")),
            size = 5,
            color = "black"
  )

Newpanel5

Newpanel5<- gg.gap(plot = Newpanel5,
                   segments = c(60,600),
                   ylim = c(-5,750),
                   tick_width = c(15,100),
                   rel_heights=c(1, 0.1, 1))

row2.1 <- plot_grid(Newpanel3, Newpanel4,
                    align = "h", nrow = 1,
                    rel_widths = c(1, 1, 1),
                    rel_heights = c(0.1, 0.1, 2))

row2 <- plot_grid(row2.1, Newpanel5,
                  align = "v", ncol = 1,
                  rel_heights = c(1,2))

# Figure5 <- plot_grid(row1 , row2.1, Newpanel5,
#                      ncol = 1,
#                      nrow = 3,
#                      align = "v",
#                      rel_heights = c(1, 1, 1.5))


### Panel 6 (Biomass Rugosity) ######
Newpanel6 <- biomassData %>%
  filter(General.Reef.Type.x == "Designed reefs") %>%
  tidyr::pivot_longer(
    cols = c(Rugosity.x),
    names_to = "Factor",
    values_to = "Level"
  ) %>%
  group_by(Factor, Level, General.Reef.Type.x) %>%
  summarise(
    biomass_offset = biomass_offset,
    n = n_distinct(ID),
    N_effect_size = length(biomass.CohenD),
    Reef_Type = General.Reef.Type.x,
    t_crit = qt(0.975, df = n - 1), 
    transformed_mean = mean(biomass.CohenD, na.rm = TRUE),
    transformed_sd = sd(biomass.CohenD, na.rm = TRUE),
    transformed_se = transformed_sd / sqrt(n),
    transformed_margin = t_crit * transformed_se,
    transformed_ci_lower = transformed_mean - transformed_margin,
    transformed_ci_upper = transformed_mean + transformed_margin,
    
    mean = exp(transformed_mean) - biomass_offset,
    ci_lower = exp(transformed_ci_lower) - biomass_offset,
    ci_upper = exp(transformed_ci_upper) - biomass_offset,
    
    shape_type = ifelse(
      is.na(ci_lower) | is.nan(ci_lower),  # Check if ci_lower is NA or NaN
      16,                                # If true (it is missing), set to 17
      ifelse(ci_lower > 0, 16, 1)         # Otherwise, run the original logic
    ),
    .groups = "drop"
  ) %>%
  distinct(Factor, Level, Reef_Type, .keep_all = TRUE) %>%
  
  ggplot(aes(x = Level, y = mean, color = Level)) +
  geom_point(aes(shape = as.factor(shape_type)), size = 5, stroke = 1) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Low rugosity" = "#73DFFF",
                                "High rugosity" = "#73DFFF")) +
  scale_shape_manual(values = c("16" = 16, "1" = 1),
                     labels = c("Low rugosity" = "Low Rugosity", "High rugosity" = "High Rugosity")) +  # map numeric to shape
  scale_x_discrete(limits = c("Low rugosity", "High rugosity"),
                   labels = scales::label_wrap(5)(c("Low Rugosity", "High Rugosity"))) +
  ylim(c(-15, 60)) +
  labs(
    x = NULL,
    y = "Standard Effect Size:\n Biomass",
    color = "Rugosity",
    shape = "Significance"
  ) +
  theme(
    legend.position = c("none"),
    legend.background = element_rect(fill = NA),
    legend.key = element_rect(fill = NA),
    panel.background = element_rect(fill = "white"),
    panel.grid = element_line(color = "grey95"),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    axis.title = element_text(size = 23),
    axis.text = element_text(size = 23),
    plot.title = element_text(hjust = 0.5)
  )+
  geom_text(aes(y = ci_upper + 2.5,
                label = paste0("n = ", n, " (", N_effect_size, ")", sep = "")),
            size = 5,
            color = "black"
  )

### Panel 7 (Biomass Vertical Relief) ######
Newpanel7 <- biomassData %>%
  filter(General.Reef.Type.x == "Designed reefs") %>%
  tidyr::pivot_longer(
    cols = c(Vertical.relief.x),
    names_to = "Factor",
    values_to = "Level"
  ) %>%
  group_by(Factor, Level, General.Reef.Type.x) %>%
  summarise(
    biomass_offset = biomass_offset,
    n = n_distinct(ID),
    N_effect_size = length(biomass.CohenD),
    Reef_Type = General.Reef.Type.x,
    t_crit = qt(0.975, df = n - 1), 
    transformed_mean = mean(biomass.CohenD, na.rm = TRUE),
    transformed_sd = sd(biomass.CohenD, na.rm = TRUE),
    transformed_se = transformed_sd / sqrt(n),
    transformed_margin = t_crit * transformed_se,
    transformed_ci_lower = transformed_mean - transformed_margin,
    transformed_ci_upper = transformed_mean + transformed_margin,
    
    mean = exp(transformed_mean) - biomass_offset,
    ci_lower = exp(transformed_ci_lower) - biomass_offset,
    ci_upper = exp(transformed_ci_upper) - biomass_offset,
    
    shape_type = ifelse(
      is.na(ci_lower) | is.nan(ci_lower),  # Check if ci_lower is NA or NaN
      16,                                # If true (it is missing), set to 17
      ifelse(ci_lower > 0, 16, 1)         # Otherwise, run the original logic
    ),
    .groups = "drop"
  ) %>%
  distinct(Factor, Level, Reef_Type, .keep_all = TRUE) %>%
  
  ggplot(aes(x = Level, y = mean, color = Level)) +
  geom_point(aes(shape = as.factor(shape_type)), size = 5, stroke = 1) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Low vertical relief" = "#73DFFF",
                                "High vertical relief" = "#73DFFF")) +
  scale_shape_manual(values = c("16" = 16, "1" = 1),
                     labels = c("Low vertical relief" = "Low Vertical Relief", "High vertical relief" = "High Vertical Relief")) +  # map numeric to shape
  scale_x_discrete(limits = c("Low vertical relief", "High vertical relief"),
                   labels = scales::label_wrap(5)(c("Low Vertical Relief", "High Vertical Relief"))) +
  ylim(c(-15, 35)) +
  labs(
    x = NULL,
    y = NULL,
    color = "Vertical Relief",
    shape = "Significance"
  ) +
  theme(
    legend.position = c("none"),
    legend.background = element_rect(fill = NA),
    legend.key = element_rect(fill = NA),
    panel.background = element_rect(fill = "white"),
    panel.grid = element_line(color = "grey95"),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    axis.title = element_text(size = 23),
    axis.text = element_text(size = 23),
    plot.title = element_text(hjust = 0.5)
  )+
  geom_text(aes(y = ci_upper + 1.5,
                label = paste0("n = ", n, " (", N_effect_size, ")", sep = "")),
            size = 5,
            color = "black"
  )

row2 <- plot_grid(Newpanel6, Newpanel7, align = "h", nrow = 1)

Figure5 <- plot_grid(row1, row2, Newpanel4,
                     ncol = 1)


ggsave(filename = "NewFigure5.png",
       path = "Figures/",
       plot = Figure5,
       width = 40,
       height = 35,
       units = "cm",
       bg = "white")

