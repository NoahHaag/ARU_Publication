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

# Functions ####
## Standard error
SE <- function(vec){
  sd(vec, na.rm= T) / sqrt(sum(!is.na(vec)))
}

# Load data ####
reefData <- read.csv("Artificial reefs (9-25-25).csv",  fileEncoding="latin1")
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
                         "Specific.Reef.Type",
                         "Substrate.type",
                         "Experimental.group",
                         "Averageable",
                         "Fish.Abundance",
                         "Fish.Biomass",
                         "Fish.Species.Richness")]

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

## Biomass
biom_mean <- with(experimental3,
                  ave(Fish.Biomass, ID, Rugosity, Vertical.relief, Openings,
                      FUN = function(x) if (all(is.na(x))) NA else mean(x, na.rm = TRUE)))
experimental3$Fish.Biomass <- ifelse(is.na(biom_mean),
                                     experimental3$Fish.Biomass,
                                     biom_mean)

## Species richness
rich_mean <- with(experimental3,
                  ave(Fish.Species.Richness, ID, Rugosity, Vertical.relief, Openings,
                      FUN = function(x) if (all(is.na(x))) NA else mean(x, na.rm = TRUE)))
experimental3$Fish.Species.Richness <- ifelse(is.na(rich_mean),
                                              experimental3$Fish.Species.Richness,
                                              rich_mean)


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

control3$Fish.Biomass <- with(control3,
                              ave(Fish.Biomass, ID, FUN = function(x) mean(x, na.rm = TRUE))
)

control3$Fish.Species.Richness <- with(control3,
                                       ave(Fish.Species.Richness, ID, FUN = function(x) mean(x, na.rm = TRUE))
)

control3 <- distinct(control3,
                     across(-c(Authors, Latitude, Longitude, Specific.Reef.Type)),
                     .keep_all = TRUE)

## each paper has one control value, and multiple experimental values
effect.df.better3 <- merge(experimental3, control3, by.x = c("ID"), by.y = c("ID"))

## calculate log response ratio for papers without multi-month reporting
effect.df.better3$LRR.abundance <-log(((effect.df.better3$Fish.Abundance.x) +1) / ((effect.df.better3$Fish.Abundance.y)+1))
effect.df.better3$LRR.Richness <-log(((effect.df.better3$Fish.Species.Richness.x) + 1) / ((effect.df.better3$Fish.Species.Richness.y)+1))
effect.df.better3$LRR.Biomass <- log(((effect.df.better3$Fish.Biomass.x)+1) / ((effect.df.better3$Fish.Biomass.y)+1))

## clean up any Inf issues
effect.df.better3$LRR.abundance[is.na(effect.df.better3$LRR.abundance) | effect.df.better3$LRR.abundance=="Inf"] = NA
effect.df.better3$LRR.Richness[is.na(effect.df.better3$LRR.Richness) | effect.df.better3$LRR.Richness=="Inf"] = NA
effect.df.better3$LRR.Biomass[is.na(effect.df.better3$LRR.Biomass) | effect.df.better3$LRR.Biomass=="Inf"] = NA

write.csv(effect.df.better3, "LRR_data_frame.csv")
# Create data frames for summary statistics and model selection ####

# find means for each papers LRR
IDs <- data.frame(unique(effect.df.better3$ID))
colnames(IDs) <- "ID"

# mean metric by paper 
abundanceData <- subset(effect.df.better3, !is.na(effect.df.better3$LRR.abundance))
biomassData <- subset(effect.df.better3, !is.na(effect.df.better3$LRR.Biomass))
richnessData <- subset(effect.df.better3, !is.na(effect.df.better3$LRR.Richness))

# abundance metrics for Supplemental Table 1A
abundanceData %>%
  summarise(
    n = n_distinct(ID),
    t_crit = qt(0.975, df = n - 1),
    mean = mean(LRR.abundance, na.rm = TRUE),
    sd = sd(LRR.abundance, na.rm = TRUE),
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
    mean = mean(LRR.Biomass, na.rm = TRUE),
    sd = sd(LRR.Biomass, na.rm = TRUE),
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
    mean = mean(LRR.Richness, na.rm = TRUE),
    sd = sd(LRR.Richness, na.rm = TRUE),
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
    mean = mean(LRR.abundance, na.rm = TRUE),
    sd = sd(LRR.abundance, na.rm = TRUE),
    se = SE(LRR.abundance),
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
    mean = mean(LRR.Biomass, na.rm = TRUE),
    sd = sd(LRR.Biomass, na.rm = TRUE),
    se = SE(LRR.Biomass),
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
    mean = mean(LRR.Richness, na.rm = TRUE),
    sd = sd(LRR.Richness, na.rm = TRUE),
    se = SE(LRR.Richness),
    margin = t_crit * se,
    ci_lower = mean - margin,
    ci_upper = mean + margin,
    .groups = "drop"
  )

richness_summary

# Summary statistics for latitude zones 
bioExplain <- data.frame((effect.df.better3$Rugosity.x), effect.df.better3$Vertical.relief.x, effect.df.better3$Openings.x, (effect.df.better3$ID))
colnames(bioExplain) <- c("Rugosity", "Vertical Relief", "Openings","ID")
es.LatZone_Lat <- data.frame(ID = bioExplain$ID,
                             LatZone = septData5$LatZone,
                             meanAbundance=abundanceData$LRR.abundance[pmatch(bioExplain$ID, abundanceData$ID)],
                             meanBiomass=biomassData$LRR.Biomass[pmatch(bioExplain$ID, biomassData$ID)],
                             meanRichness=richnessData$LRR.Richness[pmatch(bioExplain$ID, richnessData$ID)])

es.LatZone_Lat <- filter(es.LatZone_Lat, !is.na(meanAbundance | meanBiomass | meanRichness))

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

es.LatZone_Basin <- data.frame(ID = bioExplain$ID,
                               # LatZone = septData5$LatZone,
                               OceanBasin = septData5$Ocean.Basin.x,
                               meanAbundance=abundanceData$LRR.abundance[pmatch(bioExplain$ID, abundanceData$ID)],
                               meanBiomass=biomassData$LRR.Biomass[pmatch(bioExplain$ID, biomassData$ID)],
                               meanRichness=richnessData$LRR.Richness[pmatch(bioExplain$ID, richnessData$ID)])

es.LatZone_Basin <- filter(es.LatZone_Basin, !is.na(meanAbundance | meanBiomass | meanRichness))

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

latModelAbund <- list()
latModelAbund[[1]] <- lmer(LRR.abundance ~ 1 + (1|ID), data = septData5)
latModelAbund[[2]] <- lmer(LRR.abundance ~ Ocean.Basin.x + (1|ID), data = septData5)
latModelAbund[[3]] <- lmer(LRR.abundance ~ LatZone + (1|ID), data = septData5)
latModelAbund[[4]] <- lmer(LRR.abundance ~ Ocean.Basin.x + LatZone + (1|ID), data = septData5)

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
latModelBio <- list()
latModelBio[[1]] <- lmer((LRR.Biomass) ~ 1 + (1|ID), data = septData5)
latModelBio[[2]] <- lmer((LRR.Biomass) ~ Ocean.Basin.x + (1|ID), data = septData5)
latModelBio[[3]] <- lmer((LRR.Biomass) ~ LatZone + (1|ID), data = septData5)
latModelBio[[4]] <- lmer((LRR.Biomass) ~ Ocean.Basin.x + LatZone + (1|ID), data = septData5)

latModelBio_AIC<-AICctab(latModelBio,
        weights=TRUE,
        delta=TRUE,
        base = T,
        sort = F)
latModelBio_AIC
shapiro.test(resid(latModelBio[[which.max(latModelBio_AIC$weight)]]))
qqnorm(resid(latModelBio[[which.max(latModelBio_AIC$weight)]]))
qqline(resid(latModelBio[[which.max(latModelBio_AIC$weight)]]))

# Latitude zone biomass without subtropical and subarctic, removed for limited data
biomassModelData <- subset(septData5, c(septData5$LatZone != "subtropical" & septData5$LatZone != "subarctic"))
latModelBio2 <- list()
latModelBio2[[1]] <- lmer((LRR.Biomass) ~ 1 + (1|ID), data = biomassModelData)
latModelBio2[[2]] <- lmer((LRR.Biomass) ~ Ocean.Basin.x + (1|ID), data = biomassModelData)
latModelBio2[[3]] <- lmer((LRR.Biomass) ~ LatZone + (1|ID), data = biomassModelData)
latModelBio2[[4]] <- lmer((LRR.Biomass) ~ Ocean.Basin.x + LatZone + (1|ID), data = biomassModelData)

latModelBio2_AIC<-AICctab(latModelBio2,
        weights=TRUE,
        delta=TRUE,
        base = T,
        sort = F)
latModelBio2_AIC
shapiro.test(resid(latModelBio2[[which.max(latModelBio2_AIC$weight)]]))
qqnorm(resid(latModelBio2[[which.max(latModelBio2_AIC$weight)]]))
qqline(resid(latModelBio2[[which.max(latModelBio2_AIC$weight)]]))

# latitude zone and ocean basin with full dataset
latModelRich <- list()
latModelRich[[1]] <- lmer((LRR.Richness) ~ 1 + (1|ID), data = septData5)
latModelRich[[2]] <- lmer((LRR.Richness) ~ Ocean.Basin.x + (1|ID), data = septData5)
latModelRich[[3]] <- lmer((LRR.Richness) ~ LatZone + (1|ID), data = septData5)
latModelRich[[4]] <- lmer((LRR.Richness) ~ Ocean.Basin.x + LatZone + (1|ID), data = septData5)

latModelRich_AIC<- AICctab(latModelRich,
        weights=TRUE,
        delta=TRUE,
        base = T,
        sort = F)
latModelRich_AIC
shapiro.test(resid(latModelRich[[which.max(latModelRich_AIC$weight)]]))
qqnorm(resid(latModelRich[[which.max(latModelRich_AIC$weight)]]))
qqline(resid(latModelRich[[which.max(latModelRich_AIC$weight)]]))

# latitude zone and ocean basin with subtropical and Indian ocean removed for lack of data
richnessModelData <- subset(septData5, c(septData5$LatZone != "subtropical" & septData5$Ocean.Basin.x != "Indian Ocean"))
latModelRich2 <- list()
latModelRich2[[1]] <- lmer((LRR.Richness) ~ 1 + (1|ID), data = richnessModelData)
latModelRich2[[2]] <- lmer((LRR.Richness) ~ Ocean.Basin.x + (1|ID), data = richnessModelData)
latModelRich2[[3]] <- lmer((LRR.Richness) ~ LatZone + (1|ID), data = richnessModelData)
latModelRich2[[4]] <- lmer((LRR.Richness) ~ Ocean.Basin.x + LatZone + (1|ID), data = richnessModelData)

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

ecoModelAbund <- list()
ecoModelAbund[[1]] <- lmer(LRR.abundance~1 + (1|ID), data = abundanceData)
##### Single factor
ecoModelAbund[[2]] <- lmer(LRR.abundance ~ Rugosity.x + (1|ID), data = abundanceData)
ecoModelAbund[[3]] <- lmer(LRR.abundance ~ Vertical.relief.x + (1|ID), data = abundanceData)
ecoModelAbund[[4]] <- lmer(LRR.abundance ~ Openings.x + (1|ID), data = abundanceData)
##### Double factor
ecoModelAbund[[5]] <- lmer(LRR.abundance ~ Rugosity.x + Vertical.relief.x + (1|ID), data = abundanceData)
ecoModelAbund[[6]] <- lmer(LRR.abundance ~ Rugosity.x + Openings.x + (1|ID), data = abundanceData)
ecoModelAbund[[7]] <- lmer(LRR.abundance ~ Vertical.relief.x + Openings.x + (1|ID), data = abundanceData)
##### Triple factor
ecoModelAbund[[8]] <- lmer(LRR.abundance ~ Rugosity.x + Vertical.relief.x + Openings.x + (1|ID), data = abundanceData)

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

ecoModelBio <- list()
ecoModelBio[[1]] <- lmer(LRR.Biomass~1 + (1|ID), data = biomassData)
##### Single factor
ecoModelBio[[2]] <- lmer(LRR.Biomass ~ Rugosity.x + (1|ID), data = biomassData)
ecoModelBio[[3]] <- lmer(LRR.Biomass ~ Vertical.relief.x + (1|ID), data = biomassData)
ecoModelBio[[4]] <- lmer(LRR.Biomass ~ Openings.x + (1|ID), data = biomassData)
##### Double factor
ecoModelBio[[5]] <- lmer(LRR.Biomass ~ Rugosity.x + Vertical.relief.x + (1|ID), data = biomassData)
ecoModelBio[[6]] <- lmer(LRR.Biomass ~ Rugosity.x + Openings.x + (1|ID), data = biomassData)
ecoModelBio[[7]] <- lmer(LRR.Biomass ~ Vertical.relief.x + Openings.x + (1|ID), data = biomassData)
##### Triple factor
ecoModelBio[[8]] <- lmer(LRR.Biomass ~ Rugosity.x + Vertical.relief.x + Openings.x + (1|ID), data = biomassData)

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

ecoModelRich <- list()
ecoModelRich[[1]] <- lmer(LRR.Richness~1 + (1|ID), data = richnessData)
##### Single factor
ecoModelRich[[2]] <- lmer(LRR.Richness ~ Rugosity.x + (1|ID), data = richnessData)
ecoModelRich[[3]] <- lmer(LRR.Richness ~ Vertical.relief.x + (1|ID), data = richnessData)
ecoModelRich[[4]] <- lmer(LRR.Richness ~ Openings.x + (1|ID), data = richnessData)
##### Double factor
ecoModelRich[[5]] <- lmer(LRR.Richness ~ Rugosity.x + Vertical.relief.x + (1|ID), data = richnessData)
ecoModelRich[[6]] <- lmer(LRR.Richness ~ Rugosity.x + Openings.x + (1|ID), data = richnessData)
ecoModelRich[[7]] <- lmer(LRR.Richness ~ Vertical.relief.x + Openings.x + (1|ID), data = richnessData)
##### Triple factor
ecoModelRich[[8]] <- lmer(LRR.Richness ~ Rugosity.x + Vertical.relief.x + Openings.x + (1|ID), data = richnessData)

ecoModelRich_AIC<- AICtab(ecoModelRich,
       weights=TRUE,
       delta=TRUE,
       base = T, 
       sort = F)
ecoModelRich_AIC
shapiro.test(resid(ecoModelRich[[which.max(ecoModelRich_AIC$weight)]]))
qqnorm(resid(ecoModelRich[[which.max(ecoModelRich_AIC$weight)]]))
qqline(resid(ecoModelRich[[which.max(ecoModelRich_AIC$weight)]]))

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

Figure3 <- ggplot()+
  geom_point(data = abundanceData,
             aes(x = "Abundance (23)",
                 y = mean(abundanceData$LRR.abundance, na.rm = TRUE),
                 color = 'Abundance'),
             size = 5,
             shape = 15)+
  geom_errorbar(data = abundanceData,
                aes(x = "Abundance (23)" ,
                    ymin = mean(abundanceData$LRR.abundance, na.rm = TRUE) - (qt(0.975, df=(length(unique(abundanceData$ID))-1))*SE(abundanceData$LRR.abundance)),
                    ymax = mean(abundanceData$LRR.abundance, na.rm = TRUE) + (qt(0.975, df=(length(unique(abundanceData$ID))-1))*SE(abundanceData$LRR.abundance)),
                    color = 'Abundance'),
                size = 1,
                width = 0
  )+
  geom_point(data = biomassData,
             aes(x = "Biomass (12)",
                 y = mean(biomassData$LRR.Biomass, na.rm = TRUE),
                 color = 'Biomass'),
             size = 5, 
             shape = 16, 
             stroke = 1)+
  geom_errorbar(data = biomassData,
                aes(x = "Biomass (12)" ,
                    ymin = mean(biomassData$LRR.Biomass, na.rm = TRUE) - (qt(0.975, df=(length(unique(biomassData$ID))-1))*SE(biomassData$LRR.Biomass)),
                    ymax = mean(biomassData$LRR.Biomass, na.rm = TRUE) + (qt(0.975, df=(length(unique(biomassData$ID))-1))*SE(biomassData$LRR.Biomass)),
                    color = 'Biomass'),
                size = 1,
                width = 0
  )+
  geom_point(data = richnessData,
             aes(x = "Richness (21)",
                 y = mean(richnessData$LRR.Richness, na.rm = TRUE),
                 color = 'Richness'),
             size = 5, 
             shape = 17, 
             stroke = 1)+
  geom_errorbar(data = richnessData,
                aes(x = "Richness (21)" ,
                    ymin = mean(richnessData$LRR.Richness, na.rm = TRUE) - (qt(0.975, df=(length(unique(richnessData$ID))-1))*SE(richnessData$LRR.Richness)),
                    ymax = mean(richnessData$LRR.Richness, na.rm = TRUE) + (qt(0.975, df=(length(unique(richnessData$ID))-1))*SE(richnessData$LRR.Richness)),
                    color = 'Richness'),
                size = 1,
                width = 0 
  )+
  scale_color_manual(values = c('#73DFFF','#e0e002', '#38A800'))+
  scale_x_discrete(limits = c(paste("Abundance (", length(unique(abundanceData$ID)),")", sep = ""), paste("Biomass (", length(unique(biomassData$ID)), ")", sep = ""),paste("Richness (", length(unique(richnessData$ID)),")", sep = "")))+
  theme(legend.position = "none")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  xlab(NULL)+
  ylab("Log response ratio")+
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid = element_line(color = "grey95"),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    axis.title = element_text(size = 23),
    axis.text = element_text(size = 23))

Figure3

## Figure 4 ####

level_order <- c("Tropical", "Subtropical", "Temperate", "Subarctic")
es.LatZone_Lat$LatZone <- factor(es.LatZone_Lat$LatZone, ordered = T, level_order)

### Figure 4A ####
Figure4A <- ggplot(data = es.LatZone_Lat)+
  geom_point(data = subset(es.LatZone_Lat, !is.na(es.LatZone_Lat$meanBiomass) & es.LatZone_Lat$LatZone == "Tropical"),
             aes(x = "Tropical",
                 y = mean(meanBiomass),
                 color = 'Biomass'),
             size = 5,
             position = position_nudge(x = 0), 
             shape = 16, 
             stroke = 1)+
  geom_errorbar(data = subset(es.LatZone_Lat, !is.na(es.LatZone_Lat$meanBiomass) & es.LatZone_Lat$LatZone == "Tropical"),
                aes(x = "Tropical" ,
                    ymin = mean(meanBiomass, na.rm = T) - (qt(0.975, df=(length(unique(ID))-1))*SE(meanBiomass)),
                    ymax = mean(meanBiomass, na.rm = T) + (qt(0.975, df=(length(unique(ID))-1))*SE(meanBiomass)),
                    color = 'Biomass'),
                size = 1,
                width = 0,
                position = position_nudge(x = 0)
  )+
  geom_point(data = subset(es.LatZone_Lat, !is.na(es.LatZone_Lat$meanAbundance) & es.LatZone_Lat$LatZone == "Tropical"),
             aes(x = "Tropical",
                 y = mean(meanAbundance),
                 color = 'Abundance'),
             size = 5,
             shape = 15,
             position = position_nudge(x = -0.25))+
  geom_errorbar(data = subset(es.LatZone_Lat, !is.na(es.LatZone_Lat$meanAbundance) & es.LatZone_Lat$LatZone == "Tropical"),
                aes(x = "Tropical" ,
                    ymin = mean(meanAbundance, na.rm = T) - (qt(0.975, df=(length(unique(ID))-1))*SE(meanAbundance)),
                    ymax = mean(meanAbundance, na.rm = T) + (qt(0.975, df=(length(unique(ID))-1))*SE(meanAbundance)),
                    color = 'Abundance'),
                size = 1,
                width = 0,
                position = position_nudge(x =-0.25)
  )+
  geom_point(data = subset(es.LatZone_Lat, !is.na(es.LatZone_Lat$meanRichness) & es.LatZone_Lat$LatZone == "Tropical"),
             aes(x = "Tropical",
                 y = mean(meanRichness),
                 color = 'Richness'),
             size = 4,
             position = position_nudge(x = 0.25), 
             shape = 17, 
             stroke = 1)+
  geom_errorbar(data = subset(es.LatZone_Lat, !is.na(es.LatZone_Lat$meanRichness) & es.LatZone_Lat$LatZone == "Tropical"),
                aes(x = "Tropical" ,
                    ymin = mean(meanRichness, na.rm = T) - (qt(0.975, df=(length(unique(ID))-1))*SE(meanRichness)),
                    ymax = mean(meanRichness, na.rm = T) + (qt(0.975, df=(length(unique(ID))-1))*SE(meanRichness)),
                    color = 'Richness'),
                size = 1,
                width = 0,
                position = position_nudge(x =0.25)
  )+
  geom_point(data = subset(es.LatZone_Lat, !is.na(es.LatZone_Lat$meanBiomass) & es.LatZone_Lat$LatZone == "Subtropical"),
             aes(x = "Subtropical",
                 y = mean(meanBiomass),
                 color = 'Biomass'),
             size = 5,
             position = position_nudge(x = 0), 
             shape = 16, 
             stroke = 1)+
  geom_errorbar(data = subset(es.LatZone_Lat, !is.na(es.LatZone_Lat$meanBiomass) & es.LatZone_Lat$LatZone == "Subtropical"),
                aes(x = "Subtropical" ,
                    ymin = mean(meanBiomass, na.rm = T) - (qt(0.975, df=(length(unique(ID))-1))*SE(meanBiomass)),
                    ymax = mean(meanBiomass, na.rm = T) + (qt(0.975, df=(length(unique(ID))-1))*SE(meanBiomass)),
                    color = 'Biomass'),
                size = 1,
                width = 0,
                position = position_nudge(x = 0)
  )+
  geom_point(data = subset(es.LatZone_Lat, !is.na(es.LatZone_Lat$meanAbundance) & es.LatZone_Lat$LatZone == "Subtropical"),
             aes(x = "Subtropical",
                 y = mean(meanAbundance),
                 color = 'Abundance'),
             size = 5,
             position = position_nudge(x = -0.25), 
             shape = 0, 
             stroke = 1)+
  geom_errorbar(data = subset(es.LatZone_Lat, !is.na(es.LatZone_Lat$meanAbundance) & es.LatZone_Lat$LatZone == "Subtropical"),
                aes(x = "Subtropical" ,
                    ymin = mean(meanAbundance, na.rm = T) - (qt(0.975, df=(length(unique(ID))-1))*SE(meanAbundance)),
                    ymax = mean(meanAbundance, na.rm = T) + (qt(0.975, df=(length(unique(ID))-1))*SE(meanAbundance)),
                    color = 'Abundance'),
                size = 1,
                width = 0,
                position = position_nudge(x = -0.25)
  )+
  geom_point(data = subset(es.LatZone_Lat, !is.na(es.LatZone_Lat$meanRichness) & es.LatZone_Lat$LatZone == "Subtropical"),
             aes(x = "Subtropical",
                 y = mean(meanRichness),
                 color = 'Richness'),
             size = 4,
             position = position_nudge(x = 0.25), 
             shape = 17, 
             stroke = 1)+
  geom_errorbar(data = subset(es.LatZone_Lat, !is.na(es.LatZone_Lat$meanRichness) & es.LatZone_Lat$LatZone == "Subtropical"),
                aes(x = "Subtropical" ,
                    ymin = mean(meanRichness, na.rm = T) - (qt(0.975, df=(length(unique(ID))-1))*SE(meanRichness)),
                    ymax = mean(meanRichness, na.rm = T) + (qt(0.975, df=(length(unique(ID))-1))*SE(meanRichness)),
                    color = 'Richness'),
                size = 1,
                width = 0,
                position = position_nudge(x =0.25)
  )+
  geom_point(data = subset(es.LatZone_Lat, !is.na(es.LatZone_Lat$meanBiomass) & es.LatZone_Lat$LatZone == "Temperate"),
             aes(x = "Temperate",
                 y = mean(meanBiomass),
                 color = 'Biomass'),
             size = 5,
             position = position_nudge(x = 0), 
             shape = 1, 
             stroke = 1)+
  geom_errorbar(data = subset(es.LatZone_Lat, !is.na(es.LatZone_Lat$meanBiomass) & es.LatZone_Lat$LatZone == "Temperate"),
                aes(x = "Temperate" ,
                    ymin = mean(meanBiomass, na.rm = T) - (qt(0.975, df=(length(unique(ID))-1))*SE(meanBiomass)),
                    ymax = mean(meanBiomass, na.rm = T) + (qt(0.975, df=(length(unique(ID))-1))*SE(meanBiomass)),
                    color = 'Biomass'),
                size = 1,
                width = 0,
                position = position_nudge(x = 0)
  )+
  geom_point(data = subset(es.LatZone_Lat, !is.na(es.LatZone_Lat$meanAbundance) & es.LatZone_Lat$LatZone == "Temperate"),
             aes(x = "Temperate",
                 y = mean(meanAbundance),
                 color = 'Abundance'),
             size = 5,
             position = position_nudge(x = -0.25), 
             shape = 0, 
             stroke = 1)+
  geom_errorbar(data = subset(es.LatZone_Lat, !is.na(es.LatZone_Lat$meanAbundance) & es.LatZone_Lat$LatZone == "Temperate"),
                aes(x = "Temperate" ,
                    ymin = mean(meanAbundance, na.rm = T) - (qt(0.975, df=(length(unique(ID))-1))*SE(meanAbundance)),
                    ymax = mean(meanAbundance, na.rm = T) + (qt(0.975, df=(length(unique(ID))-1))*SE(meanAbundance)),
                    color = 'Abundance'),
                size = 1,
                width = 0,
                position = position_nudge(x = -0.25)
  )+
  geom_point(data = subset(es.LatZone_Lat, !is.na(es.LatZone_Lat$meanRichness) & es.LatZone_Lat$LatZone == "Temperate"),
             aes(x = "Temperate",
                 y = mean(meanRichness),
                 color = 'Richness'),
             size = 4,
             position = position_nudge(x = 0.25), 
             shape = 17, 
             stroke = 1)+
  geom_errorbar(data = subset(es.LatZone_Lat, !is.na(es.LatZone_Lat$meanRichness) & es.LatZone_Lat$LatZone == "Temperate"),
                aes(x = "Temperate" ,
                    ymin = mean(meanRichness, na.rm = T) - (qt(0.975, df=(length(unique(ID))-1))*SE(meanRichness)),
                    ymax = mean(meanRichness, na.rm = T) + (qt(0.975, df=(length(unique(ID))-1))*SE(meanRichness)),
                    color = 'Richness'),
                size = 1,
                width = 0,
                position = position_nudge(x =0.25)
  )+
  geom_point(data = subset(es.LatZone_Lat, !is.na(es.LatZone_Lat$meanBiomass) & es.LatZone_Lat$LatZone == "Subarctic"),
             aes(x = "Subarctic",
                 y = mean(meanBiomass),
                 color = 'Biomass'),
             size = 5,
             position = position_nudge(x = 0), 
             shape = 1, 
             stroke = 1)+
  geom_errorbar(data = subset(es.LatZone_Lat, !is.na(es.LatZone_Lat$meanBiomass) & es.LatZone_Lat$LatZone == "Subarctic"),
                aes(x = "Subarctic" ,
                    ymin = mean(meanBiomass, na.rm = T) - (qt(0.975, df=(length(unique(ID))-1))*SE(meanBiomass)),
                    ymax = mean(meanBiomass, na.rm = T) + (qt(0.975, df=(length(unique(ID))-1))*SE(meanBiomass)),
                    color = 'Biomass'),
                size = 1,
                width = 0,
                position = position_nudge(x = 0)
  )+
  geom_point(data = subset(es.LatZone_Lat, !is.na(es.LatZone_Lat$meanAbundance) & es.LatZone_Lat$LatZone == "Subarctic"),
             aes(x = "Subarctic",
                 y = mean(meanAbundance),
                 color = 'Abundance'),
             size = 5,
             position = position_nudge(x = -0.25), 
             shape = 0, 
             stroke = 1)+
  geom_errorbar(data = subset(es.LatZone_Lat, !is.na(es.LatZone_Lat$meanAbundance) & es.LatZone_Lat$LatZone == "Subarctic"),
                aes(x = "Subarctic" ,
                    ymin = mean(meanAbundance, na.rm = T) - (qt(0.975, df=(length(unique(ID))-1))*SE(meanAbundance)),
                    ymax = mean(meanAbundance, na.rm = T) + (qt(0.975, df=(length(unique(ID))-1))*SE(meanAbundance)),
                    color = 'Abundance'),
                size = 1,
                width = 0,
                position = position_nudge(x =-0.25 )
  )+
  geom_point(data = subset(es.LatZone_Lat, !is.na(es.LatZone_Lat$meanRichness) & es.LatZone_Lat$LatZone == "Subarctic"),
             aes(x = "Subarctic",
                 y = mean(meanRichness),
                 color = 'Richness'),
             size = 4,
             position = position_nudge(x = 0.25), 
             shape = 2, 
             stroke = 1)+
  geom_errorbar(data = subset(es.LatZone_Lat, !is.na(es.LatZone_Lat$meanRichness) & es.LatZone_Lat$LatZone == "Subarctic"),
                aes(x = "Subarctic" ,
                    ymin = mean(meanRichness, na.rm = T) - (qt(0.975, df=(length(unique(ID))-1))*SE(meanRichness)),
                    ymax = mean(meanRichness, na.rm = T) + (qt(0.975, df=(length(unique(ID))-1))*SE(meanRichness)),
                    color = 'Richness'),
                size = 1,
                width = 0,
                position = position_nudge(x =0.25)
  )+
  scale_color_manual(values = c('#e0e002','#73DFFF', '#38A800'),
                     breaks=c("Abundance", "Biomass", "Richness"))+
  scale_x_discrete(limits = c("Tropical", "Subtropical", "Temperate", "Subarctic"),
                   guide = guide_axis(n.dodge = 2))+
  ylab("Log response ratio")+
  xlab(NULL)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  # ylim(c(-2.5,3.5))+
  guides(color = guide_legend(
    ncol = 1,
    keyheight = 1,
    keywidth = 1,
    reverse = F,
    title = "Metric",
    override.aes = list(shape = c(15,16,17), linetype = "blank")
  ))+
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid = element_line(color = "grey95"),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    axis.title = element_text(size = 23),
    axis.text = element_text(size = 23))+
  geom_text(data = subset(es.LatZone_Lat, !is.na(es.LatZone_Lat$meanBiomass) & es.LatZone_Lat$LatZone == "Tropical"),
            aes(x="Tropical",
                y= mean(meanBiomass, na.rm = T) + (qt(0.975, df=(length(unique(ID))-1))*SE(meanBiomass)) + 0.25,
                label= paste("n =", length(unique(ID)))),
            position = position_nudge(x = 0),
            size = 5)+
  geom_text(data = subset(es.LatZone_Lat, !is.na(es.LatZone_Lat$meanAbundance) & es.LatZone_Lat$LatZone == "Tropical"),
            aes(x="Tropical",
                y= mean(meanAbundance, na.rm = T) + (qt(0.975, df=(length(unique(ID))-1))*SE(meanAbundance)) + 0.25,
                label= paste("n =", length(unique(ID)))),
            position = position_nudge(x = -0.25),
            size = 5)+
  geom_text(data = subset(es.LatZone_Lat, !is.na(es.LatZone_Lat$meanRichness) & es.LatZone_Lat$LatZone == "Tropical"),
            aes(x="Tropical",
                y= mean(meanRichness, na.rm = T) + (qt(0.975, df=(length(unique(ID))-1))*SE(meanRichness)) + 0.25,
                label= paste("n =", length(unique(ID)))),
            position = position_nudge(x =0.25),
            size = 5)+
  geom_text(data = subset(es.LatZone_Lat, !is.na(es.LatZone_Lat$meanBiomass) & es.LatZone_Lat$LatZone == "Subtropical"),
            aes(x="Subtropical",
                y= mean(meanBiomass, na.rm = T) + 0.5,
                label= paste("n = 1")),
            position = position_nudge(x = 0),
            size = 5)+
  geom_text(data = subset(es.LatZone_Lat, !is.na(es.LatZone_Lat$meanAbundance) & es.LatZone_Lat$LatZone == "Subtropical"),
            aes(x="Subtropical",
                y= mean(meanAbundance, na.rm = T) + (qt(0.975, df=(length(unique(ID))-1))*SE(meanAbundance)) + 0.25,
                label= paste("n =", length(unique(ID)))),
            position = position_nudge(x = -0.25),
            size = 5)+
  geom_text(data = subset(es.LatZone_Lat, !is.na(es.LatZone_Lat$meanRichness) & es.LatZone_Lat$LatZone == "Subtropical"),
            aes(x="Subtropical",
                y= mean(meanRichness, na.rm = T) + 0.5,
                label= paste("n = 1")),
            position = position_nudge(x =0.25),
            size = 5)+
  geom_text(data = subset(es.LatZone_Lat, !is.na(es.LatZone_Lat$meanBiomass) & es.LatZone_Lat$LatZone == "Temperate"),
            aes(x="Temperate",
                y= mean(meanBiomass, na.rm = T) + (qt(0.975, df=(length(unique(ID))-1))*SE(meanBiomass)) + 0.25,
                label= paste("n =", length(unique(ID)))),
            position = position_nudge(x = 0),
            size = 5)+
  geom_text(data = subset(es.LatZone_Lat, !is.na(es.LatZone_Lat$meanAbundance) & es.LatZone_Lat$LatZone == "Temperate"),
            aes(x="Temperate",
                y= mean(meanAbundance, na.rm = T) + (qt(0.975, df=(length(unique(ID))-1))*SE(meanAbundance)) + 0.25,
                label= paste("n =", length(unique(ID)))),
            position = position_nudge(x = -0.25),
            size = 5)+
  geom_text(data = subset(es.LatZone_Lat, !is.na(es.LatZone_Lat$meanRichness) & es.LatZone_Lat$LatZone == "Temperate"),
            aes(x="Temperate",
                y= mean(meanRichness, na.rm = T) + (qt(0.975, df=(length(unique(ID))-1))*SE(meanRichness)) + 0.25,
                label= paste("n =", length(unique(ID)))),
            position = position_nudge(x =0.25),
            size = 5)+
  geom_text(data = subset(es.LatZone_Lat, !is.na(es.LatZone_Lat$meanBiomass) & es.LatZone_Lat$LatZone == "Subarctic"),
            aes(x="Subarctic",
                y= mean(meanBiomass, na.rm = T) + (qt(0.975, df=(length(unique(ID))-1))*SE(meanBiomass)) + 0.25,
                label= paste("n =", length(unique(ID)))),
            position = position_nudge(x = 0),
            size = 5)+
  geom_text(data = subset(es.LatZone_Lat, !is.na(es.LatZone_Lat$meanAbundance) & es.LatZone_Lat$LatZone == "Subarctic"),
            aes(x="Subarctic",
                y= mean(meanAbundance, na.rm = T) + (qt(0.975, df=(length(unique(ID))-1))*SE(meanAbundance)) + 0.75,
                label= paste("n =", length(unique(ID)))),
            position = position_nudge(x = -0.25),
            size = 5)+
  geom_text(data = subset(es.LatZone_Lat, !is.na(es.LatZone_Lat$meanRichness) & es.LatZone_Lat$LatZone == "Subarctic"),
            aes(x="Subarctic",
                y= mean(meanRichness, na.rm = T) + (qt(0.975, df=(length(unique(ID))-1))*SE(meanRichness)) + 0.25,
                label= paste("n =", length(unique(ID)))),
            position = position_nudge(x =0.25),
            size = 5)

Figure4A

Figure4A<- plot_grid(Figure4A+theme(legend.position = c(0.9,0.9)),
                     ncol = 1)
Figure4A

### Figure 4B #####

Figure4B <- ggplot(data = es.LatZone_Basin)+
  geom_point(data = subset(es.LatZone_Basin, !is.na(es.LatZone_Basin$meanAbundance) & es.LatZone_Basin$OceanBasin == "Atlantic"),
             aes(x = "Atlantic Ocean",
                 y = mean(meanAbundance),
                 color = 'Abundance'),
             size = 5,
             position = position_nudge(x = -0.25), 
             shape = 0, 
             stroke = 1)+
  geom_errorbar(data = subset(es.LatZone_Basin, !is.na(es.LatZone_Basin$meanAbundance) & es.LatZone_Basin$OceanBasin == "Atlantic"),
                aes(x = "Atlantic Ocean" ,
                    ymin = mean(meanAbundance, na.rm = T) - (qt(0.975, df=(length(unique(ID))-1))*SE(meanAbundance)),
                    ymax = mean(meanAbundance, na.rm = T) + (qt(0.975, df=(length(unique(ID))-1))*SE(meanAbundance)),
                    color = 'Abundance'),
                size = 1,
                width = 0,
                position = position_nudge(x = -0.25)
  )+
  geom_point(data = subset(es.LatZone_Basin, !is.na(es.LatZone_Basin$meanBiomass) & es.LatZone_Basin$OceanBasin == "Atlantic"),
             aes(x = "Atlantic Ocean",
                 y = mean(meanBiomass),
                 color = 'Biomass'),
             size = 5,
             position = position_nudge(x = 0), 
             shape = 1, 
             stroke = 1)+
  geom_errorbar(data = subset(es.LatZone_Basin, !is.na(es.LatZone_Basin$meanBiomass) & es.LatZone_Basin$OceanBasin == "Atlantic"),
                aes(x = "Atlantic Ocean" ,
                    ymin = mean(meanBiomass, na.rm = T) - (qt(0.975, df=(length(unique(ID))-1))*SE(meanBiomass)),
                    ymax = mean(meanBiomass, na.rm = T) + (qt(0.975, df=(length(unique(ID))-1))*SE(meanBiomass)),
                    color = 'Biomass'),
                size = 1,
                width = 0,
                position = position_nudge(x = 0)
  )+
  geom_point(data = subset(es.LatZone_Basin, !is.na(es.LatZone_Basin$meanRichness) & es.LatZone_Basin$OceanBasin == "Atlantic"),
             aes(x = "Atlantic Ocean",
                 y = mean(meanRichness),
                 color = 'Richness'),
             size = 5,
             position = position_nudge(x =0.25), 
             shape = 17,
             stroke = 1)+
  geom_errorbar(data = subset(es.LatZone_Basin, !is.na(es.LatZone_Basin$meanRichness) & es.LatZone_Basin$OceanBasin == "Atlantic"),
                aes(x = "Atlantic Ocean" ,
                    ymin = mean(meanRichness, na.rm = T) - (qt(0.975, df=(length(unique(ID))-1))*SE(meanRichness)),
                    ymax = mean(meanRichness, na.rm = T) + (qt(0.975, df=(length(unique(ID))-1))*SE(meanRichness)),
                    color = 'Richness'),
                size = 1,
                width = 0,
                position = position_nudge(x =0.25)
  )+
  geom_point(data = subset(es.LatZone_Basin, !is.na(es.LatZone_Basin$meanAbundance) & es.LatZone_Basin$OceanBasin == "Indian Ocean"),
             aes(x = "Indian Ocean",
                 y = mean(meanAbundance),
                 color = 'Abundance'),
             size = 5,
             position = position_nudge(x = -0.25), 
             shape = 0, 
             stroke = 1)+
  geom_errorbar(data = subset(es.LatZone_Basin, !is.na(es.LatZone_Basin$meanAbundance) & es.LatZone_Basin$OceanBasin == "Indian Ocean"),
                aes(x = "Indian Ocean" ,
                    ymin = mean(meanAbundance, na.rm = T) - (qt(0.975, df=(length(unique(ID))-1))*SE(meanAbundance)),
                    ymax = mean(meanAbundance, na.rm = T) + (qt(0.975, df=(length(unique(ID))-1))*SE(meanAbundance)),
                    color = 'Abundance'),
                size = 1,
                width = 0,
                position = position_nudge(x = -0.25)
  )+
  geom_point(data = subset(es.LatZone_Basin, !is.na(es.LatZone_Basin$meanBiomass) & es.LatZone_Basin$OceanBasin == "Indian Ocean"),
             aes(x = "Indian Ocean",
                 y = mean(meanBiomass),
                 color = 'Biomass'),
             size = 5,
             position = position_nudge(x = 0), 
             shape = 1, 
             stroke = 1)+
  geom_errorbar(data = subset(es.LatZone_Basin, !is.na(es.LatZone_Basin$meanBiomass) & es.LatZone_Basin$OceanBasin == "Indian Ocean"),
                aes(x = "Indian Ocean" ,
                    ymin = mean(meanBiomass, na.rm = T) - (qt(0.975, df=(length(unique(ID))-1))*SE(meanBiomass)),
                    ymax = mean(meanBiomass, na.rm = T) + (qt(0.975, df=(length(unique(ID))-1))*SE(meanBiomass)),
                    color = 'Biomass'),
                size = 1,
                width = 0,
                position = position_nudge(x = 0)
  )+
  geom_point(data = subset(es.LatZone_Basin, !is.na(es.LatZone_Basin$meanRichness) & es.LatZone_Basin$OceanBasin == "Indian Ocean"),
             aes(x = "Indian Ocean",
                 y = mean(meanRichness),
                 color = 'Richness'),
             shape = 17,
             size = 5,
             position = position_nudge(x =0.25))+
  geom_errorbar(data = subset(es.LatZone_Basin, !is.na(es.LatZone_Basin$meanRichness) & es.LatZone_Basin$OceanBasin == "Indian Ocean"),
                aes(x = "Indian Ocean" ,
                    ymin = mean(meanRichness, na.rm = T) - (qt(0.975, df=(length(unique(ID))-1))*SE(meanRichness)),
                    ymax = mean(meanRichness, na.rm = T) + (qt(0.975, df=(length(unique(ID))-1))*SE(meanRichness)),
                    color = 'Richness'),
                size = 1,
                width = 0,
                position = position_nudge(x =0.25)
  )+
  geom_point(data = subset(es.LatZone_Basin, !is.na(es.LatZone_Basin$meanAbundance) & es.LatZone_Basin$OceanBasin == "Pacific"),
             aes(x = "Pacific Ocean",
                 y = mean(meanAbundance),
                 color = 'Abundance'),
             size = 5,
             shape = 15,
             position = position_nudge(x = -0.25))+
  geom_errorbar(data = subset(es.LatZone_Basin, !is.na(es.LatZone_Basin$meanAbundance) & es.LatZone_Basin$OceanBasin == "Pacific"),
                aes(x = "Pacific Ocean" ,
                    ymin = mean(meanAbundance, na.rm = T) - (qt(0.975, df=(length(unique(ID))-1))*SE(meanAbundance)),
                    ymax = mean(meanAbundance, na.rm = T) + (qt(0.975, df=(length(unique(ID))-1))*SE(meanAbundance)),
                    color = 'Abundance'),
                size = 1,
                width = 0,
                position = position_nudge(x = -0.25)
  )+
  geom_point(data = subset(es.LatZone_Basin, !is.na(es.LatZone_Basin$meanBiomass) & es.LatZone_Basin$OceanBasin == "Pacific"),
             aes(x = "Pacific Ocean",
                 y = mean(meanBiomass),
                 color = 'Biomass'),
             size = 5,
             position = position_nudge(x = 0), 
             shape = 1, 
             stroke = 1)+
  geom_errorbar(data = subset(es.LatZone_Basin, !is.na(es.LatZone_Basin$meanBiomass) & es.LatZone_Basin$OceanBasin == "Pacific"),
                aes(x = "Pacific Ocean" ,
                    ymin = mean(meanBiomass, na.rm = T) - (qt(0.975, df=(length(unique(ID))-1))*SE(meanBiomass)),
                    ymax = mean(meanBiomass, na.rm = T) + (qt(0.975, df=(length(unique(ID))-1))*SE(meanBiomass)),
                    color = 'Biomass'),
                size = 1,
                width = 0,
                position = position_nudge(x = 0)
  )+
  geom_point(data = subset(es.LatZone_Basin, !is.na(es.LatZone_Basin$meanRichness) & es.LatZone_Basin$OceanBasin == "Pacific"),
             aes(x = "Pacific Ocean",
                 y = mean(meanRichness),
                 color = 'Richness'),
             size = 5,
             shape = 17,
             position = position_nudge(x =0.25))+
  geom_errorbar(data = subset(es.LatZone_Basin, !is.na(es.LatZone_Basin$meanRichness) & es.LatZone_Basin$OceanBasin == "Pacific"),
                aes(x = "Pacific Ocean" ,
                    ymin = mean(meanRichness, na.rm = T) - (qt(0.975, df=(length(unique(ID))-1))*SE(meanRichness)),
                    ymax = mean(meanRichness, na.rm = T) + (qt(0.975, df=(length(unique(ID))-1))*SE(meanRichness)),
                    color = 'Richness'),
                size = 1,
                width = 0,
                position = position_nudge(x =0.25)
  )+
  scale_color_manual(values = c('#e0e002','#73DFFF', '#38A800'),
                     breaks=c("Abundance", "Biomass", "Richness"))+
  scale_x_discrete(limits = c("Atlantic Ocean", "Indian Ocean", "Pacific Ocean"),
                   guide = guide_axis(n.dodge = 2))+
  ylab("Log response ratio")+
  xlab(NULL)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  # ylim(c(-1,4))+
  guides(color = guide_legend(
    ncol = 1,
    keyheight = 1,
    keywidth = 1,
    reverse = F,
    title = "Metric",
    override.aes = list(shape = c(15,16,17), linetype = "blank")
  ))+
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid = element_line(color = "grey95"),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    axis.title = element_text(size = 23),
    axis.text = element_text(size = 23))+
  geom_text(data = subset(es.LatZone_Basin, !is.na(es.LatZone_Basin$meanBiomass) & es.LatZone_Basin$OceanBasin == "Atlantic"),
            aes(x="Atlantic Ocean",
                y= mean(meanBiomass, na.rm = T) + (qt(0.975, df=(length(unique(ID))-1))*SE(meanBiomass)) + 1,
                label= paste("n =", length(unique(ID)))),
            position = position_nudge(x = 0),
            size = 5)+
  geom_text(data = subset(es.LatZone_Basin, !is.na(es.LatZone_Basin$meanAbundance) & es.LatZone_Basin$OceanBasin == "Atlantic"),
            aes(x="Atlantic Ocean",
                y= mean(meanAbundance, na.rm = T) + (qt(0.975, df=(length(unique(ID))-1))*SE(meanAbundance)) + 1,
                label= paste("n =", length(unique(ID)))),
            position = position_nudge(x = -0.25),
            size = 5)+
  geom_text(data = subset(es.LatZone_Basin, !is.na(es.LatZone_Basin$meanRichness) & es.LatZone_Basin$OceanBasin == "Atlantic"),
            aes(x="Atlantic Ocean",
                y= mean(meanRichness, na.rm = T) + (qt(0.975, df=(length(unique(ID))-1))*SE(meanRichness)) + 1,
                label= paste("n =", length(unique(ID)))),
            position = position_nudge(x =0.25),
            size = 5)+
  geom_text(data = subset(es.LatZone_Basin, !is.na(es.LatZone_Basin$meanBiomass) & es.LatZone_Basin$OceanBasin == "Indian Ocean"),
            aes(x="Indian Ocean",
                y= mean(meanBiomass, na.rm = T) + (qt(0.975, df=(length(unique(ID))-1))*SE(meanBiomass)) + 1,
                label= paste("n =", length(unique(ID)))),
            position = position_nudge(x = 0),
            size = 5)+
  geom_text(data = subset(es.LatZone_Basin, !is.na(es.LatZone_Basin$meanAbundance) & es.LatZone_Basin$OceanBasin == "Indian Ocean"),
            aes(x="Indian Ocean",
                y= mean(meanAbundance, na.rm = T) + (qt(0.975, df=(length(unique(ID))-1))*SE(meanAbundance)) + 1,
                label= paste("n =", length(unique(ID)))),
            position = position_nudge(x = -0.25),
            size = 5)+
  geom_text(data = subset(es.LatZone_Basin, !is.na(es.LatZone_Basin$meanRichness) & es.LatZone_Basin$OceanBasin == "Indian Ocean"),
            aes(x="Indian Ocean",
                y= mean(meanRichness, na.rm = T) + 1,
                label= paste("n = 1")),
            position = position_nudge(x =0.25),
            size = 5)+
  geom_text(data = subset(es.LatZone_Basin, !is.na(es.LatZone_Basin$meanBiomass) & es.LatZone_Basin$OceanBasin == "Pacific"),
            aes(x="Pacific Ocean",
                y= mean(meanBiomass, na.rm = T) + (qt(0.975, df=(length(unique(ID))-1))*SE(meanBiomass)) + 1,
                label= paste("n =", length(unique(ID)))),
            position = position_nudge(x = 0),
            size = 5)+
  geom_text(data = subset(es.LatZone_Basin, !is.na(es.LatZone_Basin$meanAbundance) & es.LatZone_Basin$OceanBasin == "Pacific"),
            aes(x="Pacific Ocean",
                y= mean(meanAbundance, na.rm = T) + (qt(0.975, df=(length(unique(ID))-1))*SE(meanAbundance)) + 1,
                label= paste("n =", length(unique(ID)))),
            position = position_nudge(x = -0.25),
            size = 5)+
  geom_text(data = subset(es.LatZone_Basin, !is.na(es.LatZone_Basin$meanRichness) & es.LatZone_Basin$OceanBasin == "Pacific"),
            aes(x="Pacific Ocean",
                y= mean(meanRichness, na.rm = T) + (qt(0.975, df=(length(unique(ID))-1))*SE(meanRichness)) + 1,
                label= paste("n =", length(unique(ID)))),
            position = position_nudge(x =0.25),
            size = 5)

Figure4B

Figure4B <- plot_grid(Figure4B + theme(legend.position = "none"),
                      ncol = 1)
Figure4B

Figure4AB <- Figure4A / Figure4B
Figure4AB + plot_annotation(tag_levels = 'A')

## Figure 5 ########
### Abundance/Rugosity #########

abundance_summary <- abundance_summary %>%
  mutate(
    shape_type = ifelse(ci_lower > 0, 15, 0)
  )

panel1 <- abundance_summary %>%
  filter(Factor == "Rugosity.x") %>%
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
  ylim(-1, 3.5) +
  labs(
    x = NULL,
    y = "Log response ratio:\n abundance",
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
  
  geom_text(aes(y = ci_upper + 0.25,
                label = paste0("n = ", n)),
    size = 5,
    color = "black"
  )

panel1

### Abundance/Vertical Relief #########

panel2 <- abundance_summary %>%
  filter(Factor == "Vertical.relief.x") %>%
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
  ylim(-1, 3.5) +
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
  
  geom_text(aes(y = ci_upper + 0.25,
                label = paste0("n = ", n)),
            size = 5,
            color = "black"
  )

panel2

### Abundance/Porosity #########

panel3 <- abundance_summary %>%
  filter(Factor == "Openings.x") %>%
  ggplot(aes(x = Level, y = mean, color = Level)) +
  geom_point(aes(shape = as.factor(shape_type)), size = 5, stroke = 1) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("No openings" = "#e0e002",
                                "With openings" = "#e0e002")) +
  scale_shape_manual(values = c("15" = 15, "0" = 0),
                     labels = c("No openings" = "No Openings", "With openings" = "With Openings")) +  # map numeric to shape
  scale_x_discrete(limits = c("No openings", "With openings"),
                   labels = scales::label_wrap(5)(c("No Openings", "With Openings"))) +
  ylim(-1, 3.5) +
  labs(
    x = NULL,
    y = NULL,
    color = "Openings",
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
  
  geom_text(aes(y = ci_upper + 0.25,
                label = paste0("n = ", n)),
            size = 5,
            color = "black"
  )

panel3

row1 <- plot_grid(
  panel1+ theme(axis.text.x=element_blank()),
  panel2+ theme(axis.text.x=element_blank()),
  panel3+ theme(axis.text.x=element_blank()),
  ncol = 3,
  nrow = 1)

### Biomass/Rugosity #########

biomass_summary <- biomass_summary %>%
  mutate(
    shape_type = ifelse(ci_lower > 0, 16, 1)
  )

panel4 <- biomass_summary %>%
  filter(Factor == "Rugosity.x") %>%
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
  ylim(-5,12) +
  labs(
    x = NULL,
    y = "Log response ratio:\n biomass",
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
  geom_text(aes(y = ci_upper + 0.75,
                label = paste0("n = ", n)),
            size = 5,
            color = "black"
  )

panel4

### Biomass/Vertical Relief #########

panel5 <- biomass_summary %>%
  filter(Factor == "Vertical.relief.x") %>%
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
  ylim(-5,12) +
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
  geom_text(aes(y = ci_upper + 0.75,
                label = paste0("n = ", n)),
            size = 5,
            color = "black"
  )

panel5

### Biomass/Porosity #########

panel6 <- biomass_summary %>%
  filter(Factor == "Openings.x") %>%
  ggplot(aes(x = Level, y = mean, color = Level)) +
  geom_point(aes(shape = as.factor(shape_type)), size = 5, stroke = 1) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("No openings" = "#73DFFF",
                                "With openings" = "#73DFFF")) +
  scale_shape_manual(values = c("16" = 16, "1" = 1),
                     labels = c("No openings" = "No Openings", "With openings" = "With Openings")) +  # map numeric to shape
  scale_x_discrete(limits = c("No openings", "With openings"),
                   labels = scales::label_wrap(5)(c("No Openings", "With Openings"))) +
  ylim(-5,12) +
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
  geom_text(aes(y = ci_upper + 0.75,
                label = paste0("n = ", n)),
            size = 5,
            color = "black"
  )

panel6

row2<- plot_grid(
  panel4+ theme(axis.text.x=element_blank()),
  panel5+ theme(axis.text.x=element_blank()),
  panel6+ theme(axis.text.x=element_blank()),
  ncol = 3,
  nrow = 1)

### Richness/Rugosity #########

richness_summary <- richness_summary %>%
  mutate(
    shape_type = ifelse(ci_lower > 0, 17, 2)
  )

panel7 <- richness_summary %>%
  filter(Factor == "Rugosity.x") %>%
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
  ylim(-1, 2) +
  labs(
    x = NULL,
    y = "Log response ratio:\n richness",
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
  geom_text(aes(y = ci_upper + 0.25,
                label = paste0("n = ", n)),
            size = 5,
            color = "black"
  )

panel7

### Richness/Vertical Relief #########

panel8 <- richness_summary %>%
  filter(Factor == "Vertical.relief.x") %>%
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
  ylim(-1, 2) +
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
  geom_text(aes(y = ci_upper + 0.25,
                label = paste0("n = ", n)),
            size = 5,
            color = "black"
  )

panel8

### Richness/Porosity #########

panel9 <- richness_summary %>%
  filter(Factor == "Openings.x") %>%
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
  ylim(-1,2) +
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
  geom_text(aes(y = ci_upper + 0.25,
                label = paste0("n = ", n)),
            size = 5,
            color = "black"
  )

panel9

row3 <- plot_grid(panel7,
                  panel8,
                  panel9,
                  nrow = 1)

plot_grid(panel1 + theme(axis.text.x=element_blank()),
          panel2 + theme(axis.text.x=element_blank()),
          panel3 + theme(axis.text.x=element_blank()),
          panel4 + theme(axis.text.x=element_blank()),
          panel5 + theme(axis.text.x=element_blank()),
          panel6 + theme(axis.text.x=element_blank()),
          panel7,
          panel8,
          panel9,
          ncol = 3,
          nrow = 3, 
          labels = "AUTO") 


