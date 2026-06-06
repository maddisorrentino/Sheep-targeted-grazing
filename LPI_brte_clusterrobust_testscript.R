#### GLM and visuals of grazing treatment effects on fungroup cover ####
### Code by MGS, KAH, ASW 

library(tidyverse)
library(lme4)
library(marginaleffects)
library(performance)
library(glmmTMB)
library(emmeans)
library(ggeffects)
library(ggstats)
library(GLMMadaptive)
library(ggstance)
library(broom.mixed)
library(gridExtra)
library(dplyr)
library(plotrix)


## set working directory
#change to your directory as needed 

#Maddi:
setwd("/Users/madelynsorrentino/Documents/GitHub/Sheep-targeted-grazing")

## load data for plots across years 
all_lpi_2022 <- read.csv("percent_cover_all_plots_with_species_info_2022.csv")

head(all_lpi_2022)
str(all_lpi_2022)

all_lpi_2023 <- read.csv("percent_cover_all_plots_with_species_info_2023.csv")

head(all_lpi_2023)
str(all_lpi_2023)

all_lpi_2024 <- read.csv("percent_cover_all_plots_with_species_info_2024.csv")

head(all_lpi_2024)
str(all_lpi_2024)


## select columns needed from both years and stack into one dataframe
all_2022 <- all_lpi_2022 %>%
  dplyr::select("Site","Plot","Trt","Year","species","count","functional_group", "native_status") %>%
  as.data.frame()
head(all_2022)
str(all_2022)


all_2023 <- all_lpi_2023 %>%
  dplyr::select("Site","Plot","Trt","Year","species","count","functional_group", "native_status") %>%
  as.data.frame()
head(all_2023)
str(all_2023)

all_2024 <- all_lpi_2024 %>%
  dplyr::select("Site","Plot","Trt","Year","species","count","functional_group", "native_status") %>%
  as.data.frame()
head(all_2024)
str(all_2024)


#combine data and stack them 
all_lpi <- rbind(all_2022, all_2023, all_2024)
head(all_lpi)
tail(all_lpi)


## remove seed bank plots, dead shrubs, dead forbs, dead grasses
all_lpi2 <- filter(all_lpi, Plot < 33 & species != "DS" & species != "DF" & species != "DG")


## add factors for fall, spring, summer grazing
all_lpi3 <- all_lpi2 %>%
  transform(Fall =
              ifelse(Trt == "FA" | Trt == "FASP", 1, 0)) %>%
  transform(Spring =
              ifelse(Trt == "SP" | Trt == "FASP", 1, 0)) %>%
  transform(Summer = 
              ifelse(Trt == "SU",1,0))
head(all_lpi3)

#make each treatment column and year into factors
all_lpi3$Fall <- as.factor(all_lpi3$Fall)
all_lpi3$Spring <- as.factor(all_lpi3$Spring) 
all_lpi3$Summer <- as.factor(all_lpi3$Summer)
all_lpi3$Year <- as.factor(all_lpi3$Year) 


#make new treatment column to add in ungrazed plots in 2022 and 2023
all_lpi3$Grazed <- 1

#make "Grazed" 0 for plots 29, 30, 31, 32, 6 in 2022 and 2023, NOT 2024
all_lpi4 <- all_lpi3 %>% 
  mutate(Grazed = replace(Grazed, Plot == 29 & Year != 2024, 0)) %>%
  mutate(Grazed = replace(Grazed, Plot == 30 & Year!= 2024, 0)) %>% 
  mutate(Grazed = replace(Grazed, Plot == 31 & Year != 2024, 0)) %>% 
  mutate(Grazed = replace(Grazed, Plot == 32 & Year != 2024, 0)) %>%
  mutate(Grazed = replace(Grazed, Plot == 6 & Year != 2024, 0))

#make a new functional group column "functional_group_new"

all_lpi4a <- all_lpi4 %>%
  transform(functional_group_new = 
              ifelse(functional_group == "forb" & native_status == "native", "native forb",
                     ifelse(functional_group == "forb" & native_status == "introduced", "introduced forb", functional_group)))

## create new column for amount of "not" hits of 150 - this is the number of times along our transects we did NOT hit a fungroup, and group by fungroup
all_lpi5 <- all_lpi4a %>%
  group_by(Plot,Year,functional_group_new, Trt, Grazed) %>%
  dplyr::summarise(functional_group_cover = sum(count)) %>%
  mutate(nohit = 150-functional_group_cover) %>%
  as.data.frame()

#add treatment factors back? idk why this got rid of them - need to work on this code
all_lpi6 <- all_lpi5 %>%
  transform(Fall =
              ifelse(Trt == "FA" | Trt == "FASP", 1, 0)) %>%
  transform(Spring =
              ifelse(Trt == "SP" | Trt == "FASP", 1, 0)) %>%
  transform(Summer = 
              ifelse(Trt == "SU",1,0))
head(all_lpi6)

#make year, seasons into factors
all_lpi6$Fall <- as.factor(all_lpi6$Fall)
all_lpi6$Spring <- as.factor(all_lpi6$Spring) 
all_lpi6$Summer <- as.factor(all_lpi6$Summer)
all_lpi6$Year <- as.factor(all_lpi6$Year) 

#reorder so summer is first
all_lpi6$Trt <- factor(all_lpi6$Trt, levels=c("SU","FA", "SP", "FASP"))
str(all_lpi6)

##### BRTE cover mod #####

lpi_brte <- all_lpi6 %>%
  filter(functional_group_new == "cheatgrass")

#make binomial glm comparing each treatment to ungrazed control for the forb functional group 
#use "nohit" column as the number of "failures"
#model with new structure - year interacts with each treatment - now we have 5: fall, spring, fasp, summer, ungrazed
brte_cover_mod <- glmer(cbind(functional_group_cover, nohit) ~ Fall*Spring*Year + Grazed + (1|Plot), data = lpi_brte, family="binomial")
#did not converge?

summary(brte_cover_mod)

#rerun model without grazed factor for comparison
brte_cover_mod2 <- glmer(cbind(functional_group_cover, nohit) ~ Fall*Spring*Year + (1|Plot), data = lpi_brte, family="binomial")

summary(brte_cover_mod2)
performance(brte_cover_mod2)

##### marginal effects code for BRTE model #####

brte_new <- datagrid(model=brte_cover_mod2,
                     Fall=c(0, 1), 
                     Spring=c(0, 1))
## Row 1 of this grid = controls (0,0)
## Row 2 of this grid = treatment2, but not trt1 (0,1)
## Row 3 = treatment 1, but not 2 (1, 0)
## Row 4 = treatment 1 and treatment 2 (1,1)....

brte_compare <- comparisons(
  brte_cover_mod2,
  variables = "Year",
  # above = calculate the difference in the effect of change over time for each combo in the data grid created below...
  newdata = brte_new,
  ## These, in the resulting contrasts, will repeat twice, one for each post-treatment year
  comparison = "difference", ## Calculate contrast as a difference
  type = "response",
  vcov = TRUE) 

## old revpairwise option was depricated, this approximates the same workflow 
brte_hyp <- hypotheses(brte_compare, hypothesis = ~ -pairwise)
brte_compare2 <- brte_hyp[c(1:3, 23:25),] ## select pairs we care about 

## We can rename our 
brte_compare2$treateff <- c("Spring", "Fall", "Fall*Spring",
                            "Spring", "Fall", "Fall*Spring")

brte_compare2$year <- c(rep("2023", 3), rep("2024", 3))

#reorder treatments for plotting
brte_compare2$treateff <- factor(brte_compare2$treateff, levels = c("Fall", "Spring", "Fall*Spring"))


brte_me_plot <- ggplot(brte_compare2, aes(x=year, y=estimate*100, color=treateff)) +
  geom_pointrange(aes(ymin=conf.low*100, ymax=conf.high*100), position= position_dodge2(
    width = 0.5), size = 0.5, show.legend = FALSE) + 
  geom_hline(yintercept=0, linetype="dashed") +
  scale_color_manual(values = colors_usfs_022025, name = "Grazing treatment", labels = c("Fall", "Spring", "Fall x Spring")) +
  ylab("Estimated effect of targeted grazing on cheatgrass cover (%)\ncompared to summer control plots")+
  xlab("Year") +
  theme_light() +
  theme(strip.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 10, color = "black", face = "bold"), # changes size of axis title
        axis.text = element_text(size = 10, color = "black"), # changes size of axis text
        legend.title = element_text(size=10, color = "black", face = "bold"), # changes size of legend title
        legend.text = element_text(size=10, color = "black"),
        panel.spacing = unit(0.5, "cm"), 
        aspect.ratio = 1.5)

brte_me_plot

##### test of cluster robust standard errors #####
library(clubSandwich)
library(lmtest)

X <- model.matrix(brte_cover_mod2)
cl <- lpi_brte$Year
V_year <- sandwich::vcovCL(
  brte_cover_mod2,
  cluster = cl,
  type = "HC0")

vcov_year <- vcovCR(brte_cover_mod2, cluster = lpi_brte$Year, type = "CR2")

coef_test(brte_cover_mod2, vcov = vcov_year)

vcov_year <- vcovCR(
  brte_cover_mod2,
  cluster = lpi_brte$Year,
  type = "CR2")

coef_test(brte_cover_mod2, vcov = vcov_year)

## none of these seem to work with glmmtmb models