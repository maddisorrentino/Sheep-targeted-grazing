# Mineral nitrogen data cleaning and modeling 
# This R script contains code to clean mineral n data for use in models and visualizations
# Code by MGS, KAH, ASW

#install tidyverse
library(tidyverse)
library(ggplot2)
library(lme4)
library(glmmTMB)
library(performance)
library(marginaleffects)
library(sjPlot)
library(gridExtra)
library(plotrix)
library(ggeffects)

setwd("/Users/madelynsorrentino/Documents/GitHub/Sheep-targeted-grazing")

######################################################
########### Bring in mineral n datasheet #############
######################################################

## 2022 and 2023 data
nmin <- read.csv("Nmin_2022_2023.csv")

#first, remove columns we will not be using
#per notes from Arden Engel who analyzed the mineral n (ardenengel@boisestate.edu), standard curve 1 was the best out of the three standard curves run, so for analyses, we will use the mineral n numbers generated in standard curve 1 
#want to keep columns for site, plot, year, treatment, sample id, mg_NO3_per_g_stcurv1, mg_NH4_per_g_stcurv1 
#important to note there are lots of other columns of data used to come up with the final value in mg/g - refer back to raw data file to look at these numbers! 

nmin2 <- nmin %>%
  select(c("Site","Plot","Treatment","Year","mg_NO3_per_g_stcurv1","mg_NH4_per_g_stcurv1")) 

str(nmin2)
head(nmin2)

## bring in 2024 data
nmin_2024 <- read.csv("Nmin_2024.csv")

## select data we are using to match that of the 2022 and 2023 spreadsheet
## read notes above; applies to this dataset as well

nmin_2024_2 <- nmin_2024 %>%
  select(c("site", "plot", "treatment", "year", "ug_g_NH4", "ug_g_NO3"))

## rename variables to make consistent across both datasheets

nmin3 <- nmin2 %>%
  rename("nh4_mg_g" = "mg_NH4_per_g_stcurv1", "no3_mg_g" = "mg_NO3_per_g_stcurv1")

nmin_2024_3 <- nmin_2024_2 %>%
  rename("Site" = "site", "Plot" = "plot", "Treatment" = "treatment", "Year" = "year", "nh4_mg_g" = "ug_g_NH4", "no3_mg_g" = "ug_g_NO3")

## filter out values with NA year from 2024 dataframe

nmin_2024_3 <- nmin_2024_3 %>%
  filter(Year == 2024)

## make year, trt into factors
nmin_2024_3$Year <- as.factor(nmin_2024_3$Year)
nmin_2024_3$Treatment <- as.factor(nmin_2024_3$Treatment)

## combine 2022-2023 and 2024 dataframes

nmin4 <- rbind(nmin3, nmin_2024_3)

## check structure of combined dataframe
str(nmin4)

nmin4$Year <- as.factor(nmin4$Year)
nmin4$Treatment <- as.factor(nmin4$Treatment)
nmin4$no3_mg_g <- as.numeric(nmin4$no3_mg_g)
nmin4$nh4_mg_g <- as.numeric(nmin4$nh4_mg_g)


## aggregate data by plot level for analyses 

nmin5 <- nmin4 %>%
  group_by(Plot, Year) %>%
  mutate(NO3_plot = mean(no3_mg_g, na.rm=TRUE),
         NH4_plot = mean(nh4_mg_g, na.rm=TRUE)) %>%
  distinct(Site, Plot, Treatment, Year, NO3_plot, NH4_plot)


## add summer, spring, fall as factors
nmin6 <- nmin5 %>%
  transform(Fall =
              ifelse(Treatment == "FA" | Treatment == "FASP", 1, 0)) %>%
  transform(Spring =
              ifelse(Treatment == "SP" | Treatment == "FASP", 1, 0)) %>%
  transform(Summer = 
              ifelse(Treatment == "SU",1,0))

## make ug per g into mg per kg for analysis 
nmin6$no3_plot_mgkg <- nmin6$NO3_plot*1000
nmin6$nh4_plot_mgkg <- nmin6$NH4_plot*1000

## raw data mean all plots 2022
nmin_2022 <- nmin6 %>%
  filter(Year == 2022)

mean(nmin_2022$nh4_plot_mgkg)
std.error(nmin_2022$nh4_plot_mgkg)

mean(nmin_2022$no3_plot_mgkg)
std.error(nmin_2022$no3_plot_mgkg)

nmin6$Treatment <- factor(nmin6$Treatment, levels = c("SU", "FA", "SP", "FASP"))


#################################################
################## glm of no3  ##################
#################################################

##### gamma model #####
NO3_mod <- glmmTMB(no3_plot_mgkg~Fall*Spring*Year + (1|Plot), data=nmin6,
                   family=Gamma(link="log"))
summary(NO3_mod)
performance(NO3_mod)

## check residuals 
pearson_res_no3 <- residuals(NO3_mod, type = "pearson")
plot(fitted(NO3_mod), pearson_res_no3,
     xlab = "Fitted Values", ylab = "Pearson Residuals",
     main = "Pearson Residuals vs. Fitted Values")
abline(h = 0, col = "red")

##### test using brms bayesian model to help diagnose #####
library(brms)

no3_test <- brm(no3_plot_mgkg~Fall*Spring*Year + (1|Plot), family = lognormal(link = "identity"), data=nmin6)
summary(no3_test)
pp_check(no3_test)

##### test model with log normal distribution #####
NO3_mod_lognormal <- glmmTMB(no3_plot_mgkg~Fall*Spring*Year + (1|Plot), data=nmin6,
                             family=lognormal(link = "identity"))
summary(NO3_mod_lognormal) ## INSANE effect sizes?? Likely due to massive variance in 2023 and much larger numbers?
performance(NO3_mod_lognormal)

## check residuals
## note cannot calculate person residuals for this type of model 
plot(
  fitted(NO3_mod_lognormal),
  residuals(NO3_mod_lognormal, type = "response"),
  xlab = "Fitted values",
  ylab = "Response residuals")
abline(h = 0, lty = 2)

## another way
library(DHARMa)
sim_res <- simulateResiduals(
  fittedModel = NO3_mod_lognormal,
  n = 1000)

plot(sim_res)

## marginal effects code 

no3_new <- datagrid(model=NO3_mod,
                    Fall=c(0, 1), 
                    Spring=c(0, 1))
## Row 1 of this grid = controls (0,0)
## Row 2 of this grid = treatment2, but not trt1 (0,1)
## Row 3 = treatment 1, but not 2 (1, 0)
## Row 4 = treatment 1 and treatment 2 (1,1)....

no3_compare <- comparisons(
  NO3_mod,
  variables = "Year",
  # above = calculate the difference in the effect of change over time for each combo in the data grid created below...
  newdata = no3_new,
  ## These, in the resulting contrasts, will repeat twice, one for each post-treatment year
  comparison = "difference", ## Calculate contrast as a difference
  type = "response",
  vcov = TRUE) ## (once you move to binomial glm, you are going to want to specify that you'd like it on scale of response too, I think).

## old version where we specify hypothesis = "revpairwise" was deprecated, so now we need to pass our comparisons to the hypotheses function. "B" is now "row", this information should still apply! 

no3_hyp <- hypotheses(no3_compare, hypothesis = ~ -pairwise)

no3_compare2 <- no3_hyp[c(1:3, 23:25),] ## select the pairs we want 

## We can rename our 
no3_compare2$treateff <- c("Spring", "Fall", "Fall*Spring",
                           "Spring", "Fall", "Fall*Spring")

no3_compare2$year <- c(rep("2023", 3), rep("2024", 3))

#reorder treatments for plotting
no3_compare2$treateff <- factor(no3_compare2$treateff, levels = c("Fall", "Spring", "Fall*Spring"))

#colors for treatments
three_colors <- c("#e9c46a","#2a9d8f","#264653")
colors_usfs_022025 <- c("#71467f", "#7e8133", "#0d8182")
four_usfs_022025 <- c("#db9d27","#71467f", "#7e8133", "#0d8182")

no3_me_plot <- ggplot(no3_compare2, aes(x=year, y=estimate, color=treateff)) +
  geom_pointrange(aes(ymin=conf.low, ymax=conf.high), position= position_dodge2(
    width = 0.5), size = 0.5) + 
  geom_hline(yintercept=0, linetype="dashed") +
  scale_color_manual(values = colors_usfs_022025, name = "Grazing treatment", labels = c("Fall", "Spring", "Fall+Spring")) +
  ylab("Estimated effect of targeted grazing on nitrate (mg/kg)\n(controlling for time/group differences)")+
  xlab("Year") +
  theme_bw() +
  theme(strip.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, color = "black", face = "bold"), # changes size of axis title
        axis.text = element_text(size = 12, color = "black"), # changes size of axis text
        legend.title = element_text(size=12, color = "black", face = "bold"), # changes size of legend title
        legend.text = element_text(size=12, color = "black")) # changes size of legend text )

no3_me_plot

#####################################
## other ways to model these data ##
#####################################

hist(nmin6$no3_plot_mgkg)

## positive, right skewed data - typical for a gamma glm, which does not assume equal variance 
## Box–Cox alternative to Gamma GLMM for NO3

library(MASS)
library(forecast)

##### Estimate Box–Cox lambda (Gaussian, fixed effects only) #####
bc_no3 <- boxcox(
  lm(no3_plot_mgkg ~ Fall * Spring * Year, data = nmin6),
  plotit = FALSE
)

lambda_no3 <- bc_no3$x[which.max(bc_no3$y)]

## Apply Box–Cox transformation
nmin6$no3_bc <- BoxCox(nmin6$no3_plot_mgkg, lambda = -0.2)

## Fit Gaussian mixed model on transformed response
NO3_mod_bc <- glmmTMB(
  no3_bc ~ Fall * Spring * Year + (1 | Plot),
  data = nmin6,
  family = gaussian())

## Model summary
summary(NO3_mod_bc)
performance(NO3_mod_bc)


##### final box cox workflow  #####
bc_no3 <- boxcox(
  lm(no3_plot_mgkg ~ Fall * Spring * Year, data = nmin6),
  plotit = TRUE)

lambda_no3 <- bc_no3$x[which.max(bc_no3$y)]
lambda_no3


## Apply Box–Cox transformation
nmin6$no3_bc <- BoxCox(nmin6$no3_plot_mgkg, lambda = lambda_no3)


## check if this made it more normal
hist(nmin6$no3_bc)
shapiro.test(nmin6$no3_bc)
## still not perfect, but a lot better


## now fit a model using the box-cox transformed data
## https://rvlenth.github.io/emmeans/reference/make.tran.html

# Fit a model using an oddball transformation:
bctran <- make.tran("boxcox", lambda_no3)

NO3_mod_bc <- with(bctran, 
                   glmmTMB(
                     linkfun(no3_plot_mgkg) ~ Fall * Spring * Year + (1 | Plot),
                     data = nmin6,
                     family = gaussian())) ## how to apply bias adjustment to this model?
summary(NO3_mod_bc)
performance(NO3_mod_bc)

sim <- simulate(NO3_mod_bc, nsim = 1000)

# back-transform each simulation
sim_back <- lapply(sim, bctran$linkinv)

mean_pred <- Reduce("+", sim_back) / length(sim_back)

resid_bc  <- residuals(NO3_mod_bc, type = "response")
fitted_bc <- fitted(NO3_mod_bc)
ggplot(data.frame(fitted_bc, resid_bc),
       aes(x = fitted_bc, y = resid_bc)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "Fitted values (Box–Cox scale)",
    y = "Residuals (Box–Cox scale)",
    title = "Residuals vs Fitted"
  ) +
  theme_bw()


# Obtain back-transformed LS means:    
em_no3 <- emmeans(NO3_mod_bc,
                  ~ Fall * Spring * Year,
                  type="response")

em_no3 <- emmeans(
  NO3_mod_bc,
  ~ Fall * Spring * Year
) %>%
  regrid(transform = "response",
         bias.adjust = TRUE) ## apply bias adjustment 

## define control 
trt_vs_summer <- contrast(
  em_no3,
  method = list(
    "Fall"   = c(-1, 1, 0, 0),
    "Spring" = c(-1, 0, 1, 0),
    "Fall + Spring"   = c(-1, 0, 0, 1)),
  by = "Year")

## calculate "contrasts" compared to control 
did <- contrast(
  trt_vs_summer,
  method = "revpairwise",
  by = "contrast")

did_sum <- summary(did, infer = c(TRUE, TRUE))

diddf <- as.data.frame(did_sum) %>%
  filter(contrast1 != "Year2024 - Year2023") %>%
  separate(contrast1,
           into = c("year_left", "year_right"),
           sep = " - ") %>%
  mutate(
    Year = gsub("Year", "", year_left))

## try plotting back transformed values?
colors <- c("#71467f", "#7e8133", "#0d8182")

no3_bc_plot <- ggplot(diddf, aes(x=Year, y=estimate, color=contrast)) +
  geom_pointrange(aes(ymin=lower.CL, ymax=upper.CL), position= position_dodge2(
    width = 0.5), size = 0.5, show.legend = FALSE) + 
  geom_hline(yintercept=0, linetype="dashed") +
  scale_color_manual(values = colors, name = "Grazing treatment", labels = c("Fall", "Spring", "Fall+Spring")) +
  ylab("")+
  xlab("") +
  theme_bw() +
  theme(axis.title = element_text(size = 10, color = "black", face = "bold"),
        axis.text = element_text(size = 10, color = "black"), 
        legend.title = element_text(size=10, color = "black", face = "bold"),
        legend.text = element_text(size=10, color = "black"),
        aspect.ratio = 2/1) 

no3_bc_plot

