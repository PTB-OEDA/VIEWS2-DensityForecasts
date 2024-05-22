# Builds off of the Brandt-ViEWS2.Rmd document to do forecasts for 2022+
# 
# Patrick T. Brandt
#
#
# 20240514 : Initial version from combined-rolled=2022.R
#            Generates forecasts for months 503+ which cover 2023-2024.
#
#

library(arrow)
library(dplyr)

#### Read data ####
# Read in the various training datasets
cmbn_df <- read_parquet("~/VIEWS2/shared_competition_data/cm_features.parquet")
actuals <- read_parquet("~/VIEWS2/shared_competition_data/cm_actuals_2022.parquet")  

#### Simple plots as checks ####

# Simple summary plot to see that we have things correct...
ms <- sort(unique(cmbn_df$month_id))

roll.ged.sb <- matrix(NA, nrow=length(ms), ncol=3)
sbdata <- vector(mode="list", length=length(ms))

for(i in 1: length(ms))
{
  sb <- cmbn_df[cmbn_df$month_id==ms[i],4] 
  sbdata[[i]] <- as.vector(unlist(sb))
  roll.ged.sb[i,] <- c(ms[i], mean(sb$ged_sb), sd(sb$ged_sb))
}

rm(sb)

plot(roll.ged.sb[,1:2], type="l", 
     xlab = "Month", ylab = "Mean GED SB events")
plot(roll.ged.sb[,c(1,3)], type="l", xlab = "Month", 
     ylab = "Standard Deviation of GED SB events")

roll.ged.sb <- ts(roll.ged.sb, start=c(1990,1), freq=12)
colnames(roll.ged.sb) <- c("month_id", "Mean SB", "Std.Dev. SB")

plot(roll.ged.sb[,2:3], lwd=2, col=1:2, main="", cex.lab=0.9, cex.axis=0.7)

#### Packages and data formatting ####

# Packages for estimation
library(mgcv)
library(glmmTMB)

# Packages require a df not a tbl!!!!  Should be prior data
dt1 <- as.data.frame(cmbn_df[,c("month_id", "country_id", "ged_sb", "ged_sb_tlag_1")])
dt1$country_id <- as.factor(dt1$country_id)
rm(cmbn_df)

#### Truncation functions ####
# Truncate the NA overflows and the large values at upper.tol value
# This automates so we can look at sensitivity later

trunctop <- function (x, upper.tol=1e9)
{
  x <- replace(x, is.na(x), upper.tol)  # NA and NaN cases
  x <- replace(x, x>upper.tol, upper.tol)  # Truncate the top
  return(x)
}


#### Best models as selected from Berlin and past AIC results ####
# See details: https://github.com/PTB-OEDA/VIEWS2-DensityForecasts

#### Tweedie GAM with local effects ####

tw.local <- bam(ged_sb ~ ged_sb_tlag_1 + s(month_id, country_id, bs="fs"),
                data=dt1,
                family = tw(),                 
                discrete = TRUE,
                control = gam.control(trace=TRUE))

# Plot effects
plot(tw.local, page=1)
abline(h=0, col="red", lwd=2)

#### GLMMs ####

# Need time as a factor for the next commands...
dt1$month_id <- factor(dt1$month_id)

options(glmmTMB.cores=4)
psnar1.glmmtmb <- glmmTMB(ged_sb ~ ar1(month_id + 0|country_id), family = poisson,
                          data=dt1)

nbar1.glmmtmb <- glmmTMB(ged_sb ~ ar1(month_id + 0|country_id), family = nbinom1(),
                         data=dt1)

twar1.glmmtmb <- glmmTMB(ged_sb ~ ar1(month_id + 0|country_id), family = tweedie(),
                         data=dt1)

#### AIC Fit Statistics ####
library(bbmle)
AICtab(tw.local, psnar1.glmmtmb, nbar1.glmmtmb, twar1.glmmtmb)

#### Predictions ####

# TW local
tw.local.p <- cbind(dt1[,1:2],
                    predict.bam(tw.local, type = "response"))

colnames(tw.local.p) <- c("month_id", "country_id", "TW")
 
# GLMM preds
glmm.preds <- cbind(dt1[,1:2],
                    predict(psnar1.glmmtmb, type = "response"),
                    predict(nbar1.glmmtmb, type = "response"),
                    predict(twar1.glmmtmb, type = "response"))
colnames(glmm.preds) <- c( "month_id", "country_id", "P", "NB", "TW")

# Get last obs for comparison of an "in-sample"
lastobs <- dt1[dt1$month_id=="502",]
local.last <- tw.local.p[tw.local.p$month_id=="502",]
glmm.last <- glmm.preds[glmm.preds$month_id==502,]

# Make draws from the relevant modal predictions

# Local preds in-sample pdf
N <- 1000;    # Number of draws
n <- dim(lastobs)[1]
k <-  1;      #  number of forecasts

set.seed(1234)

# Local GAM Tweedie
local.TW.fc <- t(replicate(N, rTweedie(trunctop(local.last$TW), p=1.581)))

# GLMM preds in-sample pdf
glmm.P.fc <- sapply(1:n, function(i) {rpois(N, glmm.last$P[i])})
glmm.NB.fc <- sapply(1:n, function(i) {rnbinom(N, size=summary(nbar1.glmmtmb)$sigma, 
                                               mu=glmm.last$NB[i])})
glmm.TW.fc <- t(replicate(N, rTweedie(trunctop(glmm.last$TW), p=1.36)))

# Below is a vastly simplified version of rhe combined-rolled-* code in the earlier
# analyses.

nsteps <- 27 # Number of months ahead!
local.TW <- vector(mode = "list", length=nsteps)

k1 <- 502  # Starting time period = max(dt1$month_id)

for(i in 1:nsteps)
{
  # Local Tweedie model
  old.tw <- apply(local.TW.fc, 2, mean)
  tmp.tw <- predict(tw.local, data.frame(country_id=local.last$country_id,
                                         month_id = rep((k1+i), 191),
                                         ged_sb_tlag_1 = old.tw),
                    type = "response")
  local.TW.fc <- t(replicate(N, rTweedie(trunctop(tmp.tw), p=1.581)))

  local.TW[[i]] <- local.TW.fc
}


# glmms -- these are hard to do, since they do not easily admit more time

glmm.P <- glmm.NB <- glmm.TW <- vector(mode = "list", length=nsteps)

for(i in 1:nsteps)
{
  # glmm models
  old.p <- apply(glmm.P.fc, 2, mean)
  tmp.p <- predict(psnar1.glmmtmb, newdata=data.frame(country_id=glmm.last$country_id,
                                                      month_id = as.factor(rep((k1+i), 191)),
                                                      ged_sb_tlag_1 = old.p),
                   type = "response", allow.new.levels=TRUE)
  glmm.P.fc <- sapply(1:n, function(i) {rpois(N, tmp.p[i])})

  old.nb <- apply(glmm.NB.fc, 2, mean)
  tmp.nb <- predict(nbar1.glmmtmb, data.frame(country_id=glmm.last$country_id,
                                              month_id = as.factor(rep((k1+i), 191)),
                                              ged_sb_tlag_1 = old.nb),
                    type = "response", allow.new.levels=TRUE)

  glmm.NB.fc <- sapply(1:n, function(i) {rnbinom(N, size=summary(nbar1.glmmtmb)$sigma, 
                                                 mu=tmp.nb[i])})

  old.tw <- apply(glmm.TW.fc, 2, mean)
  tmp.tw <- predict(twar1.glmmtmb, data.frame(country_id=glmm.last$country_id,
                                              month_id = as.factor(rep((k1+i), 191)),
                                              ged_sb_tlag_1 = old.tw),
                    type = "response", allow.new.levels=TRUE)
  glmm.TW.fc <- t(replicate(N, rTweedie(trunctop(tmp.tw), p=1.581)))

  glmm.P[[i]] <- glmm.P.fc
  glmm.NB[[i]] <- glmm.NB.fc
  glmm.TW[[i]] <- glmm.TW.fc
}

# Little bit of clean up
rm(old.p, old.nb, old.tw, 
   local.TW.fc, 
   glmm.P.fc, glmm.NB.fc, glmm.TW.fc)

#### Forecasts setup ####
# Now make some comparisons against the 2022 data for months
table(actuals$month_id)

# compared to
range(as.numeric(levels(dt1$month_id)))

# f = forecast list
# ds = data frame of the actuals
# freq = 12; data frequency or number of forecast periods
# N = number of forecast samples per unit time
# mname = model name in quotes

convertforc <- function(f, ds, N, freq = 12, model = "M")
{
  ds <- as.data.frame(ds)
  maxt <- max(ds$month_id)
  mint <- min(ds$month_id)
  nc <- dim(f[[1]])[2]
  all <- data.frame("country_id" = rep(rep(ds[ds$month_id==maxt,]$country_id, each=N), freq), 
                    "prediction" = unlist(lapply(f, "as.vector")),
                    "sample" = 1:N,
                    "month_id" = rep(mint:maxt, each=N*nc),
                    "model" = as.factor(model))
  
  # Merge on the actuals as "true_value" here
  all <- merge(all, data.frame("month_id" = ds[,1], 
                               "country_id" = ds[,2], 
                               "true_value" = ds[,3]))
  return(all)
}

# Convert year of forecasts for 2022

# 2022 alignment:
# Remember, we are not given the last 3 months of 2021, so this needs to drop
# those from the forecasting, even though we generated them and more

t1 <- convertforc(local.TW[4:15], actuals, N, 12, model="Tweedie GAM")
t2 <- convertforc(glmm.P[4:15], actuals, N, 12, model="Poisson GLMM")
t3 <- convertforc(glmm.NB[4:15], actuals, N, 12, model="Neg Binom GLMM")
t4 <- convertforc(glmm.TW[4:15], actuals, N, 12, model="Tweedie GLMM")

forcs22 <- rbind(t1,t2,t3,t4)
rm(t1,t2,t3,t4); gc()

# Add a year to each for later: prolly should have the dates
# better organized above.
forcs22$year <- 2022
forcs22$month <- forcs22$month_id-504
  
#### Forecast Performance ####
library(scoringutils)
library(magrittr)
library(dplyr)

# Adjust here for the latest version of the scoringutils package and add
# country names
country.labels <- read.csv("~/VIEWS2/matching_tables/countries.csv")[,c(1,2,5)]
country.labels$country_id <- country.labels$id
forcs22 <- merge(forcs22, country.labels[,2:4])

colnames(forcs22) <- c("country_id", "month_id", "predicted",
                       "sample_id", "model", "observed", "year", "month", 
                       "location", "isoab")

scored.out <- score(as_forecast(forcs22, 
                                forecast_unit = c("model", "month", "isoab")))

crps22 <- scored.out %>% summarise_scores(by=c("model", "month"), na.rm=TRUE) %>% 
  select(c("model", "month", "crps"))

# CRPS by model for 2022
by(crps22[,3], crps22[,1], sum)/12

# Now, look just at the 3 GLMM models...

crps.rel <- scored.out %>% summarise_scores(by=c("model", "month"), 
                                                  relative_skill=TRUE,
                                                  baseline = "Tweedie GLMM") 
# %>% 
#   select(c("model", "month", "scaled_rel_skill"))

# Look at the measures....
tmp <- add_relative_skill(scored.out, by=c("model", "month"), 
                   metric="crps", baseline="Tweedie GLMM")

tmp1 <- tmp %>% summarise_scores(by=c("model", "month"))

cbind(tmp1[,1:2], round(tmp1[,c(5,10,11)],3))


#### SD Ratios ####
model.ses <- scored.out %>% 
  summarise_scores(by = c("model", "month"), na.rm=TRUE) %>%
  group_by(model) %>% select(model, se_mean, month) %>% arrange(model, month)

sd.model <- matrix(sqrt(model.ses$se_mean), ncol=4)
colnames(sd.model) <- unique(model.ses$model)

sd.obs <- as.vector(by(actuals$ged_sb, actuals$month_id, sd))

sd.model <- ts(sd.model, start=505)
sd.obs <- ts(sd.obs, start=505) 

plot(sd.model[,2:4]/sd.obs, plot.type = "single", col=1:3, lwd=2,
     xlab = "Month", ylab = expression(Model[SD]/Data[SD]),
     main="Standard deviation ratios plot")

abline(h=1, lty=2)
legend(506,3, legend = colnames(sd.model[,2:4]), col=1:3, lty=1, lwd=2)

#### Taylor diagram ####

library(openair)

# Need to do some date coversions to work smartly with the data here
# since the package objects expect real date variables
# tmp <- forcs22
# tmp <- tmp[tmp$model!="NB" & tmp$model!="ZINB",]
# tmp$dt <- as.factor(tmp$year)

forcs22$date <- ISOdate(forcs22$year, forcs22$month, 1)

# Normalized Taylor diagram
pdf("Taylors-2022.pdf")
TaylorDiagram(forcs22[forcs22$model!="Tweedie GAM",], 
              obs="observed", mod="predicted", 
              group=c("model"), type = "month",
              key.title = "Model",
              normalise=TRUE,
              main = "Normalized Taylor Diagram")

# UN-Normalized Taylor diagram

TaylorDiagram(forcs22[forcs22$model!="Tweedie GAM",], 
              obs="observed", mod="predicted", 
              group=c("model"), type = "month",
              key.title = "Model",
              normalise=FALSE, 
              main = "Unnormalized Taylor Diagram")
dev.off()


# Do the above without month 515
pdf("Taylors-2022-515.pdf")
TaylorDiagram(forcs22[forcs22$model!="Tweedie GAM" & forcs22$month_id!=515,], 
              obs="observed", mod="predicted", 
              group=c("model"), type = "month",
              key.title = "Model",
              normalise=TRUE,
              main = "Normalized Taylor Diagram")

# UN-Normalized Taylor diagram

TaylorDiagram(forcs22[forcs22$model!="Tweedie GAM"& forcs22$month_id!=515,], 
              obs="observed", mod="predicted", 
              group=c("model"), type = "month",
              key.title = "Model",
              normalise=FALSE, 
              main = "Unnormalized Taylor Diagram")
dev.off()

#### ggHoriz plot ####
library(ggHoriPlot)
library(ggthemes)

# Can see more options for how to use this method for plotting here
# https://rivasiker.github.io/ggHoriPlot/articles/examples.html

forcs22 %>% filter(model=="Poisson GLMM" & country_id < 120) %>% ggplot() + 
  geom_horizon(aes(x=date, y=predicted, horizonscale=4), rm.outliers=T) + 
  scale_fill_hcl(palette = 'BluGrn') +
  facet_grid(location~.) + 
  theme_bw() +
  theme(
    panel.spacing.y=unit(0, "lines"),
    strip.text.y = element_text(size = 4, angle = 0),
    legend.position = 'none',
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.border = element_blank(),panel.grid = element_blank()
  ) 

forcs22 %>% filter(model=="Poisson GLMM" & country_id < 40) %>% ggplot() + 
  geom_horizon(aes(x=date, y=predicted, horizonscale=4, 
                   origin=0, fill = ..Cutpoints..), rm.outliers=F) + 
  scale_fill_hcl(palette = 'RdBu') +
  facet_grid(location~.) + 
  theme_bw() +
  theme(
    panel.spacing.y=unit(0, "lines"),
    strip.text.y = element_text(size = 4, angle = 0),
    legend.position = 'none',
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.border = element_blank(),panel.grid = element_blank()
  ) 

#### Fan plot of the full set of forecasts for each model ####

# For details see https://guyabel.github.io/fanplot/index.html
# There is also another example / version here:
# https://surveillance.r-forge.r-project.org/pkgdown/reference/fanplot.html

library(fanplot)

# This is an aggregate forecast of the distribution over the whole globe for 2022
# Notice here we want the data with time in columns

pdf(file = "2023-2024-cm-models.pdf", width=6, height=4)
par(mfrow=c(1,4))
tmp <- matrix(as.numeric(unlist(lapply(glmm.P, as.vector))), ncol=27)
plot(NULL, type = "n", xlim = c(1, 27), ylim = c(0,300), ylab = "GLMM P")
fan(tmp, style="spaghetti", type="interval", med.ln=TRUE)

tmp <- matrix(as.numeric(unlist(lapply(glmm.NB, as.vector))), ncol=27)
plot(NULL, type = "n", xlim = c(1, 27), ylim = c(0,300), ylab = "GLMM NB")
fan(tmp, style="spaghetti", type="interval")

tmp <- matrix(as.numeric(unlist(lapply(glmm.TW, as.vector))), ncol=27)
plot(NULL, type = "n", xlim = c(1, 27), ylim = c(0,300), ylab = "GLMM TW")
fan(tmp, style="spaghetti", type="interval")

tmp <- matrix(as.numeric(unlist(lapply(local.TW, as.vector))), ncol=27)
plot(NULL, type = "n", xlim = c(1, 27), ylim = c(0,300), ylab = "GAM TW")
fan(tmp, style="spaghetti", type="interval")

# Boxplot version
par(mfrow=c(1,4))
tmp <- matrix(as.numeric(unlist(lapply(glmm.P, as.vector))), ncol=27)
plot(NULL, type = "n", xlim = c(1, 27), ylim = c(0,500), ylab = "GLMM P")
fan(tmp, style="boxfan", type="interval", med.ln=TRUE)

tmp <- matrix(as.numeric(unlist(lapply(glmm.NB, as.vector))), ncol=27)
plot(NULL, type = "n", xlim = c(1, 27), ylim = c(0,500), ylab = "GLMM NB")
fan(tmp, style="boxfan", type="interval")

tmp <- matrix(as.numeric(unlist(lapply(glmm.TW, as.vector))), ncol=27)
plot(NULL, type = "n", xlim = c(1, 27), ylim = c(0,500), ylab = "GLMM TW")
fan(tmp, style="boxfan", type="interval")

tmp <- matrix(as.numeric(unlist(lapply(local.TW, as.vector))), ncol=27)
plot(NULL, type = "n", xlim = c(1, 27), ylim = c(0,500), ylab = "GAM TW")
fan(tmp, style="boxfan", type="interval")

dev.off()

# save.image("combined-rolled-2022-local.RData")
