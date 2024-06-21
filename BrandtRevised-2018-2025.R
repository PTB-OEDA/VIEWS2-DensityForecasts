# Brandt Forecasts for 2018-2025 inclusive 
# 
# Patrick T. Brandt
#
# 20240603 : Initial version from combined-rolled-2022-local.R
#            Generates forecasts for all months in 2018++
#            Uses updated datasets provided May 30, 2024
# 20240606 : Revisions to complete GAM and GLMM forecasts over new data
#            and periods.
#            Scoring...

library(arrow)
library(dplyr)

#### Read data ####

# Read in the various training datasets: These are the updated ones supplied in
# May 2024

# Featuresd data
cmbn_df <- read_parquet("~/VIEWS2/shared_competition_data/cm_features.parquet")

# Check the actuals, since we have slightly changed up the start / end dates
# from the earlier analysis, per [here](https://www.dropbox.com/scl/fo/a3t4yp1rfi2ymo2725mqy/ADCI_EufZ0qVWPrGNU2fcOY?e=1&preview=readme.txt&rlkey=who0s3nukkytdcbs4hypjkmch&st=mnffer15&dl=0)

a18 <- read_parquet("~/VIEWS2/shared_competition_data/actuals/cm/window=Y2018/cm_actuals_2018.parquet")
a19 <- read_parquet("~/VIEWS2/shared_competition_data/actuals/cm/window=Y2019/cm_actuals_2019.parquet")
a20 <- read_parquet("~/VIEWS2/shared_competition_data/actuals/cm/window=Y2020/cm_actuals_2020.parquet")
a21 <- read_parquet("~/VIEWS2/shared_competition_data/actuals/cm/window=Y2021/cm_actuals_2021.parquet")
a22 <- read_parquet("~/VIEWS2/shared_competition_data/actuals/cm/window=Y2022/cm_actuals_2022.parquet")
a23 <- read_parquet("~/VIEWS2/shared_competition_data/actuals/cm/window=Y2023/cm_actuals_2023.parquet")
a24 <- read_parquet("~/VIEWS2/shared_competition_data/actuals/cm/window=Y2024/cm_actuals_2024.parquet")

# Ranges of month_ids by training and forecast periods
timestamps <- as.data.frame(cbind(2018:2024, 
              rbind(range(a18$month_id),
                    range(a19$month_id),
                    range(a20$month_id),
                    range(a21$month_id),
                    range(a22$month_id),
                    range(a23$month_id),
                    range(a24$month_id))))

colnames(timestamps) <- c("ForecastYear", "start", "end")

cat("These are the time and month stamps for the training data\n")
print(timestamps)

# Here's an indexing of all the dates used in the analysis for later to 
# make things easier via window
datemap <- ts(1:532, start=c(1980,1), freq=12)

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

# Packages require a df not a tbl!!!!  Should be prior data all training sets here.
# Here we just select the variables.  Subsets by rows will be below.
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



#### GAMs with local effects ####

## traingam.fit: fits a single GAM specification for a dataset
#                intention is to use it in parallel

traingam.fit <- function(i, family="tw()") {
  
  # Fit the main gam model on a given node for a time sample through index
  # equal to "i"
  
  out <- bam(ged_sb ~ ged_sb_tlag_1 + s(month_id, country_id, bs="fs"),
             data=dt1, subset = (month_id<i), 
             family = family,   discrete = TRUE)
  return(list(i,out))
  }

##### Parallel compute version ####
# Setup
library(parallel)
nodes <- nrow(timestamps)+1  # Number of cores you need 
cl <- makeCluster(spec=nodes, type="SOCK")

# Pass the libraries, functions, and data to the cluster nodes
clusterEvalQ(cl, {library(mgcv); library(glmmTMB)})
clusterExport(cl, "traingam.fit")
clusterExport(cl, "dt1")

# Set the start time list over which the training samples will be made
st <- c(timestamps$start-2, 533)

cat("These are the start months for each sample:\n")
print(st)

# Fit all the models over the appropriate sub-samples of training data
# This is done in calls to the nodes for each sample and model density.  

system.time(tw.gam <- clusterApply(cl, st, traingam.fit))
system.time(p.gam <- clusterApply(cl, st, traingam.fit, family="poisson()"))
system.time(nb.gam <- clusterApply(cl, st, traingam.fit, family="nb()"))

#### GLMMs ####

## trainglmm.fit: fits a single GLMM specification for a dataset
#                intention is to use it in parallel

trainglmm.fit <- function(i, distr=tweedie()) {

  # This differs from above since glmm does not take a "subset =" arg
  dt2 <- subset(dt1, month_id<i)
  # # Need time as a factor for the next commands...
  dt2$month_id <- factor(dt2$month_id)
    
  # Fit the main gam model on a given node for a time sample through index
  # equal to "i"
  out <- glmmTMB(ged_sb ~ ar1(month_id + 0|country_id), 
                 family = distr, data=dt2)
  return(list(i, out))
}

##### Parallel compute version ####
system.time(tw.glmm <- clusterApply(cl, st, trainglmm.fit, distr=tweedie()))
system.time(p.glmm <- clusterApply(cl, st, trainglmm.fit, distr=poisson))
system.time(nb.glmm <- clusterApply(cl, st, trainglmm.fit, distr=nbinom1()))

# Shut down the local compute cluster
stopCluster(cl)

# #### AIC Fit Statistics ####

p.glmm.AIC <- sapply(1:7, function(i) {AIC(p.glmm[[i]][2][[1]])})
tw.glmm.AIC <- sapply(1:7, function(i) {AIC(tw.glmm[[i]][2][[1]])})
nb.glmm.AIC <- sapply(1:7, function(i) {AIC(nb.glmm[[i]][2][[1]])})

p.gam.AIC <- sapply(1:7, function(i) {AIC(p.gam[[i]][2][[1]])})
tw.gam.AIC <- sapply(1:7, function(i) {AIC(tw.gam[[i]][2][[1]])})
nb.gam.AIC <- sapply(1:7, function(i) {AIC(nb.gam[[i]][2][[1]])})

# Table of the fits for each sub-sample of the data
tmp <- cbind(timestamps[,1]-1, 
        p.gam.AIC, nb.gam.AIC, tw.gam.AIC,
        p.glmm.AIC, nb.glmm.AIC, tw.glmm.AIC)

colnames(tmp) <- c("Training End",
                   "Poisson GAM", "NegBin GAM", "Tweedie GAM",
                   "Poisson GLMM", "NegBin GLMM", "Tweedie GLMM")

print(round(tmp,0))

# Check that TW GLMM is best
print(tmp[,2:6] > tmp[,7])

library(xtable)

print(xtable(tmp, digits = 0))

#### Predictions ####

# Estimation costs in this component are low.  So there is little use
# of parallel code or multiple cores

# Uses one of the outputs from each of the above gam list objection
# created by traingam.fit() to generate predictions for k periods 
# into the future.

predictgam <- function(tmp, N, nstep=12)
  {
  k <- tmp[[1]][1]    # Index of the data / months
  fit <- tmp[2][[1]]  # fit GAM model

  # Storage / output
  local.out <- vector(mode = "list", length=nstep)
  
  # Setup constants
  lp <- which(fit$model$month_id==(k-1)) # last period indices
  n <- length(lp)                        # cases / countries
  local.last <- fit$model[lp,]           # additional data
  
  # Predict the last period of the sample -- same for all GAM
  last.period <- predict(fit, type="response")[lp]
  
  # Model specific updates and simulations for nsteps each
  
  # Poisson
  if(fit$family$family=="poisson")
  {
    local.fc <- sapply(1:n, function(i) {rpois(N, last.period[i])})
    
    for(i in 1:nstep) {
    old.p <- apply(local.fc, 2, mean)
    tmp.p <- predict(fit, data.frame(country_id=local.last$country_id,
                                           month_id = rep((k+i), n),
                                           ged_sb_tlag_1 = old.p),
                     type = "response")
    local.fc <- sapply(1:n, function(i) {rpois(N, tmp.p[i])})
    local.out[[i]] <- local.fc
    }
  }

  # Tweedie
  if(startsWith(fit$family$family, "Tweedie")==TRUE)
  {
    local.fc <- t(replicate(N, rTweedie(trunctop(last.period),
                                        p=fit$family$getTheta(TRUE))))

    for(i in 1:nstep) {
      old.p <- apply(local.fc, 2, mean)
      rm(local.fc)
      tmp.p <- predict(fit, data.frame(country_id=local.last$country_id,
                                       month_id = rep((k+i), n),
                                       ged_sb_tlag_1 = old.p),
                       type = "response")
      local.fc <- t(replicate(N, rTweedie(trunctop(tmp.p), 
                                          p=fit$family$getTheta(TRUE))))
      local.out[[i]] <- local.fc
      gc(); gc()
    }
  }

  # Negbin
  if(startsWith(fit$family$family, "Negative Binomial")==TRUE)
  {
    theta <- exp(fit$family$getTheta())
    local.fc <- sapply(1:n, function(i) {rnbinom(N, size=theta, 
                                                 mu=last.period[i])})

    for(i in 1:nstep) {
      old.p <- apply(local.fc, 2, mean)
      tmp.p <- predict(fit, data.frame(country_id=local.last$country_id,
                                       month_id = rep((k+i), n),
                                       ged_sb_tlag_1 = old.p),
                       type = "response")
      local.fc <- sapply(1:n, function(i) {rnbinom(N, size=theta, mu=tmp.p[i])})
      local.out[[i]] <- local.fc
    }
  }
  
  # Return the correctly predicted and simulated object
  return(local.out)
}

set.seed(213234)
# Now using basic R apply functions (no multiple CPUs needed here)
system.time(p.gam.forecasts <- lapply(p.gam[1:7],FUN=predictgam, N=1000, nstep=12))
system.time(nb.gam.forecasts <- lapply(nb.gam[1:7],FUN=predictgam, N=1000, nstep=12))
system.time(tw.gam.forecasts <- lapply(tw.gam[1:7],FUN=predictgam, N=1000, nstep=12))

# True future predictions for 2024:5-2025:6
set.seed(1234)

p.gam.future <- predictgam(p.gam[[8]], N=1000, 14)
nb.gam.future <- predictgam(nb.gam[[8]], N=1000, 14)
tw.gam.future <- predictgam(tw.gam[[8]], N=1000, 14)

# Uses one of the outputs from each of the above glmm list objects
# created by trainglmm.fit() to generate predictions for k periods 
# into the future.

predictglmm <- function(tmp, N, nstep=12)
{
  k <- tmp[[1]][1]    # Index of the data / months
  fit <- tmp[2][[1]]  # fit GLMM model
  
  # Storage / output
  local.out <- vector(mode = "list", length=nstep)
  
  # Setup constants
  lp <- which(fit$frame$month_id==(k-1)) # last period indices
  n <- length(lp)                        # cases / countries
  local.last <- fit$frame[lp,]           # additional data
  
  # Predict the last period of the sample -- same for all GAM
  last.period <- predict(fit, type="response")[lp]
  
  # Model specific updates and simulations for nsteps each
  
  # Poisson
  if(fit$modelInfo$family$family=="poisson")
  {
    local.fc <- sapply(1:n, function(i) {rpois(N, last.period[i])})
    
    for(i in 1:nstep) {
      old.p <- apply(local.fc, 2, mean)
      tmp.p <- predict(fit, 
                       newdata=data.frame(country_id=local.last$country_id,
                                       month_id = as.factor(rep((k+i), n)),
                                       ged_sb_tlag_1 = old.p),
                       type = "response", allow.new.levels=TRUE)
      local.fc <- sapply(1:n, function(i) {rpois(N, tmp.p[i])})
      local.out[[i]] <- local.fc
    }
  }
  
  # Tweedie
  if(fit$modelInfo$family$family=="tweedie")
  {
    local.fc <- t(replicate(N, rTweedie(last.period, p=family_params(fit))))
    
    for(i in 1:nstep) {
      old.p <- apply(local.fc, 2, mean)
      tmp.p <- predict(fit, 
                       newdata=data.frame(country_id=local.last$country_id,
                                       month_id = as.factor(rep((k+i), n)),
                                       ged_sb_tlag_1 = old.p),
                       type = "response", allow.new.levels=TRUE)
      local.fc <- t(replicate(N, rTweedie(tmp.p, p=family_params(fit))))
      local.out[[i]] <- local.fc
      gc(); gc()
    }
  }
  
  # Negbin
  if(fit$modelInfo$family$family=="nbinom1")
  {
    theta <- sigma(fit)
    local.fc <- sapply(1:n, function(i) {rnbinom(N, size=theta, 
                                                 mu=last.period[i])})

    for(i in 1:nstep) {
      old.p <- apply(local.fc, 2, mean)
      tmp.p <- predict(fit, newdata=data.frame(country_id=local.last$country_id,
                                       month_id = as.factor(rep((k+i), n)),
                                       ged_sb_tlag_1 = old.p),
                       type = "response", allow.new.levels=TRUE)
      local.fc <- sapply(1:n, function(i) {rnbinom(N, size=theta, mu=tmp.p[i])})
      local.out[[i]] <- local.fc
    }
  }
  
  # Return the correctly predicted and simulated object
  return(local.out)
}


# Now forecast the GLMM over the samples
# Again, no parallelization here, but it could be done using parLapply()
set.seed(98776)
system.time(p.glmm.forecasts <- lapply(p.glmm[1:7], FUN=predictglmm, N=1000, nstep=12))
system.time(nb.glmm.forecasts <- lapply(nb.glmm[1:7], FUN=predictglmm, N=1000, nstep=12))
system.time(tw.glmm.forecasts <- lapply(tw.glmm[1:7], FUN=predictglmm, N=1000, nstep=12))

set.seed(129856)
p.glmm.future <- predictglmm(p.glmm[[8]], N=1000, 14)
nb.glmm.future <- predictglmm(nb.glmm[[8]], N=1000, 14)
tw.glmm.future <- predictglmm(tw.glmm[[8]], N=1000, 14)

# Save the forecasts and the models
save.image("BrandtRevised.RData")

# 
load("BrandtRevised.RData")

# Drop models since they are not needed in what follows and they are
# holding up a lot of RAM
rm(nb.gam, p.gam, tw.gam, p.glmm, nb.glmm, tw.glmm)
gc()

#### Forecasts formatting setup ####

# Now make some comparisons against the actuals

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
                    "draw" = 1:N,
                    "month_id" = rep(mint:maxt, each=N*nc),
                    "model" = as.factor(model))

  # Merge on the actuals as "true_value" here
  all <- merge(all, data.frame("month_id" = ds$month_id,
                               "country_id" = ds$country_id,
                               "true_value" = ds$outcome))
  return(all)
}

actuals <- list(a18,a19,a20,a21,a22,a23,a24)
rm(a18,a19,a20,a21,a22,a23,a24)

# Constructs stacked forecasts with the actuals for comparisons
stackedforecasts<-function(f, a=actuals, mname="M")
{
  for(i in 1:length(actuals))
  {
    tmp <- convertforc(f[[i]], actuals[[i]],
                       N=1000, 12, model=mname)
    tmp$year <- 2017+i
    tmp$month <- tmp$month_id-min(tmp$month_id)+1
    if(i==1) {tout <- tmp } else { tout <- rbind(tout, tmp)}
  }
  rm(tmp)
  return(tout)
}


forcs <- rbind(stackedforecasts(p.gam.forecasts, actuals, "Poisson GAM"),
               stackedforecasts(nb.gam.forecasts, actuals, "Neg Binom GAM"),
               stackedforecasts(tw.gam.forecasts, actuals, "Tweedie GAM"),
               stackedforecasts(p.glmm.forecasts, actuals, "Poisson GLMM"),
               stackedforecasts(nb.glmm.forecasts, actuals, "Neg Binom GLMM"),
               stackedforecasts(tw.glmm.forecasts, actuals, "Tweedie GLMM"))


# Topcode forecasts with truncation and get correct naming for this step
forcs$prediction <- trunctop(forcs$prediction)
forcs$outcome <- forcs$prediction

rm(nb.gam.forecasts, nb.glmm.forecasts,
    p.gam.forecasts, p.glmm.forecasts,
    tw.gam.forecasts, tw.glmm.forecasts)

#### Write out forecasts for VIEWS submission ####

# This is vastly simplified from the 2023 version, since the above
# functions have done all of the list / array coversions needed
# for the data setups need for the various down-stream applications of
# forecast submission and post-processing

# Borrows from the earlier forc-org.R script in the github....
# write_parquet()....
setwd("~/VIEWS2/brandt_VIEWS2024/")

# Make a vector of names for each of the models
fnames <- names(table(forcs$model))
lower.fnames <- gsub(" ", "_", tolower(fnames))

##### 2018-2023 Forecasts ####

# Loop over the file names and the outputs needed
for(i in 1:length(fnames))
{
  for(j in 2018:2023) {

    fpathname <- paste(lower.fnames[i], "/cm/window=Y", j, "/",
                  lower.fnames[i], "_", j, ".parquet", sep="")
    # print(fpathname)
    write_parquet(subset(forcs, year==j & model==fnames[i],
                         select=c(month_id, country_id, draw, outcome)),
                  sink = fpathname)
  }
}

##### 2024-2025 Forecasts ####

# Now write out the true future forecasts for 2024:5-2025:6....
#
# "Jul 2024 - Jun 2025 using training data up until Apr 2024 
# (window=Y2024) (true future)"
#
# Data stop       April 2024, 532
# Forecast start   July 2024, 535
# Forecast end.    June 2025, 546
N <- 1000
ds <- dt1[dt1$month_id == max(dt1$month_id),]
nc <- nrow(ds)
mint <- max(ds$month_id) + 1
k <- 14 # Number of periods predicted
print(maxt <- mint + k - 1)
all <- data.frame("country_id" = rep(rep(ds$country_id, each=N), k),
                  "draw" = 1:N,
                  "month_id" = rep(mint:maxt, each=N*nc))

# Write out the results for the true future data

# GAM
tmp <- cbind(all,outcome=unlist(lapply(p.gam.future, "as.vector")))
write_parquet(tmp[tmp$month_id>534,], 
              sink="poisson_gam/cm/window=Y2024/poisson_gam_2024.parquet")

tmp <- cbind(all,outcome=unlist(lapply(nb.gam.future, "as.vector")))
write_parquet(tmp[tmp$month_id>534,], 
              sink="neg_binom_gam/cm/window=Y2024/neg_binom_gam_2024.parquet")

tmp <- cbind(all,outcome=unlist(lapply(tw.gam.future, "as.vector")))
write_parquet(tmp[tmp$month_id>534,], 
              sink="tweedie_gam/cm/window=Y2024/tweedie_gam_2024.parquet")
# GLMM
tmp <- cbind(all,outcome=unlist(lapply(p.glmm.future, "as.vector")))
write_parquet(tmp[tmp$month_id>534,], 
              sink="poisson_glmm/cm/window=Y2024/poisson_glmm_2024.parquet")

tmp <- cbind(all,outcome=unlist(lapply(nb.glmm.future, "as.vector")))
write_parquet(tmp[tmp$month_id>534,], 
              sink="neg_binom_glmm/cm/window=Y2024/neg_binom_glmm_2024.parquet")

tmp <- cbind(all,outcome=unlist(lapply(tw.glmm.future, "as.vector")))
write_parquet(tmp[tmp$month_id>534,], 
              sink="tweedie_glmm/cm/window=Y2024/tweedie_glmm_2024.parquet")

rm(tmp,all)

rm(p.gam.future,nb.gam.future,tw.gam.future,
   p.glmm.future,nb.glmm.future,tw.glmm.future)



#### Forecast Performance ####
library(scoringutils)
library(magrittr)
library(dplyr)

# Keep only results through 2023 for now in scoring and drop redundant column
forcs <- forcs[,-9]
forcs <- forcs %>% filter(year< 2024)

# Adjust here for the latest version of the scoringutils package and add
# country names
country.labels <- read.csv("~/VIEWS2/matching_tables/countries.csv")[,c(1,2,5)]
country.labels$country_id <- country.labels$id
forcs <- merge(forcs, country.labels[,2:4])

colnames(forcs) <- c("country_id", "month_id", "predicted",
                       "sample_id", "model", "observed", "year", "month",
                       "location", "isoab")

# Need to do some date conversions to work smartly with the data here
# since the package objects expect real date variables for the later Taylor
# diagrams

forcs$date <- ISOdate(forcs$year, forcs$month, 1)

# Merge in the other forecasts for the baselines
#
# Produced by baseline-views-forecasts.R
#
# Include the VIEWS benchmarks

load("~/Documents/GitHub/VIEWS2-DensityForecasts/viewsforcs.RData")

# Keep only through 2023 for now...
viewsforcs <- viewsforcs %>% filter(year < 2024)

# Make model a factor so they can be merged.
viewsforcs$model <- as.factor(viewsforcs$model)

# Merge up the baseline with the ones produced here
identical(colnames(viewsforcs), colnames(forcs))
forcs <- rbind(forcs,viewsforcs)
rm(viewsforcs)
gc()

# Save forecasts setup through here....
save.image("stackedforecasts-2018-2025.RData")

# Now do the forecast scorings
scored.out <- score(as_forecast(forcs,
                                forecast_unit = c("model", "month_id", "isoab")))

crps <- scored.out %>% summarise_scores(by=c("model", "month_id"), na.rm=TRUE) %>%
  select(c("model", "month_id", "crps"))

# Make a time series of the CRPS for each model / forecast period
crps.ts <- score(as_forecast(forcs,
                                forecast_unit = c("model", "date", "isoab"))) %>%
  select(c("model", "date", "crps")) %>% arrange(model, date)

tmp <- crps.ts %>% summarise_scores(by=c("model", "date"))

mnames <- unique(tmp$model)
crps.ts1 <- matrix(tmp$crps, ncol=length(mnames))
colnames(crps.ts1) <- mnames

crps.ts1 <- ts(crps.ts1, start=c(2018,1), freq=12)

# These are the annualized CRPS to compare to their metrics

# So this is all the years and models.
xtable(round(aggregate.ts(crps.ts1, nfrequency = 1,FUN = mean), 2))



# But this belies the changes in the data over time that make this hard:

plot(window(crps.ts1[,4:9], end=c(2020,6)), plot.type="s", lty=1:6)
plot(window(crps.ts1[,4:9], start=c(2020,7)), plot.type="s", lty=1:6)

# Fancier plot
tmp <- window(cbind(roll.ged.sb[,2], crps.ts1[,4:8]),start=c(2017,1))
colnames(tmp) <- c(colnames(roll.ged.sb)[2], 
                   paste(colnames(crps.ts1)[4:8], "CRPS"))

pdf(file="~/VIEWS2/crpsbymonth.pdf", width=6, height=3.5)
par(mfrow=c(1,2))
plot(window(tmp, start=c(2018,1), end=c(2020,6)), 
     ylim=c(0,700), plot.type="s",
     cex.lab=0.8, cex.axis=0.8, main="", col=1:7, lwd=2, lty=1:7, ylab="")

plot(window(tmp, start=c(2020,7), end=c(2023,12)), 
     ylim=c(0,700), plot.type="s",
     cex.lab=0.8, cex.axis=0.8, main="", col=1:7, lwd=2, lty=1:7, ylab="")
dev.off()

# # crps.rel <- scored.out %>% summarise_scores(by=c("model", "month_id"),
# #                                                   relative_skill=TRUE,
# #                                                   baseline = "Tweedie GLMM")
# # 
# # # Look at the measures....
# # tmp <- add_relative_skill(scored.out, by=c("model", "month_id"),
# #                    metric="crps", baseline="Tweedie GLMM")
# # 
# # tmp1 <- tmp %>% summarise_scores(by=c("model", "month_id"))
# # 
# # cbind(tmp1[,1:2], round(tmp1[,c(5,10,11)],3))
# 
# 
#### SD Ratios ####
model.ses <- scored.out %>%
  summarise_scores(by = c("model", "month_id"), na.rm=TRUE) %>%
  group_by(model) %>% select(model, se_mean, month_id) %>% arrange(model, month_id)

# Organize and label columns with model names
mdls <- unique(model.ses$model)
sd.model <- matrix(sqrt(model.ses$se_mean), ncol=length(mdls))
colnames(sd.model) <- mdls

# Get the relevant time series setup
sd.model <- ts(sd.model, start=c(2018,1), freq=12)
sd.obs <- window(roll.ged.sb[,3], start=c(2018,1), end=c(2023,12))

# Plot the SDs
par(mfcol=c(3,2))
plot(window(sd.model[,4:9]/sd.obs, start=c(2018,1), end=c(2018,12)),
     plot.type = "single", col=1:6, lwd=2,
     xlab = "Month", ylab = expression(Model[SD]/Data[SD]),
     main="Standard deviation ratios, 2018")
abline(h=1, lty=2)
plot(window(sd.model[,4:9]/sd.obs, start=c(2019,1), end=c(2019,12)),
     plot.type = "single", col=1:6, lwd=2,
     xlab = "Month", ylab = expression(Model[SD]/Data[SD]),
     main="Standard deviation ratios, 2019")
abline(h=1, lty=2)

plot(window(sd.model[,4:9]/sd.obs, start=c(2020,1), end=c(2020,12)),
     plot.type = "single", col=1:6, lwd=2,
     xlab = "Month", ylab = expression(Model[SD]/Data[SD]),
     main="Standard deviation ratios, 2020")
abline(h=1, lty=2)
plot(window(sd.model[,4:9]/sd.obs, start=c(2021,1), end=c(2021,12)),
     plot.type = "single", col=1:6, lwd=2,
     xlab = "Month", ylab = expression(Model[SD]/Data[SD]),
     main="Standard deviation ratios, 2021")
abline(h=1, lty=2)

plot(window(sd.model[,4:9]/sd.obs, start=c(2022,1), end=c(2022,12)),
     plot.type = "single", col=1:6, lwd=2,
     xlab = "Month", ylab = expression(Model[SD]/Data[SD]),
     main="Standard deviation ratios, 2022")
abline(h=1, lty=2)
plot(window(sd.model[,4:9]/sd.obs, start=c(2023,1), end=c(2023,12)),
     plot.type = "single", col=1:6, lwd=2,
     xlab = "Month", ylab = expression(Model[SD]/Data[SD]),
     main="Standard deviation ratios, 2023")
abline(h=1, lty=2)




# Now can add / contrast the prediction errors with the data and its
# observed volatility...
pdf(file="~/VIEWS2/sdratios.pdf", width=6, height=6)
par(mfcol=c(2,2), mar=c(4,4,1,1))
plot(window(roll.ged.sb[,2], start=c(2018,1), end=c(2022,12)), ylab="Mean Counts")
plot(window(sd.model[,4:9]/sd.obs, start=c(2018,1), end=c(2022,12)),
     plot.type = "single", col=1:6, lwd=2,
     xlab = "Month", ylab = expression(Model[SD]/Data[SD]))
abline(h=1, lty=2)
plot(window(roll.ged.sb[,2], start=c(2023,1), end=c(2023,12)), ylab="Mean Counts")
plot(window(sd.model[,4:9]/sd.obs, start=c(2023,1), end=c(2023,12)),
     plot.type = "single", col=1:6, lwd=2,
     xlab = "Month", ylab = expression(Model[SD]/Data[SD]))
abline(h=1, lty=2)
dev.off()

#### Taylor diagram ####

library(openair)

# Normalized Taylor diagrams

TaylorDiagram(forcs[forcs$model!="Poisson GAM" &
                      forcs$model!="Neg Binom GAM" &
                      forcs$model!="Tweedie GAM",],
              obs="observed", mod="predicted",
              group=c("model"), type = "year",
              key.title = "Model",
              normalise=TRUE,
              main = "Normalized Taylor Diagram")

# UN-Normalized Taylor diagram for GLMMs

TaylorDiagram(forcs[forcs$model!="Poisson GAM" &
                      forcs$model!="Neg Binom GAM" &
                      forcs$model!="Tweedie GAM",],
              obs="observed", mod="predicted",
              group=c("model"), type = "year",
              key.title = "Model",
              normalise=FALSE,
              main = "Unnormalized Taylor Diagram")


# Drop 2023
TaylorDiagram(forcs[forcs$model!="Poisson GAM" &
                      forcs$model!="Neg Binom GAM" &
                      forcs$model!="Tweedie GAM" & 
                      forcs$year < 2023,],
              obs="observed", mod="predicted",
              group=c("model"), type = "year",
              key.title = "Model",
              normalise=TRUE,
              main = "Normalized Taylor Diagram")

# # Experimenting...
TaylorDiagram(forcs[forcs$model!="Poisson GAM" &
                      forcs$model!="Neg Binom GAM" &
                      forcs$model!="Tweedie GAM",],
              obs="observed", mod="predicted",
              group=c("model", "country_id"), type = "year",
              key.title = "Model",
              normalise=TRUE,
              main = "")

TaylorDiagram(forcs[forcs$model!="Poisson GAM" &
                      forcs$model!="Neg Binom GAM" &
                      forcs$model!="Tweedie GAM",],
              obs="observed", mod="predicted",
              group=c("model", "year"))

# TaylorDiagram(forcs[forcs$model!="Poisson GAM" & 
#                       forcs$model!="Neg Binom GAM" & 
#                       forcs$model!="Tweedie GAM" & 
#                       forcs$model!="Poisson GLMM" & 
#                       forcs$year<2023,],
#               obs="observed", mod="predicted",
#               group=c("year", "country_id"), type = "model",
#               key.title = "Model",
#               normalise=TRUE,
#               main = "")
# 

scored.out[-grep("GAM", scored.out$model),] %>% 
  ggplot() + geom_point(aes(x=se_mean, y=bias, group=model)) + 
  facet_wrap(~model) + coord_cartesian(xlim = c(0, 3000))

# #### ggHoriz plot ####
# library(ggHoriPlot)
# library(ggthemes)
# 
# # Can see more options for how to use this method for plotting here
# # https://rivasiker.github.io/ggHoriPlot/articles/examples.html
# 
# forcs22 %>% filter(model=="Poisson GLMM" & country_id < 120) %>% ggplot() + 
#   geom_horizon(aes(x=date, y=predicted, horizonscale=4), rm.outliers=T) + 
#   scale_fill_hcl(palette = 'BluGrn') +
#   facet_grid(location~.) + 
#   theme_bw() +
#   theme(
#     panel.spacing.y=unit(0, "lines"),
#     strip.text.y = element_text(size = 4, angle = 0),
#     legend.position = 'none',
#     axis.text.y = element_blank(),
#     axis.title.y = element_blank(),
#     axis.ticks.y = element_blank(),
#     panel.border = element_blank(),panel.grid = element_blank()
#   ) 
# 
# forcs22 %>% filter(model=="Poisson GLMM" & country_id < 40) %>% ggplot() + 
#   geom_horizon(aes(x=date, y=predicted, horizonscale=4, 
#                    origin=0, fill = ..Cutpoints..), rm.outliers=F) + 
#   scale_fill_hcl(palette = 'RdBu') +
#   facet_grid(location~.) + 
#   theme_bw() +
#   theme(
#     panel.spacing.y=unit(0, "lines"),
#     strip.text.y = element_text(size = 4, angle = 0),
#     legend.position = 'none',
#     axis.text.y = element_blank(),
#     axis.title.y = element_blank(),
#     axis.ticks.y = element_blank(),
#     panel.border = element_blank(),panel.grid = element_blank()
#   ) 
# 

#### Checks on actuals ####
actuals.given <- merge(rbind(data.frame(actuals[[3]], year=2019, 
                                        month=actuals[[3]]$month_id - min(actuals[[3]]$month_id) + 1),
                             data.frame(actuals[[4]], year=2020, 
                                        month=actuals[[4]]$month_id - min(actuals[[4]]$month_id) + 1),
                             data.frame(actuals[[5]], year=2021,
                                        month=actuals[[5]]$month_id - min(actuals[[5]]$month_id) + 1),
                             data.frame(actuals[[6]], year=2022,
                                        month=actuals[[6]]$month_id - min(actuals[[6]]$month_id) + 1),
                             data.frame(actuals[[7]], year=2023,
                                        month=actuals[[7]]$month_id - min(actuals[[7]]$month_id) + 1)),
                       country.labels[,2:4])

actuals.given$date <- ISOdate(year = actuals.given$year, month = actuals.given$month, day=1)

quantile(actuals.given$outcome, probs = seq(0.85,0.9999, by=0.005))

# So for post 2020 what are the cases in the highest 0.005 of all?
out <- actuals.given %>% filter(outcome>1940) %>% arrange(date)


print(xtable(cbind(out[,c(2,6)],format(out[,8], "%B %Y"))),
      file="worstsince2018.tex")

par(mfcol=c(2,2), mar=c(3,2,3,1))
plot(roll.ged.sb[,2], type="l", 
     xlab = "Month", ylab = "Mean GED SB events",
     main = "Mean by month")
abline(v=2021, lty=2, lwd=2)

plot(roll.ged.sb[,3], type="l", xlab = "Month", 
     ylab = "Standard Deviation of GED SB events",
     main = "Std Dev. by month")
abline(v=2021, lty=2, lwd=2)

plot(ts(unlist(lapply(sbdata, max)), start=c(1990,1), freq=12), 
     ylab="", xlab="", cex.lab=0.8, 
     main = "Maximums by month")
abline(v=2021, lty=2, lwd=2)
plot(ts(matrix(unlist(lapply(sbdata, quantile, probs=c(0.95,0.99))), ncol=2, byrow=TRUE), start=c(1990,1), freq=12),
     plot.type="s", col=1:2, lty=1:2, lwd=2, ylab="", xlab="",
     cex.axis=0.5,
     main = "95th and 99th percentiles" )
abline(v=2021, lty=2, lwd=2)


#### Fan plot of the full set of forecasts for each model ####

# For details see https://guyabel.github.io/fanplot/index.html
# There is also another example / version here:
# https://surveillance.r-forge.r-project.org/pkgdown/reference/fanplot.html

# library(fanplot)
# 
# # This is an aggregate forecast of the distribution over the whole globe for 2022
# # Notice here we want the data with time in columns
# 
# pdf(file = "2023-2024-cm-models.pdf", width=6, height=4)
# par(mfrow=c(1,4))
# tmp <- matrix(as.numeric(unlist(lapply(glmm.P, as.vector))), ncol=27)
# plot(NULL, type = "n", xlim = c(1, 27), ylim = c(0,300), ylab = "GLMM P")
# fan(tmp, style="spaghetti", type="interval", med.ln=TRUE)
# 
# tmp <- matrix(as.numeric(unlist(lapply(glmm.NB, as.vector))), ncol=27)
# plot(NULL, type = "n", xlim = c(1, 27), ylim = c(0,300), ylab = "GLMM NB")
# fan(tmp, style="spaghetti", type="interval")
# 
# tmp <- matrix(as.numeric(unlist(lapply(glmm.TW, as.vector))), ncol=27)
# plot(NULL, type = "n", xlim = c(1, 27), ylim = c(0,300), ylab = "GLMM TW")
# fan(tmp, style="spaghetti", type="interval")
# 
# tmp <- matrix(as.numeric(unlist(lapply(local.TW, as.vector))), ncol=27)
# plot(NULL, type = "n", xlim = c(1, 27), ylim = c(0,300), ylab = "GAM TW")
# fan(tmp, style="spaghetti", type="interval")
# 
# # Boxplot version
# par(mfrow=c(1,4))
# tmp <- matrix(as.numeric(unlist(lapply(glmm.P, as.vector))), ncol=27)
# plot(NULL, type = "n", xlim = c(1, 27), ylim = c(0,500), ylab = "GLMM P")
# fan(tmp, style="boxfan", type="interval", med.ln=TRUE)
# 
# tmp <- matrix(as.numeric(unlist(lapply(glmm.NB, as.vector))), ncol=27)
# plot(NULL, type = "n", xlim = c(1, 27), ylim = c(0,500), ylab = "GLMM NB")
# fan(tmp, style="boxfan", type="interval")
# 
# tmp <- matrix(as.numeric(unlist(lapply(glmm.TW, as.vector))), ncol=27)
# plot(NULL, type = "n", xlim = c(1, 27), ylim = c(0,500), ylab = "GLMM TW")
# fan(tmp, style="boxfan", type="interval")
# 
# tmp <- matrix(as.numeric(unlist(lapply(local.TW, as.vector))), ncol=27)
# plot(NULL, type = "n", xlim = c(1, 27), ylim = c(0,500), ylab = "GAM TW")
# fan(tmp, style="boxfan", type="interval")
# 
# dev.off()
# 


#### Deprecated / test code ####

# Now unit test on one sample / forecast setup at a time....

# # Poisson test
# tmp <- p.gam[[1]]
# test <- predictgam(tmp, 10 ,nstep = 6)
# print(str(test))
# 
# # Tweedie test
# tmp <- tw.gam[[1]]
# test <- predictgam(tmp, 100 ,nstep = 12)
# print(str(test))
# 
# # Negbin
# tmp <- nb.gam[[1]]
# test <- predictgam(tmp, 20 ,nstep = 7)
# print(str(test))

# Now unit test on one sample / forecast setup at a time....

# # Poisson test
# tmp <- p.glmm[[1]]
# system.time(test <- predictglmm(tmp, 10 ,nstep = 6))
# print(str(test))
# lapply(test, mean, na.rm=TRUE)
# 
# # Tweedie test
# tmp <- tw.glmm[[1]]
# system.time(test <- predictglmm(tmp, 100 ,nstep = 12))
# print(str(test))
# 
# # Negbin
# tmp <- nb.glmm[[1]]
# system.time(test <- predictglmm(tmp, 20 ,nstep = 7))
# print(str(test))
# lapply(test, sd, na.rm=TRUE)


# # Convert forecasts  as a test
# 
# for(i in 1:length(actuals))
# {
#   tmp <- convertforc(nb.glmm.forecasts[[i]], actuals[[i]], 
#                      N=1000, 12, model="Neg Binom GLMM")
#   tmp$year <- 2017+i
#   tmp$month <- tmp$month_id-min(tmp$month_id)+1
#   
#   if(i==1) {tout <- tmp } else { tout <- rbind(tout, tmp)}
# }
# rm(tmp)
