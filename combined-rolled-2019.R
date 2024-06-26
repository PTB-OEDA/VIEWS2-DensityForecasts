# Builds off of the Brandt-ViEWS2.Rmd document to do forecasts for 2019-2021

# Read in the various training sets
library(arrow)

df2017 <- read_parquet("shared_competition_data/cm_features_to_oct2017.parquet")
df2018 <- read_parquet("shared_competition_data/cm_features_to_oct2018.parquet") 
df2019 <- read_parquet("shared_competition_data/cm_features_to_oct2019.parquet")
df2020 <- read_parquet("shared_competition_data/cm_features_to_oct2020.parquet")

library(dplyr)
cmbn_df2018 <- bind_rows(df2017, df2018)
cmbn_df2019 <- bind_rows(cmbn_df2018, df2019)
cmbn_df2020 <- bind_rows(cmbn_df2019, df2020)


##################################################################
# CHANGE BELOW HERE TO ADJUST YEAR FOR SAMPLES
##################################################################

# Packages require a df not a tbl!!!!
dt1 <- as.data.frame(cmbn_df2018[,c("month_id", "country_id", "ged_sb", "ged_sb_tlag_1")])
dt1$country_id <- as.factor(dt1$country_id)

# For estimation
library(mgcv)
library(glmmTMB)

# Fit models over periods

psn.local <- bam(ged_sb ~ ged_sb_tlag_1 + s(month_id, country_id, bs="fs"),
                 data=dt1,
                 family = "poisson", discrete=TRUE)

psn.local.te <- bam(ged_sb ~ ged_sb_tlag_1 + te(month_id, country_id, bs="fs"),
                    data=dt1,
                    family = "poisson", discrete=TRUE)

nb.local <- bam(ged_sb ~ ged_sb_tlag_1 + s(month_id, country_id, bs="fs"),
                data=dt1,
                family = nb(), 
                discrete = TRUE,
                control = gam.control(trace=TRUE))

nb.local.te <- bam(ged_sb ~ ged_sb_tlag_1 + te(month_id, country_id, bs="fs"),
                   data=dt1,
                   family = nb(),                 
                   discrete = TRUE,
                   control = gam.control(trace=TRUE))

tw.local <- bam(ged_sb ~ ged_sb_tlag_1 + s(month_id, country_id, bs="fs"),
                data=dt1,
                family = tw(),                 
                discrete = TRUE,
                control = gam.control(trace=TRUE))

tw.local.te <- bam(ged_sb ~ ged_sb_tlag_1 + te(month_id, country_id, bs="fs"),
                   data=dt1,
                   family = tw(),                 
                   discrete = TRUE,
                   control = gam.control(trace=TRUE))

# GLMMs
# Need time as a factor for the next commands...
dt1$month_id <- factor(dt1$month_id)

psnar1.glmmtmb <- glmmTMB(ged_sb ~ ar1(month_id + 0|country_id), family = poisson,
                          data=dt1, control=glmmTMBControl(parallel = 4))

nbar1.glmmtmb <- glmmTMB(ged_sb ~ ar1(month_id + 0|country_id), family = nbinom1(),
                         data=dt1)

twar1.glmmtmb <- glmmTMB(ged_sb ~ ar1(month_id + 0|country_id), family = tweedie(),
                         data=dt1)

cat("----------------------------------------------------------------------\n")
timestamp()
cat("----------------------------------------------------------------------\n")

# Predictions in-sample
local.preds <- cbind(dt1[,1:2],
                     predict.bam(psn.local, type = "response"),
                     predict.bam(nb.local, type = "response"),
                     predict.bam(tw.local, type = "response"))
colnames(local.preds) <- c("month_id", "country_id", "P", "NB", "TW")

# Tensor preds
tensor.preds <- cbind(dt1[,1:2],
                      predict(psn.local.te, type = "response"),
                      predict(nb.local.te, type = "response"),
                      predict(tw.local.te, type = "response"))
colnames(tensor.preds) <- c("month_id", "country_id", "P", "NB", "TW")

# GLMM preds
glmm.preds <- cbind(dt1[,1:2],
                    predict(psnar1.glmmtmb, type = "response"),
                    predict(nbar1.glmmtmb, type = "response"),
                    #                    predict(cmpar1.glmmtmb, type = "response"),
                    predict(twar1.glmmtmb, type = "response"))
colnames(glmm.preds) <- c( "month_id", "country_id", "P", "NB", "TW")


# Get last obs for comparison of an "in-sample"
lastobs <- dt1[dt1$month_id=="466",]
local.last <- local.preds[local.preds$month_id=="466",]
tensor.last <- tensor.preds[tensor.preds$month_id==466,]
glmm.last <- glmm.preds[glmm.preds$month_id==466,]

# Make draws from the relevant modal predictions

# Local preds in-sample pdf
N <- 1000;    # Number of draws
n <- dim(lastobs)[1]
k <-  1;      #  number of forecasts

set.seed(1234)
local.P.fc <- sapply(1:n, function(i) {rpois(N, local.last$P[i])})
theta <- exp(nb.local$family$getTheta())
local.NB.fc <- sapply(1:n, function(i) {rnbinom(N, size=theta, mu=local.last$NB[i])})
local.TW.fc <- t(replicate(N, rTweedie(local.last$TW, p=1.581)))

# Tensor preds in-sample pdf
tensor.P.fc <- sapply(1:n, function(i) {rpois(N, tensor.last$P[i])})
theta <- exp(nb.local.te$family$getTheta())
tensor.NB.fc <- sapply(1:n, function(i) {rnbinom(N, size=theta, mu=tensor.last$NB[i])})
tensor.TW.fc <- t(replicate(N, rTweedie(tensor.last$TW,
                                        p=tw.local.te$family$getTheta(TRUE))))

# glmm preds in-sample pdf
glmm.P.fc <- sapply(1:n, function(i) {rpois(N, glmm.last$P[i])})
glmm.NB.fc <- sapply(1:n, function(i) {rnbinom(N, size=34, mu=glmm.last$NB[i])})
glmm.TW.fc <- t(replicate(N, rTweedie(glmm.last$TW, p=1.36)))


local.P.2019 <- local.NB.2019 <- local.TW.2019 <- vector(mode = "list", length=12)
 
k1 <- 466
timestamp()
for(i in 1:12)
{
  # Local models
  old.p <- apply(local.P.fc, 2, mean)
  tmp.p <- predict(psn.local, data.frame(country_id=local.last$country_id,
                                         month_id = rep((k1+i), 191),
                                         ged_sb_tlag_1 = old.p),
                   type = "response")
  local.P.fc <- sapply(1:n, function(i) {rpois(N, tmp.p[i])})

  old.nb <- apply(local.NB.fc, 2, mean)
  tmp.nb <- predict(nb.local, data.frame(country_id=local.last$country_id,
                                         month_id = rep((k1+i), 191),
                                         ged_sb_tlag_1 = old.nb),
                    type = "response")
  theta <- exp(nb.local$family$getTheta())
  local.NB.fc <- sapply(1:n, function(i) {rnbinom(N, size=theta, mu=tmp.nb[i])})

  old.tw <- apply(local.TW.fc, 2, mean)
  tmp.tw <- predict(tw.local, data.frame(country_id=local.last$country_id,
                                         month_id = rep((k1+i), 191),
                                         ged_sb_tlag_1 = old.tw),
                    type = "response")
  local.TW.fc <- t(replicate(N, rTweedie(tmp.tw, p=1.581)))

  local.P.2019[[i]] <- local.P.fc
  local.NB.2019[[i]] <- local.NB.fc
  local.TW.2019[[i]] <- local.TW.fc
}
timestamp()

# Now the same for the other models...

tensor.P.2019 <- tensor.NB.2019 <- tensor.TW.2019 <- vector(mode = "list", length=12)

k1 <- 466
timestamp()
for(i in 1:12)
{
  # tensor models
  old.p <- apply(tensor.P.fc, 2, mean)
  tmp.p <- predict(psn.local.te, data.frame(country_id=tensor.last$country_id,
                                            month_id = rep((k1+i), 191),
                                            ged_sb_tlag_1 = old.p),
                   type = "response")
  tensor.P.fc <- sapply(1:n, function(i) {rpois(N, tmp.p[i])})

  old.nb <- apply(tensor.NB.fc, 2, mean)
  tmp.nb <- predict(nb.local.te, data.frame(country_id=tensor.last$country_id,
                                            month_id = rep((k1+i), 191),
                                            ged_sb_tlag_1 = old.nb),
                    type = "response")
  theta <- exp(nb.local.te$family$getTheta())
  tensor.NB.fc <- sapply(1:n, function(i) {rnbinom(N, size=theta, mu=tmp.nb[i])})

  old.tw <- apply(tensor.TW.fc, 2, mean)
  tmp.tw <- predict(tw.local.te, data.frame(country_id=tensor.last$country_id,
                                            month_id = rep((k1+i), 191),
                                            ged_sb_tlag_1 = old.tw),
                    type = "response")
  tensor.TW.fc <- t(replicate(N, rTweedie(tmp.tw, p=1.581)))

  tensor.P.2019[[i]] <- tensor.P.fc
  tensor.NB.2019[[i]] <- tensor.NB.fc
  tensor.TW.2019[[i]] <- tensor.TW.fc
}


# glmm -- these are hard to do, since they do not easily admit more time
glmm.P.2019 <- glmm.NB.2019 <- glmm.TW.2019 <- vector(mode = "list", length=12)

k1 <- 466
timestamp()
for(i in 1:12)
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
  #  theta <- exp(nbar1.glmmtmb$family$getTheta())
  glmm.NB.fc <- sapply(1:n, function(i) {rnbinom(N, size=34, mu=tmp.nb[i])})

  old.tw <- apply(glmm.TW.fc, 2, mean)
  tmp.tw <- predict(twar1.glmmtmb, data.frame(country_id=glmm.last$country_id,
                                              month_id = as.factor(rep((k1+i), 191)),
                                              ged_sb_tlag_1 = old.tw),
                    type = "response", allow.new.levels=TRUE)
  glmm.TW.fc <- t(replicate(N, rTweedie(tmp.tw, p=1.581)))

  glmm.P.2019[[i]] <- glmm.P.fc
  glmm.NB.2019[[i]] <- glmm.NB.fc
  glmm.TW.2019[[i]] <- glmm.TW.fc
}
timestamp()

save.image("ViEWS2-2019.RData")
