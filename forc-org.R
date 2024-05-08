# forc-org.R  : Forecast organization and output for later uses
#
# Patrick T. Brandt
#
# 20230926 : Reorganize and output forecasts for views
# 20231002 : Corrected country and months to be correct integers (decodes the factors)
# 20231004 : Add top codings for NA and extreme values trunctop()
# 20240508 : Adds 2022 forecasts output to the results.

library(arrow)
source("reorg.R")
setwd("~/VIEWS2/")

# 2018 models / forecasts / data
load("ViEWS2-prelim.RData")

# Get final forecast targets so we can id the final countries
dfa2021 <- read_parquet("shared_competition_data/cm_actuals_2021.parquet")

#################################
# Set country numbers correctly
#################################

# Do a quick review of the correct ones
cids <- scan("viewscountryid.txt", sep=",")

cout <- cbind(as.numeric(lastobs$series), 
              labels(lastobs$series), 
              unique(dfa2021$country_id),
              cids)

# Keep only the objects we absolutely need
cnum <- as.integer(cout[,4])

# 2018 forecast data
local.poisson.2018 <- reorg(local.P.2018, cnum=cnum , k1=456, k2=12)
local.negbin.2018 <- reorg(local.NB.2018, cnum=cnum, k1=456, k2=12)
local.tweedie.2018 <- reorg(local.TW.2018, cnum=cnum, k1=456, k2=12)

tensor.poisson.2018 <- reorg(tensor.P.2018, cnum=cnum, k1=456, k2=12)
tensor.negbin.2018 <- reorg(tensor.NB.2018, cnum=cnum, k1=456, k2=12)
tensor.tweedie.2018 <- reorg(tensor.TW.2018, cnum=cnum, k1=456, k2=12)

glmm.poisson.2018 <- reorg(glmm.P.2018, cnum=cnum, k1=456, k2=12)
glmm.negbin.2018 <- reorg(glmm.NB.2018, cnum=cnum, k1=456, k2=12)
glmm.tweedie.2018 <- reorg(glmm.TW.2018, cnum=cnum, k1=456, k2=12)

# Code to inspect overflow or NA cases
table(local.poisson.2018[(is.na(local.poisson.2018$outcome)),1:2])
table(local.negbin.2018[(is.na(local.negbin.2018$outcome)),1:2])


# Truncate the NA overflows and the large values at upper.tol value
# This automates so we can look at sensitivity later
trunctop <- function (x, upper.tol=1e9)
{
  x <- replace(x, is.na(x), upper.tol)  # NA and NaN cases
  x <- replace(x, x>upper.tol, upper.tol)  # Truncate the top
  return(x)
}

# Deal with the NA and upper limits
local.poisson.2018$outcome <- trunctop(local.poisson.2018$outcome)
local.negbin.2018$outcome <- trunctop(local.negbin.2018$outcome)
local.tweedie.2018$outcome <- trunctop(local.tweedie.2018$outcome)
tensor.poisson.2018$outcome <- trunctop(tensor.poisson.2018$outcome)
tensor.negbin.2018$outcome <- trunctop(tensor.negbin.2018$outcome)
tensor.tweedie.2018$outcome <- trunctop(tensor.tweedie.2018$outcome)
glmm.poisson.2018$outcome <- trunctop(glmm.poisson.2018$outcome)
glmm.negbin.2018$outcome <- trunctop(glmm.negbin.2018$outcome)
glmm.tweedie.2018$outcome <- trunctop(glmm.tweedie.2018$outcome)

# Check again
table(local.poisson.2018[(is.na(local.poisson.2018$outcome)),1:2])
table(local.negbin.2018[(is.na(local.negbin.2018$outcome)),1:2])


# Write the outputs
write_parquet(local.poisson.2018, 
              "Brandt_VIEWS2023/poisson_gamlocal/cm/test_window_2018/local_poisson_2018.parquet")

write_parquet(local.negbin.2018, 
              "Brandt_VIEWS2023/negbin_gamlocal/cm/test_window_2018/local_negbin_2018.parquet")

write_parquet(local.tweedie.2018, 
              "Brandt_VIEWS2023/tweedie_gamlocal/cm/test_window_2018/local_tweedie_2018.parquet")

write_parquet(tensor.poisson.2018, 
              "Brandt_VIEWS2023/poisson_gamtensor/cm/test_window_2018/tensor_poisson_2018.parquet")

write_parquet(tensor.negbin.2018, 
              "Brandt_VIEWS2023/negbin_gamtensor/cm/test_window_2018/tensor_negbin_2018.parquet")

write_parquet(tensor.tweedie.2018, 
              "Brandt_VIEWS2023/tweedie_gamtensor/cm/test_window_2018/tensor_tweedie_2018.parquet")

write_parquet(glmm.poisson.2018, 
              "Brandt_VIEWS2023/poisson_glmm/cm/test_window_2018/glmm_poisson_2018.parquet")

write_parquet(glmm.negbin.2018, 
              "Brandt_VIEWS2023/negbin_glmm/cm/test_window_2018/glmm_negbin_2018.parquet")

write_parquet(glmm.tweedie.2018, 
              "Brandt_VIEWS2023/tweedie_glmm/cm/test_window_2018/glmm_tweedie_2018.parquet")

##############
# 2019
##############

load("ViEWS2-2019.RData")

local.poisson.2019 <- reorg(local.P.2019, cnum=cnum , k1=468, k2=12)
local.negbin.2019 <- reorg(local.NB.2019, cnum=cnum, k1=468, k2=12)
local.tweedie.2019 <- reorg(local.TW.2019, cnum=cnum, k1=468, k2=12)

tensor.poisson.2019 <- reorg(tensor.P.2019, cnum=cnum, k1=468, k2=12)
tensor.negbin.2019 <- reorg(tensor.NB.2019, cnum=cnum, k1=468, k2=12)
tensor.tweedie.2019 <- reorg(tensor.TW.2019, cnum=cnum, k1=468, k2=12)

glmm.poisson.2019 <- reorg(glmm.P.2019, cnum=cnum, k1=468, k2=12)
glmm.negbin.2019 <- reorg(glmm.NB.2019, cnum=cnum, k1=468, k2=12)
glmm.tweedie.2019 <- reorg(glmm.TW.2019, cnum=cnum, k1=468, k2=12)

# Deal with the NA and upper limits
local.poisson.2019$outcome <- trunctop(local.poisson.2019$outcome)
local.negbin.2019$outcome <- trunctop(local.negbin.2019$outcome)
local.tweedie.2019$outcome <- trunctop(local.tweedie.2019$outcome)
tensor.poisson.2019$outcome <- trunctop(tensor.poisson.2019$outcome)
tensor.negbin.2019$outcome <- trunctop(tensor.negbin.2019$outcome)
tensor.tweedie.2019$outcome <- trunctop(tensor.tweedie.2019$outcome)
glmm.poisson.2019$outcome <- trunctop(glmm.poisson.2019$outcome)
glmm.negbin.2019$outcome <- trunctop(glmm.negbin.2019$outcome)
glmm.tweedie.2019$outcome <- trunctop(glmm.tweedie.2019$outcome)

# Write the outputs
write_parquet(local.poisson.2019, 
              "Brandt_VIEWS2023/poisson_gamlocal/cm/test_window_2019/local_poisson_2019.parquet")

write_parquet(local.negbin.2019, 
              "Brandt_VIEWS2023/negbin_gamlocal/cm/test_window_2019/local_negbin_2019.parquet")

write_parquet(local.tweedie.2019, 
              "Brandt_VIEWS2023/tweedie_gamlocal/cm/test_window_2019/local_tweedie_2019.parquet")

write_parquet(tensor.poisson.2019, 
              "Brandt_VIEWS2023/poisson_gamtensor/cm/test_window_2019/tensor_poisson_2019.parquet")

write_parquet(tensor.negbin.2019, 
              "Brandt_VIEWS2023/negbin_gamtensor/cm/test_window_2019/tensor_negbin_2019.parquet")

write_parquet(tensor.tweedie.2019, 
              "Brandt_VIEWS2023/tweedie_gamtensor/cm/test_window_2019/tensor_tweedie_2019.parquet")

write_parquet(glmm.poisson.2019, 
              "Brandt_VIEWS2023/poisson_glmm/cm/test_window_2019/glmm_poisson_2019.parquet")

write_parquet(glmm.negbin.2019, 
              "Brandt_VIEWS2023/negbin_glmm/cm/test_window_2019/glmm_negbin_2019.parquet")

write_parquet(glmm.tweedie.2019, 
              "Brandt_VIEWS2023/tweedie_glmm/cm/test_window_2019/glmm_tweedie_2019.parquet")

#####################
# 2020
#####################
load("ViEWS2-2020.RData")

local.poisson.2020 <- reorg(local.P.2020, cnum=cnum , k1=480, k2=12)
local.negbin.2020 <- reorg(local.NB.2020, cnum=cnum, k1=480, k2=12)
local.tweedie.2020 <- reorg(local.TW.2020, cnum=cnum, k1=480, k2=12)

tensor.poisson.2020 <- reorg(tensor.P.2020, cnum=cnum, k1=480, k2=12)
tensor.negbin.2020 <- reorg(tensor.NB.2020, cnum=cnum, k1=480, k2=12)
tensor.tweedie.2020 <- reorg(tensor.TW.2020, cnum=cnum, k1=480, k2=12)

glmm.poisson.2020 <- reorg(glmm.P.2020, cnum=cnum, k1=480, k2=12)
glmm.negbin.2020 <- reorg(glmm.NB.2020, cnum=cnum, k1=480, k2=12)
glmm.tweedie.2020 <- reorg(glmm.TW.2020, cnum=cnum, k1=480, k2=12)

# Deal with the NA and upper limits
local.poisson.2020$outcome <- trunctop(local.poisson.2020$outcome)
local.negbin.2020$outcome <- trunctop(local.negbin.2020$outcome)
local.tweedie.2020$outcome <- trunctop(local.tweedie.2020$outcome)
tensor.poisson.2020$outcome <- trunctop(tensor.poisson.2020$outcome)
tensor.negbin.2020$outcome <- trunctop(tensor.negbin.2020$outcome)
tensor.tweedie.2020$outcome <- trunctop(tensor.tweedie.2020$outcome)
glmm.poisson.2020$outcome <- trunctop(glmm.poisson.2020$outcome)
glmm.negbin.2020$outcome <- trunctop(glmm.negbin.2020$outcome)
glmm.tweedie.2020$outcome <- trunctop(glmm.tweedie.2020$outcome)


# Write the outputs
write_parquet(local.poisson.2020, 
              "Brandt_VIEWS2023/poisson_gamlocal/cm/test_window_2020/local_poisson_2020.parquet")

write_parquet(local.negbin.2020, 
              "Brandt_VIEWS2023/negbin_gamlocal/cm/test_window_2020/local_negbin_2020.parquet")

write_parquet(local.tweedie.2020, 
              "Brandt_VIEWS2023/tweedie_gamlocal/cm/test_window_2020/local_tweedie_2020.parquet")

write_parquet(tensor.poisson.2020, 
              "Brandt_VIEWS2023/poisson_gamtensor/cm/test_window_2020/tensor_poisson_2020.parquet")

write_parquet(tensor.negbin.2020, 
              "Brandt_VIEWS2023/negbin_gamtensor/cm/test_window_2020/tensor_negbin_2020.parquet")

write_parquet(tensor.tweedie.2020, 
              "Brandt_VIEWS2023/tweedie_gamtensor/cm/test_window_2020/tensor_tweedie_2020.parquet")

write_parquet(glmm.poisson.2020, 
              "Brandt_VIEWS2023/poisson_glmm/cm/test_window_2020/glmm_poisson_2020.parquet")

write_parquet(glmm.negbin.2020, 
              "Brandt_VIEWS2023/negbin_glmm/cm/test_window_2020/glmm_negbin_2020.parquet")

write_parquet(glmm.tweedie.2020, 
              "Brandt_VIEWS2023/tweedie_glmm/cm/test_window_2020/glmm_tweedie_2020.parquet")


#####################
# 2021
#####################
load("ViEWS2-2021.RData")

local.poisson.2021 <- reorg(local.P.2021, cnum=cnum , k1=492, k2=12)
local.negbin.2021 <- reorg(local.NB.2021, cnum=cnum, k1=492, k2=12)
local.tweedie.2021 <- reorg(local.TW.2021, cnum=cnum, k1=492, k2=12)

tensor.poisson.2021 <- reorg(tensor.P.2021, cnum=cnum, k1=492, k2=12)
tensor.negbin.2021 <- reorg(tensor.NB.2021, cnum=cnum, k1=492, k2=12)
tensor.tweedie.2021 <- reorg(tensor.TW.2021, cnum=cnum, k1=492, k2=12)

glmm.poisson.2021 <- reorg(glmm.P.2021, cnum=cnum, k1=492, k2=12)
glmm.negbin.2021 <- reorg(glmm.NB.2021, cnum=cnum, k1=492, k2=12)
glmm.tweedie.2021 <- reorg(glmm.TW.2021, cnum=cnum, k1=492, k2=12)

# Deal with the NA and upper limits
local.poisson.2021$outcome <- trunctop(local.poisson.2021$outcome)
local.negbin.2021$outcome <- trunctop(local.negbin.2021$outcome)
local.tweedie.2021$outcome <- trunctop(local.tweedie.2021$outcome)
tensor.poisson.2021$outcome <- trunctop(tensor.poisson.2021$outcome)
tensor.negbin.2021$outcome <- trunctop(tensor.negbin.2021$outcome)
tensor.tweedie.2021$outcome <- trunctop(tensor.tweedie.2021$outcome)
glmm.poisson.2021$outcome <- trunctop(glmm.poisson.2021$outcome)
glmm.negbin.2021$outcome <- trunctop(glmm.negbin.2021$outcome)
glmm.tweedie.2021$outcome <- trunctop(glmm.tweedie.2021$outcome)


# Write the outputs
write_parquet(local.poisson.2021, 
              "Brandt_VIEWS2023/poisson_gamlocal/cm/test_window_2021/local_poisson_2021.parquet")

write_parquet(local.negbin.2021, 
              "Brandt_VIEWS2023/negbin_gamlocal/cm/test_window_2021/local_negbin_2021.parquet")

write_parquet(local.tweedie.2021, 
              "Brandt_VIEWS2023/tweedie_gamlocal/cm/test_window_2021/local_tweedie_2021.parquet")

write_parquet(tensor.poisson.2021, 
              "Brandt_VIEWS2023/poisson_gamtensor/cm/test_window_2021/tensor_poisson_2021.parquet")

write_parquet(tensor.negbin.2021, 
              "Brandt_VIEWS2023/negbin_gamtensor/cm/test_window_2021/tensor_negbin_2021.parquet")

write_parquet(tensor.tweedie.2021, 
              "Brandt_VIEWS2023/tweedie_gamtensor/cm/test_window_2021/tensor_tweedie_2021.parquet")

write_parquet(glmm.poisson.2021, 
              "Brandt_VIEWS2023/poisson_glmm/cm/test_window_2021/glmm_poisson_2021.parquet")

write_parquet(glmm.negbin.2021, 
              "Brandt_VIEWS2023/negbin_glmm/cm/test_window_2021/glmm_negbin_2021.parquet")

write_parquet(glmm.tweedie.2021, 
              "Brandt_VIEWS2023/tweedie_glmm/cm/test_window_2021/glmm_tweedie_2021.parquet")

#####################
# 2022
#####################
load("ViEWS2-2022.RData")

local.poisson.2022 <- reorg(local.P.2022, cnum=cnum , k1=492, k2=12)
local.negbin.2022 <- reorg(local.NB.2022, cnum=cnum, k1=492, k2=12)
local.tweedie.2022 <- reorg(local.TW.2022, cnum=cnum, k1=492, k2=12)

tensor.poisson.2022 <- reorg(tensor.P.2022, cnum=cnum, k1=492, k2=12)
tensor.negbin.2022 <- reorg(tensor.NB.2022, cnum=cnum, k1=492, k2=12)
tensor.tweedie.2022 <- reorg(tensor.TW.2022, cnum=cnum, k1=492, k2=12)

glmm.poisson.2022 <- reorg(glmm.P.2022, cnum=cnum, k1=492, k2=12)
glmm.negbin.2022 <- reorg(glmm.NB.2022, cnum=cnum, k1=492, k2=12)
glmm.tweedie.2022 <- reorg(glmm.TW.2022, cnum=cnum, k1=492, k2=12)

# Deal with the NA and upper limits
local.poisson.2022$outcome <- trunctop(local.poisson.2022$outcome)
local.negbin.2022$outcome <- trunctop(local.negbin.2022$outcome)
local.tweedie.2022$outcome <- trunctop(local.tweedie.2022$outcome)
tensor.poisson.2022$outcome <- trunctop(tensor.poisson.2022$outcome)
tensor.negbin.2022$outcome <- trunctop(tensor.negbin.2022$outcome)
tensor.tweedie.2022$outcome <- trunctop(tensor.tweedie.2022$outcome)
glmm.poisson.2022$outcome <- trunctop(glmm.poisson.2022$outcome)
glmm.negbin.2022$outcome <- trunctop(glmm.negbin.2022$outcome)
glmm.tweedie.2022$outcome <- trunctop(glmm.tweedie.2022$outcome)


# Write the outputs
write_parquet(local.poisson.2022, 
              "Brandt_VIEWS2023/poisson_gamlocal/cm/test_window_2022/local_poisson_2022.parquet")

write_parquet(local.negbin.2022, 
              "Brandt_VIEWS2023/negbin_gamlocal/cm/test_window_2022/local_negbin_2022.parquet")

write_parquet(local.tweedie.2022, 
              "Brandt_VIEWS2023/tweedie_gamlocal/cm/test_window_2022/local_tweedie_2022.parquet")

write_parquet(tensor.poisson.2022, 
              "Brandt_VIEWS2023/poisson_gamtensor/cm/test_window_2022/tensor_poisson_2022.parquet")

write_parquet(tensor.negbin.2022, 
              "Brandt_VIEWS2023/negbin_gamtensor/cm/test_window_2022/tensor_negbin_2022.parquet")

write_parquet(tensor.tweedie.2022, 
              "Brandt_VIEWS2023/tweedie_gamtensor/cm/test_window_2022/tensor_tweedie_2022.parquet")

write_parquet(glmm.poisson.2022, 
              "Brandt_VIEWS2023/poisson_glmm/cm/test_window_2022/glmm_poisson_2022.parquet")

write_parquet(glmm.negbin.2022, 
              "Brandt_VIEWS2023/negbin_glmm/cm/test_window_2022/glmm_negbin_2022.parquet")

write_parquet(glmm.tweedie.2022, 
              "Brandt_VIEWS2023/tweedie_glmm/cm/test_window_2022/glmm_tweedie_2022.parquet")


# Final misc checks to make sure trunctop did what it should -- these 
# should all be empty
table(local.negbin.2019[(is.na(local.negbin.2019$outcome)),1:2])
table(local.negbin.2021[(is.na(local.negbin.2021$outcome)),1:2])
table(tensor.poisson.2021[(is.na(tensor.poisson.2021$outcome)),1:2])
table(tensor.negbin.2019[(is.na(tensor.negbin.2019$outcome)),1:2])
table(tensor.negbin.2021[(is.na(tensor.negbin.2021$outcome)),1:2])

###########################################################
###########################################################
# Get the objects we need for further analysis and 
# dump all of the rest to free some memory 
# Lower storage!
###########################################################
###########################################################

# Keep the things only in the list we want!
ls18 <- ls()[grep("2018", ls())]
ls19 <- ls()[grep("2019", ls())]
ls20 <- ls()[grep("2020", ls())]
ls21 <- ls()[grep("2021", ls())]
ls22 <- ls()[grep("2022", ls())]

rm(list=setdiff(ls(), 
                c("reorg","countries", "cout", 
                  ls18, ls19, ls20, ls21, ls22)))

# Write out the final forecast results for later use and comparisons
save.image("ForecastDensities_2018-2022.RData")

