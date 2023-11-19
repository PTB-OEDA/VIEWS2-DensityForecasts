---
output:
  pdf_document: default
  html_document: default
---
# VIEWS2-DensityForecasts
Patrick T. Brandt 

## VIEWS 2 Density Forecasts

These are the files used to generate the models for the [2023/24 VIEWS Prediction Challenge](https://viewsforecasting.org/prediction-competition-2/).  These setup and the data and replicate the analysis presented there on Bayesian Density Forecasts. You access the data used directly from [here](https://www.dropbox.com/sh/yxk5w04p2e1xtqk/AACU2k5EUOuEeMq2kZ3gpZZwa?dl=0) since they are large not included in this repo.

The files here have only lightly been edited from what was used for the forecasts in October 2023

1.  Main file for the data setup and estimation in R is the `Brandt-VIEWS2-Estimation.Rmd`.  This sets up the data and does the basic model exploration and selection via AIC and some in-sample CRPS computations.  This generates the forecast output for 2018.

2.  Forecasts for the 2019-2021 periods are handled in separate R scripts:
    1.  `combined-rolled-2019.R` for 2019
    2.  `combined-rolled-2020.R` for 2020
    3.  `combined-rolled-2021.R` for 2021

3.  Organization and reporting of the output for what was submitted to the competition is in `reorg.R` and `forc-ord.r` scripts (to be called in this order).  Final output comes in a file named `ForecastDensities_2018-2021.RData`. 

The forecasts from these methods are not included here.  They are easy to generate and can be done by running the above on a *single* CPU core in less than 24 hours.  So this is not some crazy computational job, nor has it really been optimized (the time could be greatly trimmed with some smart / easy parallelization of the models estimation and forecast generation). But the goal is to create simple, easy to follow replication code here so others can build on it.


Questions to Patrick Brandt at pbrandt_at_utdallas_dot_edu.

