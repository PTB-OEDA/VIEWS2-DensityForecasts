# reorg.R script
# Transforms the lists of forecast lists by model and period to be 
# in the format VIEWS 2.0 accepts
#
# Patrick T. Brandt 
#
# 20230926 : Revised version as a stand alone file
# 20231002 : Added type check on the output for the country_id to force integer.
# 
# Now format the output for the submission for 2018

# Quick function to do that setup for ViEWS, per instructions
# cm-level files should have these column names: 
# "month_id", "country_id", "draw", "outcome"

# x = forecast list array made above
# cnum = factors for the countries
# k1 = starting time period
# k2 = integer sequence of forecast periods
reorg <- function(x, cnum=lastobs[,1], k1=454, k2=length(x))
{ 
  # Get constants
  tmp <- dim(x[[1]])
  draws <- tmp[1]; cty <- as.integer(tmp[2])
  drawvec <- as.integer(1:draws)
  
  # Do the first one
  xtmp <- cbind(rep(k1+1, draws*cty),    # month_id
                rep(cnum, each=draws),   # country_id
                rep(drawvec, cty),       # draw idx
                matrix(x[[1]], ncol=1))  # forecast value
  
  # Now iterate over the periods
  for(i in 2:k2)
  {
    x1 <- cbind(rep(k1+i, draws*cty),    # month_id
                rep(cnum, each=draws),   # country_id
                rep(drawvec, cty),       # draw idx
                matrix(x[[i]], ncol=1))  # forecast value
    xtmp <- rbind(xtmp, x1)
  }
  
  colnames(xtmp) <- c("month_id", "country_id",
                      "draw", "outcome")
  return(as.data.frame(xtmp))
}

# Objects we need in the output:
#
# local.*.2018
# tensor.*.2018
# glmm.*.2018
