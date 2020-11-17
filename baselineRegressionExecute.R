# clear memory; set working directory
rm(list = ls())
setwd("~/Box/Jeff Work/Jeff work/ximing/codeR")
########################################################################################
# choose parameters for model;
########################################################################################
# which crop should be used
# 	crop 1: corn
# 	crop 2: soybeans
cropList = c(1,2)

# which trimester of growing season should be used (0 is entire growing season)
trimesterList = c(0,2)

# which model specification should be used
# 	model 0: levels with no control
# 	model 1: levels after taking out soil effects - soil model
# 	model 2: levels after taking out soil and weather effects - weatherSoil model
# 	model 3: levels after taking out time effects (trends and year FE) - time model
# 	model 4: levels after taking out soil, weather, time effects (trends and year FE) - weatherSoilTime model
# 	model 11: shocks - taking out time effects (trends and year FE) as well county fixed effects - timeCountyFE model
# 	model 12: shocks - taking out weather, time effects (trends and year FE) as well county fixed effects - weatherTimeCountyFE model
modelList = c(0,1,2,3,4,11,12)

# which subset of the data should be used
# 	subset 0: all observations
# 	subset 1: hotter third of observations
# 	subset 2: colder third of observations
# 	subset 11: first third of years
# 	subset 12: second third of years
# 	subset 13: last third of years
subsetList = c(0,1,2,11,12,13)

########################################################################################
# execute code for given model choices;
########################################################################################
source('BaselineRegression_version2018.r')
