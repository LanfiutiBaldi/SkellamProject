#Code to obtain all figures as in the paper

library(tidyverse)

#FIG 1
load("Output/Output 1961-2000 Age 0+.RData")
plot_initial(dt)

#FIG 2
load("Output/Output 1961-2000 Age 0+.RData")
source("Functions.R")

accuracy_heatmap(check_for, "RMSE")

#FIG 3
load("Output/Output 1961-2000 Age 40+.RData")
source("Functions.R")

accuracy_heatmap(check_for, "RMSE")

#FIG 4
load("Output/Output 1961-2000 Age 0+.RData")
source("Functions.R")

plot_agesfitting(check_tot)

#FIG 5
load("Output/Output 1961-2000 Age 40+ (A).RData")
source("Functions.R")

plot_agesfitting(check_tot)



