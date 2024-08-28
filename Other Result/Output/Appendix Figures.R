#Code to obtain all figures as in the appendix

library(tidyverse)

load("Output/Output 1971-2000 Age 0+.RData")
source("Functions.R")

accuracy_heatmap(check_for, "RMSE")
plot_agesfitting(check_tot%>% 
                   filter(Year>1970))

load("Output/Output 1981-2000 Age 0+.RData")
source("Functions.R")

accuracy_heatmap(check_for, "RMSE")
plot_agesfitting(check_tot%>% 
                   filter(Year>1980))

load("Output/Output 1971-2000 Age 40+.RData")
source("Functions.R")

accuracy_heatmap(check_for, "RMSE")
plot_agesfitting(check_tot %>% 
                   filter(Year>1970))

load("Output/Output 1981-2000 Age 40+.RData")
source("Functions.R")

accuracy_heatmap(check_for, "RMSE")
plot_agesfitting(check_tot%>% 
                   filter(Year>1980))

