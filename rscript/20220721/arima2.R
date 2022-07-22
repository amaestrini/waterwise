#
#```{r zoo001,echo=FALSE,eval=TRUE,message=FALSE,warning=FALSE,print=FALSE,results="hold"}
library(zoo)
library(dygraphs)
library(stringr)
library(dplyr)
library(reshape2)
library(hydroGOF)
library(lubridate)
library(RMAWGEN)
#library(quantreg)
#library(forecast)

selected_global_model <- "MPI-M-MPI-ESM-LR"
station0 <- "WWRLUN"
## INPUT PANEL SELECTION
####
it_station <- c("WWTRES","WWRLUN","T0099","T0083","T0236","SMICH","T0090","T0169")
out50_file <- '/home/ecor/activity/2022/fbk/local/waterwise/data/climate_historical_rcp_%s_with_observations.rds' %>%
  sprintf(paste(it_station,collapse="_"))

## SAVE THE DATA
out50 <- readRDS(out50_file)
###


###
