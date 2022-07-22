#
#```{r zoo001,echo=FALSE,eval=TRUE,message=FALSE,warning=FALSE,print=FALSE,results="hold"}
library(zoo)
library(dygraphs)
library(stringr)
library(dplyr)
library(reshape2)
library(hydroGOF)
library(lubridate)
library(quantreg)


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
rzoo_code <- function(
  data_in=out50,
  station0="WWRLUN",
  in_selected_global_model="MPI-M-MPI-ESM-LR",
  varsel="temp"
)  {

  outf <- list()


  variables <- c("tasmin","tasmax","pr")


  out70 <- data_in %>% filter(experiment %in% c("historical","rcp85","observation"),gcm %in% c(in_selected_global_model,"observation"))  %>% dplyr::filter(station==station0)

  ##%>% filter(model %in% n_models)
  ttttz <- out70 %>% filter(variable %in% variables,!is.na(time)) %>% select(value,variable,model,time)


  tasminz <- ttttz %>% filter(variable=="tasmin") %>% select(-variable) %>% dcast(time  ~ model,fun.aggregate=mean)
  time <- tasminz$time
  tasminz <- tasminz %>% select(-time) %>% as.zoo()
  index(tasminz) <- time
  ###

  tasmaxz <- ttttz %>% filter(variable=="tasmax") %>% select(-variable) %>% dcast(time  ~ model,fun.aggregate=mean)
  time <- tasmaxz$time
  tasmaxz <- tasmaxz %>% select(-time) %>% as.zoo()
  index(tasmaxz) <- time
  ####

  prz <- ttttz %>% filter(variable=="pr") %>% select(-variable) %>% dcast(time  ~ model,fun.aggregate=mean)
  time <- prz$time
  prz <- prz %>% select(-time) %>% as.zoo()
  index(prz) <- time

  ## ZOO OBJECTS
  if (varsel=="temp") {
    tempz <- (tasminz+tasmaxz)/2
    outf$main <- "Monthly averaged (max+min)/2 temperature vs time (%s) at %s" %>% sprintf(in_selected_global_model,station0)
    outf$ylab <- "Temperature [deg C]"
    outf$stepPlot <- FALSE
    outf$zoov <- tempz
    #dygraph(zoov,main=main,xlab="Time",ylab=ylab) %>% dyRangeSelector() %>% dySeries("observation", stepPlot = stepPlot, color = "blue")

  }

  if (varsel=="tasmin") {
    outf$main <- "Monthly averaged Daily Minium temperature vs time (%s) at %s" %>% sprintf(in_selected_global_model,station0)
    outf$ylab <- "Temperature [deg C]"
    outf$stepPlot <- FALSE
    outf$zoov <- tasminz
    #dygraph(zoov,main=main,xlab="Time",ylab=ylab) %>% dyRangeSelector() %>% dySeries("observation", stepPlot = stepPlot, color = "blue")

  }

  if (varsel=="tasmax") {
    outf$main <- "Monthly averaged Daily Maximum temperature vs time (%s) at %s" %>% sprintf(in_selected_global_model,station0)
    outf$ylab <- "Temperature [deg C]"
    outf$stepPlot <- FALSE
    outf$zoov <- tasmaxz
    #dygraph(zoov,main=main,xlab="Time",ylab=ylab) %>% dyRangeSelector() %>% dySeries("observation", stepPlot = stepPlot, color = "blue")

  }

  if (varsel=="diff_temp") {
    difftempz <- (tasmaxz-tasminz)
    outf$main <- "Monthly averaged diurnal thermal range vs time (%s) at %s" %>% sprintf(in_selected_global_model,station0)
    outf$ylab <- "Temperature diff. [deg C]"
    outf$stepPlot <- FALSE
    outf$zoov <- difftempz
    #dygraph(zoov,main=main,xlab="Time",ylab=ylab) %>% dyRangeSelector() %>% dySeries("observation", stepPlot = stepPlot, color = "blue")

  }

  if (varsel=="pr") {
    outf$main <- "Monthly averaged daily precipitation vs time (%s) at %s" %>% sprintf(in_selected_global_model,station0)
    outf$ylab <- "Daily precipitation [mm/day]"
    outf$zoov <- prz
    outf$stepPlot <- TRUE
    # dygraph(zoov,main=main,xlab="Time",ylab=ylab) %>% dyRangeSelector() %>% dySeries("observation", stepPlot = stepPlot,color = "blue")
  }


  ### GOF
  outf$gof <- outf$zoov %>% as.data.frame() %>%  apply(MARGIN=2,gof,obs=outf$zoov$observation,na.rm=TRUE,simplify=FALSE) %>% do.call(what="cbind") %>% as.data.frame()
  names(outf$gof) <- names(outf$zoov)
  ####

  #### for bokeh module
  outf$ttt0 <- as.data.frame(outf$zoov)
  outf$ttt0$month <- month(as.Date(index(outf$zoov)))
  outf$ttt0$year <- year(as.Date(index(outf$zoov)))
  outf$ttt0$time <- as.Date(index(outf$zoov))
  outf$ttt0 <- melt(outf$ttt0,id=rev(c("year","month","time","observation")))
  names(outf$ttt0)[names(outf$ttt0)=="variable"] <- "model"
  outf$ttt <- outf$ttt0 %>% filter(!is.na(observation))


  #### MODEL
  ##outf$dowscaling_model <- glm(observation ~ value+model+factor(month),data=outf$ttt)
  ##outf$dowscaling_model <- glm(observation ~ value+model,data=outf$ttt)
  outf$dowscaling_model <- rq(observation ~ value+model,data=outf$ttt)
  ####
  outf$ttt0$predicted <- predict(outf$dowscaling_model,newdata=outf$ttt0)
  outf$ttt <- outf$ttt0 %>% filter(!is.na(observation))

  ####
  outf$zoov_modeled <- outf$ttt0 %>% select(model,time,predicted) %>% dcast(time  ~ model,fun.aggregate=mean)
  time <- outf$zoov_modeled$time
  outf$zoov_modeled  <- outf$zoov_modeled %>% select(-time) %>% as.zoo()
  index(outf$zoov_modeled) <- time
  outf$zoov_modeled$observation <- outf$zoov$observation
  outf$main2 <- sprintf("%s at %s (gcm: %s)",varsel,station0,in_selected_global_model)
  varsel="temp"
  ### GOF
  outf$gof_modeled <- outf$zoov_modeled %>% as.data.frame() %>%  apply(MARGIN=2,gof,obs=outf$zoov_modeled$observation,na.rm=TRUE,simplify=FALSE) %>% do.call(what="cbind") %>% as.data.frame()
  names(outf$gof_modeled) <- names(outf$zoov_modeled)
  ####

  return(outf)
}


##outf <<- rzoo_code()

##})


rr <<- rzoo_code()

test_temp  <- rr$ttt %>% filter(model %in% "MPI-M-MPI-ESM-LR_CLMcom-CCLM4-8-17")
itesting <- 1:120

# out <- arima(test_temp$observation,xreg=test_temp$predicted)
# mar <- arima(test_temp$observation[itesting],xreg=test_temp$predicted[itesting],order=c(12,0,0))
#
# test_temp$predicted2 <- predict(mar,newxreg=test_temp$predicted)$pred %>% as.numeric()

##########
mar2 <- arima(test_temp$observation[itesting],order=c(3,0,3),xreg=test_temp$value[itesting])
###pp <- predict(mar2,n.ahead = nrow(test_temp))
# > as.character(unique(rr$ttt$model))
# [1] "MPI-M-MPI-ESM-LR_CLMcom-CCLM4-8-17"            "MPI-M-MPI-ESM-LR_CLMcom-ETH-COSMO-crCLIM-v1-1"
# [3] "MPI-M-MPI-ESM-LR_DMI-HIRHAM5"                  "MPI-M-MPI-ESM-LR_KNMI-RACMO22E"
# [5] "MPI-M-MPI-ESM-LR_MOHC-HadREM3-GA7-05"          "MPI-M-MPI-ESM-LR_MPI-CSC-REMO2009"
# [7] "MPI-M-MPI-ESM-LR_SMHI-RCA4"                    "MPI-M-MPI-ESM-LR_UHOH-WRF361H"

test_temp0  <- rr$ttt0 %>% filter(model %in% "MPI-M-MPI-ESM-LR_MOHC-HadREM3-GA7-05" ) ##"MPI-M-MPI-ESM-LR_DMI-HIRHAM5") ##"MPI-M-MPI-ESM-LR_CLMcom-CCLM4-8-17")
pp <- predict(mar2,n.ahead = nrow(test_temp0),newxreg=test_temp0$value)
test_temp0$predict2 <- pp$pred
test_temp0$se <- pp$se
k <- qnorm(0.95)
test_temp0$vmax <-  test_temp0$predict2+k*test_temp0$se
test_temp0$vmin <-  test_temp0$predict2-k*test_temp0$se

itime <- test_temp0$time
mmodel <- test_temp$model[1]
test_temp0z <- test_temp0 %>% select(-time,-model) %>% select(observation,predict2,vmin,vmax,value) %>% as.zoo()
index(test_temp0z) <- itime

rmsev <- rmse(obs=test_temp0z$observation,sim=test_temp0z$predict2)
label <- sprintf("%s (RMSE=%f)(dowscaled)",mmodel,rmsev)
dy0 <- dygraph(test_temp0z, main = "Dygraph title") %>%
  dySeries(c("vmin", "predict2", "vmax"), label = label) %>%
 dySeries(c("value"), label = as.character(mmodel)) %>%
  dyRangeSelector()
####





stop("HERE")
##########
library(ggplot2)

test_temp2 <- test_temp0 %>% select(time,predict2,observation,vmin,vmax) %>% melt(id=c("time","vmax","vmin")) ###,"se")
gg <- ggplot()+geom_line(aes(x=time,y=value,group=variable,color=variable),data=test_temp2)+theme_bw()+geom_ribbon(aes(ymin=test_temp2$vmin, ymax=test_temp2$vmax,x=test_temp2$time),fill="grey70",alpha=0.5)
##+geom_line(aes(x=time,y=predict2),data=test_temp)
##gg <- ggplot()+geom_line(aes(x=time,y=observation),data=test_temp)+theme_bw()+geom_line(aes(x=time,y=predict2),data=test_temp)

### GUARDARE QUI:
#https://stackoverflow.com/questions/63929512/dygraphs-in-r-plot-ribbon-and-mean-line-of-different-groups
## SEARCH FOR ARIMA ::

out_temp0 <- test_temp2 %>% select(time,variable,value) %>% dcast(time ~ variable)
itime <- out_temp0$time
out_temp0 <- out_temp %>% select(-time) %>% as.zoo()
