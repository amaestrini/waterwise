#
#```{r zoo001,echo=FALSE,eval=TRUE,message=FALSE,warning=FALSE,print=FALSE,results="hold"}
library(zoo)
library(dygraphs)
library(stringr)
library(dplyr)
library(reshape2)
library(hydroGOF)
library(lubridate)
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
rzoo_code <- function(
  dataset=out50,
  station0="WWRLUN",
  global_climate_model="MPI-M-MPI-ESM-LR",
  ncalibration=10*12,
  icalibration=1:ncalibration,
  varsel="temp"
)  {

  outf <- list()


  variables <- c("tasmin","tasmax","pr")


  out70 <- dataset %>% filter(experiment %in% c("historical","rcp85","observation"),gcm %in% c(global_climate_model,"observation"))  %>% dplyr::filter(station==station0)

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
    outf$main <- "Monthly averaged (max+min)/2 temperature vs time (%s) at %s" %>% sprintf(global_climate_model,station0)
    outf$ylab <- "Temperature [deg C]"
    outf$stepPlot <- FALSE
    outf$zoov <- tempz
    #dygraph(zoov,main=main,xlab="Time",ylab=ylab) %>% dyRangeSelector() %>% dySeries("observation", stepPlot = stepPlot, color = "blue")

  }

  if (varsel=="tasmin") {
    outf$main <- "Monthly averaged Daily Minium temperature vs time (%s) at %s" %>% sprintf(global_climate_model,station0)
    outf$ylab <- "Temperature [deg C]"
    outf$stepPlot <- FALSE
    outf$zoov <- tasminz
    #dygraph(zoov,main=main,xlab="Time",ylab=ylab) %>% dyRangeSelector() %>% dySeries("observation", stepPlot = stepPlot, color = "blue")

  }

  if (varsel=="tasmax") {
    outf$main <- "Monthly averaged Daily Maximum temperature vs time (%s) at %s" %>% sprintf(global_climate_model,station0)
    outf$ylab <- "Temperature [deg C]"
    outf$stepPlot <- FALSE
    outf$zoov <- tasmaxz
    #dygraph(zoov,main=main,xlab="Time",ylab=ylab) %>% dyRangeSelector() %>% dySeries("observation", stepPlot = stepPlot, color = "blue")

  }

  if (varsel=="diff_temp") {
    difftempz <- (tasmaxz-tasminz)
    outf$main <- "Monthly averaged diurnal thermal range vs time (%s) at %s" %>% sprintf(global_climate_model,station0)
    outf$ylab <- "Temperature diff. [deg C]"
    outf$stepPlot <- FALSE
    outf$zoov <- difftempz
    #dygraph(zoov,main=main,xlab="Time",ylab=ylab) %>% dyRangeSelector() %>% dySeries("observation", stepPlot = stepPlot, color = "blue")

  }

  if (varsel=="pr") {
    outf$main <- "Monthly averaged daily precipitation vs time (%s) at %s" %>% sprintf(global_climate_model,station0)
    outf$ylab <- "Daily precipitation [mm/day]"
    outf$zoov <- prz
    outf$stepPlot <- TRUE
    # dygraph(zoov,main=main,xlab="Time",ylab=ylab) %>% dyRangeSelector() %>% dySeries("observation", stepPlot = stepPlot,color = "blue")
  }

  outf$downscaling_models <- list()
  nns <- names(outf$zoov)[names(outf$zoov)!="observation"]
  iobs <- which(!is.na(outf$zoov$observation))
  iobs <<- iobs
  icalibration <<- icalibration

  zzoov <<- outf$zoov
  for (it in nns) {
  print(it)
    x <<- as.numeric(outf$zoov$observation[iobs][icalibration])
    climate_model_value <<- as.numeric(outf$zoov[,it][iobs][icalibration])
    outf$downscaling_models[[it]] <- arima(x=x,order=c(3,1,1),xreg=climate_model_value)
    ## GAUSSIANIZATION

    ##
    ##
    ##
    ###https://otexts.com/fpp2/seasonal-arima.html
    mar <<-  outf$downscaling_models[[it]]
    climate_model_value_all <- as.numeric(outf$zoov[,it])
    pred <<- predict(outf$downscaling_models[[it]],newxreg=climate_model_value_all) ##,n.ahead=2) ##,n.ahead = nrow(outf$zoov))
    preds <- cbind(pred=as.numeric(pred$pred),se=as.numeric(pred$se)) %>% as.zoo()
    preds <<- preds
    index(preds) <- index(outf$zoov)
    names(preds) <- sprintf("%s_%s",it,names(preds)) ##c(sprintf("%s_pred",it),sprintf("%s_se",it))
    outf$zoov <- cbind(outf$zoov,preds)
 #   outf$zoov <- cbind(outf$zoov,pred=as.zoo(as.numeric(pred$pred)))
 #    outf$zoov <- cbind(outf$zoov,se=as.zoo(as.numeric(pred$se)))
#    names(outf$zoov)[names(outf$zoov)=="pred"] <- sprintf("%s_pred",it)
    #names(outf$zoov)[names(outf$zoov)=="se"] <- sprintf("%s_se",it)
  #  outf$zoov[,sprintf("%s_pred",it)] <- as.numeric(NA)
  #  outf$zoov[,sprintf("%s_pred",it)] <- as.numeric(pred$pred)
  #  outf$zoov[,sprintf("%s_se",it)] <- as.numeric(pred$se)
  }

  ###

  return(outf)

}

### VERIFICARE::
(ecdf(climate_model_value)(climate_model_value) %>% quantile(x=climate_model_value))-climate_model_value
### verificare questo:
vv <- (ecdf(climate_model_value)(climate_model_value) %>% quantile(x=climate_model_value,type=3))-climate_model_value

rr <- rzoo_code()

##

#
itn <- "MPI-M-MPI-ESM-LR_CLMcom-ETH-COSMO-crCLIM-v1-1" ##"MPI-M-MPI-ESM-LR_CLMcom-CCLM4-8-17"
nn <- c("observation",itn,sprintf("%s_%s",itn,c("pred","se")))
k <- qnorm(0.95)
zoo2 <- rr$zoov[,nn]
zoo2$vmax <- as.numeric(zoo2[,sprintf("%s_pred",itn)])+k*as.numeric(zoo2[,sprintf("%s_se",itn)])
zoo2$vmin <- as.numeric(zoo2[,sprintf("%s_pred",itn)])-k*as.numeric(zoo2[,sprintf("%s_se",itn)])
zoo2 <- zoo2[,names(zoo2)!=sprintf("%s_se",itn)]

mmodel <- itn
rmsev <- rmse(obs=zoo2$observation,sim=zoo2[,sprintf("%s_pred",itn)])
label <- sprintf("%s (RMSE=%f)(dowsncaled)",mmodel,rmsev)
dy0 <- dygraph(zoo2, main = "Dygraph title") %>%
  dySeries(c("vmin",sprintf("%s_pred",itn) , "vmax"), label = label) %>%
  dySeries(c(itn), label = as.character(mmodel)) %>%
  dyRangeSelector()

stop("QUI")
####




#
#   ### GOF
#   outf$gof <- outf$zoov %>% as.data.frame() %>%  apply(MARGIN=2,gof,obs=outf$zoov$observation,na.rm=TRUE,simplify=FALSE) %>% do.call(what="cbind") %>% as.data.frame()
#   names(outf$gof) <- names(outf$zoov)
#   ####
#
#   #### for bokeh module
#   outf$ttt0 <- as.data.frame(outf$zoov)
#   outf$ttt0$month <- month(as.Date(index(outf$zoov)))
#   outf$ttt0$year <- year(as.Date(index(outf$zoov)))
#   outf$ttt0$time <- as.Date(index(outf$zoov))
#   outf$ttt0 <- melt(outf$ttt0,id=rev(c("year","month","time","observation")))
#   names(outf$ttt0)[names(outf$ttt0)=="variable"] <- "model"
#   outf$ttt <- outf$ttt0 %>% filter(!is.na(observation))
#
#
#   #### MODEL
#   ##outf$dowscaling_model <- glm(observation ~ value+model+factor(month),data=outf$ttt)
#   ##outf$dowscaling_model <- glm(observation ~ value+model,data=outf$ttt)
#   outf$dowscaling_model <- rq(observation ~ value+model,data=outf$ttt)
#   ####
#   outf$ttt0$predicted <- predict(outf$dowscaling_model,newdata=outf$ttt0)
#   outf$ttt <- outf$ttt0 %>% filter(!is.na(observation))
#
#   ####
#   outf$zoov_modeled <- outf$ttt0 %>% select(model,time,predicted) %>% dcast(time  ~ model,fun.aggregate=mean)
#   time <- outf$zoov_modeled$time
#   outf$zoov_modeled  <- outf$zoov_modeled %>% select(-time) %>% as.zoo()
#   index(outf$zoov_modeled) <- time
#   outf$zoov_modeled$observation <- outf$zoov$observation
#   outf$main2 <- sprintf("%s at %s (gcm: %s)",varsel,station0,global_climate_model)
#   varsel="temp"
#   ### GOF
#   outf$gof_modeled <- outf$zoov_modeled %>% as.data.frame() %>%  apply(MARGIN=2,gof,obs=outf$zoov_modeled$observation,na.rm=TRUE,simplify=FALSE) %>% do.call(what="cbind") %>% as.data.frame()
#   names(outf$gof_modeled) <- names(outf$zoov_modeled)
#   ####
#
#   return(outf)
# }


##outf <<- rzoo_code()

##})








##
stop("HERE")
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
