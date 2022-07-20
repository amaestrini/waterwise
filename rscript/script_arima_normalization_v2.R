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
rzoo_code <- function(
  dataset=out50,
  station0="WWRLUN",
  global_climate_model="MPI-M-MPI-ESM-LR",
  ncalibration=10*12,
  icalibration=1:ncalibration,
  varsel="tasmax",
  order=c(3,1,1), ## order input for arima
  add_dygraph=TRUE,
  ...

)  {

  outf <- list()


  variables <- c("tasmin","tasmax","pr")


  out70 <- dataset %>% filter(experiment %in% c("historical","rcp85","observation"),gcm %in% c(global_climate_model,"observation"))  %>% dplyr::filter(station==station0)

  ##%>% filter(model %in% n_models)
  ttttz <- out70 %>% filter(variable %in% variables,!is.na(time)) %>% dplyr::select(value,variable,model,time)


  tasminz <- ttttz %>% filter(variable=="tasmin") %>% dplyr::select(-variable) %>% dcast(time  ~ model,fun.aggregate=mean)
  time <- tasminz$time
  tasminz <- tasminz %>% dplyr::select(-time) %>% as.zoo()
  index(tasminz) <- time
  ###

  tasmaxz <- ttttz %>% filter(variable=="tasmax") %>% dplyr::select(-variable) %>% dcast(time  ~ model,fun.aggregate=mean)
  time <- tasmaxz$time
  tasmaxz <- tasmaxz %>% dplyr::select(-time) %>% as.zoo()
  index(tasmaxz) <- time
  ####

  prz <- ttttz %>% filter(variable=="pr") %>% dplyr::select(-variable) %>% dcast(time  ~ model,fun.aggregate=mean)
  time <- prz$time
  prz <- prz %>% dplyr::select(-time) %>% as.zoo()
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
  outf$rcm_names <- nns
  iobs <- which(!is.na(outf$zoov$observation))
  iobs <- iobs
  icalibration <- icalibration

  zzoov <- outf$zoov
  for (it in nns) {
  print(it)
  ###  x0- <- as.numeric(outf$zoov$observation[iobs][icalibration])

    deltaclimate <- as.numeric(outf$zoov$observation[iobs][icalibration])####-as.numeric(outf$zoov[,it][iobs][icalibration])
    climate_model_value_all <- as.numeric(outf$zoov[,it])
    climate_model_value <-  as.numeric(outf$zoov[,it][iobs][icalibration])
  ###


    x <- normalizeGaussian(x=deltaclimate) %>% ts(frequency=12) ####ecdf(deltaclimate)(deltaclimate)*() %>% qnorm()
    xreg <- normalizeGaussian(x=climate_model_value,data=climate_model_value_all) %>% ts(frequency=12)#### ecdf(climate_model_value_all)(climate_model_value) %>% qnorm()
    newxreg <- normalizeGaussian(x=climate_model_value_all,data=climate_model_value_all) %>% ts(frequency=12)####ecdf(climate_model_value_all)(climate_model_value_all) %>% qnorm()
  ###
    outf$downscaling_models[[it]] <- arima(x=x,order=order,xreg=xreg,...)
    attr(outf$downscaling_models[[it]],"x") <- x
    attr(outf$downscaling_models[[it]],"xreg") <- xreg
    attr(outf$downscaling_models[[it]],"newxreg") <- newxreg

  ###
    mar <-  outf$downscaling_models[[it]]
    k <- qnorm(0.95)
    pred <- predict(outf$downscaling_models[[it]],newxreg=newxreg) ##,n.ahead=2) ##,n.ahead = nrow(outf$zoov))
    preds <- cbind(min=as.numeric(pred$pred)-k*as.numeric(pred$se),pred=as.numeric(pred$pred),max=as.numeric(pred$pred)+k*as.numeric(pred$se))
    climate_preds <- preds
    for (c in 1:ncol(preds)) {

      climate_preds[,c] <-  normalizeGaussian(x=preds[,c],data=deltaclimate,inverse=TRUE)###+as.numeric(outf$zoov[,it])  ###quantile(x=climate_model_value,probs=pnorm(preds[,c]),type=3)
    }
    ####
    ###ouf$zoo


    climate_preds <- as.data.frame(climate_preds) %>% as.zoo()
    names(climate_preds) <- c("min","pred","max")
    index(climate_preds) <- index(outf$zoov)
    print(names(climate_preds))
    print(it)
    names(climate_preds) <- sprintf("%s_%s",it,names(climate_preds)) ##c(sprintf("%s_pred",it),sprintf("%s_se",it))
    print(names(climate_preds))
    print("c")
    outf$zoov <- cbind(outf$zoov,climate_preds)

  }
  ### DYGRAPH
  if (add_dygraph==TRUE) {

    outf$dy <- list()
    for (itn in outf$rcm_names) {
      ##itn <- "MPI-M-MPI-ESM-LR_CLMcom-CCLM4-8-17"  ####nns[1] ### "MPI-M-MPI-ESM-LR_CLMcom-ETH-COSMO-crCLIM-v1-1" ##"MPI-M-MPI-ESM-LR_CLMcom-CCLM4-8-17"
      nn <- c("observation",itn,sprintf("%s_%s",itn,c("min","pred","max")))

      print(itn)
      print(names(outf$zoov))
      zoo2 <- outf$zoov[,nn]
      mmodel <- itn
      rmsev <- rmse(obs=zoo2$observation,sim=zoo2[,sprintf("%s_pred",itn)])
      label <- sprintf("%s (RMSE=%f)(dowsncaled)",mmodel,rmsev)
      main <- outf$main
      outf$dy[[itn]] <- dygraph(zoo2, main = main) %>%
        dySeries(sprintf("%s_%s",itn,c("min","pred","max")), label = label,color="black") %>%
        dySeries(c(itn), label = as.character(mmodel)) %>%
        dySeries("observation",color="blue") %>% dyRangeSelector()

    }


  }

  return(outf)

}

#####

station0="WWTRES"
selected_global_climate_model=selected_global_model
varsel="tasmin"
order=c(3,1,0)

###
rr <- rzoo_code(station0=station0,global_climate_model=selected_global_model,varsel=varsel,order=order)



stop("HERE")
itn <- "MPI-M-MPI-ESM-LR_CLMcom-CCLM4-8-17"  ####nns[1] ### "MPI-M-MPI-ESM-LR_CLMcom-ETH-COSMO-crCLIM-v1-1" ##"MPI-M-MPI-ESM-LR_CLMcom-CCLM4-8-17"
nn <- c("observation",itn,sprintf("%s_%s",itn,c("min","pred","max")))
zoo2 <- rr$zoov[,nn]






#zoo2$vmax <- as.numeric(zoo2[,sprintf("%s_pred",itn)])+k*as.numeric(zoo2[,sprintf("%s_se",itn)])
#zoo2$vmin <- as.numeric(zoo2[,sprintf("%s_pred",itn)])-k*as.numeric(zoo2[,sprintf("%s_se",itn)])
zoo2 <- zoo2[,names(zoo2)!=sprintf("%s_se",itn)]

mmodel <- itn
rmsev <- rmse(obs=zoo2$observation,sim=zoo2[,sprintf("%s_pred",itn)])
label <- sprintf("%s (RMSE=%f)(dowsncaled)",mmodel,rmsev)
main <-rr$main
dy00 <- dygraph(zoo2, main = main) %>%
  dySeries(sprintf("%s_%s",itn,c("min","pred","max")), label = label,color="black") %>%
  dySeries(c(itn), label = as.character(mmodel)) %>%
  dySeries("observation",color="blue") %>% dyRangeSelector()

#
# ###
# zoo2 <- rr$zoov
# main <- rr$main
# dy0 <- dygraph(zoo2, main = main)
#
# for (itn in rr$rcm_names) {
#   print(itn)
#   mmodel <- itn
#   rmsev <- rmse(obs=zoo2$observation,sim=zoo2[,sprintf("%s_pred",itn)])
#   label <- sprintf("%s (RMSE=%f)(dowsncaled)",mmodel,rmsev)
#   dy0 <- dy0 %>% dySeries(sprintf("%s_%s",itn,c("min","pred","max")), label = label)
#   dy0 <- dy0 %>% dySeries(c(itn), label = as.character(mmodel))
# }
# dy0 <- dy0 %>% dySeries("observation",color="blue") %>% dyRangeSelector()
#
#
