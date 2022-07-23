##
##
## Author : Emanuele Cordano
## Date: 2022 07 20
## This function creates a list containg a "zoo" object (time series)("rzoo") for each variable and each station and each global climate model
##
##
##
#
# > source("~/activity/2022/fbk/local/waterwise/rmd/functions.R")
# > rr <- rzoo_code()
# > class(rr)
# [1] "list"
# > names(rr)
# [1] "main"               "ylab"               "stepPlot"           "zoov"               "variable"
# [6] "station"            "downscaling_models" "rcm_names"          "dy"                 "scplot"
# [11] "gof"
# >
NULL
#' It generates a list with downscaled time series for the selected station from the values of the related Climate Projection model grid.
#'
#' @param dataset dataset containing input time series (time domain - monthly);
#' @param station id name of the station
#' @param global_climate_model acronym of the global climate model;
#' @param iclalibaration record used for dowscaling model (ARIMA) calibration
#' @param varsel selected vartiable, i.e one of \code{c("pr","tasmax","tasmin")}
#' @param order,... arguments for ARIMA model. See \code{\link{arima}} and \code{\link{predict.Arima}}
#' @param add_dygraph,add_scatteplot logical, graphic options.
#'
#' @author Emanuele Cordano
#' @date 2022 07 20
#'
#'
rzoo_code <- function(
  dataset=out50,
  station0="WWRLUN",
  global_climate_model="MPI-M-MPI-ESM-LR",
  ncalibration=12*12, ## 12 years##10*12, ## 10 years
  icalibration=1:ncalibration,
  varsel="tasmax",
  order=c(3,1,1), ## order input for arima
  add_dygraph=TRUE,
  add_scatterplot=TRUE,
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
  outf$variable <- varsel
  outf$station <- station0
  outf$downscaling_models <- list()
  nns <- names(outf$zoov)[names(outf$zoov)!="observation"]
  outf$rcm_names <- nns
  iobs <- which(!is.na(outf$zoov$observation))
  iobs <- iobs
  icalibration <- icalibration

  zzoov <- outf$zoov
  if (varsel %in% c("tasmin","tasmax")) for (it in nns) {

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

    names(climate_preds) <- sprintf("%s_%s",it,names(climate_preds)) ##c(sprintf("%s_pred",it),sprintf("%s_se",it))

    outf$zoov <- cbind(outf$zoov,climate_preds)

  }
  ### DYGRAPH
  if (add_dygraph==TRUE) {

    outf$dy <- list()
    for (itn in outf$rcm_names) {
      ##itn <- "MPI-M-MPI-ESM-LR_CLMcom-CCLM4-8-17"  ####nns[1] ### "MPI-M-MPI-ESM-LR_CLMcom-ETH-COSMO-crCLIM-v1-1" ##"MPI-M-MPI-ESM-LR_CLMcom-CCLM4-8-17"
      if (varsel %in% c("tasmin","tasmax")) {
        nn <- c("observation",itn,sprintf("%s_%s",itn,c("min","pred","max")))
      } else {

        nn <- c("observation",itn)
      }

      zoo2 <- outf$zoov[,nn]
      mmodel <- itn
      if (varsel %in% c("tasmin","tasmax")) {
        rmsev <- rmse(obs=zoo2$observation,sim=zoo2[,sprintf("%s_pred",itn)])
      } else {

        rmsev <- rmse(obs=zoo2$observation,sim=zoo2[,itn])

      }
      main <- outf$main
      outf$dy[[itn]] <- dygraph(zoo2, main = main) %>%
        ##dySeries(sprintf("%s_%s",itn,c("min","pred","max")), label = label,color="black") %>%
        dySeries(c(itn), label = as.character(mmodel),color="green") %>%
        dySeries("observation",color="blue") %>% dyRangeSelector()
        if (varsel %in% c("tasmin","tasmax")) {
          label=sprintf("%s (RMSE=%f",it,rmsev)
          outf$dy[[itn]] <- outf$dy[[itn]] %>% dySeries(sprintf("%s_%s",itn,c("min","pred","max")), label = label,color="black")

        }
    }

  }

  if (!(varsel %in% c("tasmin","tasmax"))) add_scatterplot <- FALSE

  if (add_scatterplot==TRUE) {
    outf$scplot <- list()

    for (itn in outf$rcm_names) {


      main <- sprintf("%s at %s - model %s",outf$varsel,outf$station,itn)
      dftemp <- outf$zoov[,nn]
      names(dftemp) <- c("observation","grid_value","predicted_min","predicted_value","predicted_max")
      time_df <-  (index(dftemp))
      dftemp <- as.data.frame(dftemp)
      dftemp$time <- time_df
      dftemp <- dftemp %>% dplyr::filter(!is.na(observation)) %>% melt(id=c("observation","time","predicted_min","predicted_max"))
      gg <- ggplot(data=dftemp)+geom_point(aes(x=observation,y=value,group=variable,color=variable))+theme_bw()
      gg <- gg+geom_ribbon(aes(ymin=predicted_min,ymax=predicted_max,x=observation),fill="grey",alpha=0.4)+geom_abline(slope=1,intercept = 0)
      gg <- gg+ylab("model outcome")+ggtitle(main)+theme(aspect.ratio = 1)
      outf$scplot[[itn]] <- gg
    }
  }

  ###

  outf$gof <- list()

  for (itn in c("observation",outf$rcm_names)) {
    if (itn!="observation" & (varsel %in% c("tasmin","tasmax"))) {
      itn_pred <- sprintf("%s_%s",itn,"pred")
    } else {
      itn_pred <- itn
    }
    outf$gof[[itn]] <- gof(obs=outf$zoov$observation,sim=outf$zoov[,itn_pred])
  }
  nnn <- names(outf$gof)
  nngof <- c("RMSE","MAE","KGE")
  outf$gof <- as.data.frame(outf$gof)[nngof,] %>% t() %>% as.data.frame()
  outf$gof$model <- nnn
  outf$gof$station <- outf$station
  outf$gof$variable <-  outf$variable
  outf$gof <- outf$gof[,c("station","variable","model",nngof)]

  return(outf)

}

