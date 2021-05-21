### R version 4.0.3
## code for generating Maxent predictions at different model settings
# requires output of native_regions.R and occurrence_records.R

library(spThin)
library(ENMeval)
library(dismo)
library(sp)

rem_correlations <-  function(predictors, backgroundpoints, selectEVs, correlationthreshold){
  # variables that are correlated with a correlation coefficient higher than the 'correlationthreshold' will be  excluded
  
  selected_predictors<-NULL
  for(p in 1:nrow(selectEVs)){
    
    if(length(selected_predictors)==0){
      selected_predictors <- raster::subset(predictors, grep(selectEVs$variable[p], names(predictors), value = T))
    } else {
      selected_predictors_add <- raster::subset(predictors, grep(selectEVs$variable[p], names(predictors), value = T))
      selected_predictors<-stack(selected_predictors, selected_predictors_add)
      selected_predictors_add<-NULL
    }
  }
  
  
  # The first variablename is expected to be the variable with LOWEST contribution to the maxent model
  backgroundvalues <- data.frame(extract(selected_predictors, backgroundpoints))
  backgroundvalues<-na.exclude(backgroundvalues)
  # calculate correlation coefficients
  correlation.coefficients <- sqrt(cor(backgroundvalues, method = "spearman")^2)
  tt1<-list()
  for(p in 1:nrow(selectEVs)){
    tt<-which(correlation.coefficients[p,] >= correlationthreshold)
    tt<-tt[which(tt!=p)]
    
    if(length(tt)>0){
      for (l in 1:length(tt)) {
        if((as.numeric(tt[l])>p)==TRUE){
          tt_add<-as.numeric(p)
          tt1<-c(tt1, tt_add)
        }
      }
      
    }
    
  }
  
  tt1<-as.numeric(tt1)
  tt1<-unique(tt1)
  
  if(length(tt1)>0){
    selectEVs<-selectEVs[-tt1,]
  }
  
  selected_predictors<-NULL
  for(p in 1:nrow(selectEVs)){
    
    if(length(selected_predictors)==0){
      selected_predictors <- raster::subset(predictors, grep(selectEVs$variable[p], names(predictors), value = T))
    } else {
      selected_predictors_add <- raster::subset(predictors, grep(selectEVs$variable[p], names(predictors), value = T))
      selected_predictors<-stack(selected_predictors, selected_predictors_add)
      selected_predictors_add<-NULL
    }
  }
  
  
  return(selected_predictors)
}

xToCoords <- function(x, longLat = NULL, sp = TRUE) {
  
  if (class(x) == 'data.frame' | class(x) == 'matrix') {
    if (is.null(longLat) & ncol(x) == 2) longLat <- 1:2
    x <- x[ , longLat]
    crs <- getCRS('wgs84', asCRS=TRUE)
  } else {
    crs <- sp::CRS(raster::projection(x))
    x <- sp::coordinates(x)
  }
  
  if (sp) x <- sp::SpatialPoints(x, proj4string=crs)
  x
  
}

elimCellDups <- function(x, r, longLat = NULL, priority = NULL) {
  
  # check CRS of points and raster
  if (class(x) %in% c('SpatialPoints', 'SpatialPointsDataFrame')) {
    if (is.na(raster::projection(r))) {
      warning('Raster will be assumed to have same coordinate reference system as points.', .immediate=TRUE)
      raster::projection(r) <- raster::projection(x)
    } else if (raster::projection(r) != raster::projection(x)) {
      stop('Raster and points do not have the same coordinate reference system.')
    }
  }
  
  # get coordinates
  xy <- xToCoords(x, longLat, sp=FALSE)
  
  # get cell numbers for each point and adjoin with data frame
  cellNum <- raster::cellFromXY(r, xy)
  
  # remember original row names
  index <- 1:nrow(xy)
  
  # define priority
  if (is.null(priority)) priority <- 1:nrow(xy)
  
  # index of points to remove
  removeThese <- integer()
  
  # remove redundant points in each cell
  uniCells <- unique(cellNum)
  for (thisCell in uniCells) {
    
    # if more than one point per cell
    if (sum(cellNum == thisCell) > 1) {
      
      thisRow <- index[thisCell == cellNum]
      thisPriority <- priority[thisCell == cellNum]
      
      thisRow <- thisRow[order(thisPriority)]
      removeThese <- c(removeThese, thisRow[2:length(thisRow)])
      
    }
    
  }
  
  # remove redundant points
  if (length(removeThese) > 0) {
    x <- if (class(x) %in% c('data.frame', 'matrix', 'SpatialPointsDataFrame')) {
      x[-removeThese, ]
    } else {
      x[-removeThese]
    }
  }
  
  x
  
}


load("data/environmental_predictors")
load("data/fishnet")

#create points to sample background data at all non-NA cells
background_points <- spsample(as(environmental_predictors, 'SpatialPolygons'), n=length(environmental_predictors[[1]][!is.na(environmental_predictors[[1]])]), type="regular")
background_points <- spTransform(background_points,CRS("+proj=longlat +datum=WGS84"))

generate_maxent_prediction <- function(species, occurrence_records, native_regions, 
                                       environmental_predictors, parallel, ncores, 
                                       background_points, fishnet){
 
  wgsrpd_regional_list<-native_regions[[3]]
  wgsrpd_country_list<-native_regions[[2]]
  nat_reg<-native_regions[[1]]
  
  nat_reg_mask <- fishnet[!is.na(sp::over(fishnet, sp::geometry(nat_reg))), ] 
  crop_environmental_predictors <- raster::crop(environmental_predictors, extent(nat_reg_mask))
  crop_environmental_predictors <- raster::mask(crop_environmental_predictors, nat_reg_mask)
  
  # select background points in native regions
  background_points_nat_reg <- background_points[!is.na(sp::over(background_points, sp::geometry(nat_reg_mask))), ] 
  
  # raw occurrence data for Model 1
  occ_points_raw<-occurrence_records
  
  # presence cells for Model 2
  ras<-crop_environmental_predictors[[1]]
  presence_cells<-occ_points_raw
  projection(presence_cells)<-projection(ras)
  projection(ras)<-projection(presence_cells)
  presence_cells<-elimCellDups(presence_cells, r=ras)
  
  # thinned presence cells for Model 3
  thinned_presence_cells <- presence_cells 
  thinned_presence_cells <- raster::rasterize(SpatialPoints(coordinates(thinned_presence_cells)), 
                                              crop_environmental_predictors)
  thinned_presence_cells[!is.na(thinned_presence_cells)]<-1
  thinned_presence_cells<-raster::rasterToPoints(thinned_presence_cells)
  sds<-NULL
  
  # spatial thinning at the distance of approx. 2 cells
  sds<-suppressWarnings(spThin::thin.algorithm(data.frame(coordinates(thinned_presence_cells)),55.66*1.99999, rep=1))
  thinned_presence_cells<-SpatialPoints(data.frame(x=sds[[1]][1], y=sds[[1]][2]), 
                                        proj4string=CRS("+proj=longlat +datum=WGS84"))
  thinned_presence_cells<-thinned_presence_cells
  
  
  ########################
  if(length(occ_points_raw)>=5){
    
    
    data<-NULL
    model0<-NULL
    selectEVs<-NULL
    eval<-NULL
    
    if(length(occ_points_raw)>=25){
      eval <- suppressWarnings(ENMevaluate(coordinates(occ_points_raw), crop_environmental_predictors, 
                                           coordinates(background_points_nat_reg), method='block', RMvalues=c(1), 
                                           fc=c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"), parallel=parallel,
                                           numCores = ncores, bin.output=TRUE, algorithm = 'maxent.jar'))
    } else {
      eval <- suppressWarnings(ENMevaluate(coordinates(occ_points_raw), crop_environmental_predictors,
                                           coordinates(background_points_nat_reg), method = 'jackknife', RMvalues = c(1),
                                           fc=c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"), parallel=parallel,
                                           numCores = ncores, bin.output = TRUE, algorithm = 'maxent.jar'))
    }
    
    
    
    if(length(eval)>0){
      
      model0_eval_results<-eval@results
      
      # select the best model based on AICc if this metric is available for more than 50% of generated models
      if((nrow(subset(eval@results, !is.na(eval@results$delta.AICc)))/nrow(eval@results))>0.5){
        model0_results <- eval@results[which(eval@results$delta.AICc==0),]
        model0_results <- model0_results[order(model0_results$avg.test.AUC, decreasing = TRUE),]
        model0_results <- model0_results[1,]
        model0 <- eval@models[[which(eval@results$settings==model0_results$settings)[1]]]
      } else {
        model0_results <- eval@results[which(eval@results$avg.test.AUC==max(eval@results$avg.test.AUC)),][1,]
        model0 <- eval@models[[which(eval@results$settings==model0_results$settings)[1]]]
      }
      
      
      prediction_model0 <- dismo::predict(model0, crop_environmental_predictors, args = c("outputformat=cloglog"), na.rm = TRUE)

      presence_data_raw <- data.frame(raster::extract(crop_environmental_predictors, occ_points_raw, fun=mean, na.rm=TRUE), presence=1)
      presence_data_raw<-na.exclude(presence_data_raw)
      background_data <- data.frame(raster::extract(crop_environmental_predictors, background_points_nat_reg, fun=mean, na.rm=TRUE), presence=0)
      
      evaluate_model0 <- dismo::evaluate(presence_data_raw, background_data, model0)
      threshold_model0 <- dismo::threshold(evaluate_model0)
      
      selectEVs <-  var.importance(model0)
      
      selectEVs<-selectEVs[order(selectEVs$permutation.importance),]
      selectEVs<-subset(selectEVs, selectEVs$permutation.importance>=0)
      selected_environmental_predictors<-NULL
      selected_environmental_predictors<-rem_correlations(crop_environmental_predictors, 
                                                          background_points_nat_reg, selectEVs, 
                                                          correlationthreshold = 0.7)
      
      
      
      eval<-NULL
      if(length(occ_points_raw)>=25){
        eval <- suppressWarnings( ENMevaluate(coordinates(occ_points_raw), selected_environmental_predictors, 
                                              coordinates(background_points_nat_reg), method='block', RMvalues=c(1,2,3,5,10),
                                              fc=c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"), parallel=parallel,
                                              numCores = ncores, bin.output=TRUE,algorithm = 'maxent.jar'))
      } else {
        eval <- suppressWarnings( ENMevaluate(coordinates(occ_points_raw), selected_environmental_predictors, 
                                              coordinates(background_points_nat_reg), method='jackknife', RMvalues=c(1,2,3,5,10),
                                              fc=c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"), parallel=parallel, numCores = ncores, 
                                              bin.output=TRUE,algorithm = 'maxent.jar'))
      }
      
      
      
      model1_results<-NULL
      model1<-NULL
      prediction_model1<-NULL
      evaluate_model1<-NULL
      threshold_model1<-NULL
      rawEVs<-NULL
      model1_eval_results<-NULL
      
      if(length(eval)>0){
        
        model1_eval_results<-eval@results
        
        if((nrow(subset(eval@results, !is.na(eval@results$delta.AICc)))/nrow(eval@results))>0.5){
          model1_results <- eval@results[which(eval@results$delta.AICc==0),]
          model1_results<-model1_results[order(model1_results$avg.test.AUC, decreasing = TRUE),]
          model1_results <- model1_results[1,]
          model1 <- eval@models[[which(eval@results$settings==model1_results$settings)[1]]]
        } else {
          model1_results <- eval@results[which(eval@results$avg.test.AUC==max(eval@results$avg.test.AUC)),][1,]
          model1 <- eval@models[[which(eval@results$settings==model1_results$settings)[1]]]
        }
        
        prediction_model1 <- dismo::predict(model1, crop_environmental_predictors, args = c("outputformat=cloglog"), na.rm = TRUE)
        
        evaluate_model1 <- dismo::evaluate(presence_data_raw, background_data, model1)
        threshold_model1 <- dismo::threshold(evaluate_model1)
        
        rawEVs <-  var.importance(model1)
        
      }
      
      
      model2_results<-NULL
      model2<-NULL
      prediction_model2<-NULL
      evaluate_model2<-NULL
      threshold_model2<-NULL
      cellEVs<-NULL
      model2_eval_results<-NULL
      
      if(length(presence_cells)>=5){
        
        eval<-NULL
        if(length(presence_cells)>=25){
          eval <- suppressWarnings( ENMevaluate(coordinates(presence_cells), selected_environmental_predictors, 
                                                coordinates(background_points_nat_reg), method='block', RMvalues=c(1,2,3,5,10),
                                                fc=c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"), parallel=parallel, numCores = ncores,
                                                bin.output=TRUE, algorithm = 'maxent.jar'))
        } else {
          eval <- suppressWarnings( ENMevaluate(coordinates(presence_cells), selected_environmental_predictors, 
                                                coordinates(background_points_nat_reg), method='jackknife', RMvalues=c(1,2,3,5,10),
                                                fc=c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"), parallel=parallel, numCores = ncores, 
                                                bin.output=TRUE, algorithm = 'maxent.jar'))
        }
        
        
        
        if(length(eval)>0){
          
          model2_eval_results<-eval@results
          
          if((nrow(subset(eval@results, !is.na(eval@results$delta.AICc)))/nrow(eval@results))>0.5){
            model2_results <- eval@results[which(eval@results$delta.AICc==0),]
            model2_results<-model2_results[order(model2_results$avg.test.AUC, decreasing = TRUE),]
            model2_results <- model2_results[1,]
            model2 <- eval@models[[which(eval@results$settings==model2_results$settings)[1]]]
          } else {
            model2_results <- eval@results[which(eval@results$avg.test.AUC==max(eval@results$avg.test.AUC)),][1,]
            model2 <- eval@models[[which(eval@results$settings==model2_results$settings)[1]]]
          }
          
          presence_data_cells <- data.frame(raster::extract(crop_environmental_predictors, presence_cells, fun=mean, na.rm=TRUE), presence=1)
          presence_data_cells <- na.exclude(presence_data_cells)
          
          prediction_model2 <- dismo::predict(model2, crop_environmental_predictors, args = c("outputformat=cloglog"), na.rm = TRUE)
          
          evaluate_model2 <- dismo::evaluate(presence_data_cells, background_data, model2)
          threshold_model2 <- dismo::threshold(evaluate_model2)
          cellEVs <-  var.importance(model2)
          
        }
        
      }
      
      model3_results<-NULL
      model3<-NULL
      prediction_model3<-NULL
      evaluate_model3<-NULL
      threshold_model3<-NULL
      thinEVs<-NULL
      model3_eval_results<-NULL
      
      if(length(thinned_presence_cells)>=5){
        
        eval<-NULL
        if(length(thinned_presence_cells)>=25){
          eval <- suppressWarnings( ENMevaluate(coordinates(thinned_presence_cells), selected_environmental_predictors,
                                                coordinates(background_points_nat_reg), method='block', RMvalues=c(1,2,3,5,10),
                                                fc=c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"), parallel=parallel, numCores = ncores,
                                                bin.output=TRUE,algorithm = 'maxent.jar'))
        } else {
          eval <- suppressWarnings( ENMevaluate(coordinates(thinned_presence_cells), selected_environmental_predictors,
                                                coordinates(background_points_nat_reg), method='jackknife', RMvalues=c(1,2,3,5,10),
                                                fc=c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"), parallel=parallel, numCores = ncores,
                                                bin.output=TRUE,algorithm = 'maxent.jar'))
        }
        
        
    
        
        if(length(eval)>0){
          
          model3_eval_results<-eval@results
          
          if((nrow(subset(eval@results, !is.na(eval@results$delta.AICc)))/nrow(eval@results))>0.5){
            model3_results <- eval@results[which(eval@results$delta.AICc==0),]
            model3_results<-model3_results[order(model3_results$avg.test.AUC, decreasing = TRUE),]
            model3_results <- model3_results[1,]
            model3 <- eval@models[[which(eval@results$settings==model3_results$settings)[1]]]
          } else {
            model3_results <- eval@results[which(eval@results$avg.test.AUC==max(eval@results$avg.test.AUC)),][1,]
            model3 <- eval@models[[which(eval@results$settings==model3_results$settings)[1]]]
          }
          
          
          presence_data_thinned <- data.frame(raster::extract(crop_environmental_predictors, thinned_presence_cells, fun=mean, na.rm=TRUE), presence=1)
          presence_data_thinned <- na.exclude(presence_data_thinned)
          
          training_data_raw <- rbind(presence_data_raw, background_data)
          training_data_raw <- na.exclude(training_data_raw)
          
          prediction_model3 <- dismo::predict(model3, crop_environmental_predictors, args = c("outputformat=cloglog"), na.rm = TRUE)
          evaluate_model3 <- dismo::evaluate(presence_data_thinned, background_data, model3)
          threshold_model3 <- dismo::threshold(evaluate_model3)
          thinEVs <-  var.importance(model3)
          
          
          
        }} 
      
      eval_results<-list(model0_eval_results,model1_eval_results,model2_eval_results,model3_eval_results)
      names(eval_results)<-c("Model0","Model1","Model2","Model3")
      model_results<-list(model0_results,model1_results,model2_results,model3_results)
      names(model_results)<-c("Model0","Model1","Model2","Model3")
      models<-list(model0,model1,model2,model3)
      names(models)<-c("Model0","Model1","Model2","Model3")
      EVs<-list(selectEVs,rawEVs,cellEVs,thinEVs)
      names(EVs)<-c("Model0","Model1","Model2","Model3")
      predictions<-list(prediction_model0,prediction_model1,prediction_model2,prediction_model3)
      names(predictions)<-c("Model0","Model1","Model2","Model3")
      evaluate<-list(evaluate_model0,evaluate_model1,evaluate_model2,evaluate_model3)
      names(evaluate)<-c("Model0","Model1","Model2","Model3")
      thresholds<-list(threshold_model0,threshold_model1,threshold_model2,threshold_model3)
      names(thresholds)<-c("Model0","Model1","Model2","Model3")
      
      data_list<-list(eval_results, model_results, models, EVs, predictions, evaluate, thresholds)
      names(data_list)<-c("eval_results", "model_results", "models", "EVs", "predictions", "evaluate", "thresholds")
      return(data_list)
      
    }}
}


## usage
maxent_prediction<-generate_maxent_prediction("Amomum pterocarpum", occurrence_records, 
                                              native_regions, environmental_predictors, 
                                              parallel = TRUE, ncores = 4, background_points, fishnet)
summary(maxent_prediction)

plot(maxent_prediction$predictions$Model1)

str(maxent_prediction)
