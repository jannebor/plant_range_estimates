### R version 4.0.3
## code for evaluating generated Maxent predictions
# requires output of native_regions.R, occurrence_records.R, and maxent_predictions.R

library(PRROC)
library(ROCR)
library(raster)
library(sp)
library(spThin)
library(stringr)

evaluate_maxent_prediction <- function(maxent_prediction, native_regions, 
                                     occurrence_records, environmental_predictors,
                                     fishnet, reference_poly=NULL){

wgsrpd_regional_list<-native_regions[[3]]
wgsrpd_country_list<-native_regions[[2]]
nat_reg<-native_regions[[1]]

nat_reg_mask <- fishnet[!is.na(sp::over(fishnet, sp::geometry(nat_reg))), ] 
crop_environmental_predictors<-crop(environmental_predictors, nat_reg_mask)
crop_environmental_predictors<- mask(crop_environmental_predictors, nat_reg_mask)


# raw occurrence data for Model 1
occ_points_raw<-occurrence_records

# presence cells for Model 2
ras<-crop_environmental_predictors[[1]]
presence_cells<-occ_points_Model1
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
sds<-suppressWarnings(thin.algorithm(data.frame(coordinates(thinned_presence_cells)),55.66*1.99999, rep=1))
thinned_presence_cells<-SpatialPoints(data.frame(x=sds[[1]][1], y=sds[[1]][2]), 
                                      proj4string=CRS("+proj=longlat +datum=WGS84"))
thinned_presence_cells<-thinned_presence_cells

occ_ras<-rasterize(occ_points_raw, maxent_prediction$predictions$Model0, fun="count")      
occ_ras<-occ_ras$ID
mor_gbif_points<-raster::Moran(occ_ras)
gbif_points_Model1<-length(occ_points_raw)
gbif_points_Model2<-length(presence_cells)
gbif_points_Model3<-length(thinned_presence_cells)


if(exists("reference_poly")){
  if(length(reference_poly)>0){
    reference_poly<-subset(reference_poly, reference_poly$ORIGIN==1)
    reference_poly<-subset(reference_poly, reference_poly$PRESENCE<3)
  }
} else {
  reference_poly <- NULL
}



datatypes<-c("Model 0", "Model 1", "Model 2", "Model 3")
# retrieving environmental variables
for (d in 1:length(datatypes)) {
  if(length(maxent_prediction$models[[d]])>0){
    EVs<-var.importance(maxent_prediction$models[[d]])
    EVs<-EVs[order(EVs$permutation.importance),]
    EVs<-subset(EVs, EVs$permutation.importance>0)
    E<-paste(EVs$variable) 
    assign(paste("variables_n_",str_replace(datatypes[d]," ","_"),sep=""), length(E))
    E<-str_remove(E, "_LC_300")
    E<-str_replace(E, "_bio10_","_BIO")
    E<-str_replace(E,"_BIO0","_BIO")
    E<-str_replace(E, "ESACCI","ESA_CCI")
    assign(paste("variables_",str_replace(datatypes[d]," ","_"),sep=""),paste(E, collapse = "; "))
    
    pred<-maxent_prediction$predictions[[d]]
    
    if(length(pred)>0){
      ### performance based on presence cells
      if(length(presence_cells)>0){
        
        reference<-presence_cells
        reference<-rasterize(reference, pred, fun="count")
        reference<-reference$ID
        reference[!is.na(reference)]<-1
        reference[is.na(reference)]<-0
        reference<-mask(reference, pred)
        reference<-raster::as.matrix(reference)
        reference<-as.vector(reference)
        reference<-na.exclude(reference)
        
        
        if(length(unique(reference))>1){
          l2<-as.numeric(reference)
          
          pred_num<-raster::as.matrix(pred)
          pred_num<-as.vector(pred_num)
          pred_num<-na.exclude(pred_num)
          pred_num<-as.numeric(pred_num)
          
          dat<-data.frame(score1=pred_num, label=l2)
          
          pred1 <- suppressWarnings(prediction(dat$score1,dat$label))
          auc_ROCR <- suppressWarnings(performance(pred1, measure = "auc"))
          f1_score <- suppressWarnings(performance(pred1, measure = "f"))
          pr <- pr.curve(scores.class0=dat[dat$label=="1",]$score1,
                         scores.class1=dat[dat$label=="0",]$score1,
                         curve=F)
          assign(paste("Cellauc_",str_replace(datatypes[d]," ","_"),sep=""),auc_ROCR@y.values[[1]])
          assign(paste("Cellaucpr_",str_replace(datatypes[d]," ","_"),sep=""),pr$auc.integral)
          assign(paste("Cellf1_score_",str_replace(datatypes[d]," ","_"),sep=""),max(na.exclude(as.data.frame(f1_score@y.values))))
          
          
        }
      } else {
        assign(paste("Cellauc_",str_replace(datatypes[d]," ","_"),sep=""),NA)
        assign(paste("Cellaucpr_",str_replace(datatypes[d]," ","_"),sep=""),NA)
        assign(paste("Cellf1_score_",str_replace(datatypes[d]," ","_"),sep=""),NA)
      }
      
      
      ### performance based on reference_poly input
      if(exists("reference_poly")){
        if(length(reference_poly)>0){
          
          reference<-reference_poly
          reference<-rasterize(reference, pred)
          reference[!is.na(reference)]<-1
          reference[is.na(reference)]<-0
          reference<-mask(reference, pred)
          reference<-raster::as.matrix(reference)
          reference<-as.vector(reference)
          reference<-na.exclude(reference)
          
          
          if(length(unique(reference))>1){
            l2<-as.numeric(reference)
            
            pred_num<-raster::as.matrix(pred)
            pred_num<-as.vector(pred_num)
            pred_num<-na.exclude(pred_num)
            pred_num<-as.numeric(pred_num)
            
            dat<-data.frame(score1=pred_num, label=l2)
            
            pred1 <- suppressWarnings(prediction(dat$score1,dat$label))
            auc_ROCR <- suppressWarnings(performance(pred1, measure = "auc"))
            f1_score <- suppressWarnings(performance(pred1, measure = "f"))
            pr <- pr.curve(scores.class0=dat[dat$label=="1",]$score1,
                           scores.class1=dat[dat$label=="0",]$score1,
                           curve=F)
            assign(paste("REFauc_",str_replace(datatypes[d]," ","_"),sep=""),auc_ROCR@y.values[[1]])
            assign(paste("REFaucpr_",str_replace(datatypes[d]," ","_"),sep=""),pr$auc.integral)
            assign(paste("REFf1_score_",str_replace(datatypes[d]," ","_"),sep=""),max(na.exclude(as.data.frame(f1_score@y.values))))
            
            
          }
        } 
      } else {
          assign(paste("REFauc_",str_replace(datatypes[d]," ","_"),sep=""),NA)
          assign(paste("REFaucpr_",str_replace(datatypes[d]," ","_"),sep=""),NA)
          assign(paste("REFf1_score_",str_replace(datatypes[d]," ","_"),sep=""),NA)
        }
      
    }
    
    
  } else {
    assign(paste("variables_n_",str_replace(datatypes[d]," ","_"),sep=""), NA)
    assign(paste("variables_",str_replace(datatypes[d]," ","_"),sep=""), NA)
    assign(paste("Cellauc_",str_replace(datatypes[d]," ","_"),sep=""),NA)
    assign(paste("Cellaucpr_",str_replace(datatypes[d]," ","_"),sep=""),NA)
    assign(paste("Cellf1_score_",str_replace(datatypes[d]," ","_"),sep=""),NA)
    assign(paste("REFauc_",str_replace(datatypes[d]," ","_"),sep=""),NA)
    assign(paste("REFaucpr_",str_replace(datatypes[d]," ","_"),sep=""),NA)
    assign(paste("REFf1_score_",str_replace(datatypes[d]," ","_"),sep=""),NA)
  }
}



data<-data.frame(species=(rep(maxent_prediction$species,4)))

type<-data.frame(type=c("Model 0","Model 1","Model 2", "Model 3"))

n_points<-rbind(gbif_points_Model1, gbif_points_Model1, gbif_points_Model2, gbif_points_Model3)

gbif_raw<-rbind(gbif_points_Model1, gbif_points_Model1, gbif_points_Model1, gbif_points_Model1)

moransI<-rbind(mor_gbif_points,mor_gbif_points,mor_gbif_points,mor_gbif_points)

### OOB performance
oob_performance<-rbind(maxent_prediction$model_results$Model0[1:16], maxent_prediction$model_results$Model1[1:16])
for(d in 3:length(datatypes)){
  if(!is.null(maxent_prediction$model_results[[d]])){
    oob_performance<-rbind(oob_performance,maxent_prediction$model_results[[d]][1:16])
  } else {
    oob_performance[d,]<-NA
      }
}


#maxent_prediction$model_results$Model2[1:16], maxent_prediction$model_results$Model3[1:16])

variables<-data.frame(n_variables=rbind(variables_n_Model_0, variables_n_Model_1, variables_n_Model_2, variables_n_Model_3),
                      variables=rbind(variables_Model_0, variables_Model_1, variables_Model_2, variables_Model_3))

Cell_metrics<-data.frame(AUC_Cell=c(Cellauc_Model_0, Cellauc_Model_1, Cellauc_Model_2, Cellauc_Model_3),
                    AUCPR_Cell=c(Cellaucpr_Model_0, Cellaucpr_Model_1, Cellaucpr_Model_2, Cellaucpr_Model_3),
                    F_score_Cell=c(Cellf1_score_Model_0, Cellf1_score_Model_1, Cellf1_score_Model_2, Cellf1_score_Model_3))

REF_metrics<-data.frame(AUC_REF=c(REFauc_Model_0, REFauc_Model_1, REFauc_Model_2, REFauc_Model_3),
                         AUCPR_REF=c(REFaucpr_Model_0, REFaucpr_Model_1, REFaucpr_Model_2, REFaucpr_Model_3),
                         F_score_REF=c(REFf1_score_Model_0, REFf1_score_Model_1, REFf1_score_Model_2, REFf1_score_Model_3))


### Thresholds
thresholds<-rbind(maxent_prediction$thresholds$Model0,maxent_prediction$thresholds$Model1)
for(d in 3:length(datatypes)){
  if(!is.null(maxent_prediction$model_results[[d]])){
    thresholds<-rbind(thresholds,maxent_prediction$thresholds[[d]])
  } else {
    thresholds[d,]<-NA
  }
}



colnames(thresholds)<-c("CutOff_kappa","CutOff_spec_sens","CutOff_no_omission",
                        "CutOff_prevalence","CutOff_equal_sens_spec","CutOff_sensitivity")

df<-cbind(data, type, gbif_raw, n_points, moransI, oob_performance, variables, Cell_metrics, REF_metrics, thresholds)
row.names(df)<-datatypes
df<-df[2:4,]
return(df)


}


# usage
eval_df <- evaluate_maxent_prediction(maxent_prediction, native_regions, occurrence_records, environmental_predictors, fishnet)

eval_df
