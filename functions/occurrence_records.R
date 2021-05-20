### R version 4.0.3
## code for retrieving occurrence records from GBIF
# requires output of native_regions.R


library(rgbif)


get_occurrence_records<-function(species, native_regions){
  
  wgsrpd_regional_list<-native_regions[[3]]
  wgsrpd_country_list<-native_regions[[2]]
  nat_reg<-native_regions[[1]]
  
  tryCatch({
    ppow<-get_pow(species, accepted = TRUE, rows = 1, messages=FALSE)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  tryCatch({
    ppow_data<-pow_lookup(ppow[1])
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  key<-name_backbone(name=paste(ppow_data$meta$name))$speciesKey
  
  if (length(key)!=0){
    if (length(wgsrpd_country_list)>0){
      
      #create empty DataFrame
      occ_points_raw <- na.exclude(data.frame(x=NA,y=NA,country=NA))
      
      for (c in 1:length(wgsrpd_country_list)){
        occ<-NULL
        # retrieve data from native countries, without geospatial issues and only records from 2000-2020
        tryCatch({
          occ<-occ_search(taxonKey=key, country = paste(wgsrpd_country_list[c]), year="2000,2020", fields="all", hasCoordinate = T, hasGeospatialIssue = F,limit=100000)
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
        
        if (!is.null(occ$data)){
          occ_add <- data.frame(x=occ$data$decimalLongitude,y=occ$data$decimalLatitude,country=wgsrpd_country_list[c])
          occ_points_raw<-rbind(occ_points_raw,occ_add)
          
        }
        
      }
        
      # remove duplicates 
      occ_points_raw <- occ_points_raw[!duplicated(occ_points_raw),]
      
      occ_points_raw <- SpatialPointsDataFrame(data.frame(x=occ_points_raw$x, y=occ_points_raw$y), occ_points_raw,proj4string=CRS("+proj=longlat +datum=WGS84"))
      # remove non-native observations
      occ_points_raw <- occ_points_raw[!is.na(sp::over(occ_points_raw, sp::geometry(nat_reg))), ] 
      
      # if only few records available, remove temporal filter
      if(length(occ_points_raw)<25){
        #create empty DataFrame
        occ_points_raw <- na.exclude(data.frame(x=NA,y=NA,country=NA))
        
        for (c in 1:length(wgsrpd_country_list)){
          occ<-NULL
          tryCatch({
            occ<-occ_search(taxonKey=key, country = paste(wgsrpd_country_list[c]), fields="all", hasCoordinate = T, hasGeospatialIssue = F,limit=100000)
          }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
          
          if (!is.null(occ$data)){
            occ_add <- data.frame(x=occ$data$decimalLongitude,y=occ$data$decimalLatitude,country=wgsrpd_country_list[c])
            occ_points_raw<-rbind(occ_points_raw,occ_add)
            
          }
          
        }
        
        # remove duplicates 
        occ_points_raw <- occ_points_raw[!duplicated(occ_points_raw),]
        
        occ_points_raw <- SpatialPointsDataFrame(data.frame(x=occ_points_raw$x, y=occ_points_raw$y), occ_points_raw,proj4string=CRS("+proj=longlat +datum=WGS84"))
        # remove non-native observations
        occ_points_raw <- occ_points_raw[!is.na(sp::over(occ_points_raw, sp::geometry(nat_reg))), ] 
        
      }
      
      
      return(occ_points_raw)
    }
  
}











}


# usage
occurrence_records<-get_occurrence_records("Amomum pterocarpum", native_regions)
plot(occurrence_records)
