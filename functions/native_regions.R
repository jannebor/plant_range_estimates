### R version 4.0.3
## code for retrieving occurrence records from GBIF

library(taxize)
library(rvest)
library(stringr)
library(rgeos)

load("data/tdwg_level4")
load("data/tdwg_conversion_scheme")

# web scraper function for data retrievel from plants of the world online
pow_wgsrpd <- function(species, type){
  ppow <- NULL
  
  while(length(ppow)<1){
    t0<-proc.time()
    tryCatch({
      ppow<-get_pow(species, accepted = TRUE, rows = 1, messages=FALSE)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
    Sys.sleep(2)
    t<-t0-proc.time()
    
    if(abs(t[3])>30){
      break
    }
    
  }
  
  if(!is.na(ppow[1])){
    
    ppow_data<-pow_lookup(ppow[1])
    suppressWarnings(url<-read_html(paste("http://plantsoftheworldonline.org/taxon/",ppow[1],sep="")))
    
    selector.type<-"#distribution-listing > h3:nth-child(1)"
    fype<-rvest::html_nodes(x=url, css=selector.type) %>%
      html_text(trim=TRUE)
    if((length(grep(tolower(type),tolower(fype)))>0)){
      selector.name<-"#distribution-listing > p:nth-child(2)"
      
      fnames<-html_nodes(x=url, css=selector.name) %>%
        html_text(trim=TRUE)
    } else {
      selector.type<-"#distribution-listing > h3:nth-child(3)"
      fype<-html_nodes(x=url, css=selector.type) %>%
        html_text(trim=TRUE)
      if((length(grep(tolower(type),tolower(fype)))>0)){
        selector.name<-"#distribution-listing > p:nth-child(4)"
        
        fnames<-html_nodes(x=url, css=selector.name) %>%
          html_text(trim=TRUE)
      } else { selector.type<-"#distribution-listing > h3:nth-child(5)"
      fype<-html_nodes(x=url, css=selector.type) %>%
        html_text(trim=TRUE)
      if((length(grep(tolower(type),tolower(fype)))>0)){
        selector.name<-"#distribution-listing > p:nth-child(6)"
        
        fnames<-html_nodes(x=url, css=selector.name) %>%
          html_text(trim=TRUE)
      }
      }
    }
    
    if(length(fnames)>0){
      
      #distribution-listing > p:nth-child(2)
      fnames<-gsub("\n","", fnames)
      fnames<-gsub("\r","", fnames)
      fnames<-gsub(" ","", fnames)
      fnames<-str_split(fnames, pattern=",")
      fnames<-unlist(fnames)
      
      if(length(fnames)>0){
        for (t in 1:length(fnames)) {
          
          if (fnames[t]=="Panam???"){
            fnames[t]<-"Panama"
            
          }
          if (fnames[t]=="NorthwestTerritorie"){
            fnames[t]<-"NorthwestTerritori"
            
          }
          
          
          
        }
        
        return(fnames)
        
      }
    }
  }
}


###
# takes the output of pow_wgsrpd as input
# format defines to destination format during conversion:
# one of: "Continent", "Continent code", "Sub continent", "Sub continent code",
# "Region", "isocode5", "Country", "isocode2"
wgsrpd_conversion <-function(wgsrpd_regions, format){
  
  if(length(wgsrpd_regions)>0){
    for (t in 1:length(wgsrpd_regions)) {
      
      if (wgsrpd_regions[t]=="Panam???"){
        wgsrpd_regions[t]<-"Panama"
      }
    }
    
    sub_df<-subset(tdwg_conversion_scheme, tdwg_conversion_scheme$Country==wgsrpd_regions[1])
    if(nrow(sub_df)==0){
      sub_df<-subset(tdwg_conversion_scheme, tdwg_conversion_scheme$Region==wgsrpd_regions[1])
      if(nrow(sub_df)==0){
        sub_df<-subset(tdwg_conversion_scheme, tdwg_conversion_scheme$Sub_cont==wgsrpd_regions[1])
        if(nrow(sub_df)==0){
          sub_df<-subset(tdwg_conversion_scheme, tdwg_conversion_scheme$Continent==wgsrpd_regions[1])
        }
      }
    }
    
    
    
    
    if(length(wgsrpd_regions)>1){
      for (l in 2:length(wgsrpd_regions)){
        sub_add<-NA
        sub_add<-subset(tdwg_conversion_scheme, tdwg_conversion_scheme$Country==wgsrpd_regions[l])
        
        if(nrow(sub_add)<1){
          
          sub_add<-subset(tdwg_conversion_scheme, tdwg_conversion_scheme$Region==wgsrpd_regions[l])
          
          if(nrow(sub_add)<1){
            sub_add<-subset(tdwg_conversion_scheme, tdwg_conversion_scheme$Sub_cont==wgsrpd_regions[l])
            
            if(nrow(sub_add)<1){
              sub_add<-subset(tdwg_conversion_scheme, tdwg_conversion_scheme$Continent==wgsrpd_regions[l])
            }
          }
        }
        
        sub_df<-rbind(sub_df, sub_add)
        
      }
      
    }
    
    if(nrow(sub_df)>0){
      
      for (k in 1:nrow(sub_df)){
        
        if (sub_df$ISO[k]=="Na"){
          sub_df$ISO[k]<-"NA"
          
        }
      }
    }
    
    if(tolower(format)=="continent code"){
      wgsrpd_country_list<-unique(sub_df$Cont_code)
    }
    if(tolower(format)=="continent"){
      wgsrpd_country_list<-unique(sub_df$Continent)
    }
    if(tolower(format)=="sub continent code"){
      wgsrpd_country_list<-unique(sub_df$Sub_cont_code)
    }
    if(tolower(format)=="sub continent"){
      wgsrpd_country_list<-unique(sub_df$Sub_cont)
    }
    if(tolower(format)=="region"){
      wgsrpd_country_list<-unique(sub_df$Region)
    }
    if(format=="isocode5"){
      wgsrpd_country_list<-unique(sub_df$X5Letter)
    }
    if(tolower(format)=="country"){
      wgsrpd_country_list<-unique(sub_df$Country)
    }
    if(format=="isocode2"){
      wgsrpd_country_list<-unique(sub_df$ISO)
    }
    
    
    return(wgsrpd_country_list)
    
  }
}

get_native_regions<-function(species){
   
  fnames<-pow_wgsrpd(species, type="Native")
  
  if (length(fnames)>0){
    wgsrpd_country_list<-wgsrpd_conversion(fnames, format = "isocode2")
    wgsrpd_regional_list<-wgsrpd_conversion(fnames, format = "isocode5")
    
    nat_reg<-NULL
    for(y in 1:length(wgsrpd_regional_list)){
      if(length(nat_reg)==0){
        nat_reg<-subset(tdwg_level4, tdwg_level4$Level4_cod==wgsrpd_regional_list[y])
        nat_reg <- aggregate(nat_reg)
      } else {
        t_add<-subset(tdwg_level4, tdwg_level4$Level4_cod==wgsrpd_regional_list[y])
        t_add <- aggregate(t_add)
        nat_reg<-gUnion(t_add, nat_reg)
      }
    }
    
    return(list(nat_reg, wgsrpd_country_list, wgsrpd_regional_list))
    
  }
  
  
  
  
  
}

# usage:
native_regions<-get_native_regions("Amomum pterocarpum")
plot(native_regions[[1]])
print(native_regions[[2]])
print(native_regions[[3]])


