library(stringr)
library(taxize)
library(rvest)

#species name
#type one of "Native", "Introduced", "Uncertain"
pow_wgsrpd <- function(species, type){
  ppow <- NULL
  fnames <- NULL
  
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
    fype<-html_nodes(x=url, css=selector.type) %>%
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
    } else {
      message(paste("no regions reported as", type))
    }
  }
}

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
    
    tdwg<-read.csv("https://raw.githubusercontent.com/jannebor/supporting_files/master/tdwg.csv")
    tdwg$ISO<-as.character(tdwg$ISO)
    
    tdwg$Continent<-gsub(" ","", tdwg$Continent)
    tdwg$Continent<-substr(tdwg$Continent, 1, 18)
    
    tdwg$Sub_cont<-gsub(" ","", tdwg$Sub_cont)
    tdwg$Sub_cont<-substr(tdwg$Sub_cont, 1, 18)
    
    tdwg$Region<-gsub(" ","", tdwg$Region)
    tdwg$Region<-substr(tdwg$Region, 1, 18)
    
    tdwg$Country<-gsub(" ","", tdwg$Country)
    tdwg$Country<-substr(tdwg$Country, 1, 18)
    
    
    sub_df<-subset(tdwg, tdwg$Country==wgsrpd_regions[1])
    if(nrow(sub_df)==0){
      sub_df<-subset(tdwg, tdwg$Region==wgsrpd_regions[1])
      if(nrow(sub_df)==0){
        sub_df<-subset(tdwg, tdwg$Sub_cont==wgsrpd_regions[1])
        if(nrow(sub_df)==0){
          sub_df<-subset(tdwg, tdwg$Continent==wgsrpd_regions[1])
        }
      }
    }
    
    
    
    
    if(length(wgsrpd_regions)>1){
      for (l in 2:length(wgsrpd_regions)){
        sub_add<-NA
        sub_add<-subset(tdwg, tdwg$Country==wgsrpd_regions[l])
        
        if(nrow(sub_add)<1){
          
          sub_add<-subset(tdwg, tdwg$Region==wgsrpd_regions[l])
          
          if(nrow(sub_add)<1){
            sub_add<-subset(tdwg, tdwg$Sub_cont==wgsrpd_regions[l])
            
            if(nrow(sub_add)<1){
              sub_add<-subset(tdwg, tdwg$Continent==wgsrpd_regions[l])
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
      country_list<-unique(sub_df$Cont_code)
    }
    if(tolower(format)=="continent"){
      country_list<-unique(sub_df$Continent)
    }
    if(tolower(format)=="sub continent code"){
      country_list<-unique(sub_df$Sub_cont_code)
    }
    if(tolower(format)=="sub continent"){
      country_list<-unique(sub_df$Sub_cont)
    }
    if(tolower(format)=="region"){
      country_list<-unique(sub_df$Region)
    }
    if(format=="isocode5"){
      country_list<-unique(sub_df$X5Letter)
    }
    if(tolower(format)=="country"){
      country_list<-unique(sub_df$Country)
    }
    if(format=="isocode2"){
      country_list<-unique(sub_df$ISO)
    }
    
    return(country_list)
    
  }
}


# examples
#type needs to be one one of "Native", "Introduced", "Uncertain"
wgsrpd <- pow_wgsrpd("Ranunculus glacialis", type="Introduced")

wgsrpd <- pow_wgsrpd("Ranunculus glacialis", type="Native")

# taking the output of pow_wgsrpd as input as well as a country format:
# one of: "Continent", "Continent code", "Sub continent", "Sub continent code",
# "Region", "isocode5", "Country", "isocode2"
countrylist <- wgsrpd_conversion(wgsrpd, format="isocode5")
