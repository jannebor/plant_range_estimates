## Native range estimates for red-listed vascular plant species

Distribution data for individual species (e.g. range maps) is central for understanding both global biodiversity patterns and associated anthropogenic impacts. Our aim is to provide a comprehensive set of spatial data for a large number of, so far unmapped, vascular plant species, supporting the development of future assessments of global biodiversity impacts and distributions.


![Range map examples - Maxent predictions transformed into potential binary range maps at different levels of confidence using different cut-off thresholds](demo/slideshow.gif)


### Dataset

The dataset consists of native regions for 47,675 species retrieved from [Plants of the World online](http://plantsoftheworldonline.org/), density of available native occurrence records for 30,906 species retrieved from [GBIF](https://www.gbif.org/), and standardized, large-scale [Maxent](https://biodiversityinformatics.amnh.org/open_source/maxent/) range estimates for 27,208 species, highlighting environmentally suitable areas within species' native regions. The data can be explored using our  [data explorer](https://plant-ranges.indecol.no/).

The full dataset can be downloaded in 30 arc minutes spatial resolution (approx. 56 km at the equator).

The folder basic contains two netCDF files (default output and raw output) that assemble the best performing Maxent prediction (<i>varname</i>: Maxent prediction) for each species selected based on the highest harmonic mean between AUC and AUC<sub>PR</sub>, along with number of occurrence records per cell (<i>varname</i>: Presence cells) and rasterized native WGSRPD-regions (<i>varname</i>: Native region). 

The folder advanced contains two netCDF files (default output and raw output). These files contain number of occurrence records per cell (<i>varname</i>: Presence cells), rasterized native WGSRPD-regions (<i>varname</i>: Native region), and one Maxent prediction for each occurrence data type (<i>varname</i>: Model 1, Model 2 or Model 3). 
NOTE that <i>varname</i> "Maxent prediction" is not applicable

Each band in the netCDF files assembles the mentioned variables for one species. The corresponding bands can be looked up in the metadata (i.e. speciesID).

The metadata can be used to filter models based on species, performance or desired data types, to look up the relevant study extent for masking individual predictions, or to select different cut-off thresholds for generating binary range maps. For information about the given thresholds visit <https://www.rdocumentation.org/packages/dismo/versions/1.1-4/topics/threshold>.

------------------------------------------------------------------------

### Examples

Species selection based on metadata:

```{r setup, include = FALSE}
library(raster)
library(ncdf4)

path <- "./doi_10.5061_dryad.qbzkh18h9__v3"

dataset <- "basic"

format <- "default"

# metadata of the suggested (basic) dataset (only the best prediction for each species)
metadata<-read.csv(paste(path,"/",dataset,"/","metadata_",format,".csv",sep=""))

# for a given species, e.g.:
species_name <- "Amomum pterocarpum"

# speciesID
speciesID <- metadata$speciesID[which(metadata$scientificname==species_name)][1]

# type of provided data
metadata$data[which(metadata$scientificname==species_name)][1]

# and, if a Maxent prediction has been generated, the best performing model
# can be looked up in the metadata
best_model <- metadata$model[which(metadata$scientificname==species_name)][1]


# each band of the netCDF file assembles data for one species
# variables represent the different data types within each band

ras <- raster(paste(path,"/",dataset,"/","range_data_",format,".nc",sep=""), varname = "Native region", band = speciesID)

# varname needs to be one of "Native region", "Presence cells", "Maxent prediction", or
# in advanced datasets a specific Model needs to be specified instead of "Maxent prediction": 
# "Model 1", "Model 2", "Model 3"

plot(ras)
```
All bands in the netCDF file have the same spatial extent, to crop to species-specific extents:

```{r setup, include = FALSE}
ras<-crop(ras, extent(metadata$extent.xmin[which(metadata$speciesID==speciesID)][1],
                      metadata$extent.xmax[which(metadata$speciesID==speciesID)][1],
                      metadata$extent.ymin[which(metadata$speciesID==speciesID)][1],
                      metadata$extent.ymax[which(metadata$speciesID==speciesID)][1]))
plot(ras)
```



the above described steps summarized in a function:

```{r setup, include = FALSE}
read_data<-function(species_name, variables, path, dataset, format){
  
  metadata<-read.csv(paste(path,"/",dataset,"/","metadata_",format,".csv",sep=""))
  
  speciesID <- metadata$speciesID[which(metadata$scientificname==species_name)][1]
  
  ras<-stack()
  
  for(v in variables){
    if(dataset=="advanced"){
      if(v == "Maxent prediction"){
        metadata_suggested<-read.csv(paste(path,"/basic/","metadata_",format,".csv",sep=""))
        v <- metadata_suggested$model[which(metadata_suggested$scientificname==species_name)][1]
      }
    }
    
    ras<-stack(ras,raster(paste(path,"/",dataset,"/","range_data_",format,".nc",sep=""), varname = v, band = speciesID))
  }
  
  names(ras)<-variables
  
  return(crop(ras, extent(metadata$extent.xmin[which(metadata$speciesID==speciesID)][1],
                          metadata$extent.xmax[which(metadata$speciesID==speciesID)][1],
                          metadata$extent.ymin[which(metadata$speciesID==speciesID)][1],
                          metadata$extent.ymax[which(metadata$speciesID==speciesID)][1]))
  )
}
```


show Native regions, Occurrence data and best performing Maxent prediction (cloglog transformed) for <i>Amomum pterocarpum</i>

```{r setup, include = FALSE}
variables <- c("Native region","Presence cells","Maxent prediction")

ras<-read_data(species_name = "Amomum pterocarpum", variables = variables, path = path, dataset = "basic", format = "default")

plot(ras)
```




### Metadata

#### General species information

speciesID, scientificname, redlistcategory, rank, family, order, class, subdivision, division, and data (type of provided data one of Native region, Occurrence records or Maxent prediction)


#### Original extent of native regions, e.g. for cropping bands of the netCDF file:

extent.xmin: Minimum x-value of the original study extent

extent.xmax: Maximum x-value of the original study extent

extent.ymin: Minimum y-value of the original study extent

extent.ymax: Maximum y-value of the original study extent


#### Occurrence data

moransI: Moran's Index based on the raw occurrence data

raw.points: Number of raw occurrence points for the given species


#### Maxent settings

n.points: Number of occurrence points used for training the given model

datatype: Data used for training the model, i.e. raw occurrence data, presence cells or thinned presence cells

model: Model 1, Model 2 or Model 3

settings: Best performing settings for the given model. Feature transformations out of linear, quadratic, hinge, product and threshold. Regularization multipliers out of 1,2,3,5 and 10.

features: Selected feature transformations out of linear, quadratic, hinge, product and threshold

rm: Selected regularization multiplier

AICc: Corrected Akaike Information Criterion

delta.AICc: Difference to best AICc in the same run

parameters: Number of parameters considered in the given model

n.variables: Number of variables considered in the given model

variables: Variables considered in the given model


#### Maxent performance

train.AUC: obtained AUC during model calibration

avg.test.AUC: Average test AUC during model calibration

var.test.AUC: Variance in testing AUC during model calibration

avg.diff.AUC: Average difference in testing AUC during model calibration

var.diff.AUC: Variance in testing AUC during model calibration

avg.test.orMTP: Average proportion of testing localities with Maxent output values lower than the value associated with the training locality with the lowest value

var.test.orMTP: Variance in proportion of testing localities with Maxent output values lower than the value associated with the training locality with the lowest value

avg.test.or10pct: Average proportion of testing localities with Maxent output values lower than the value associated with the value that excludes the 10 percent of training localities with the lowest predicted suitability

var.test.or10pct: Variance in proportion of testing localities with Maxent output values lower than the value associated with the value that excludes the 10 percent of training localities with the lowest predicted suitability

cell.AUC: AUC values based on comparison to presence cells only

cell.AUCPR: AUC<sub>PR</sub> values based on comparison to presence cells only

cell.maxF1: Maximum F<sub>1</sub>-score based on comparison to presence cells only

ref.AUC: AUC values based on comparison to expert-based reference range map

ref.AUCPR: AUC<sub>PR</sub> values based on comparison to expert-based reference range map

ref.maxF1: Maximum F<sub>1</sub>-score based on comparison to expert-based reference range map


#### Potential thresholds for converting Maxent range estimates into binary range maps

cutoff.kappa: Threshold that maximizes kappa

cutoff.spec.sens: Threshold at which sum of specificity and sensitivity is highest

cutoff.no.omission: Threshold allowing no omission of occurrence data

cutoff.prevalence: Threshold where observed and modeled prevalence are closest

cutoff.equal.sens.spec: Threshold where sensitivity and specificity are equal

cutoff.sensitivity: Threshold that maximizes sensitivity
