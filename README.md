## Native range estimates for red-listed vascular plant species

Distribution data for individual species (e.g. range maps) is central for understanding both global biodiversity patterns and associated anthropogenic impacts. Our aim is to provide a comprehensive set of spatial data for a large number of, so far unmapped, vascular plant species, supporting the development of future assessments of global biodiversity impacts and distributions.


![Range map examples - Maxent predictions transformed into potential binary range maps at different levels of confidence using different cut-off thresholds](demo/slideshow.gif)


### Dataset

The dataset consists of native regions for 47,675 species retrieved from [Plants of the World online](http://plantsoftheworldonline.org/), density of available native occurrence records for 30,906 species retrieved from [GBIF](https://www.gbif.org/), and standardized, large-scale [Maxent](https://biodiversityinformatics.amnh.org/open_source/maxent/) range estimates for 27,208 species, highlighting environmentally suitable areas within species' native regions. The data can be explored using our  [data explorer](https://plant-ranges.indecol.no/).

The full dataset can be downloaded in 30 arc minutes spatial resolution (approx. 56 km at the equator).

The suggested dataset “range_data.nc” contains the gathered and generated data allocated to three variables that can be called by specifying a varname: rasterized native WGSPRD-regions (varname: Native region) the number of occurrence records per cell (varname: Presence cells), and the best performing Maxent prediction (cloglog transformed, varname: Maxent prediction) for each species selected based on the highest harmonic mean between AUC and AUC<sub>PR</sub>.

“complete_dataset.zip” contains all generated netCDF files in cloglog and raw output format. The files “range_data_full.nc” and “range_data_full_rawOutput.nc” contain one Maxent prediction for each occurrence data type (i.e. varname: Model 1, Model 2 or Model 3), while “range_data.nc” and “range_data_rawOutput.nc” contain only the best performing Maxent prediction for each species (varname: Maxent prediction), selected based on the highest harmonic mean between AUC and AUCPR (see section Technical Validation). Number of occurrence records per cell (varname: Presence cells) and rasterized native WGSRPD-regions (varname: Native region) are provided in all netCDF files.

Each band in the netCDF files assembles the mentioned variables for one species. The corresponding bands can be looked up in the metadata (i.e. speciesID).

The metadata can be used to filter models based on species, performance or desired data types, to look up the relevant study extent for masking individual predictions, or to select different cut-off thresholds for generating binary range maps. For information about the given thresholds visit <https://www.rdocumentation.org/packages/dismo/versions/1.1-4/topics/threshold>.

------------------------------------------------------------------------

### Examples

Species selection based on metadata:

```{r setup, include = FALSE}
library(raster)

# metadata of the suggested dataset (only the best prediction for each species)
metadata<-read.csv(".../metadata.csv")

# for a given species, e.g.:
species_name <- "Amomum pterocarpum"

# speciesID
speciesID <- metadata$speciesID[which(metadata$scientificname==species_name)][1]
# type of provided data
metadata$data[which(metadata$scientificname==species_name)][1]
# and, if a Maxent prediction has been generated, the best performing model
best_model <- metadata$model[which(metadata$scientificname==species_name)][1]
# can be looked up in the metadata

# each band of the netCDF file assembles data for one species
# variables represent the different data types within each band
variable <- paste(best_model)

ras <- raster("C:/Users/janbor/Documents/GitLab/plant_ranges/app/demo/range_data.nc", varname = variable, band = speciesID)
# varname needs to be one of "Native region", "Presence cells", or a specific Model: "Model 1", "Model 2", "Model 3"

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

read_data<-function(variables, speciesID, filelocation){
  
  ras<-stack()
  
  for(v in variables){
  ras<-stack(ras,raster(filelocation, varname = v, band = speciesID))
  }
  
  return(crop(ras, extent(metadata$extent.xmin[which(metadata$speciesID==speciesID)][1],
                          metadata$extent.xmax[which(metadata$speciesID==speciesID)][1],
                          metadata$extent.ymin[which(metadata$speciesID==speciesID)][1],
                          metadata$extent.ymax[which(metadata$speciesID==speciesID)][1]))
  )
}

species_name <- "Amomum pterocarpum"
speciesID <- metadata$speciesID[which(metadata$scientificname==species_name)][1]
best_model <- metadata$model[which(metadata$scientificname==species_name)][1]
variables <- c("Native region","Presence cells",paste(best_model))

ras<-read_data(variables = variables, speciesID = speciesID, filelocation = ".../range_data.nc")
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
