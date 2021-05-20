## Native range estimates for red-listed vascular plant species

Distribution data for individual species (e.g. range maps) is central for understanding both global biodiversity patterns and associated anthropogenic impacts. Our aim is to provide a comprehensive set of spatial data for a large number of, so far unmapped, vascular plant species, supporting the development of future assessments of global biodiversity impacts and distributions.


![Range map examples - Maxent predictions transformed into potential binary range maps at different levels of confidence using different cut-off thresholds](demo/slideshow.gif)


### Dataset

The dataset consists of native regions for 47,675 species retrieved from [Plants of the World online](http://plantsoftheworldonline.org/), density of available native occurrence records for 30,906 species retrieved from [GBIF](https://www.gbif.org/), and standardized, large-scale [Maxent](https://biodiversityinformatics.amnh.org/open_source/maxent/) range estimates for 27,208 species, highlighting environmentally suitable areas within species' native regions. The data can be explored using our  [data explorer](https://plant-ranges.indecol.no/).

The full dataset can be downloaded [here](https://www.dropbox.com/sh/meeb6ru84778k94/AADrdCleHeMujip60C7EuMH1a?dl=1) at 30 arc minutes spatial resolution (approx. 56 km at the equator).

["range_data.nc""](https://www.dropbox.com/s/vqoiep3y8703yp5/range_estimates.nc?dl=1) contains all generated predictions based on the different models (i.e. raw data, presence cells and thinned presence cells), along with density of occurrence records and rasterized native regions.

["range_estimates_suggested.nc"](https://www.dropbox.com/s/9az5yoq4ayx0139/range_estimates_suggested.nc?dl=1) contains only the best performing prediction for each species, selected based on AUC and AUC<sub>PR</sub>.

The [metadata](https://www.dropbox.com/s/ktf1pk6hsk62d80/metadata_full.csv?dl=1) can be used to select appropriate cut-off thresholds for generating binary range maps, filter models based on species, performance or desired data types, and to look up the relevant study extent for masking individual predictions.

Underlying Maxent models are available as R objects upon request.

For information about the given thresholds visit <https://www.rdocumentation.org/packages/dismo/versions/1.1-4/topics/threshold>.

------------------------------------------------------------------------

### Examples

Species selection based on metadata:

```{r setup, include = FALSE}
library(raster)

metadata<-read.csv(".../metadata.csv")
#or
metadata_full<-read.csv(".../metadata_full.csv")

species_name <- "Amomum pterocarpum"

###speciesIDs can be looked up in the metadata
speciesID <- metadata$speciesID[which(metadata$scientificname==species_name)][1]

pred_ras <- raster(".../range_estimates_suggested.nc", band = speciesID)
#or
pred_ras <- raster(".../range_estimates.nc", varname = "Model 3", band = speciesID)
###varname needs to be "Model 1", "Model 2" or "Model 3"; "Occurrence records"; "Native region"

plot(pred_ras)


```

Cropping the predictions to original extent:

```{r setup, include = FALSE}

pred_ras_cropped<-crop(pred_ras, extent(metadata$extent.xmin[which(metadata$scientificname==species_name)][1],
                      metadata$extent.xmax[which(metadata$scientificname==species_name)][1],
                      metadata$extent.ymin[which(metadata$scientificname==species_name)][1],
                      metadata$extent.ymax[which(metadata$scientificname==species_name)][1]))

plot(pred_ras_cropped)

```

### Metadata

#### General species information

speciesID, scientificname, redlistcategory, rank, family, order, class, subdivision, division, and data (type of provided data out of Native region, Occurrence records or Maxent prediction)


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
