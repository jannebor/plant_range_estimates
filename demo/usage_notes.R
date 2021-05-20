library(raster)
library(viridis)
library(smoothr)

setwd("X:/indecol/USERS/JanB/Data")

metadata<-read.csv("C:/Users/janbor/Desktop/OneDrive - NTNU/Data/iucn_gbif/Idun/dataset/metadata.csv")



#random species from manuscript
species<-c("Cedrus libani","Laburnum anagyroides","Amomum pterocarpum",
            "Megistostegium nodulosum","Eucalyptus elliptica","Mammillaria grahamii",
            "Memecylon elegantulum","Siparuna conica","Psidium salutare",
            "Pyracantha angustifolia","Terminalia macrostachya","Trisetaria dufourei")

species_name<-species[3]

species_name<-"Cornutia odorata"
species_name<-"Pinus sylvestris"
###speciesIDs can be looked up in the metadata
speciesID <- metadata$speciesID[which(metadata$scientificname==species_name)][1]
#loading occ points
load(paste("iucn_gbif/Ranges/Model4/",species_name,"/",species_name,"_occ_points2", sep=""))

#a random species
speciesID <- metadata$speciesID[which(metadata$scientificname==metadata$scientificname[runif(1, min=1, max=nrow(metadata))])][1]
#loading occ points
load(paste("iucn_gbif/Ranges/Model4/",metadata$scientificname[which(metadata$speciesID==speciesID)][1],"/",metadata$scientificname[which(metadata$speciesID==speciesID)][1],"_occ_points2", sep=""))

pred_ras <- raster("C:/Users/janbor/Desktop/OneDrive - NTNU/Data/iucn_gbif/Idun/dataset/range_estimates_suggested.nc", band = speciesID)
pred_ras_cropped<-crop(pred_ras, extent(metadata$extent.xmin[which(metadata$speciesID==speciesID)][1],
                                        metadata$extent.xmax[which(metadata$speciesID==speciesID)][1],
                                        metadata$extent.ymin[which(metadata$speciesID==speciesID)][1],
                                        metadata$extent.ymax[which(metadata$speciesID==speciesID)][1]))
par(mfrow=c(1,2))
#plotting prediction
plot(pred_ras_cropped, col=viridis(512, direction = -1),
     main=paste("AUC = ", metadata$cell.AUC[which(metadata$speciesID==speciesID)]))
plot(occ_points, pch=21, col="black", bg="white", add=T)


#plotting range map
presence_cells<-rasterize(occ_points, pred_ras_cropped, fun="count")
presence_cells<-presence_cells$ID
presence_cells[!is.na(presence_cells)]<-1
presence_cells <- rasterToPolygons(presence_cells,dissolve=T)
presence_cells<-smoothr::smooth(presence_cells, method="ksmooth", smoothness=2)



range_map_pos<-pred_ras_cropped
range_map_pos[range_map_pos<metadata$cutoff.prevalence[which(metadata$speciesID==speciesID)][1]]<-NA
range_map_pos[range_map_pos>=metadata$cutoff.prevalence[which(metadata$speciesID==speciesID)][1]]<-1
range_map_pos <- rasterToPolygons(range_map_pos,dissolve=T)
range_map_pos<-smoothr::smooth(range_map_pos, method="ksmooth", smoothness=2)
range_map_pro<-pred_ras_cropped
range_map_pro[range_map_pro<metadata$cutoff.no.omission[which(metadata$speciesID==speciesID)][1]]<-NA
range_map_pro[range_map_pro>=metadata$cutoff.no.omission[which(metadata$speciesID==speciesID)][1]]<-1
range_map_pro <- rasterToPolygons(range_map_pro,dissolve=T)
range_map_pro<-smoothr::smooth(range_map_pro, method="ksmooth", smoothness=2)


plot(pred_ras_cropped, col="grey")
plot(range_map_pos, col=viridis(3)[1], add=T)
plot(range_map_pro, col=viridis(3)[2], add=T)
plot(presence_cells, col=viridis(3)[3], add=T)

par(mfrow=c(1,1))


presence_cells<-rasterize(occ_points, pred_ras_cropped, fun="count")
presence_cells<-presence_cells$ID
presence_cells[!is.na(presence_cells)]<-1
presence_cells <- rasterToPolygons(presence_cells,dissolve=T)
presence_cells<-smoothr::smooth(presence_cells, method="ksmooth", smoothness=2)



range_map_pos<-pred_ras_cropped
range_map_pos[range_map_pos<metadata$cutoff.no.omission[which(metadata$speciesID==speciesID)][1]]<-NA
range_map_pos[range_map_pos>=metadata$cutoff.no.omission[which(metadata$speciesID==speciesID)][1]]<-1
range_map_pos <- rasterToPolygons(range_map_pos,dissolve=T)
range_map_pos<-smoothr::smooth(range_map_pos, method="ksmooth", smoothness=2)
range_map_pro<-pred_ras_cropped
range_map_pro[range_map_pro<metadata$cutoff.sensitivity[which(metadata$speciesID==speciesID)][1]]<-NA
range_map_pro[range_map_pro>=metadata$cutoff.sensitivity[which(metadata$speciesID==speciesID)][1]]<-1
range_map_pro <- rasterToPolygons(range_map_pro,dissolve=T)
range_map_pro<-smoothr::smooth(range_map_pro, method="ksmooth", smoothness=2)


plot(pred_ras_cropped, col="grey")
plot(range_map_pos, col=viridis(3)[1], add=T)
plot(range_map_pro, col=viridis(3)[2], add=T)
plot(presence_cells, col=viridis(3)[3], add=T)







