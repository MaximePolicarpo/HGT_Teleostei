###### Import libraries -------


set.seed(2712)

rm(list = ls())


library("dplyr")
library(tidyr)
library(ape)
library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library(gt)
library("phytools")
library(adephylo)
library(viridis)
library(ggnewscale)
library(robis)
library(ggmap)
library(mregions)
library(geojson)
library(sp)
library(rgbif)

###### Initate a world map with ggmap  -------


world <- map_data("world")

worldmap <- 
  ggplot(world, aes(x=long, y=lat)) +
  geom_polygon(aes(group=group, alpha=0.1)) +
  scale_y_continuous(breaks = (-2:2) * 30) +
  scale_x_continuous(breaks = (-4:4) * 45) +
  theme(
    panel.background = element_rect(fill = "#F0F8FF"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  coord_equal() +
  theme(legend.position = "none")

###### Lets look at our species distribution -------

my_occs.Clupea_harengus <- occurrence(scientificname = "Clupea harengus")
my_occs.Hypomesus_transpacificus <- occurrence(scientificname = "Hypomesus transpacificus")
my_occs.Osmerus_eperlanus <- occurrence(scientificname = "Osmerus eperlanus")


pdf("Map_Osmeriformes_Clupea.pdf", width = 10, height = 7)

worldmap + 
  geom_point(data = my_occs.Clupea_harengus, 
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "#CC79A7", shape = 21, alpha = 2/3) + 
  geom_point(data = my_occs.Hypomesus_transpacificus, 
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "#009E73", shape = 21, alpha = 2/3) +
  geom_point(data = my_occs.Osmerus_eperlanus, 
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "#034E39", shape = 21, alpha = 2/3)

dev.off()


my_occs.Paramormyrops_kingsleyae <- occ_data(scientificName = "Paramormyrops kingsleyae")
my_occs.Clarias_gariepinus <- occ_data(scientificName = "Clarias gariepinus")

pdf("Map_Paramormyrops_Clarias.pdf", width = 10, height = 7)

worldmap + 
  geom_point(data = my_occs.Paramormyrops_kingsleyae$data, 
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "#45AF91", shape = 21, alpha = 2/3) +
  geom_point(data = my_occs.Clarias_gariepinus$data, 
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "#901E5D", shape = 21, alpha = 2/3) 
dev.off()



my_occs.Thalassophryne_megalops <- occurrence(scientificname = "Thalassophryne megalops")
my_occs.Thalassophryne_amazonica <- occ_data(scientificName = "Thalassophryne amazonica")$data
my_occs.Thalassophryne_maculosa <- occurrence(scientificname = "Thalassophryne maculosa")
my_occs.Thalassophryne_punctata <- occurrence(scientificname = "Thalassophryne punctata")
my_occs.Opsanus_beta <- occurrence(scientificname = "Opsanus beta")


pdf("Map_Thalassophryne_Clupea.pdf", width = 10, height = 7)

worldmap + 
  geom_point(data = my_occs.Clupea_harengus, 
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "#CC79A7", shape = 21, alpha = 2/3) + 
  geom_point(data = my_occs.Opsanus_beta, 
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "#D55E00", shape = 21, alpha = 2/3) + 
  geom_point(data = my_occs.Thalassophryne_megalops, 
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "#E69F00", shape = 21, alpha = 2/3) +
  geom_point(data = my_occs.Thalassophryne_amazonica, 
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "#E69F00", shape = 21, alpha = 2/3) +
  geom_point(data = my_occs.Thalassophryne_maculosa, 
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "#E69F00", shape = 21, alpha = 2/3) +
  geom_point(data = my_occs.Thalassophryne_punctata, 
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "#E69F00", shape = 21, alpha = 2/3) 

dev.off()


my_occs.Chanos_chanos <- occurrence(scientificname = "Chanos chanos")
my_occs.Neoarius_graeffei <- occurrence(scientificname = "Neoarius graeffei")


pdf("Map_Chanos_Neoarius.pdf", width = 10, height = 7)

worldmap + 
  geom_point(data = my_occs.Chanos_chanos, 
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "#1B9A76", shape = 21, alpha = 2/3) + 
  geom_point(data = my_occs.Neoarius_graeffei, 
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "#901E5D", shape = 21, alpha = 2/3) 

dev.off()

my_occs.Megalops_cyprinoides <- occurrence(scientificname = "Megalops cyprinoides")
my_occs.Centroberyx_gerrardi <- occurrence(scientificname = "Centroberyx gerrardi")


pdf("Map_Megalops_Beryx.pdf", width = 10, height = 7)

worldmap + 
  geom_point(data = my_occs.Megalops_cyprinoides, 
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "#1B9A76", shape = 21, alpha = 2/3) + 
  geom_point(data = my_occs.Centroberyx_gerrardi, 
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "#E69F00", shape = 21, alpha = 2/3) 

dev.off()



my_occs.Clupea_harengus <- occurrence(scientificname = "Clupea harengus")
my_occs.Alosa_alosa <- occurrence(scientificname = "Alosa alosa")
my_occs.Sardina_pilchardus <- occurrence(scientificname = "Sardina pilchardus")
my_occs.Caranx_melampygus <- occurrence(scientificname = "Caranx melampygus")
my_occs.Seriola_lalandi <- occurrence(scientificname = "Seriola lalandi")
my_occs.Lates_calcarifer <- occurrence(scientificname = "Lates calcarifer")


pdf("Map_Clupeiformes_Caranx.pdf", width = 10, height = 7)

worldmap + 
  geom_point(data = my_occs.Clupea_harengus, 
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "#CC79A7", shape = 21, alpha = 2/3) +
  geom_point(data = my_occs.Alosa_alosa, 
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "#901E5D", shape = 21, alpha = 2/3) +
  geom_point(data = my_occs.Sardina_pilchardus, 
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "#DC9FC1", shape = 21, alpha = 2/3) +
  geom_point(data = my_occs.Caranx_melampygus, 
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "#0979B7", shape = 21, alpha = 2/3) +
  geom_point(data = my_occs.Seriola_lalandi, 
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "#62C7FF", shape = 21, alpha = 2/3) +
  geom_point(data = my_occs.Lates_calcarifer, 
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "darkblue", shape = 21, alpha = 2/3) 


dev.off()




my_occs.Thunnus_maccoyi <- occ_data(scientificName = "Thunnus maccoyi")$data #no data in obis nor gbif
my_occs.Scomber_scombrus <-  occurrence(scientificname = "Scomber scombrus")
my_occs.Thunnus_atlanticus <-  occurrence(scientificname = "Thunnus atlanticus")


pdf("Map_Thunnus_Scomber.pdf", width = 10, height = 7)

worldmap + 
  geom_point(data = my_occs.Clupea_harengus, 
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "#CC79A7", shape = 21, alpha = 2/3) +
  geom_point(data = my_occs.Alosa_alosa, 
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "#901E5D", shape = 21, alpha = 2/3) +
  geom_point(data = my_occs.Sardina_pilchardus, 
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "#DC9FC1", shape = 21, alpha = 2/3)  +
  geom_point(data = my_occs.Thunnus_atlanticus, 
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "#009E73", shape = 21, alpha = 2/3) +
  geom_point(data = my_occs.Scomber_scombrus, 
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "#034E39", shape = 21, alpha = 2/3)

dev.off()

my_occs.Pangasius_djambal <-  occ_data(scientificName = "Pangasius djambal")$data
my_occs.Pangasianodon_hypophthalmus <-  occ_data(scientificName = "Pangasianodon hypophthalmus")$data
my_occs.Onychostoma_macrolepis <-  occ_data(scientificName = "Onychostoma macrolepis")$data
my_occs.Labeo_rohita <-  occ_data(scientificName = "Labeo_rohita")$data


pdf("Map_Pangas_Cypriniforme.pdf", width = 10, height = 7)

worldmap + 
  geom_point(data = my_occs.Pangasius_djambal, 
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "#CC79A7", shape = 21, alpha = 2/3) +
  geom_point(data = my_occs.Pangasianodon_hypophthalmus, 
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "#901E5D", shape = 21, alpha = 2/3) +
  geom_point(data = my_occs.Onychostoma_macrolepis, 
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "#009E73", shape = 21, alpha = 2/3) +
  geom_point(data = my_occs.Carassius_auratus, 
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "olivedrab1", shape = 21, alpha = 2/3)

dev.off()
