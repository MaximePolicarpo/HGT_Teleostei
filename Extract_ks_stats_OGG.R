##### Libraries  ---------------------------------

rm(list=ls())

set.seed(2712)

library(dplyr)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(data.table)
split_tibble <- function(tibble, column = 'col') {
  tibble %>% split(., .[,column]) %>% lapply(., function(x) x[,setdiff(names(x),column)])
}


#### Load species taxonomy  ---------------------------------

species_taxonomy <- 
  read.table("~/Horizontal_transfer_project/Genomic_data/species_taxonomy_taxid.csv",
             sep=",",
             header=FALSE)
colnames(species_taxonomy) <- 
  c("family", "order", "genus", "species", "taxid")

species_taxonomy[(species_taxonomy$species == "Amphiprion_ocellaris"),"order"] <- "Ovalentaria incertae sedis"
species_taxonomy[(species_taxonomy$species == "Stegastes_partitus"),"order"] <- "Ovalentaria incertae sedis"
species_taxonomy[(species_taxonomy$species == "Acanthochromis_polyacanthus"),"order"] <- "Ovalentaria incertae sedis"
species_taxonomy[(species_taxonomy$species == "Dicentrarchus_labrax"),"order"] <- "Eupercaria incertae sedis"
species_taxonomy[(species_taxonomy$species == "Morone_saxatilis"),"order"] <- "Eupercaria incertae sedis"
species_taxonomy[(species_taxonomy$species == "Larimichthys_crocea"),"order"] <- "Eupercaria incertae sedis sec"
species_taxonomy[(species_taxonomy$species == "Collichthys_lucidus"),"order"] <- "Eupercaria incertae sedis sec"
species_taxonomy[(species_taxonomy$species == "Nibea_albiflora"),"order"] <- "Eupercaria incertae sedis sec"
species_taxonomy[(species_taxonomy$species == "Scatophagus_argus"),"order"] <- "Eupercaria incertae sedis third"
species_taxonomy[(species_taxonomy$species == "Lates_calcarifer"),"order"] <- "Carangaria incertae sedis"
species_taxonomy[(species_taxonomy$species == "Lates_japonicus"),"order"] <- "Carangaria incertae sedis"
species_taxonomy[(species_taxonomy$species == "Toxotes_jaculatrix"),"order"] <- "Carangaria incertae sedis sec"
species_taxonomy[(species_taxonomy$species == "Parambassis_ranga"),"order"] <- "Ovalentaria incertae sedis third"
species_taxonomy[(species_taxonomy$species == "Lateolabrax_maculatus"),"order"] <- "Centrarchiformes"


toadd <- species_taxonomy %>% filter(species == "Salvelinus_namaycush") 
toadd$species <- "Salvelinus_sp_IW2_2015"
species_taxonomy <- rbind(species_taxonomy, toadd)
species_taxonomy[(species_taxonomy$species == "Salvelinus_sp_IW2_2015"),"order"] <- "Salmoniformes"

species_order <-  species_taxonomy %>% dplyr::select(species, order)


#### Import genes matching to TEs    ---------------------------------

putative_genes_TE_diamond_repbase <-
  scan("Transposable_elements_candidates.repbase.diamond.txt",
       what="character")
putative_genes_TE_blast_repbase <-
  scan("Transposable_elements_candidates.repbase.blast.txt",
       what="character")
putative_genes_TE_diamond_replib <-
  scan("Transposable_elements_candidates.repeatpeps.diamond.txt",
       what="character")
putative_genes_TE_blast_replib <-
  scan("Transposable_elements_candidates.repeatpeps.blast.txt",
       what="character")
putative_genes_TE_pfam <- 
  scan("Transposable_elements_candidates_PFAM.txt",
       what="character")



putative_TE_genes <- 
  unique(c(putative_genes_TE_diamond_repbase,putative_genes_TE_blast_repbase,
           putative_genes_TE_diamond_replib, putative_genes_TE_blast_replib,
           putative_genes_TE_pfam))

putative_OGG_TE_pfam_blast <- 
  scan("N5_OGG_Transposable_elements_candidates_PFAM_BLAST_20perc.diamond_blast_pfam.txt", what="character") 


#### Compute ks stats    ---------------------------------

args <- commandArgs(trailingOnly = TRUE)

curr_species <- args[1]

list_species <- 
  scan("all_annotated_species.txt", what="character")

non_teleost_sp <- 
  c("Polypterus_senegalus","Erpetoichthys_calabaricus","Amia_calva","Atractosteus_spatula",
    "Lepisosteus_oculatus","Acipenser_oxyrinchus_oxyrinchus","Acipenser_ruthenus","Huso_huso",
    "Polyodon_spathula")


#Look at which species to compare with
curr_order <- 
  species_taxonomy %>%
  filter(species == curr_species) %>%
  pull(order)

species_to_compare <- 
  species_taxonomy %>%
  filter(species %in% list_species) %>%
  filter(! species %in% non_teleost_sp) %>%
  filter(! order %in% curr_order) %>%
  pull(species) %>%
  unique()


#Import statistics computed on BUSCO coding sequences

myspecies_cds_df <- 
  read.table(paste("HGT_stats_Results_per_sp/", 
                   curr_species,".cds.stats",
                   sep=""),
             header=FALSE,
             sep=",")


colnames(myspecies_cds_df) <- c("OGG","seq1", "seq2", "ks", "ka", "vks", "vka", "nt_used", "codon_used", "cds_identity")

myspecies_cds_df <- 
  myspecies_cds_df %>%
  filter(ks != 9.99999900) %>% 
  filter(ks > 0) %>%
  filter(codon_used > 100) %>%
  filter(! grepl(curr_species, seq2)) 


#Remove all genes or OGG corresponding to TEs

putative_TE_genes <- gsub("---", "_", putative_TE_genes)
putative_TE_genes <- gsub("\\.", "_", putative_TE_genes)
putative_TE_genes <- gsub("-", "_", putative_TE_genes)

myspecies_cds_df <- 
  myspecies_cds_df %>%
  filter(! seq1 %in% putative_TE_genes) %>%
  filter(! seq2 %in% putative_TE_genes) %>%
  filter(! OGG %in% putative_OGG_TE_pfam_blast)



#Rename Salvelinus in the seq2 column
myspecies_cds_df$seq2 <- gsub("Salvelinus_sp_IW2-2015", "Salvelinus_sp_IW2_2015", myspecies_cds_df$seq2)
species_to_compare <- gsub("Salvelinus_sp_IW2-2015", "Salvelinus_sp_IW2_2015", species_to_compare)


#Now perform a loop to perform all pairwise comparisons

ks_stats_df <- as.data.frame(NULL)

for(second_sp in species_to_compare){
  
  myspecies_df_comp <- 
    myspecies_cds_df %>%
    filter(grepl(second_sp, seq2))
  
  
  ks_quantile_1 <- quantile(myspecies_df_comp$ks, 0.001) 
  ks_quantile_2 <- quantile(myspecies_df_comp$ks, 0.002) 
  ks_quantile_3 <- quantile(myspecies_df_comp$ks, 0.003) 
  ks_quantile_4 <- quantile(myspecies_df_comp$ks, 0.004) 
  ks_quantile_5 <- quantile(myspecies_df_comp$ks, 0.005) 
  median_ks_value <- median(myspecies_df_comp$ks)
  mean_ks_value <- mean(myspecies_df_comp$ks)
  
  curr_df <- 
    as.data.frame(
      cbind(curr_species,second_sp,ks_quantile_1, ks_quantile_2, ks_quantile_3,
            ks_quantile_4, ks_quantile_5, median_ks_value, mean_ks_value))
  
  colnames(curr_df) <- c("species_1", "species_2", "quantile_OGG_001", "quantile_OGG_002",
                         "quantile_OGG_003", "quantile_OGG_004", "quantile_OGG_005", 
                         "median_OGG","mean_OGG")
  
  ks_stats_df <- rbind(ks_stats_df, curr_df)
  
}


#Now write table :)
write.table(ks_stats_df,
            args[2],
            col.names=TRUE,
            row.names=FALSE,
            quote=FALSE,
            sep=",")

