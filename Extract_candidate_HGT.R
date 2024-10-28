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



#### Import Suspected TEs    ---------------------------------

putative_genes_TE_pfam_blast <-
  scan("Transposable_elements_candidates_PFAM_BLAST.txt",
       what="character")

putative_OGG_TE_pfam_blast <-
  scan("N5_OGG_Transposable_elements_candidates_PFAM_BLAST_30perc.txt",
       what="character")


#### Import KS quantiles (OGG + BUSCO)    ---------------------------------

#Prepare a lit of unique pairwise species combinations
list_teleost_species <- 
  scan("teleost_annotated_species.txt", what="character")
list_teleost_species <- gsub("Salvelinus_sp_IW2-2015", "Salvelinus_sp_IW2_2015", list_teleost_species)


all_comb <- apply(combn(list_teleost_species,2),2,paste,collapse=',')
pairs_df <- as.data.frame(t(sapply(all_comb, function(x) unlist(strsplit(x, ",")))))
row.names(pairs_df) <- NULL
colnames(pairs_df) <- c("species_1", "species_2")


#Import Ks statistics computed on all OGG or on BUSCO genes


OGG_ks_stats_df <- 
  read.table("All_OGG_Ks_stats.csv",
             header=TRUE,
             sep=",")
OGG_ks_stats_df <- left_join(pairs_df, OGG_ks_stats_df, by=c("species_1", "species_2"))


BUSCO_ks_stats_df <- 
  read.table("All_BUSCO_Ks_stats.csv",
             header=TRUE,
             sep=",")
BUSCO_ks_stats_df$species_2 <- gsub("Salvelinus_sp_IW2-2015", "Salvelinus_sp_IW2_2015", BUSCO_ks_stats_df$species_2 )
BUSCO_ks_stats_df <- left_join(pairs_df, BUSCO_ks_stats_df, by=c("species_1", "species_2"))


#Combine tables

OGG_BUSCO_Ks_stats_df <- 
  left_join(OGG_ks_stats_df, BUSCO_ks_stats_df, by=c("species_1", "species_2"))


#### Find candidate HGT   ---------------------------------

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


#Import statistics computed on coding sequences 

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

putative_genes_TE_pfam_blast <- gsub("---", "_", putative_genes_TE_pfam_blast)
putative_genes_TE_pfam_blast <- gsub("\\.", "_", putative_genes_TE_pfam_blast)
putative_genes_TE_pfam_blast <- gsub("-", "_", putative_genes_TE_pfam_blast)

myspecies_cds_df <- 
  myspecies_cds_df %>%
  filter(! seq1 %in% putative_genes_TE_pfam_blast) %>%
  filter(! seq2 %in% putative_genes_TE_pfam_blast) %>%
  filter(! OGG %in% putative_OGG_TE_pfam_blast)


#Rename Salvelinus in the seq2 column
myspecies_cds_df$seq2 <- gsub("Salvelinus_sp_IW2-2015", "Salvelinus_sp_IW2_2015", myspecies_cds_df$seq2)
species_to_compare <- gsub("Salvelinus_sp_IW2-2015", "Salvelinus_sp_IW2_2015", species_to_compare)

#Now perform a loop to perform all pairwise comparisons

HGT_candidate_df <- as.data.frame(NULL)

for(second_sp in species_to_compare){
  
  myspecies_df_comp <- 
    myspecies_cds_df %>%
    filter(grepl(second_sp, seq2))
  
  curr_ks_stats <- 
    OGG_BUSCO_Ks_stats_df %>%
    filter((species_1 == second_sp & species_2 == curr_species) | 
             (species_2 == second_sp & species_1 == curr_species)) 
  
  
  OGG_ks_quantile_1 <- curr_ks_stats$quantile_OGG_001
  OGG_ks_quantile_2 <- curr_ks_stats$quantile_OGG_002
  OGG_ks_quantile_3 <- curr_ks_stats$quantile_OGG_003
  OGG_ks_quantile_4 <- curr_ks_stats$quantile_OGG_004
  OGG_ks_quantile_5 <- curr_ks_stats$quantile_OGG_005
  OGG_median_ks_value <- curr_ks_stats$median_OGG
  OGG_mean_ks_value <- curr_ks_stats$mean_OGG
  

  BUSCO_ks_quantile_1 <- curr_ks_stats$quantile_busco_001
  BUSCO_ks_quantile_2 <- curr_ks_stats$quantile_busco_002
  BUSCO_ks_quantile_3 <- curr_ks_stats$quantile_busco_003
  BUSCO_ks_quantile_4 <- curr_ks_stats$quantile_busco_004
  BUSCO_ks_quantile_5 <- curr_ks_stats$quantile_busco_005
  BUSCO_median_ks_value <- curr_ks_stats$median_busco_001
  BUSCO_mean_ks_value <- curr_ks_stats$mean_busco_001
  
  

  #Find the largest quantile value (BUSCO and OGG)
  
  max_quantile_value <- 
    max(c(OGG_ks_quantile_1, OGG_ks_quantile_2, OGG_ks_quantile_3,
        OGG_ks_quantile_4, OGG_ks_quantile_5, BUSCO_ks_quantile_1,
        BUSCO_ks_quantile_2, BUSCO_ks_quantile_3, BUSCO_ks_quantile_4, 
        BUSCO_ks_quantile_5))
  
  #Filter at this "max" quantile value
  myspecies_df_comp_filt <- 
    myspecies_df_comp %>%
    filter(ks < max_quantile_value) %>%
    arrange(ks)
  
  myspecies_df_comp_filt <- 
    myspecies_df_comp_filt %>%
    mutate(species_1 = curr_species) %>%
    mutate(species_2 = second_sp)  %>%
    mutate(quantile_value_001_ks_OGG = OGG_ks_quantile_1) %>%
    mutate(quantile_value_002_ks_OGG = OGG_ks_quantile_2) %>%
    mutate(quantile_value_003_ks_OGG = OGG_ks_quantile_3) %>%
    mutate(quantile_value_004_ks_OGG = OGG_ks_quantile_4) %>%
    mutate(quantile_value_005_ks_OGG = OGG_ks_quantile_5) %>%
    mutate(median_ks_OGG = OGG_median_ks_value) %>%
    mutate(mean_ks_OGG = OGG_mean_ks_value) %>%
    mutate(quantile_value_001_ks_BUSCO = BUSCO_ks_quantile_1) %>%
    mutate(quantile_value_002_ks_BUSCO = BUSCO_ks_quantile_2) %>%
    mutate(quantile_value_003_ks_BUSCO = BUSCO_ks_quantile_3) %>%
    mutate(quantile_value_004_ks_BUSCO = BUSCO_ks_quantile_4) %>%
    mutate(quantile_value_005_ks_BUSCO = BUSCO_ks_quantile_5) %>%
    mutate(median_ks_BUSCO = BUSCO_median_ks_value) %>%
    mutate(mean_ks_BUSCO = BUSCO_mean_ks_value)
  
  HGT_candidate_df <- rbind(HGT_candidate_df, myspecies_df_comp_filt)
  
}


#Now write tables :)
write.table(HGT_candidate_df,
            args[2],
            col.names=TRUE,
            row.names=FALSE,
            quote=FALSE,
            sep=",")


