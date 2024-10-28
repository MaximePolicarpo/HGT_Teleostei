##### Libraries  ---------------------------------

rm(list=ls())

set.seed(2712)


library("caper")
library("ape")
library(dplyr)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(adephylo)
library(phylobase)
library(data.table)
library(phytools)
library("geiger")
library(phangorn)
library(castor)
split_tibble <- function(tibble, column = 'col') {
  tibble %>% split(., .[,column]) %>% lapply(., function(x) x[,setdiff(names(x),column)])
}



#### Load species taxonomy  ---------------------------------

species_taxonomy <- 
  read.table("species_taxonomy_taxid.csv",
             sep=",",
             header=FALSE)
colnames(species_taxonomy) <- 
  c("family", "order", "genus", "species", "taxid")


list_teleost_species <- 
scan("teleost_annotated_species.txt", what="character")
species_tree <- read.tree("AMAS_concatenated_alignment_BUSCO.fa.timetree.nwk")
teleost_tree <- keep.tip(species_tree, list_teleost_species)
teleost_tree$edge.length <- teleost_tree$edge.length * 1000



species_taxonomy <- 
  read.table("species_taxonomy_taxid.csv",
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




#### Import Suspected TEs    ---------------------------------

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



length(unique(c(putative_genes_TE_diamond_repbase,putative_genes_TE_blast_repbase,
                putative_genes_TE_diamond_replib, putative_genes_TE_blast_replib,
                putative_genes_TE_pfam)))


putative_TE_genes <- 
  unique(c(putative_genes_TE_diamond_repbase,putative_genes_TE_blast_repbase,
           putative_genes_TE_diamond_replib, putative_genes_TE_blast_replib,
           putative_genes_TE_pfam))




##### Strategy A -- Order method  ---------------------------------


args <- commandArgs(trailingOnly = TRUE)

#curr_species <- "Hypomesus_transpacificus"

curr_species <- args[1]


list_species <- 
  scan("all_annotated_species.txt", what="character")
non_teleost_sp <- 
  c("Polypterus_senegalus","Erpetoichthys_calabaricus","Amia_calva","Atractosteus_spatula",
    "Lepisosteus_oculatus","Acipenser_oxyrinchus_oxyrinchus","Acipenser_ruthenus","Huso_huso",
    "Polyodon_spathula")

#Define the order and species belonging to the same order ... 
curr_order <- 
  species_taxonomy %>%
  filter(species == curr_species) %>%
  pull(order)

species_same_order_df <- 
  species_taxonomy %>%
  filter(order == curr_order) %>% 
  filter(species != curr_species) %>%
  filter(species %in% list_species) 


same_order_species <- species_same_order_df %>% pull(species)



#lets continue only if there is at-least one other species in the order ... 

if(nrow(species_same_order_df) > 0){
  blastn_df <- 
    fread(paste("Blastn_results/", 
                curr_species,".blastn",
                sep=""),
          header=FALSE)
  blastn_df <- as.data.frame(blastn_df)
  
  colnames(blastn_df) <- 
    c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", 
      "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "gaps")
  
  
  #Remove putative TEs from qseqids
  
  blastn_df <- 
    blastn_df %>%
    filter(! qseqid %in% putative_TE_genes) %>%
    filter(! sseqid %in% putative_TE_genes) 
  
  #Add species to the table and remove non teleost species 
  
  blastn_df <- 
    blastn_df %>%
    mutate(qseqsp = curr_species)
  
  
  blastn_df <- 
    as.data.frame(
      blastn_df %>%
        rowwise() %>%
        mutate(sseqsp = gsub("---.*","",sseqid))
    ) %>%
    filter(! sseqsp %in% non_teleost_sp)
  
  
  #Add the target order .. 
  
  sp_order <- species_taxonomy %>% dplyr::select(species, order)
  colnames(sp_order) <- c("sseqsp", "sseqorder")
  
  blastn_df <- 
    left_join(blastn_df, 
              sp_order,
              by="sseqsp")
  

  
  
  blastn_df_bh <- 
    as.data.frame(
      blastn_df %>%
        group_by(qseqid) %>%
        filter(evalue == min(evalue)) %>%
        mutate(has_same_order = any(sseqorder == curr_order)) %>%
        filter(if_else(has_same_order, sseqorder == curr_order, TRUE)) %>%
        slice(1) %>%
        ungroup() %>%
        select(-has_same_order)
    )
  
  

  #Keep all blast results where the best match does not belong to the same order
  
  Outliers_orders_df <- 
    blastn_df_bh %>%
    filter(! sseqsp %in% same_order_species) %>%
    filter(length > 300) %>%
    filter(length > (qlen * 85 / 100)) %>%
    mutate(abs_diff_len = abs(qlen - slen)) %>%
    filter(abs_diff_len < 50) 
  
  write.table(Outliers_orders_df,
              args[2],
              col.names=TRUE,
              row.names=FALSE,
              quote=FALSE,
              sep=",")
  
  
}





##### Strategy B -- Distance method w same order species  ---------------------------------

args <- commandArgs(trailingOnly = TRUE)

#curr_species <- "Hypomesus_transpacificus"
#curr_species <- "Pangasianodon_hypophthalmus"


curr_species <- args[1]


list_species <- 
  scan("~/Horizontal_transfer_project/BetweenActino_HGT/all_annotated_species.txt", what="character")
non_teleost_sp <- 
  c("Polypterus_senegalus","Erpetoichthys_calabaricus","Amia_calva","Atractosteus_spatula",
    "Lepisosteus_oculatus","Acipenser_oxyrinchus_oxyrinchus","Acipenser_ruthenus","Huso_huso",
    "Polyodon_spathula")

#Define the order and species belonging to the same order ... 
curr_order <- 
  species_taxonomy %>%
  filter(species == curr_species) %>%
  pull(order)

species_same_order_df <- 
  species_taxonomy %>%
  filter(order == curr_order) %>% 
  filter(species != curr_species) %>%
  filter(species %in% list_species) 


same_order_species <- species_same_order_df %>% pull(species)

blastn_df <- 
  fread(paste("Blastn_results/", 
              curr_species,".blastn",
              sep=""),
        header=FALSE)
blastn_df <- as.data.frame(blastn_df)

colnames(blastn_df) <- 
  c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", 
    "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "gaps")


#Remove putative TEs from qseqids

blastn_df <- 
  blastn_df %>%
  filter(! qseqid %in% putative_TE_genes) %>%
  filter(! sseqid %in% putative_TE_genes) 

#Add species to the table and remove non teleost species 

blastn_df <- 
  blastn_df %>%
  mutate(qseqsp = curr_species)


blastn_df <- 
  as.data.frame(
    blastn_df %>%
      rowwise() %>%
      mutate(sseqsp = gsub("---.*","",sseqid))
  ) %>%
  filter(! sseqsp %in% non_teleost_sp)


#Add the target order  

sp_order <- species_taxonomy %>% dplyr::select(species, order)
colnames(sp_order) <- c("sseqsp", "sseqorder")

blastn_df <- 
  left_join(blastn_df, 
            sp_order,
            by="sseqsp")




#Add divergence times between species

blastn_df <- 
  as.data.frame(blastn_df %>%
                  rowwise() %>%
                  mutate(div_time = 
                           get_pairwise_distances(teleost_tree, 
                                                  qseqsp, 
                                                  sseqsp, 
                                                  as_edge_counts=FALSE, 
                                                  check_input=TRUE)/2))



#Keep only the best hit per query


blastn_df_bh <- 
  blastn_df %>%
  group_by(qseqid) %>%
  filter(evalue == min(evalue)) %>%
  slice(which.min(div_time))


#Add a round divergence time


blastn_df_bh <- 
  as.data.frame(blastn_df_bh %>%
                  rowwise() %>%
                  mutate(div_time_round = round(div_time, digits=0)))



total_nb_query <- nrow(blastn_df_bh)


#Extract last categories < 0.1

Div_time_01 <- 
  blastn_df_bh %>%
  group_by(div_time_round) %>%
  summarise(count = n()) %>%
  arrange(desc(div_time_round)) %>%
  mutate(cumulative_count = cumsum(count)) %>%
  mutate(cumulative_proportion = cumulative_count/total_nb_query) %>%
  filter(cumulative_proportion < 0.1) %>%
  arrange(div_time_round) %>%
  head(1) %>%
  pull(div_time_round)


#Extract Blast outlier and filter on query length ...


Blast_outlier_df_stratB <- 
  blastn_df_bh %>%
  filter(div_time_round >= Div_time_01) %>%
  filter(length > 300) %>%
  filter(length > (qlen * 85 / 100)) %>%
  mutate(abs_diff_len = abs(qlen - slen)) %>%
  filter(abs_diff_len < 50) 


write.table(Blast_outlier_df_stratB,
            args[3],
            col.names=TRUE,
            row.names=FALSE,
            quote=FALSE,
            sep=",")






##### Strategy C -- HGT index bitscore   ---------------------------------


args <- commandArgs(trailingOnly = TRUE)

#curr_species <- "Hypomesus_transpacificus"

curr_species <- args[1]


list_species <- 
  scan("~/Horizontal_transfer_project/BetweenActino_HGT/all_annotated_species.txt", what="character")
non_teleost_sp <- 
  c("Polypterus_senegalus","Erpetoichthys_calabaricus","Amia_calva","Atractosteus_spatula",
    "Lepisosteus_oculatus","Acipenser_oxyrinchus_oxyrinchus","Acipenser_ruthenus","Huso_huso",
    "Polyodon_spathula")

#Define the order and species belonging to the same order ... 
curr_order <- 
  species_taxonomy %>%
  filter(species == curr_species) %>%
  pull(order)

species_same_order_df <- 
  species_taxonomy %>%
  filter(order == curr_order) %>% 
  filter(species != curr_species) %>%
  filter(species %in% list_species) 


same_order_species <- species_same_order_df %>% pull(species)



#lets continue only if there is at-least one other species in the order ... 

if(nrow(species_same_order_df) > 0){
  blastn_df <- 
    fread(paste("Blastn_results/", 
                curr_species,".blastn",
                sep=""),
          header=FALSE)
  blastn_df <- as.data.frame(blastn_df)
  
  colnames(blastn_df) <- 
    c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", 
      "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "gaps")
  
  
  #Remove putative TEs from qseqids
  
  blastn_df <- 
    blastn_df %>%
    filter(! qseqid %in% putative_TE_genes) %>%
    filter(! sseqid %in% putative_TE_genes) 
  
  #Add species to the table and remove non teleost species 
  
  blastn_df <- 
    blastn_df %>%
    mutate(qseqsp = curr_species)
  
  
  blastn_df <- 
    as.data.frame(
      blastn_df %>%
        rowwise() %>%
        mutate(sseqsp = gsub("---.*","",sseqid))
    ) %>%
    filter(! sseqsp %in% non_teleost_sp)
  
  
  #Add the target order .. 
  
  sp_order <- species_taxonomy %>% dplyr::select(species, order)
  colnames(sp_order) <- c("sseqsp", "sseqorder")
  
  blastn_df <- 
    left_join(blastn_df, 
              sp_order,
              by="sseqsp")
  

  bh_bitscore_sameorder_df <- 
    as.data.frame(
      blastn_df %>%
        filter(length > 300) %>%
        filter(length > (qlen * 85 / 100)) %>%
        mutate(abs_diff_len = abs(qlen - slen)) %>%
        filter(abs_diff_len < 50) %>% 
        group_by(qseqid) %>%
        filter(bitscore == max(bitscore)) %>%
        slice_max(bitscore, n = 1, with_ties = FALSE) %>%
        ungroup() 
    ) %>%
    dplyr::select(qseqid, sseqid, bitscore, sseqsp, sseqorder)
  
  bh_bitscore_sameorder_df  <- bh_bitscore_sameorder_df %>% filter(sseqorder == curr_order)
  
  colnames(bh_bitscore_sameorder_df) <- 
    c("qseqid", "ingroup_sseqid", "ingroup_bitscore",
      "ingroup_sseqsp", "ingroup_sseqorder")
  
 
  #Keep only the best hit not belonging to the same order based in the bitscore
  
  bh_bitscore_Notsameorder_df <- 
    as.data.frame(
      blastn_df %>%
        rowwise() %>%
        filter(length > 300) %>%
        filter(length > (qlen * 85 / 100)) %>%
        mutate(abs_diff_len = abs(qlen - slen)) %>%
        filter(abs_diff_len < 100) %>% 
        group_by(qseqid) %>%
        filter(sseqorder != curr_order) %>%
        filter(bitscore == max(bitscore)) %>%
        slice_max(bitscore, n = 1, with_ties = FALSE) %>%
        ungroup() 
    ) %>%
    dplyr::select(qseqid, sseqid, bitscore, sseqsp, sseqorder)
  
  
  colnames(bh_bitscore_Notsameorder_df) <- 
    c("qseqid", "outgroup_sseqid", "outgroup_bitscore",
      "outgroup_sseqsp", "outgroup_sseqorder")
  
  #Merge both tables and compare bitscores to compute a HGTindex
  
  bh_bitscore_df <- 
    full_join(bh_bitscore_sameorder_df, bh_bitscore_Notsameorder_df,
              by="qseqid")
  
  bh_bitscore_df <- 
    bh_bitscore_df %>%
    mutate(HGTindex = outgroup_bitscore/ingroup_bitscore) %>%
    filter(! is.na(HGTindex))
  

  #For each species pairwise comparisons, compute a HGT index distribution
  
  pairwise_combinations <- 
    bh_bitscore_df %>%
    filter(HGTindex > 0.85) %>%
    dplyr::select(ingroup_sseqsp, outgroup_sseqsp) %>%
    distinct()
  
  nb_comp <- nrow(pairwise_combinations)
  
  print(paste(nb_comp, "pairwise species to compare..."))
  
  all_pairwise_quantile_index <- as.data.frame(NULL)
  for(curr_comp in 1:nrow(pairwise_combinations)){
    
    test_sp1 <- pairwise_combinations[curr_comp,]$ingroup_sseqsp
    test_sp2 <- pairwise_combinations[curr_comp,]$outgroup_sseqsp
    
    bh_bitscore_all_notsameorder_df <- 
      as.data.frame(
        blastn_df %>%
          rowwise() %>%
          filter(length > 300) %>%
          filter(length > (qlen * 85 / 100)) %>%
          mutate(abs_diff_len = abs(qlen - slen)) %>%
          filter(abs_diff_len < 100) %>% 
          group_by(qseqid) %>%
          filter(sseqsp == test_sp2) %>%
          filter(bitscore == max(bitscore)) %>%
          slice_max(bitscore, n = 1, with_ties = FALSE) %>%
          ungroup() 
      ) %>%
      dplyr::select(qseqid, sseqid, bitscore, sseqsp, sseqorder)
      
    distrib_HGTindex_df <- 
      full_join(bh_bitscore_sameorder_df, bh_bitscore_all_notsameorder_df,
                by="qseqid") %>%
      filter(! is.na(sseqsp)) %>%
      filter(! is.na(ingroup_sseqsp))
    
    
    distrib_HGTindex_df <- distrib_HGTindex_df %>%
      mutate(HGTindex = bitscore/ingroup_bitscore) %>%
      filter(! is.na(HGTindex))
    
      
    HGT_index_values <- distrib_HGTindex_df %>% pull(HGTindex)
    index_quantile_5perc <- quantile(HGT_index_values, 0.995)
    index_quantile_1perc <- quantile(HGT_index_values, 0.999)
    
    curr_df <- as.data.frame(cbind(test_sp1, test_sp2, index_quantile_1perc, index_quantile_5perc))
    colnames(curr_df) <- c("ingroup_sseqsp", "outgroup_sseqsp", "HGTindex_quantile01", "HGTindex_quantile05")
      
    all_pairwise_quantile_index <- rbind(all_pairwise_quantile_index, curr_df)
    
  }
    
  bh_bitscore_df_q <- 
    left_join(bh_bitscore_df, all_pairwise_quantile_index, 
              by=c("ingroup_sseqsp", "outgroup_sseqsp")) 
    
  #Remove candidates below the quantile 5% 
  
  bh_bitscore_df_q_filt <- 
    bh_bitscore_df_q %>% filter(HGTindex > 0.85) %>% filter(HGTindex > HGTindex_quantile05) 
  
  write.table(bh_bitscore_df_q_filt,
              args[4],
              col.names=TRUE,
              row.names=FALSE,
              quote=FALSE,
              sep=",")
  
  
}




