##### Libraries  ---------------------------------

rm(list=ls())

set.seed(2712)


library("caper")
library("patchwork")
library("ape")
library(dplyr)
library(ggplot2)
library(tidyverse)
library("lattice")
library(reshape2)
library(adephylo)
library(phylobase)
library(data.table)
library(phytools)
library("corrplot")
library("geiger")
library("ggpubr")
library(RColorBrewer)
library(gt)
library(scales)
library(phangorn)
library(castor)
library(VennDiagram)
library(venn)
library(ggtree)
split_tibble <- function(tibble, column = 'col') {
  tibble %>% split(., .[,column]) %>% lapply(., function(x) x[,setdiff(names(x),column)])
}




#### My functions + Colors palettes  ---------------------------------

args = commandArgs(trailingOnly=TRUE)


PGLS_pvalue <- function(pgls_rslt) {
  sum_cor <- summary(pgls_rslt)
  pvalue =  formatC(sum_cor$coefficients[8], digits = 3)
  if (pvalue == "   0"){ pvalue = 2.2e-16}
  return(pvalue)
}

PGLS_R2 <- function(pgls_rslt) {
  sum_cor <- summary(pgls_rslt)
  r2 = formatC(sum_cor$r.squared, digits = 2)
  return(r2)
}

PGLS_lambda <- function(pgls_rslt) {
  sum_cor <- summary(pgls_rslt)
  PGLS_lambda = sum_cor$param[2]
  return(PGLS_lambda)
}


GLS_function <- function(gls_rslt) {
  GLS_cc <- coef(gls_rslt)
  GLS_f <- function(x) GLS_cc[1] + GLS_cc[2]*x
  return(GLS_f)
}





#### Load species taxonomy  ---------------------------------


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



#Verify that our orders are monophyletic

species_tree <- read.tree("AMAS_concatenated_alignment_BUSCO.fa.timetree.nwk")
species_tree$tip.label[species_tree$tip.label == "Salvelinus_sp_IW2-2015"] <- "Salvelinus_sp_IW2_2015"



species_test <- "Danio_rerio"
misc_table <- as.data.frame(species_test) %>% mutate(value = 1)
colnames(misc_table) <- c("node","value")
misc_table$node <- as.numeric(misc_table$node)
node_label_corresp <- 
  left_join(species_tree, misc_table, by = 'node') %>%
  dplyr::select(node, label)



order_list <- 
  species_taxonomy %>%
  filter(species %in% species_tree$tip.label) %>%
  pull(order) %>%
  unique()

Order_mono_df <- as.data.frame(NULL)
for(curr_order in order_list){
  
  curr_sp_list <- 
    species_taxonomy %>%
    filter(species %in% species_tree$tip.label) %>%
    filter(order == curr_order) %>%
    pull(species)
  
  if(length(curr_sp_list) > 1){
    
    curr_MRCA <- getMRCA(species_tree, curr_sp_list)
    curr_desc <- Descendants(species_tree, curr_MRCA)[[1]]
    curr_desc_label <- node_label_corresp %>% filter(node %in% curr_desc) %>% pull(label)
    curr_desc_orders <- species_taxonomy %>% filter(species %in% curr_desc_label) %>% pull(order) %>% unique()
    
    if(length(curr_desc_label) > length(curr_sp_list)) {
      curr_state="paraphyletic"
      
    } else {
      
      curr_state="monophyletic"
      
    }
    
  } else {
    
    curr_state="monophyletic"
    
  }
  
  
  curr_df <- 
    as.data.frame(cbind(curr_order, curr_state))
  colnames(curr_df) <- c("order", "state")
  
  Order_mono_df <- rbind(Order_mono_df, curr_df)
  
}


Order_mono_df %>%
  group_by(state) %>%
  dplyr::summarise(n())


Order_mono_df %>%
  filter(state == "paraphyletic")

## Annotated genomes

species_annotated <- 
  c("Lateolabrax_maculatus",
    "Amia_calva",
    "Argentina_silus",
    "Centroberyx_gerrardi",
    "Aplochiton_taeniatus",
    "Galaxias_olidus",
    "Datnioides_undecimradiatus",
    "Antennarius_maculatus",
    "Lutjanus_argentimaculatus",
    "Electrona_antarctica",
    "Parambassis_ranga",
    "Borostomias_antarcticus",
    "Hoplostethus_mediterraneus",
    "Hoplostethus_atlanticus",
    "Hyperoplus_immaculatus",
    "Ammodytes_marinus",
    "Zeus_faber")

orders_annotated <- 
  species_taxonomy %>%
  filter(species %in% species_annotated) %>%
  pull(order) %>% unique()

species_taxonomy %>%
  filter(species %in% species_tree$tip.label) %>%
  filter(order %in% orders_annotated) %>%
  filter(! species %in% species_annotated)

species_taxonomy %>%
  filter(species %in% species_annotated) %>%
  filter(order == "Centrarchiformes") 



#### Load Genome assembly BUSCO table  ---------------------------------

Genome_busco_df <- 
  read.table("BUSCO.results.csv",
             sep=",",
             header=FALSE)
colnames(Genome_busco_df) <- 
  c("species", "Complete", "Complete_single", "Complete_duplicated", "Fragmented",
    "Missing", "Total", "Scaffold_nb", "Contig_nb", "Genome_size_bp", "gaps_perc", "Scaff_N50_bp",
    "Contig_N50_bp")

Genome_busco_df <- 
  Genome_busco_df %>%
  mutate(BUSCO_perc = Complete/Total) %>%
  mutate(Genome_size_Gb = Genome_size_bp/1000000000) %>%
  mutate(Scaff_N50_Mb = Scaff_N50_bp/1000000) %>%
  mutate(N50_on_Gs = Scaff_N50_Mb/Genome_size_Gb)

#Combine with taxonomy

species_table_busco <- 
  left_join(Genome_busco_df ,species_taxonomy, by="species")

#### Load Annotation BUSCO table  ---------------------------------

Annot_busco_df <- 
  read.table("BUSCO_annotations.results.csv",
             sep=",",
             header=FALSE)
colnames(Annot_busco_df) <- 
  c("species", "Complete", "Complete_single", "Complete_duplicated", "Fragmented",
    "Missing", "Total")

Annot_busco_df <- 
  Annot_busco_df %>%
  mutate(BUSCO_perc = Complete/Total)

annotated_species <- 
  Annot_busco_df %>%
  pull(species)
  
#### Erosion curve BUSCO   ---------------------------------


BUSCO_erosion_df <- c()
for (i in seq(0, 1, 0.01)){
  current_df <- 
    species_table_busco %>% 
    dplyr::filter(BUSCO_perc >= i) %>% 
    summarise(number = n()) %>%
    mutate(Busco_percentage = i)
  BUSCO_erosion_df <- rbind(BUSCO_erosion_df, current_df)
}



#Define the number and species filtered 

pdf("Erosion_BUSCO.pdf",width = 6.34,  height = 4.61)

BUSCO_erosion_df %>%
  ggplot(., aes(x=Busco_percentage, y=number)) +
  geom_point() +
  theme_classic() +
  geom_line() +
  ylab("Number of genomes") +
  xlab("Proportion of complete BUSCO genes")+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position="none")

dev.off()



BUSCO_80_nb_species <- 
  BUSCO_erosion_df %>% 
  filter(Busco_percentage == 0.8) %>%
  pull(number)


#Define species who pass the 80% BUSCO filter
species_genome_80_busco <- 
  species_table_busco %>%
  filter(BUSCO_perc >= 0.8) %>%
  pull(species)





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
  

venn.diagram(
  x = list(putative_genes_TE_diamond_repbase, putative_genes_TE_blast_repbase, 
           putative_genes_TE_diamond_replib, putative_genes_TE_blast_replib,
           putative_genes_TE_pfam),
  category.names = c("diamond_repbase" , "blast_repbase" , "diamond_replib",
                     "blast_replib", "pfam"),
  filename = 'TE_venn_diagram.png',
  output=TRUE,
  height = 4000,
  width = 4000
)



putative_TE_genes <- 
  unique(c(putative_genes_TE_diamond_repbase,putative_genes_TE_blast_repbase,
           putative_genes_TE_diamond_replib, putative_genes_TE_blast_replib,
           putative_genes_TE_pfam))

putative_OGG_TE_pfam_blast <- 
  scan("N5_OGG_Transposable_elements_candidates_PFAM_BLAST_20perc.diamond_blast_pfam.txt", what="character") 


#### Compare ks statistics OGG vs BUSCO   ---------------------------------


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

OGG_BUSCO_Ks_stats_df <- 
  OGG_BUSCO_Ks_stats_df %>%
  mutate(group_ks_median = if_else(
    median_busco_001 < median_OGG,
    "busco_lower", 
    "busco_higher"
  )) %>%
  mutate(group_ks_mean = if_else(
    mean_busco_001 < mean_OGG,
    "busco_lower", 
    "busco_higher"
  )) %>%
  mutate(group_ks_q001 = if_else(
    quantile_busco_001 < quantile_OGG_001,
    "busco_lower", 
    "busco_higher"
  )) %>%
  mutate(group_ks_q002 = if_else(
    quantile_busco_002 < quantile_OGG_002,
    "busco_lower", 
    "busco_higher"
  )) %>%
  mutate(group_ks_q003 = if_else(
    quantile_busco_003 < quantile_OGG_003,
    "busco_lower", 
    "busco_higher"
  )) %>%
  mutate(group_ks_q004 = if_else(
    quantile_busco_004 < quantile_OGG_004,
    "busco_lower", 
    "busco_higher"
  )) %>%
  mutate(group_ks_q005 = if_else(
    quantile_busco_005 < quantile_OGG_005,
    "busco_lower", 
    "busco_higher"
  ))



#Lets test the correlation ks stats 

#Medians
cor_test_result <- 
  cor.test(
    OGG_BUSCO_Ks_stats_df %>% pull(median_busco_001),
    OGG_BUSCO_Ks_stats_df %>% pull(median_OGG),
    method="pearson"
  )
estimate_p <- round(cor_test_result$estimate, 3)
pvalue_p <- cor_test_result$p.value
pvalue_p <- format.pval(pvalue_p, digits = 2)

OGG_BUSCO_Ks_stats_df %>%
  ggplot(., aes(x=median_busco_001, y=median_OGG, color=group_ks_median)) +
  geom_point() +
  theme_classic() +
  geom_smooth(method = "lm", color = "black") +
  geom_abline(slope=1, linetype="dashed", color = "gray") +
  scale_color_manual(values = c("busco_lower" = "#B18F2B", "busco_higher" = "#004D40")) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none") +
  xlab("BUSCO - median ks") +
  ylab("OGG - median ks") +
  annotate(x=0.7, y=1.8, 
           label=paste("R = ", estimate_p, "; P = ", pvalue_p), 
           geom="text", size=5) 



#Means

cor_test_result <- 
  cor.test(
    OGG_BUSCO_Ks_stats_df %>% pull(mean_busco_001),
    OGG_BUSCO_Ks_stats_df %>% pull(mean_OGG),
    method="pearson"
  )
estimate_p <- round(cor_test_result$estimate, 3)
pvalue_p <- cor_test_result$p.value
pvalue_p <- format.pval(pvalue_p, digits = 2)

OGG_BUSCO_Ks_stats_df %>%
  ggplot(., aes(x=mean_busco_001, y=mean_OGG, color=group_ks_mean)) +
  geom_point() +
  theme_classic() +
  geom_smooth(method = "lm", color = "black") +
  geom_abline(slope=1, linetype="dashed", color = "gray") +
  scale_color_manual(values = c("busco_lower" = "#B18F2B", "busco_higher" = "#004D40")) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none") +
  xlab("BUSCO - mean ks") +
  ylab("OGG - mean ks") +
  annotate(x=0.7, y=1.8, 
           label=paste("R = ", estimate_p, "; P = ", pvalue_p), 
           geom="text", size=5) 



#quantile 0.1%

cor_test_result <- 
  cor.test(
    OGG_BUSCO_Ks_stats_df %>% pull(quantile_busco_001),
    OGG_BUSCO_Ks_stats_df %>% pull(quantile_OGG_001),
    method="pearson"
  )
estimate_p <- round(cor_test_result$estimate, 3)
pvalue_p <- cor_test_result$p.value
pvalue_p <- format.pval(pvalue_p, digits = 2)

OGG_BUSCO_Ks_stats_df %>%
  ggplot(., aes(x=quantile_busco_001, y=quantile_OGG_001, color=group_ks_q001)) +
  geom_point() +
  theme_classic() +
  geom_smooth(method = "lm", color = "black") +
  geom_abline(slope=1, linetype="dashed", color = "gray") +
  scale_color_manual(values = c("busco_lower" = "#B18F2B", "busco_higher" = "#004D40")) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none") +
  xlab("BUSCO - quantile 0.1%") +
  ylab("OGG -  quantile 0.1%") +
  annotate(x=0.25, y=0.7, 
           label=paste("R = ", estimate_p, "; P = ", pvalue_p), 
           geom="text", size=5) 


#quantile 0.2%

cor_test_result <- 
  cor.test(
    OGG_BUSCO_Ks_stats_df %>% pull(quantile_busco_002),
    OGG_BUSCO_Ks_stats_df %>% pull(quantile_OGG_002),
    method="pearson"
  )
estimate_p <- round(cor_test_result$estimate, 3)
pvalue_p <- cor_test_result$p.value
pvalue_p <- format.pval(pvalue_p, digits = 2)

OGG_BUSCO_Ks_stats_df %>%
  ggplot(., aes(x=quantile_busco_002, y=quantile_OGG_002, color=group_ks_q002)) +
  geom_point() +
  theme_classic() +
  geom_smooth(method = "lm", color = "black") +
  geom_abline(slope=1, linetype="dashed", color = "gray") +
  scale_color_manual(values = c("busco_lower" = "#B18F2B", "busco_higher" = "#004D40")) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none") +
  xlab("BUSCO - quantile 0.2%") +
  ylab("OGG -  quantile 0.2%") +
  annotate(x=0.25, y=0.7, 
           label=paste("R = ", estimate_p, "; P = ", pvalue_p), 
           geom="text", size=5) 


#quantile 0.3%

cor_test_result <- 
  cor.test(
    OGG_BUSCO_Ks_stats_df %>% pull(quantile_busco_003),
    OGG_BUSCO_Ks_stats_df %>% pull(quantile_OGG_003),
    method="pearson"
  )
estimate_p <- round(cor_test_result$estimate, 3)
pvalue_p <- cor_test_result$p.value
pvalue_p <- format.pval(pvalue_p, digits = 2)

OGG_BUSCO_Ks_stats_df %>%
  ggplot(., aes(x=quantile_busco_003, y=quantile_OGG_003, color=group_ks_q003)) +
  geom_point() +
  theme_classic() +
  geom_smooth(method = "lm", color = "black") +
  geom_abline(slope=1, linetype="dashed", color = "gray") +
  scale_color_manual(values = c("busco_lower" = "#B18F2B", "busco_higher" = "#004D40")) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none") +
  xlab("BUSCO - quantile 0.3%") +
  ylab("OGG -  quantile 0.3%") +
  annotate(x=0.25, y=0.7, 
           label=paste("R = ", estimate_p, "; P = ", pvalue_p), 
           geom="text", size=5) 


#quantile 0.4%

cor_test_result <- 
  cor.test(
    OGG_BUSCO_Ks_stats_df %>% pull(quantile_busco_004),
    OGG_BUSCO_Ks_stats_df %>% pull(quantile_OGG_004),
    method="pearson"
  )
estimate_p <- round(cor_test_result$estimate, 3)
pvalue_p <- cor_test_result$p.value
pvalue_p <- format.pval(pvalue_p, digits = 2)

OGG_BUSCO_Ks_stats_df %>%
  ggplot(., aes(x=quantile_busco_004, y=quantile_OGG_004, color=group_ks_q004)) +
  geom_point() +
  theme_classic() +
  geom_smooth(method = "lm", color = "black") +
  geom_abline(slope=1, linetype="dashed", color = "gray") +
  scale_color_manual(values = c("busco_lower" = "#B18F2B", "busco_higher" = "#004D40")) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none") +
  xlab("BUSCO - quantile 0.4%") +
  ylab("OGG -  quantile 0.4%") +
  annotate(x=0.25, y=0.7, 
           label=paste("R = ", estimate_p, "; P = ", pvalue_p), 
           geom="text", size=5) 



#quantile 0.5%

cor_test_result <- 
  cor.test(
    OGG_BUSCO_Ks_stats_df %>% pull(quantile_busco_005),
    OGG_BUSCO_Ks_stats_df %>% pull(quantile_OGG_005),
    method="pearson"
  )
estimate_p <- round(cor_test_result$estimate, 3)
pvalue_p <- cor_test_result$p.value
pvalue_p <- format.pval(pvalue_p, digits = 2)

OGG_BUSCO_Ks_stats_df %>%
  ggplot(., aes(x=quantile_busco_005, y=quantile_OGG_005, color=group_ks_q005)) +
  geom_point() +
  theme_classic() +
  geom_smooth(method = "lm", color = "black") +
  geom_abline(slope=1, linetype="dashed", color = "gray") +
  scale_color_manual(values = c("busco_lower" = "#B18F2B", "busco_higher" = "#004D40")) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none") +
  xlab("BUSCO - quantile 0.5%") +
  ylab("OGG -  quantile 0.5%") +
  annotate(x=0.25, y=0.7, 
           label=paste("R = ", estimate_p, "; P = ", pvalue_p), 
           geom="text", size=5) 






#### Correlate ks statistics with divergence times -- OGG  ---------------------------------

list_teleost_species <- 
  scan("teleost_annotated_species.txt", what="character")
list_teleost_species <- gsub("Salvelinus_sp_IW2-2015", "Salvelinus_sp_IW2_2015", list_teleost_species)

#Import the teleost species tree
species_tree <- read.tree("AMAS_concatenated_alignment_BUSCO.fa.timetree.nwk")
species_tree$tip.label[species_tree$tip.label == "Salvelinus_sp_IW2-2015"] <- "Salvelinus_sp_IW2_2015"

teleost_tree <- keep.tip(species_tree, list_teleost_species)
teleost_tree$edge.length <- teleost_tree$edge.length * 1000


#Lets add divergence times to the table
OGG_BUSCO_Ks_stats_df <- 
  as.data.frame(OGG_BUSCO_Ks_stats_df %>%
                  rowwise() %>%
                  mutate(div_time = 
                           get_pairwise_distances(teleost_tree, 
                                                  species_1, 
                                                  species_2, 
                                                  as_edge_counts=FALSE, 
                                                  check_input=TRUE)/2))



#Medians
cor_test_result <- 
  cor.test(
    OGG_BUSCO_Ks_stats_df %>% pull(div_time),
    OGG_BUSCO_Ks_stats_df %>% pull(median_OGG),
    method="pearson"
  )
estimate_p <- round(cor_test_result$estimate, 3)
pvalue_p <- cor_test_result$p.value
pvalue_p <- format.pval(pvalue_p, digits = 2)

OGG_BUSCO_Ks_stats_df %>%
  ggplot(., aes(x=div_time, y=median_OGG)) +
  geom_point() +
  theme_classic() +
  geom_smooth(method = "lm", color = "black") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none") +
  xlab("Divergence time (Myr)") +
  ylab("OGG - median ks") +
  annotate(x=25, y=1.8, 
           label=paste("R = ", estimate_p, "; P = ", pvalue_p), 
           geom="text", size=5) 


#mean
cor_test_result <- 
  cor.test(
    OGG_BUSCO_Ks_stats_df %>% pull(div_time),
    OGG_BUSCO_Ks_stats_df %>% pull(mean_OGG),
    method="pearson"
  )
estimate_p <- round(cor_test_result$estimate, 3)
pvalue_p <- cor_test_result$p.value
pvalue_p <- format.pval(pvalue_p, digits = 2)

OGG_BUSCO_Ks_stats_df %>%
  ggplot(., aes(x=div_time, y=mean_OGG)) +
  geom_point() +
  theme_classic() +
  geom_smooth(method = "lm", color = "black") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none") +
  xlab("Divergence time (Myr)") +
  ylab("OGG - mean ks") #+
  #annotate(x=25, y=1.8, 
   #        label=paste("R = ", estimate_p, "; P = ", pvalue_p), 
   #        geom="text", size=5) 

#ks 0.1%
cor_test_result <- 
  cor.test(
    OGG_BUSCO_Ks_stats_df %>% pull(div_time),
    OGG_BUSCO_Ks_stats_df %>% pull(quantile_OGG_001),
    method="pearson"
  )
estimate_p <- round(cor_test_result$estimate, 3)
pvalue_p <- cor_test_result$p.value
pvalue_p <- format.pval(pvalue_p, digits = 2)

OGG_BUSCO_Ks_stats_df %>%
  ggplot(., aes(x=div_time, y=quantile_OGG_001)) +
  geom_point() +
  theme_classic() +
  geom_smooth(method = "lm", color = "black") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none") +
  xlab("Divergence time (Myr)") +
  ylab("OGG - quantile 0.1%") +
  annotate(x=25, y=1.8, 
           label=paste("R = ", estimate_p, "; P = ", pvalue_p), 
           geom="text", size=5) 

#ks 0.2%
cor_test_result <- 
  cor.test(
    OGG_BUSCO_Ks_stats_df %>% pull(div_time),
    OGG_BUSCO_Ks_stats_df %>% pull(quantile_OGG_002),
    method="pearson"
  )
estimate_p <- round(cor_test_result$estimate, 3)
pvalue_p <- cor_test_result$p.value
pvalue_p <- format.pval(pvalue_p, digits = 2)

OGG_BUSCO_Ks_stats_df %>%
  ggplot(., aes(x=div_time, y=quantile_OGG_002)) +
  geom_point() +
  theme_classic() +
  geom_smooth(method = "lm", color = "black") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none") +
  xlab("Divergence time (Myr)") +
  ylab("OGG - quantile 0.2%") +
  annotate(x=25, y=1.8, 
           label=paste("R = ", estimate_p, "; P = ", pvalue_p), 
           geom="text", size=5) 


#ks 0.3%
cor_test_result <- 
  cor.test(
    OGG_BUSCO_Ks_stats_df %>% pull(div_time),
    OGG_BUSCO_Ks_stats_df %>% pull(quantile_OGG_003),
    method="pearson"
  )
estimate_p <- round(cor_test_result$estimate, 3)
pvalue_p <- cor_test_result$p.value
pvalue_p <- format.pval(pvalue_p, digits = 2)

OGG_BUSCO_Ks_stats_df %>%
  ggplot(., aes(x=div_time, y=quantile_OGG_003)) +
  geom_point() +
  theme_classic() +
  geom_smooth(method = "lm", color = "black") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none") +
  xlab("Divergence time (Myr)") +
  ylab("OGG - quantile 0.3%") +
  annotate(x=25, y=1.8, 
           label=paste("R = ", estimate_p, "; P = ", pvalue_p), 
           geom="text", size=5) 


#ks 0.4%
cor_test_result <- 
  cor.test(
    OGG_BUSCO_Ks_stats_df %>% pull(div_time),
    OGG_BUSCO_Ks_stats_df %>% pull(quantile_OGG_004),
    method="pearson"
  )
estimate_p <- round(cor_test_result$estimate, 3)
pvalue_p <- cor_test_result$p.value
pvalue_p <- format.pval(pvalue_p, digits = 2)

OGG_BUSCO_Ks_stats_df %>%
  ggplot(., aes(x=div_time, y=quantile_OGG_004)) +
  geom_point() +
  theme_classic() +
  geom_smooth(method = "lm", color = "black") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none") +
  xlab("Divergence time (Myr)") +
  ylab("OGG - quantile 0.4%") +
  annotate(x=25, y=1.8, 
           label=paste("R = ", estimate_p, "; P = ", pvalue_p), 
           geom="text", size=5) 


#ks 0.5%
cor_test_result <- 
  cor.test(
    OGG_BUSCO_Ks_stats_df %>% pull(div_time),
    OGG_BUSCO_Ks_stats_df %>% pull(quantile_OGG_005),
    method="pearson"
  )
estimate_p <- round(cor_test_result$estimate, 3)
pvalue_p <- cor_test_result$p.value
pvalue_p <- format.pval(pvalue_p, digits = 2)

OGG_BUSCO_Ks_stats_df %>%
  ggplot(., aes(x=div_time, y=quantile_OGG_005)) +
  geom_point() +
  theme_classic() +
  geom_smooth(method = "lm", color = "black") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none") +
  xlab("Divergence time (Myr)") +
  ylab("OGG - quantile 0.5%") #+
  #annotate(x=25, y=1.8, 
    #       label=paste("R = ", estimate_p, "; P = ", pvalue_p), 
    #       geom="text", size=5) 


#### KS strategy - Import HGT candidate tables based on Ks   ---------------------------------

#Import teleost species 
list_teleost_species <- 
  scan("teleost_annotated_species.txt", what="character")
list_teleost_species <- gsub("Salvelinus_sp_IW2-2015", "Salvelinus_sp_IW2_2015", list_teleost_species)


#Import HGT table with all pairwise comparisons with a ks below the max values between
#ks 0.5% quantile based on HOG ks distributions or ks 0.5% quantile based on BUSCO distributions
#This allowed to reduce the size of this dataframe, otherwise it would be too large to load ... 

HGT_candidate_df <- 
  fread("HGT_candidates_dataframe.csv",
        header=TRUE,
        sep=",")

#Remove identitcal pairwise comparisons

all_comb <- apply(combn(list_teleost_species,2),2,paste,collapse=',')
pairs_df <- as.data.frame(t(sapply(all_comb, function(x) unlist(strsplit(x, ",")))))
row.names(pairs_df) <- NULL
colnames(pairs_df) <- c("species_1", "species_2")

HGT_candidate_df <- left_join(pairs_df, HGT_candidate_df, by=c("species_1", "species_2"))
HGT_candidate_df <- HGT_candidate_df %>% filter(! is.na(cds_identity))


#### KS strat - Add scaffold informations to HGT table   ---------------------------------

nb_gene_per_scaff_df <- 
  fread("ALL_species.inproteome.counts",
        header=FALSE,
        sep=",")
nb_gene_per_scaff_df <- as.data.frame(nb_gene_per_scaff_df)
colnames(nb_gene_per_scaff_df) <- c("species","scaffold", "gene_name", "same_scaff_nbgene")

scaffold_lengths <- 
  fread("all_scaff_length.tsv",
        header=FALSE,
        sep="\t")
scaffold_lengths <- as.data.frame(scaffold_lengths)
colnames(scaffold_lengths) <- c("species","scaffold", "scaffold_bp")

nb_gene_per_scaff_Length_df <- 
  left_join(nb_gene_per_scaff_df, scaffold_lengths, by=c("species","scaffold"))

nb_gene_per_scaff_Length_df <- 
  nb_gene_per_scaff_Length_df %>%
  dplyr::select(-species)

nb_gene_per_scaff_Length_df$gene_name <- gsub("Salvelinus_sp_IW2-2015", "Salvelinus_sp_IW2_2015", nb_gene_per_scaff_Length_df$gene_name)
nb_gene_per_scaff_Length_df$gene_name <- gsub("\\.", "_", nb_gene_per_scaff_Length_df$gene_name)
nb_gene_per_scaff_Length_df$gene_name <- gsub("-", "_", nb_gene_per_scaff_Length_df$gene_name)

colnames(nb_gene_per_scaff_Length_df) <- c("scaffold_seq1", "seq1", "same_scaff_nbgene_seq1", "scaffold_seq1_length")


#HGT_candidate_df <- 
#  left_join(HGT_candidate_df, nb_gene_per_scaff_Length_df, by="seq1")
HGT_candidate_df <- 
  left_join(HGT_candidate_df, nb_gene_per_scaff_Length_df, by="seq1")


colnames(nb_gene_per_scaff_Length_df) <- c("scaffold_seq2", "seq2", "same_scaff_nbgene_seq2", "scaffold_seq2_length")

#HGT_candidate_df <- 
#  left_join(HGT_candidate_df, nb_gene_per_scaff_Length_df, by="seq2")
HGT_candidate_df <- 
  left_join(HGT_candidate_df, nb_gene_per_scaff_Length_df, by="seq2")




#### KS strat - Add microsynteny results to HGT table  ---------------------------------

#Genes not present in MicroSynteny_scores_table.csv = no gene around, discard
microsynteny_df <- 
  fread("MicroSynteny_scores_table.csv",
        header=FALSE,
        sep=",")
microsynteny_df <- as.data.frame(microsynteny_df)

colnames(microsynteny_df) <- c("seq1", "seq2", "micro_score", "nb_nbr_seq1", "nb_nbr_seq2")

#add this to the HGT table


HGT_candidate_df <- 
  left_join(HGT_candidate_df, microsynteny_df, by=c("seq1", "seq2")) %>%
  filter(! is.na(micro_score))

#### KS strat - Filter the HGT table with ks, scaffold and microsynteny   ---------------------------------

nrow(HGT_candidate_df %>% filter(ks < quantile_value_005_ks_OGG))
length(HGT_candidate_df %>% filter(ks < quantile_value_005_ks_OGG) %>% pull(OGG) %>% unique())

nrow(HGT_candidate_df %>% filter(ks < quantile_value_005_ks_OGG) %>% filter(same_scaff_nbgene_seq1 > 1 & same_scaff_nbgene_seq2 > 1) %>% filter(micro_score == 0))
length(HGT_candidate_df %>% filter(ks < quantile_value_005_ks_OGG) %>% filter(same_scaff_nbgene_seq1 > 1 & same_scaff_nbgene_seq2 > 1) %>% filter(micro_score == 0) %>% pull(OGG) %>% unique())

#Filter the HGT table


HGT_candidate_df_filtered <- 
  HGT_candidate_df %>%
  filter(ks < quantile_value_005_ks_OGG) %>%
  filter(same_scaff_nbgene_seq1 > 1 & same_scaff_nbgene_seq2 > 1) %>%
  filter(micro_score == 0) 

HGT_candidate_df_filtered %>% pull(OGG) %>% unique()


write.table(HGT_candidate_df_filtered,
                 "HGT_candidates_dataframe.filtered.csv",
                 sep=",",
                 col.names=TRUE,
                 row.names=FALSE,
                 quote=FALSE)

write(HGT_candidate_df_filtered %>% pull(OGG) %>% unique(), 
      "good_OGG_list.txt")




#### BLASTN strat - Import and parse BLAST outliers identified   ---------------------------------

#Import outlier from blast searches
blast_outliers_df_A <- 
  fread("All_blast_outliers_table.A.csv",
             header=TRUE,
             sep=",")
blast_outliers_df_B <- 
  fread("All_blast_outliers_table.B.csv",
             header=TRUE,
             sep=",")
blast_outliers_df_C <- 
  fread("All_blast_outliers_table.C.csv",
             header=TRUE,
             sep=",")

colnames(blast_outliers_df_C) <- 
  c("qseqid", "ingroup_sseqid", "ingroup_bitscore", "ingroup_sseqsp", "ingroup_sseqorder",
    "sseqid", "outgroup_bitscore", "outgroup_sseqsp", "outgroup_sseqorder", "HGTindex", 
    "HGTindex_quantile01", "HGTindex_quantile05")

#Import the link between these outliers and their OGG name
OGG_qseq_A <- read.table("c1_genes.OGG.A.txt", header=FALSE, sep="\t")
OGG_qseq_B <- read.table("c1_genes.OGG.B.txt", header=FALSE, sep="\t")
OGG_qseq_C <- read.table("c1_genes.OGG.C.txt", header=FALSE, sep="\t")

colnames(OGG_qseq_A) <- c("qseqid", "OGG_qseq") 
colnames(OGG_qseq_B) <- c("qseqid", "OGG_qseq") 
colnames(OGG_qseq_C) <- c("qseqid", "OGG_qseq") 


OGG_sseq_A <- read.table("c2_genes.OGG.A.txt", header=FALSE, sep="\t")
OGG_sseq_B <- read.table("c2_genes.OGG.B.txt", header=FALSE, sep="\t")
OGG_sseq_C <- read.table("c2_genes.OGG.C.txt", header=FALSE, sep="\t")

colnames(OGG_sseq_A) <- c("sseqid", "OGG_sseq") 
colnames(OGG_sseq_B) <- c("sseqid", "OGG_sseq") 
colnames(OGG_sseq_C) <- c("sseqid", "OGG_sseq") 


#Add the OGG names to the tables

blast_outliers_df_A <- left_join(blast_outliers_df_A, OGG_qseq_A, by="qseqid")
blast_outliers_df_A <- left_join(blast_outliers_df_A, OGG_sseq_A, by="sseqid")

blast_outliers_df_B <- left_join(blast_outliers_df_B, OGG_qseq_B, by="qseqid")
blast_outliers_df_B <- left_join(blast_outliers_df_B, OGG_sseq_B, by="sseqid")

blast_outliers_df_C <- left_join(blast_outliers_df_C, OGG_qseq_C, by="qseqid")
blast_outliers_df_C <- left_join(blast_outliers_df_C, OGG_sseq_C, by="sseqid")

#Remove if qseq != sseq

blast_outliers_df_A <- 
  blast_outliers_df_A %>%
  filter(! is.na(OGG_qseq)) %>%
  filter(! is.na(OGG_sseq)) %>%
  filter(OGG_qseq == OGG_sseq) %>%
  mutate(OGG = OGG_qseq) %>%
  dplyr::select(- c(OGG_qseq, OGG_qseq))


blast_outliers_df_B <- 
  blast_outliers_df_B %>%
  filter(! is.na(OGG_qseq)) %>%
  filter(! is.na(OGG_sseq)) %>%
  filter(OGG_qseq == OGG_sseq) %>%
  mutate(OGG = OGG_qseq) %>%
  dplyr::select(- c(OGG_qseq, OGG_qseq))


blast_outliers_df_C <- 
  blast_outliers_df_C %>%
  filter(! is.na(OGG_qseq)) %>%
  filter(! is.na(OGG_sseq)) %>%
  filter(OGG_qseq == OGG_sseq) %>%
  mutate(OGG = OGG_qseq) %>%
  dplyr::select(- c(OGG_qseq, OGG_qseq))



#Count the number of retained pairwise comparisons and OGGs

nrow(blast_outliers_df_A)
nrow(blast_outliers_df_B)
nrow(blast_outliers_df_C)

length(blast_outliers_df_A %>% pull(OGG) %>% unique())
length(blast_outliers_df_B %>% pull(OGG) %>% unique())
length(blast_outliers_df_C %>% pull(OGG) %>% unique())

blast_outliers_OGGs_A <- blast_outliers_df_A %>% pull(OGG) %>% unique()
blast_outliers_OGGs_B <- blast_outliers_df_B %>% pull(OGG) %>% unique()
blast_outliers_OGGs_C <- blast_outliers_df_C %>% pull(OGG) %>% unique()


write.table(blast_outliers_df_A, "blast_outliers_df_filtered.A.csv", sep=",",
            col.names=TRUE, row.names=FALSE, quote=FALSE)
write.table(blast_outliers_df_B, "blast_outliers_df_filtered.B.csv", sep=",",
            col.names=TRUE, row.names=FALSE, quote=FALSE)
write.table(blast_outliers_df_C, "blast_outliers_df_filtered.C.csv", sep=",",
            col.names=TRUE, row.names=FALSE, quote=FALSE)

write(blast_outliers_OGGs_A, "blast_outliers_OGGs.A.txt")
write(blast_outliers_OGGs_B, "blast_outliers_OGGs.B.txt")
write(blast_outliers_OGGs_C, "blast_outliers_OGGs.C.txt")



#### BLASTN strat - Add microsynteny results and filter based on microsynt scores   ---------------------------------

microsynteny_blast_df <- 
  fread("MicroSynteny_scores_table.blast.csv",
        header=FALSE,
        sep=",")

microsynteny_blast_df <- as.data.frame(microsynteny_blast_df)

colnames(microsynteny_blast_df) <- c("qseqid", "sseqid", "micro_score", "nb_nbr_seq1", "nb_nbr_seq2")

#add this to the blast table - first order

blast_outliers_df_A_f <- 
  left_join(blast_outliers_df_A, microsynteny_blast_df, by=c("qseqid", "sseqid")) %>% distinct() %>%
  filter(! is.na(micro_score))

blast_outliers_df_B_f <- 
  left_join(blast_outliers_df_B, microsynteny_blast_df, by=c("qseqid", "sseqid")) %>% distinct() %>%
  filter(! is.na(micro_score))

blast_outliers_df_C_f <- 
  left_join(blast_outliers_df_C, microsynteny_blast_df, by=c("qseqid", "sseqid")) %>% distinct() %>%
  filter(! is.na(micro_score))


#add this to the blast table - second order

colnames(microsynteny_blast_df) <- c("sseqid", "qseqid", "micro_score", "nb_nbr_seq1", "nb_nbr_seq2")

blast_outliers_df_A_r <- 
  left_join(blast_outliers_df_A, microsynteny_blast_df, by=c("qseqid", "sseqid")) %>% distinct() %>%
  filter(! is.na(micro_score))

blast_outliers_df_B_r <- 
  left_join(blast_outliers_df_B, microsynteny_blast_df, by=c("qseqid", "sseqid")) %>% distinct() %>%
  filter(! is.na(micro_score))

blast_outliers_df_C_r <- 
  left_join(blast_outliers_df_C, microsynteny_blast_df, by=c("qseqid", "sseqid")) %>% distinct() %>%
  filter(! is.na(micro_score))



#Merge tables

blast_outliers_df_A_a <- rbind(blast_outliers_df_A_f, blast_outliers_df_A_r)
blast_outliers_df_B_a <- rbind(blast_outliers_df_B_f, blast_outliers_df_B_r)
blast_outliers_df_C_a <- rbind(blast_outliers_df_C_f, blast_outliers_df_C_r)

#Lets filter if there are conserved genes around candidate HGTs + not the only gene on the scaffold

blast_outliers_df_A_micro <- as.data.frame(blast_outliers_df_A_a) %>% filter(micro_score == 0) %>% filter(nb_nbr_seq1 > 1 & nb_nbr_seq2 > 1)
blast_outliers_df_B_micro <- as.data.frame(blast_outliers_df_B_a) %>% filter(micro_score == 0) %>% filter(nb_nbr_seq1 > 1 & nb_nbr_seq2 > 1)
blast_outliers_df_C_micro <- as.data.frame(blast_outliers_df_C_a) %>% filter(micro_score == 0) %>% filter(nb_nbr_seq1 > 1 & nb_nbr_seq2 > 1)

nrow(blast_outliers_df_A_micro)
nrow(blast_outliers_df_B_micro)
nrow(blast_outliers_df_C_micro)

#Write retained BLAST candidates
unique_OGG_blast_A <- blast_outliers_df_A_micro %>% pull(OGG) %>% unique() 
unique_OGG_blast_B <- blast_outliers_df_B_micro %>% pull(OGG) %>% unique() 
unique_OGG_blast_C <- blast_outliers_df_C_micro %>% pull(OGG) %>% unique() 

length(unique_OGG_blast_A)
length(unique_OGG_blast_B)
length(unique_OGG_blast_C)


write(unique_OGG_blast_A, "BLASTN_stratA_HGT_OGG.txt")
write(unique_OGG_blast_B, "BLASTN_stratB_HGT_OGG.txt")
write(unique_OGG_blast_C, "BLASTN_stratC_HGT_OGG.txt")


write.table(blast_outliers_df_A_micro, "blast_outliers_df_A_micro.tsv",
            sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
write.table(blast_outliers_df_B_micro, "blast_outliers_df_B_micro.tsv",
            sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
write.table(blast_outliers_df_C_micro, "blast_outliers_df_C_micro.tsv",
            sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)


#Here ===> use bash script to add ka ks info

#### BLASTN strat - Add Ks information   ---------------------------------



blast_outliers_df_A_micro_ks <- 
  read.table("~/Horizontal_transfer_project/Best_BLAST_strategy/blast_outliers_df_A_micro.kaks.tsv", 
           sep="\t",
           header=TRUE)
blast_outliers_df_B_micro_ks <- 
  read.table("~/Horizontal_transfer_project/Best_BLAST_strategy/blast_outliers_df_B_micro.kaks.tsv", 
             sep="\t",
             header=TRUE)
blast_outliers_df_C_micro_ks <- 
  read.table("~/Horizontal_transfer_project/Best_BLAST_strategy/blast_outliers_df_C_micro.kaks.tsv", 
             sep="\t",
             header=TRUE)
blast_outliers_df_D_micro_ks <- 
  read.table("~/Horizontal_transfer_project/Best_BLAST_strategy/blast_outliers_df_D_micro.kaks.tsv", 
             sep="\t",
             header=TRUE)
blast_outliers_df_D_micro_ks <- blast_outliers_df_D_micro_ks %>% mutate(sseqsp = outgroup_sseqsp)

#Now add ks stats at the species level 

species_ks_df1 <- 
  OGG_BUSCO_Ks_stats_df %>% 
  dplyr::select(species_1, species_2,  quantile_OGG_001, quantile_OGG_002, quantile_OGG_003,
                quantile_OGG_004, quantile_OGG_005, median_OGG, mean_OGG)
species_ks_df2 <- 
  OGG_BUSCO_Ks_stats_df %>% 
  dplyr::select(species_2, species_1, quantile_OGG_001, quantile_OGG_002, quantile_OGG_003,
                quantile_OGG_004, quantile_OGG_005, median_OGG, mean_OGG)

colnames(species_ks_df1) <- c("qseqsp", "sseqsp", "quantile_OGG_001", "quantile_OGG_002", "quantile_OGG_003",
                              "quantile_OGG_004", "quantile_OGG_005", "median_OGG", "mean_OGG")
colnames(species_ks_df2) <- colnames(species_ks_df1)

species_ks_df <- rbind(species_ks_df2, species_ks_df1)


blast_outliers_df_A_micro_ks_info <-
  left_join(blast_outliers_df_A_micro_ks, species_ks_df, by=c("qseqsp", "sseqsp")) %>%
  filter(! is.na(quantile_OGG_001)) #some comparisons disapear if same order, as no ks computed among orders

blast_outliers_df_B_micro_ks_info <-
  left_join(blast_outliers_df_B_micro_ks, species_ks_df, by=c("qseqsp", "sseqsp")) %>%
  filter(! is.na(quantile_OGG_001)) 

blast_outliers_df_C_micro_ks_info <-
  left_join(blast_outliers_df_C_micro_ks, species_ks_df, by=c("qseqsp", "sseqsp")) %>%
  filter(! is.na(quantile_OGG_001)) 

blast_outliers_df_D_micro_ks_info <-
  left_join(blast_outliers_df_D_micro_ks, species_ks_df, by=c("qseqsp", "sseqsp")) %>%
  filter(! is.na(quantile_OGG_001)) 



blast_outliers_df_A_micro_ks_info %>%
  ggplot(., aes(x=ks, y=quantile_OGG_005)) +
  geom_point() +
  theme_minimal() +
  xlab("ks") +
  ylab("ks quantile 0.5%") +
  xlim(0, 1) +
  geom_abline(slope=1, linetype="dashed", color = "gray") 


blast_outliers_df_B_micro_ks_info %>%
  ggplot(., aes(x=ks, y=quantile_OGG_005)) +
  geom_point() +
  theme_minimal() +
  xlab("ks") +
  ylab("ks quantile 0.5%") +
  xlim(0, 1) +
  geom_abline(slope=1, linetype="dashed", color = "gray") 


blast_outliers_df_C_micro_ks_info %>%
  ggplot(., aes(x=ks, y=quantile_OGG_005)) +
  geom_point() +
  theme_minimal() +
  xlab("ks") +
  ylab("ks quantile 0.5%") +
  xlim(0, 1) +
  geom_abline(slope=1, linetype="dashed", color = "gray") 



blast_outliers_df_D_micro_ks_info %>%
  filter(OGG == "N5.HOG0004457")

#### Upset plot - Nb of OGG in the different strategies   ---------------------------------

OGG_ks <- HGT_candidate_df_filtered %>% pull(OGG) %>% unique()
OGG_blast_A <- unique_OGG_blast_A
OGG_blast_B <- unique_OGG_blast_B
OGG_blast_C <- unique_OGG_blast_D # the strategy D is called C in the manuscript

all_OGG <- unique(c(OGG_ks, OGG_blast_A, OGG_blast_B, OGG_blast_C))


library(UpSetR)

input_df <- c(
  Ks = length(OGG_ks),
  Blast_A = length(OGG_blast_A),
  Blast_B = length(OGG_blast_B),
  Blast_C = length(OGG_blast_C),
  "Ks&Blast_A" = length(intersect(OGG_ks, OGG_blast_A)),
  "Ks&Blast_B" = length(intersect(OGG_ks, OGG_blast_B)),
  "Ks&Blast_C" = length(intersect(OGG_ks, OGG_blast_C)),
  "Blast_A&Blast_B" = length(intersect(OGG_blast_A, OGG_blast_B)),
  "Blast_A&Blast_C" = length(intersect(OGG_blast_A, OGG_blast_C)),
  "Blast_B&Blast_C" = length(intersect(OGG_blast_B, OGG_blast_C)),
  "Ks&Blast_A&Blast_B&Blast_C" = 
    length(intersect(intersect(intersect(OGG_ks, OGG_blast_A), OGG_blast_B), OGG_blast_C)),
  "Blast_A&Blast_B&Blast_C" = 
    length(intersect(intersect(OGG_blast_A, OGG_blast_B), OGG_blast_C)),
  "Ks&Blast_A&Blast_B" = 
    length(intersect(intersect(OGG_ks, OGG_blast_A), OGG_blast_B)),
  "Ks&Blast_A&Blast_C" = 
    length(intersect(intersect(OGG_ks, OGG_blast_A), OGG_blast_C)),
  "Ks&Blast_B&Blast_C" = 
    length(intersect(intersect(OGG_ks, OGG_blast_B), OGG_blast_C))
)




upset(fromExpression(input_df), 
      nintersects = 15, 
      nsets = 5, 
      order.by = "freq", 
      decreasing = T, 
      mb.ratio = c(0.6, 0.4),
      number.angles = 0, 
      text.scale = 1.1, 
      point.size = 2.8, 
      line.size = 1,
      empty.intersections = "on"
)



#### Figure Manuscript - Species Tree   ---------------------------------



species_tree <- read.tree("AMAS_concatenated_alignment_BUSCO.fa.timetree.nwk")
species_tree$tip.label[species_tree$tip.label == "Salvelinus_sp_IW2-2015"] <- "Salvelinus_sp_IW2_2015"


species_taxonomy_intree <- 
  species_taxonomy %>% filter(species %in% species_tree$tip.label) %>%
  distinct()
all_orders <- species_taxonomy_intree$order %>% unique()

MRCA_orders_df <- as.data.frame(NULL)
for(curr_order in all_orders){
  curr_species <- species_taxonomy_intree %>% filter(order == curr_order) %>% pull(species)
  
  if(length(curr_species) > 1){
    order_MRCA <- findMRCA(species_tree, tips=curr_species, type="node")
    
    curr_df <- as.data.frame(cbind(curr_order, order_MRCA))
    colnames(curr_df) <- c("order", "MRCA")
    MRCA_orders_df <- rbind(MRCA_orders_df, curr_df)
  }
  
}

species_test <- "Danio_rerio"
misc_table <- as.data.frame(species_test) %>% mutate(value = 1)
colnames(misc_table) <- c("node","value")
misc_table$node <- as.numeric(misc_table$node)
node_label_corresp <- 
  left_join(species_tree, misc_table, by = 'node') %>%
  dplyr::select(node, label)



#orders_colors[names(orders_colors) == "Gobiiformes"]
node_label_corresp %>% filter(grepl("Amia_calva",label))

species_taxonomy_intree %>%
  filter(grepl("Acantho",species))

species_taxonomy_intree %>%
  filter(order == "Spariformes")

table_HOG_fig <- 
  read.table("Table_HOG_Figure.tsv",
             sep="\t",
             header=TRUE)

HOG_colors <- 
  c("N5.HOG0047987" = "#E69F00",
    "N5.HOG0001647" = "#000000",
    "N5.HOG0013514" = "#56B4E9",
    "N5.HOG0018901" = "#009E73",
    "N5.HOG0008885" = "#F0E442",
    "N5.HOG0004480" = "#0072B2",
    "N5.HOG0028899" = "#D55E00",
    "N5.HOG0046344" = "#CC79A7",
    "N5.HOG0004450" = "gray",
    "N5.HOG0029633" = "#E69F00",
    "N5.HOG0010622" = "#000000",
    "N5.HOG0013500" = "#56B4E9",
    "N5.HOG0019032" = "#009E73",
    "N5.HOG0004451" = "#F0E442",
    "N5.HOG0030200" = "#0072B2",
    "N5.HOG0034670" = "#D55E00",
    "N5.HOG0019111" = "#CC79A7")
  


ggtree(species_tree, layout="inward_circular", xlim=0.5, color="gray20", size=0.2) +
  geom_hilight(node=483, fill="gray", alpha=0.2) + #Polypteriformes
  geom_hilight(node=480, fill="gray20", alpha=0.2) + #Acipenseriformes
  geom_hilight(node=479, fill="gray20", alpha=0.2) + #Semionotiformes
  #geom_hilight(node="Amia_calva", fill="gray", alpha=0.2) + #Amia_calva
  geom_hilight(node=443, fill="#A6761D", alpha=0.2) + #Cypriniformes
  geom_hilight(node=312, fill="#A6CEE3", alpha=0.2) + #Tetraodontiformes
  geom_hilight(node=267, fill="#1F78B4", alpha=0.2) + #Cichliformes
  geom_hilight(node=277, fill="#B2DF8A", alpha=0.2) + #Cyprinodontiformes
  geom_hilight(node=361, fill="#33A02C", alpha=0.2) + #Pleuronectiformes
  geom_hilight(node=263, fill="#FB9A99", alpha=0.2) + #Ovalentaria incertae sedis
  geom_hilight(node=308, fill="#E31A1C", alpha=0.2) + #Eupercaria incertae sedis sec
  geom_hilight(node=360, fill="#FDBF6F", alpha=0.2) + #Carangaria incertae sedis
  geom_hilight(node=423, fill="#FF7F00", alpha=0.2) + #Siluriformes
  geom_hilight(node=371, fill="#CAB2D6", alpha=0.2) + #Syngnathiformes
  geom_hilight(node=401, fill="#B15928", alpha=0.2) + #Salmoniformes
  geom_hilight(node=293, fill="#FBB4AE", alpha=0.2) + #Beloniformes
  geom_hilight(node=356, fill="#B3CDE3", alpha=0.2) + #Carangiformes
  geom_hilight(node=476, fill="#CCEBC5", alpha=0.2) + #Osteoglossiformes
  geom_hilight(node=392, fill="#DECBE4", alpha=0.2) + #Gadiformes
  geom_hilight(node=324, fill="#FED9A6", alpha=0.2) + #Perciformes
  geom_hilight(node=310, fill="#E5D8BD", alpha=0.2) + #Eupercaria incertae sedis
  geom_hilight(node=320, fill="#FDDAEC", alpha=0.2) + #Labriformes
  geom_hilight(node=387, fill="black", alpha=0.2) + #Gobiiformes
  geom_hilight(node=413, fill="#B3E2CD", alpha=0.2) + #Esociformes
  geom_hilight(node=472, fill="#FDCDAC", alpha=0.2) + #Anguilliformes
  geom_hilight(node=442, fill="#CBD5E8", alpha=0.2) + #Gymnotiformes
  geom_hilight(node=475, fill="#F4CAE4", alpha=0.2) + #Elopiformes
  geom_hilight(node=314, fill="#7FC97F", alpha=0.2) + #Centrarchiformes
  geom_hilight(node=439, fill="#F1E2CC", alpha=0.2) + #Characiformes
  geom_hilight(node=415, fill="#E41A1C", alpha=0.2) + #Clupeiformes
  geom_hilight(node=396, fill="#386CB0", alpha=0.2) + #Galaxiiformes
  geom_hilight(node=296, fill="#377EB8", alpha=0.2) + #Atheriniformes
  geom_hilight(node=471, fill="#FFFF33", alpha=0.2) + #Albuliformes
  geom_hilight(node=398, fill="#F781BF", alpha=0.2) + #Osmeriformes
  geom_hilight(node=384, fill="#66C2A5", alpha=0.2) + #Scombriformes
  geom_hilight(node=352, fill="#E78AC3", alpha=0.2) + #Anabantiformes
  geom_hilight(node=390, fill="#E7298A", alpha=0.2) + #Trachichthyiformes
  geom_hilight(node=319, fill="#66A61E", alpha=0.2) + #Uranoscopiformes
  geom_taxalink(data=table_HOG_fig, 
                mapping=aes(taxa1=Donor_species, 
                            taxa2=Recipient_species, 
                            color=HOG,
                            hratio=Indice,
                            linetype=Line_type), 
                size=0.4,
                ncp=10,
                offset=0.02) +
  scale_color_manual(values = HOG_colors) +
  scale_linetype_manual(values = c("dotted","solid")) +
  geom_treescale(fontsize=2, linesize=0.3, x=1, y=1) +
  theme(legend.position="none")
  




#### Codon Usage - Evaluate codon usage bias per species -- ENCprime  ---------------------------------


library(seqinr)
library(coRdon)

list_teleost_sp <- list_teleost_species
list_teleost_sp[list_teleost_sp == "Salvelinus_sp_IW2_2015"] = "Salvelinus_sp_IW2-2015"

#all_ENC_df <- as.data.frame(NULL)
mean_ENC_df <- as.data.frame(NULL)
for(species in list_teleost_sp){
  curr_CDS_name <- 
    paste("~/Horizontal_transfer_project/Genomic_data/",
          species,
          "/",
          species,
          ".cds.nostop",
          sep="") ##Change the file name accordingly. Should point out to a fasta file with all coding sequences, without any STOP codons
  
    curr_CDS <- 
      readSet(
        file=curr_CDS_name
      )
  
    curr_codon_table <- codonTable(curr_CDS)
    
    curr_ENC <- 
      ENCprime(curr_codon_table, ribosomal = FALSE,
               stop.rm = TRUE, filtering = "hard", len.threshold = 100)
    mean_ENC <- mean(curr_ENC)
    
    curr_ENC_df <- 
      as.data.frame(curr_ENC) %>%
      mutate(species = species)
    
    #all_ENC_df <- rbind(all_ENC_df, curr_ENC_df)
    mean_ENC_df <- rbind(mean_ENC_df, as.data.frame(cbind(mean_ENC, species)))

}

mean_ENC_df$mean_ENC <- as.numeric(mean_ENC_df$mean_ENC)

pdf("ENCprime.pdf",width = 8.34,  height = 4.61)

mean_ENC_df %>%
  ggplot(., aes(x=1, y=mean_ENC)) +
  geom_violin() +
  geom_jitter(position=position_jitter(0.2)) +
  ylab("ENCp") +
  geom_hline(yintercept = 35, color="#D55E00") +
  geom_hline(yintercept = 61, color="#56B4E9") +
  geom_hline(yintercept = 20, color="#56B4E9") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none") +
  ylim(20, 61)

dev.off()

mean_ENC_df %>%
  pull(mean_ENC) %>%
  min()

mean_ENC_df %>%
  pull(mean_ENC) %>%
  max()

mean_ENC_df %>%
  pull(mean_ENC) %>%
  mean()

#### Codon Usage - Evaluate codon usage bias per species -- GC3  ---------------------------------

subset_teleost_sp <- 
  c("Clupea_harengus",
    "Paramormyrops_kingsleyae",
    "Megalops_cyprinoides",
    "Danio_rerio",
    "Hypomesus_transpacificus",
    "Pangasianodon_hypophthalmus",
    "Oreochromis_niloticus")


all_GC3_df <- as.data.frame(NULL)
for(species in subset_teleost_sp){
  
  curr_CDS_name <- 
    paste("~/Horizontal_transfer_project/Genomic_data/",
          species,
          "/",
          species,
          ".cds.nostop",
          sep="") ##Change the file name accordingly. Should point out to a fasta file with all coding sequences, without any STOP codons
  
  
  curr_CDS <- 
    seqinr::read.fasta(curr_CDS_name,
                       forceDNAtolower = TRUE)
  
  
  curr_GC3_df <- as.data.frame(NULL)
  
  for (i in 1:length(curr_CDS)){
    sequence <- curr_CDS[[i]]
    seq_name <- attr(sequence, "Annot")
    seq_name <- gsub(">", "", seq_name)
    seq_name <- gsub(" .*", "", seq_name)
    
    len_sequence <- length(sequence)
    
    if(len_sequence >= 300){
      curr_gc3 <- GC3(sequence, frame = 0)
      curr_df <- as.data.frame(cbind(seq_name, curr_gc3, species))
      curr_GC3_df <- rbind(curr_GC3_df, curr_df)
    }
  }
  
  all_GC3_df <- rbind(all_GC3_df, curr_GC3_df)
}

write.table(all_GC3_df, 
            "all_RSCU_df.csv",
            sep=",",
            col.names=TRUE,
            row.names=FALSE,
            quote=FALSE)

colnames(all_GC3_df) <- c("gene_name", "GC3", "species")
all_GC3_df$GC3 <- as.numeric(all_GC3_df$GC3 )
  


all_GC3_df %>%
  ggplot(., aes(x=GC3, color=species)) +
  geom_density(linewidth = 1) +
  scale_color_manual(
    values = c(
      "Clupea_harengus" = "#CC79A7",
      "Paramormyrops_kingsleyae" = "#D55E00",
      "Megalops_cyprinoides" = "#0072B2",
      "Danio_rerio" = "#F0E442",
      "Hypomesus_transpacificus" = "#009E73",
      "Pangasianodon_hypophthalmus" = "#56B4E9",
      "Oreochromis_niloticus" = "#E69F00"
    )
  ) +
  xlab("GC3") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position="none")


#### Codon Usage - Evaluate codon usage bias per species -- RSCU  ---------------------------------

subset_teleost_sp <- 
  c("Clupea_harengus",
    "Paramormyrops_kingsleyae",
    "Megalops_cyprinoides",
    "Danio_rerio",
    "Hypomesus_transpacificus",
    "Pangasianodon_hypophthalmus",
    "Oreochromis_niloticus")

all_RSCU_df <- as.data.frame(NULL)
for(species in subset_teleost_sp){
  
  print(species)
  
  curr_CDS_name <- 
    paste("~/Horizontal_transfer_project/Genomic_data/",
          species,
          "/",
          species,
          ".cds.nostop",
          sep="") ##Change the file name accordingly. Should point out to a fasta file with all coding sequences, without any STOP codons
  
  
  curr_CDS <- 
    seqinr::read.fasta(curr_CDS_name,
                       forceDNAtolower = TRUE)
  
  
  curr_RSCU_df <- as.data.frame(NULL)
  for (i in 1:length(curr_CDS)){
    sequence <- curr_CDS[[i]]
    seq_name <- attr(sequence, "Annot")
    len_sequence <- length(sequence)
    
    if(len_sequence >= 300){
      curr_uco <- uco(sequence, index = "rscu")
      curr_df <- data.frame(codon = names(curr_uco), rscu = as.vector(curr_uco))
      curr_df_wide <- as.data.frame(curr_df %>% pivot_wider(names_from = codon, values_from = rscu))
      curr_RSCU_df <- rbind(curr_RSCU_df, curr_df_wide)
    }
    
  }
  
  mean_RSCU_df <- 
    as.data.frame(colMeans(curr_RSCU_df, na.rm = TRUE)) %>%
    rownames_to_column(., var = "codon") %>%
    mutate(species = species) 
  colnames(mean_RSCU_df) <- c("codon", "RSCU", "species")
  
  column_sd_values <- apply(curr_RSCU_df, 2, sd, na.rm = TRUE)
  column_sd_values <- as.data.frame(column_sd_values, na.rm = TRUE) %>%
    rownames_to_column(., var = "codon") %>%
    mutate(species = species) 
  
  mean_RSCU_df <- left_join(mean_RSCU_df, column_sd_values, by=c("codon", "species"))
  colnames(mean_RSCU_df) <- c("codon", "RSCU", "species", "SD")
  
  mean_RSCU_df <- 
    mean_RSCU_df %>%
    mutate(N_gene = nrow(curr_RSCU_df))
  
  
  all_RSCU_df <- rbind(all_RSCU_df, mean_RSCU_df)
  
}
  

all_RSCU_df$SE <- all_RSCU_df$SD / sqrt(all_RSCU_df$N_gene)


curr_v_codon <- head(all_RSCU_df %>% pull(codon) %>% unique(), 10)

all_RSCU_df %>%
  filter(codon %in% curr_v_codon) %>%
  ggplot(., aes(x=codon, y=RSCU, fill=species)) +
  geom_bar(stat="identity", color="black", position=position_dodge(0.5), width = 0.5) + 
  geom_errorbar(
    aes(ymin= RSCU - SE, ymax= RSCU + SE), 
    width=0.2, 
    position=position_dodge(0.5)
  ) +
  scale_fill_manual(
    values = c(
      "Clupea_harengus" = "#CC79A7",
      "Paramormyrops_kingsleyae" = "#D55E00",
      "Megalops_cyprinoides" = "#0072B2",
      "Danio_rerio" = "#F0E442",
      "Hypomesus_transpacificus" = "#009E73",
      "Pangasianodon_hypophthalmus" = "#56B4E9",
      "Oreochromis_niloticus" = "#E69F00"
    )
  ) +
  ylab("RSCU") + 
  theme_classic() +
  geom_hline(yintercept = 1, color="black", linetype = "dashed") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position="none")


all_RSCU_df_wide <- 
  as.data.frame(
    all_RSCU_df %>%
  dplyr::select(codon, RSCU, species) %>%
  pivot_wider(names_from = species, values_from = RSCU)
  )



library(rstatix)

all_RSCU_df_wide_m <- as.matrix(all_RSCU_df_wide %>% dplyr::select(-codon))
rownames(all_RSCU_df_wide_m) <- all_RSCU_df_wide$codon

cor.rscu <- cor(all_RSCU_df_wide_m, use = "pairwise.complete.obs")

cor.rscu.melted <- melt(cor.rscu)


corrplot(cor.rscu, 
         is.corr = TRUE,
         method = "circle",
         type = "lower", 
         addCoef.col = "white",
         number.cex = 0.7, 
         tl.cex = 0.7)



#### Synteny block conservations statistics  ---------------------------------


weighted_median <- function(x, w) {
  x_sorted <- sort(x)
  w_sorted <- w[order(x)]
  
  cum_weights <- cumsum(w_sorted)
  total_weight <- sum(w_sorted)
  
  median_index <- which(cum_weights >= total_weight / 2)[1]
  return(x_sorted[median_index])
}


GLS_function <- function(gls_rslt) {
  GLS_cc <- coef(gls_rslt)
  GLS_f <- function(x) GLS_cc[1] + GLS_cc[2]*x
  return(GLS_f)
}



synteny_score_sp <- 
  fread("Summary_SyntenyBlocks_Pairwise.csv",
        sep=",",
        header=FALSE)
synteny_score_sp <- as.data.frame(synteny_score_sp)
colnames(synteny_score_sp) <- c("species1", "species2", "count", "nb_common")


nrow(synteny_score_sp %>%
  dplyr::select("species1", "species2") %>%
  distinct())

refseq_species <- 
  scan("RefSeq_species_list.txt",
       what="character")

#Compute the mean per species comparisons , the median, and the number of >0

summary_common_nb_df <- 
  synteny_score_sp %>%
  group_by(species1, species2) %>%
  summarise(weighted_mean = sum(count * nb_common) / sum(count),
            weighted_median = weighted_median(nb_common, count),
            proportion_above_0 = sum(count[nb_common > 0]) / sum(count)
            ) 
  

summary_common_nb_df <- 
  as.data.frame(summary_common_nb_df %>%
                  rowwise() %>%
                  mutate(div_time = 
                           get_pairwise_distances(teleost_tree, 
                                                  species1, 
                                                  species2, 
                                                  as_edge_counts=FALSE, 
                                                  check_input=TRUE)/2))


summary_common_nb_df$divergence_time <- summary_common_nb_df$div_time

summary_common_nb_df <- 
  summary_common_nb_df %>%
  mutate(Refseq_comp = if_else(
    (species1 %in% refseq_species) & (species2 %in% refseq_species),
    "refseq",
    "other"
  ))


summary_common_nb_df %>%
  filter(species1 == "Hypomesus_transpacificus") %>%
  filter(species2 == "Clupea_harengus")






#Mean versus divergence time

mean_vs_time_lm <- lm(weighted_mean ~ divergence_time, 
                    data = summary_common_nb_df)
summary_lm <- summary(mean_vs_time_lm)
slope_lm_ks <- summary_lm$coefficients[2]
inter_lm_ks <- summary_lm$coefficients[1]
curr_function_lm <- GLS_function(summary_lm)

mean_vs_time_lm_refseq <- lm(weighted_mean ~ divergence_time, 
                      data = summary_common_nb_df %>% filter(Refseq_comp == "refseq"))
summary_lm_refseq <- summary(mean_vs_time_lm_refseq)
slope_lm_ks_refseq <- summary_lm_refseq$coefficients[2]
inter_lm_ks_refseq <- summary_lm_refseq$coefficients[1]
curr_function_lm_refseq <- GLS_function(summary_lm_refseq)

mean_vs_time_lm_Norefseq <- lm(weighted_mean ~ divergence_time, 
                             data = summary_common_nb_df %>% filter(Refseq_comp != "refseq"))
summary_lm_Norefseq <- summary(mean_vs_time_lm_Norefseq)
slope_lm_ks_Norefseq <- summary_lm_Norefseq$coefficients[2]
inter_lm_ks_Norefseq <- summary_lm_Norefseq$coefficients[1]
curr_function_lm_Norefseq <- GLS_function(summary_lm_Norefseq)




summary_common_nb_df %>%
  ggplot(., aes(x=(divergence_time), y=weighted_mean, color=Refseq_comp)) +
  geom_point() + 
  scale_color_manual(values = c("other" = "#648FFF", "refseq" = "#E69F00")) + 
  stat_function(fun = curr_function_lm, color="black") +
  stat_function(fun = curr_function_lm_refseq, color="#E69F00") +
  stat_function(fun = curr_function_lm_Norefseq, color="#648FFF") +
  theme_classic() +
  xlab("Divergence time (Mya)") +
  ylab("Mean number of syntenic genes") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none")





synteny_df_clupea_hypomesus_summary <- 
  synteny_score_sp %>%
  filter(species1 == "Hypomesus_transpacificus") %>%
  filter(species2 == "Clupea_harengus") 



synteny_df_clupea_hypomesus_summary %>%
  ggplot(., aes(x=nb_common, y=count)) +
  geom_bar(stat="identity", fill="gray", color="black") +
  theme_classic() +
  xlab("Number of common HOG around pairs of genes belonging to the same HOG") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=15),
        plot.subtitle=element_text(size=16),
        legend.position="none")




#Multiple lm model

mean_vs_time_lm_simple <- 
  lm(weighted_mean ~ divergence_time,
     data = summary_common_nb_df)
summary(mean_vs_time_lm_simple)


mean_vs_time_lm_mult <- 
  lm(weighted_mean ~ divergence_time + Refseq_comp,
     data = summary_common_nb_df)
summary(mean_vs_time_lm_mult)

AIC(mean_vs_time_lm_simple, mean_vs_time_lm_mult)

#Median versus divergence time


median_vs_time_lm <- lm(weighted_median ~ divergence_time, 
                      data = summary_common_nb_df)


summary_lm <- summary(median_vs_time_lm)
slope_lm_ks <- summary_lm$coefficients[2]
inter_lm_ks <- summary_lm$coefficients[1]
curr_function_lm <- GLS_function(summary_lm)

summary_common_nb_df %>%
  ggplot(., aes(x=(divergence_time), y=weighted_median, color=Refseq_comp)) +
  geom_point() + 
  scale_color_manual(values = c("other" = "black", "refseq" = "#E69F00")) + 
  stat_function(fun = curr_function_lm, color="black") +
  theme_classic() +
  xlab("Divergence time (Mya)") +
  ylab("Median number of syntenic genes") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none")



#Median versus divergence time


median_vs_time_lm <- lm(proportion_above_0 ~ divergence_time, 
                        data = summary_common_nb_df)


summary_lm <- summary(median_vs_time_lm)
slope_lm_ks <- summary_lm$coefficients[2]
inter_lm_ks <- summary_lm$coefficients[1]
curr_function_lm <- GLS_function(summary_lm)

summary_common_nb_df %>%
  ggplot(., aes(x=(divergence_time), y=proportion_above_0, color=Refseq_comp)) +
  geom_point() + 
  scale_color_manual(values = c("other" = "black", "refseq" = "#E69F00")) + 
  stat_function(fun = curr_function_lm, color="black") +
  theme_classic() +
  xlab("Divergence time (Mya)") +
  ylab("Proportion of syntent blocks (>0)") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none")





#### Circos plot - Hypomesus vs Clupea   ---------------------------------

library(Biostrings)
library(circlize)

#Load genomes
Clupea_harengus_genome <-
  read.table("GCF_900700415.2_Ch_v2.0.2_genomic.fna.fai", ## genome to be found on NCBI. Fai faile generated with samtools faidx genome.fa
             sep="\t")
colnames(Clupea_harengus_genome) <- 
  c("scaffold", "length", "offset", "linebases", "linewidth")
Clupea_harengus_genome <- Clupea_harengus_genome %>% mutate(start = 1)

Hypomesus_transpacificus_genome <-
  read.table("GCF_021917145.1_fHypTra1_genomic.fna.fai", ## genome to be found on NCBI. Fai faile generated with samtools faidx genome.fa
             sep="\t")
colnames(Hypomesus_transpacificus_genome) <- 
  c("scaffold", "length", "offset", "linebases", "linewidth")
Hypomesus_transpacificus_genome <- Hypomesus_transpacificus_genome %>% mutate(start = 1)

#Keep only the chromosomes

Clupea_harengus_genome_chr <- Clupea_harengus_genome %>% filter(grepl("NC_", scaffold))
Clupea_harengus_genome_chr$chr_nb <- paste("chr", seq_len(nrow(Clupea_harengus_genome_chr)))

Hypomesus_transpacificus_genome_chr <- Hypomesus_transpacificus_genome %>% filter(grepl("NC_", scaffold))
Hypomesus_transpacificus_genome_chr$chr_nb <- paste("chr", seq_len(nrow(Hypomesus_transpacificus_genome_chr)))


#Combine the genomes in one dataframe

Combined_genomes_chr <- 
  rbind(Clupea_harengus_genome_chr %>% mutate(species = "Clupea_harengus"), 
        Hypomesus_transpacificus_genome_chr %>% mutate(species = "Hypomesus_transpacificus"))

Combined_genomes_chr$factor <- paste(Combined_genomes_chr$species, Combined_genomes_chr$scaffold, sep = "_")


#Make a dataframe of links between horizontally transfered genes

HGT_positions_df <- 
  read.table("/HGT_links_Ch_Ht.tsv",
             sep="\t",
             header=TRUE,
             quote = "\"",
             comment.char = "")

HGT_positions_df <- 
  HGT_positions_df %>% 
  mutate(sector.index1 = paste0("Clupea_harengus_", sector.index1))

HGT_positions_df <- 
  HGT_positions_df %>% 
  mutate(sector.index2 = paste0("Hypomesus_transpacificus_", sector.index2))



#Initialize circos

pdf(file = "circos_plot.pdf", width = 8, height = 8)


circos.clear()

circos.par(cell.padding = c(0.02, 0, 0.02, 0), gap.degree = 0)

# Initialize the circos plot for both species with correct chromosome lengths
circos.initialize(factors = Combined_genomes_chr$factor, xlim = cbind(0, Combined_genomes_chr$length))

# Draw the first layer : Genomes

circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr_name = CELL_META$sector.index
  chr_length = CELL_META$xlim[2]
  chr_number = CELL_META$chr_nb
  
  if (grepl("Hypomesus_transpacificus", chr_name)) {
    circos.rect(0, 0, chr_length, 1, col = "gray88", border = "black")
    circos.text(chr_length / 2, 0.5, gsub("Hypomesus_transpacificus_", "", chr_name), facing = "bending.inside", niceFacing = TRUE, col = "black", cex = 0.1)
    
  } else if (grepl("Clupea_harengus", chr_name)) {
    circos.rect(0, 0, chr_length, 1, col = "gray42", border = "black")
    circos.text(chr_length / 2, 0.5, gsub("Clupea_harengus_", "", chr_name), facing = "bending.inside", niceFacing = TRUE, col = "black", cex = 0.5, h="top", bg.border=F)
  }
}, bg.border = NA, track.height = 0.1)


#Add link

for (i in 1:nrow(HGT_positions_df)) {
  circos.link(
    sector.index1 = HGT_positions_df$sector.index1[i],
    point1 = HGT_positions_df$point1[i],
    sector.index2 = HGT_positions_df$sector.index2[i],
    point2 = HGT_positions_df$point2[i],
    col = HGT_positions_df$color[i],
    lwd = HGT_positions_df$nb_cluster[i]
  )
}


dev.off()



#### dN/dS analysis of good candidates   ---------------------------------

### Selection analaysis files in the folder "Selection_analysis"

#Import dN/dS tables computed with HyPhy fitMG4

N5.HOG0047987_omega_df <- read.table("~/Horizontal_transfer_project/Candidates_analysis/N5.HOG0047987/N5.HOG0047987.dN_dS.csv", sep=",", header=FALSE)
N5.HOG0001647_omega_df <- read.table("~/Horizontal_transfer_project/Candidates_analysis/N5.HOG0001647/N5.HOG0001647.dN_dS.csv", sep=",", header=FALSE)
N5.HOG0013514_omega_df <- read.table("~/Horizontal_transfer_project/Candidates_analysis/N5.HOG0013514/N5.HOG0013514.dN_dS.csv", sep=",", header=FALSE)
N5.HOG0018901_omega_df <- read.table("~/Horizontal_transfer_project/Candidates_analysis/N5.HOG0018901/N5.HOG0018901.dN_dS.csv", sep=",", header=FALSE)
N5.HOG0008885_omega_df <- read.table("~/Horizontal_transfer_project/Candidates_analysis/N5.HOG0008885/N5.HOG0008885.dN_dS.csv", sep=",", header=FALSE)
N5.HOG0004480_omega_df <- read.table("~/Horizontal_transfer_project/Candidates_analysis/N5.HOG0004480/N5.HOG0004480.dN_dS.csv", sep=",", header=FALSE)
N5.HOG0028899_omega_df <- read.table("~/Horizontal_transfer_project/Candidates_analysis/N5.HOG0028899/N5.HOG0028899.dN_dS.csv", sep=",", header=FALSE)
N5.HOG0004450_omega_df <- read.table("~/Horizontal_transfer_project/Candidates_analysis/N5.HOG0004450/N5.HOG0004450.dN_dS.csv", sep=",", header=FALSE)
N5.HOG0029633_omega_df <- read.table("~/Horizontal_transfer_project/Candidates_analysis/N5.HOG0029633/N5.HOG0029633.dN_dS.csv", sep=",", header=FALSE)
N5.HOG0046344_omega_df <- read.table("~/Horizontal_transfer_project/Candidates_analysis/N5.HOG0029633/N5.HOG0046344.dN_dS.csv", sep=",", header=FALSE)
N5.HOG0010622_omega_df <- read.table("~/Horizontal_transfer_project/Candidates_analysis/N5.HOG0010622/N5.HOG0010622.dN_dS.csv", sep=",", header=FALSE)
N5.HOG0013500_omega_df <- read.table("~/Horizontal_transfer_project/Candidates_analysis/N5.HOG0013500/N5.HOG0013500.dN_dS.csv", sep=",", header=FALSE)
N5.HOG0019032_omega_df <- read.table("~/Horizontal_transfer_project/Candidates_analysis/N5.HOG0019032/N5.HOG0019032.dN_dS.csv", sep=",", header=FALSE)
N5.HOG0004451_omega_df <- read.table("~/Horizontal_transfer_project/Candidates_analysis/N5.HOG0004451/N5.HOG0004451.dN_dS.csv", sep=",", header=FALSE)
N5.HOG0030200_omega_df <- read.table("~/Horizontal_transfer_project/Candidates_analysis/N5.HOG0030200/N5.HOG0030200.dN_dS.csv", sep=",", header=FALSE)
N5.HOG0034670_omega_df <- read.table("~/Horizontal_transfer_project/Candidates_analysis/N5.HOG0034670/N5.HOG0034670.dN_dS.csv", sep=",", header=FALSE)
N5.HOG0019111_omega_df <- read.table("~/Horizontal_transfer_project/Candidates_analysis/N5.HOG0019111/N5.HOG0019111.dN_dS.csv", sep=",", header=FALSE)

#Rename columns
colnames(N5.HOG0047987_omega_df) <- c("species_gene_name", "LB", "MLE", "UB", "dN", "dS")
colnames(N5.HOG0001647_omega_df) <- c("species_gene_name", "LB", "MLE", "UB", "dN", "dS")
colnames(N5.HOG0013514_omega_df) <- c("species_gene_name", "LB", "MLE", "UB", "dN", "dS")
colnames(N5.HOG0018901_omega_df) <- c("species_gene_name", "LB", "MLE", "UB", "dN", "dS")
colnames(N5.HOG0008885_omega_df) <- c("species_gene_name", "LB", "MLE", "UB", "dN", "dS")
colnames(N5.HOG0004480_omega_df) <- c("species_gene_name", "LB", "MLE", "UB", "dN", "dS")
colnames(N5.HOG0028899_omega_df) <- c("species_gene_name", "LB", "MLE", "UB", "dN", "dS")
colnames(N5.HOG0004450_omega_df) <- c("species_gene_name", "LB", "MLE", "UB", "dN", "dS")
colnames(N5.HOG0029633_omega_df) <- c("species_gene_name", "LB", "MLE", "UB", "dN", "dS")
colnames(N5.HOG0046344_omega_df) <- c("species_gene_name", "LB", "MLE", "UB", "dN", "dS")
colnames(N5.HOG0010622_omega_df) <- c("species_gene_name", "LB", "MLE", "UB", "dN", "dS")
colnames(N5.HOG0013500_omega_df) <- c("species_gene_name", "LB", "MLE", "UB", "dN", "dS")
colnames(N5.HOG0019032_omega_df) <- c("species_gene_name", "LB", "MLE", "UB", "dN", "dS")
colnames(N5.HOG0004451_omega_df) <- c("species_gene_name", "LB", "MLE", "UB", "dN", "dS")
colnames(N5.HOG0030200_omega_df) <- c("species_gene_name", "LB", "MLE", "UB", "dN", "dS")
colnames(N5.HOG0034670_omega_df) <- c("species_gene_name", "LB", "MLE", "UB", "dN", "dS")
colnames(N5.HOG0019111_omega_df) <- c("species_gene_name", "LB", "MLE", "UB", "dN", "dS")

#Remove branches with a saturations of syn/non-syn mutations

N5.HOG0047987_omega_df <- N5.HOG0047987_omega_df %>%  filter(dN < 1) %>% filter(dS < 1) %>% filter(dS > 0.01) %>% filter(MLE < 5) %>% mutate(HOG = "N5.HOG0047987")
N5.HOG0001647_omega_df <- N5.HOG0001647_omega_df %>%  filter(dN < 1) %>% filter(dS < 1) %>% filter(dS > 0.01) %>% filter(MLE < 5) %>% mutate(HOG = "N5.HOG0001647")
N5.HOG0013514_omega_df <- N5.HOG0013514_omega_df %>%  filter(dN < 1) %>% filter(dS < 1) %>% filter(dS > 0.01) %>% filter(MLE < 5) %>% mutate(HOG = "N5.HOG0013514")
N5.HOG0018901_omega_df <- N5.HOG0018901_omega_df %>%  filter(dN < 1) %>% filter(dS < 1) %>% filter(dS > 0.01) %>% filter(MLE < 5) %>% mutate(HOG = "N5.HOG0018901")
N5.HOG0008885_omega_df <- N5.HOG0008885_omega_df %>%  filter(dN < 1) %>% filter(dS < 1) %>% filter(dS > 0.01) %>% filter(MLE < 5) %>% mutate(HOG = "N5.HOG0008885")
N5.HOG0004480_omega_df <- N5.HOG0004480_omega_df %>%  filter(dN < 1) %>% filter(dS < 1) %>% filter(dS > 0.01) %>% filter(MLE < 5) %>% mutate(HOG = "N5.HOG0004480")
N5.HOG0028899_omega_df <- N5.HOG0028899_omega_df %>%  filter(dN < 1) %>% filter(dS < 1) %>% filter(dS > 0.01) %>% filter(MLE < 5) %>% mutate(HOG = "N5.HOG0028899")
N5.HOG0004450_omega_df <- N5.HOG0004450_omega_df %>%  filter(dN < 1) %>% filter(dS < 1) %>% filter(dS > 0.01) %>% filter(MLE < 5) %>% mutate(HOG = "N5.HOG0004450")
N5.HOG0029633_omega_df <- N5.HOG0029633_omega_df %>%  filter(dN < 1) %>% filter(dS < 1) %>% filter(dS > 0.01) %>% filter(MLE < 5) %>% mutate(HOG = "N5.HOG0029633")
N5.HOG0046344_omega_df <- N5.HOG0046344_omega_df %>%  filter(dN < 1) %>% filter(dS < 1) %>% filter(dS > 0.01) %>% filter(MLE < 5) %>% mutate(HOG = "N5.HOG0046344")
N5.HOG0010622_omega_df <- N5.HOG0010622_omega_df %>%  filter(dN < 1) %>% filter(dS < 1) %>% filter(dS > 0.01) %>% filter(MLE < 5) %>% mutate(HOG = "N5.HOG0010622")
N5.HOG0013500_omega_df <- N5.HOG0013500_omega_df %>%  filter(dN < 1) %>% filter(dS < 1) %>% filter(dS > 0.01) %>% filter(MLE < 5) %>% mutate(HOG = "N5.HOG0013500")
N5.HOG0019032_omega_df <- N5.HOG0019032_omega_df %>%  filter(dN < 1) %>% filter(dS < 1) %>% filter(dS > 0.01) %>% filter(MLE < 5) %>% mutate(HOG = "N5.HOG0019032")
N5.HOG0004451_omega_df <- N5.HOG0004451_omega_df %>%  filter(dN < 1) %>% filter(dS < 1) %>% filter(dS > 0.01) %>% filter(MLE < 5) %>% mutate(HOG = "N5.HOG0004451")
N5.HOG0030200_omega_df <- N5.HOG0030200_omega_df %>%  filter(dN < 1) %>% filter(dS < 1) %>% filter(dS > 0.01) %>% filter(MLE < 5) %>% mutate(HOG = "N5.HOG0030200")
N5.HOG0034670_omega_df <- N5.HOG0034670_omega_df %>%  filter(dN < 1) %>% filter(dS < 1) %>% filter(dS > 0.01) %>% filter(MLE < 5) %>% mutate(HOG = "N5.HOG0034670")
N5.HOG0019111_omega_df <- N5.HOG0019111_omega_df %>%  filter(dN < 1) %>% filter(dS < 1) %>% filter(dS > 0.01) %>% filter(MLE < 5) %>% mutate(HOG = "N5.HOG0019111")

#Combine all the tables

ALL_omega_df <- 
  rbind(N5.HOG0047987_omega_df,
    N5.HOG0001647_omega_df,
    N5.HOG0013514_omega_df,
    N5.HOG0018901_omega_df,
    N5.HOG0008885_omega_df,
    N5.HOG0004480_omega_df,
    N5.HOG0028899_omega_df,
    N5.HOG0004450_omega_df,
    N5.HOG0029633_omega_df,
    N5.HOG0046344_omega_df,
    N5.HOG0010622_omega_df,
    N5.HOG0013500_omega_df,
    N5.HOG0019032_omega_df,
    N5.HOG0004451_omega_df,
    N5.HOG0030200_omega_df,
    N5.HOG0034670_omega_df,
    N5.HOG0019111_omega_df
  )


#Mark recipient branches

ALL_omega_df <- 
  ALL_omega_df %>%
  mutate(Branch_type = case_when(
    HOG == "N5.HOG0047987" & 
      species_gene_name %in% 
      c("Node105",
        "Node107",
        "Hypomesus_transpacificus_rna_XM_047033889_1",
        "Hypomesus_transpacificus_rna_XM_047020147_1",
        "Osmerus_eperlanus_rna_XM_062485952_1") ~ "Recipient_1",
    
    HOG == "N5.HOG0001647" & 
      species_gene_name %in% 
      c("Hypomesus_transpacificus_rna_XM_047032478_1") ~ "Recipient_1",
    
    HOG == "N5.HOG0013514" & 
      species_gene_name %in% 
      c("Node360",
        "Hypomesus_transpacificus_rna_XM_047045767_1",
        "Hypomesus_transpacificus_rna_XM_047045786_1",
        "Osmerus_eperlanus_rna_XM_062468699_1",
        "Node371",
        "Osmerus_eperlanus_rna_XM_062468715_1",
        "Osmerus_eperlanus_rna_XM_062468701_1") ~ "Recipient_1",
    
    HOG == "N5.HOG0018901" & 
      species_gene_name %in% 
      c("Hypomesus_transpacificus_rna_XM_047023826_1",
        "Node729",
        "Hypomesus_transpacificus_rna_XM_047025851_1",
        "Hypomesus_transpacificus_rna_XM_047044257_1",
        "Hypomesus_transpacificus_rna_XM_047043933_1") ~ "Recipient_1",
    
    
    HOG == "N5.HOG0008885" & 
      species_gene_name %in% 
      c("Hypomesus_transpacificus_rna_XM_047048755_1") ~ "Recipient_1",
    
    
    HOG == "N5.HOG0004480" & 
      species_gene_name %in% 
      c("Osmerus_eperlanus_rna_XM_062453720_1") ~ "Recipient_1",
    
    HOG == "N5.HOG0028899" & 
      species_gene_name %in% 
      c("Hypomesus_transpacificus_rna_XM_047031800_1") ~ "Recipient_1",
    
    
    HOG == "N5.HOG0004450" & 
      species_gene_name %in% 
      c("Node381",
        "Node377",
        "Lates_calcarifer_rna_XM_018688510_2",
        "Lates_calcarifer_rna_XM_018686387_2",
        "Lates_japonicus_gene_AKAME5_002034700",
        "Seriola_lalandi_dorsalis_rna_XM_023399112_1",
        "Caranx_melampygus_mRNA27263",
        "Caranx_melampygus_mRNA11584",
        "Caranx_melampygus_mRNA25461") ~ "Recipient_1",
    
    
    HOG == "N5.HOG0029633" & 
      species_gene_name %in% 
      c("Paramormyrops_kingsleyae_rna_XM_023822687_1") ~ "Recipient_1",
    
    HOG == "N5.HOG0046344" & 
      species_gene_name %in% 
      c("Paramormyrops_kingsleyae_rna_XM_023822671_1") ~ "Recipient_1",
    
    HOG == "N5.HOG0010622" & 
      species_gene_name %in% 
      c("Paramormyrops_kingsleyae_rna_XM_023831520_1") ~ "Recipient_1",
    
    
    HOG == "N5.HOG0013500" & 
      species_gene_name %in% 
      c("Node21",
        "Thalassophryne_amazonica_rna_XM_034167571_1",
        "Thalassophryne_amazonica_rna_XM_034166331_1") ~ "Recipient_1",
    
    HOG == "N5.HOG0019032" & 
      species_gene_name %in% 
      c("Node17",
        "Node19",
        "Neoarius_graeffei_rna_XM_060938829_1",
        "Neoarius_graeffei_rna_XM_060938827_1",
        "Neoarius_graeffei_rna_XM_060939486_1") ~ "Recipient_1",
    
    
    HOG == "N5.HOG0004451" & 
      species_gene_name %in% 
      c("Centroberyx_gerrardi_g23848_t1") ~ "Recipient_1",
    
    
    HOG == "N5.HOG0030200" & 
      species_gene_name %in% 
      c("Node14",
        "Node15",
        "Pangasianodon_hypophthalmus_rna_XM_026945431_3",
        "Pangasius_djambal_PDJAM_T00103720",
        "Pangasius_djambal_PDJAM_T00012360") ~ "Recipient_1",
    
    
    HOG == "N5.HOG0034670" & 
      species_gene_name %in% 
      c("Pangasianodon_gigas_PGIGA_T00260940",
        "Pangasianodon_gigas_PGIGA_T00102290",
        "Node73",
        "Pangasianodon_hypophthalmus_rna_XM_053237209_1",
        "Pangasianodon_hypophthalmus_rna_XM_053236680_1",
        "Node80",
        "Pangasianodon_hypophthalmus_rna_XM_053236679_1",
        "Node79",
        "Pangasianodon_hypophthalmus_rna_XM_053237136_1",
        "Node78",
        "Pangasianodon_gigas_PGIGA_T00260840",
        "Pangasius_djambal_PDJAM_T00258610",
        "Node85",
        "Node77",
        "Pangasianodon_gigas_PGIGA_T00260930",
        "Node76",
        "Node72") ~ "Recipient_1",
    
    
    HOG == "N5.HOG0019111" & 
      species_gene_name %in% 
      c("Hypomesus_transpacificus_rna_XM_047049897_1") ~ "Recipient_1",
    
    HOG == "N5.HOG0019111" & 
      species_gene_name %in% 
      c("Clupea_harengus_rna_XM_031577667_2",
        "Node1109",
        "Clupea_harengus_rna_XM_031577666_1",
        "Node1108",
        "Clupea_harengus_rna_XM_031578285_2",
        "Node1107",
        "Alosa_alosa_rna_XM_048231371_1",
        "Alosa_sapidissima_rna_XM_042075309_1",
        "Node1118",
        "Sardina_pilchardus_rna_XM_062535392_1",
        "Node1117") ~ "Recipient_2"
    
  ))

ALL_omega_df <- 
  ALL_omega_df %>%
  mutate(Branch_type = ifelse(is.na(Branch_type), "Other", Branch_type))


#Mark branches with positive or relaxed selection detected

ALL_omega_df <- 
  ALL_omega_df %>%
  mutate(Selection_type = case_when(
    HOG == "N5.HOG0047987" & 
      species_gene_name %in% 
      c("Osmerus_eperlanus_rna_XM_062485952_1",
        "Hypomesus_transpacificus_rna_XM_047033889_1") ~ "Diversifying",
    
    HOG == "N5.HOG0001647" & 
      species_gene_name %in% 
      c("Hypomesus_transpacificus_rna_XM_047032478_1") ~ "Diversifying",
    
    HOG == "N5.HOG0013514" & 
      species_gene_name %in% 
      c("Hypomesus_transpacificus_rna_XM_047045786_1",
        "Node371") ~ "Diversifying",


    HOG == "N5.HOG0004480" & 
      species_gene_name %in% 
      c("Osmerus_eperlanus_rna_XM_062453720_1") ~ "Diversifying",
    
    HOG == "N5.HOG0028899" & 
      species_gene_name %in% 
      c("Hypomesus_transpacificus_rna_XM_047031800_1") ~ "Diversifying",
    
    
    HOG == "N5.HOG0004450" & 
      species_gene_name %in% 
      c("Node377") ~ "Diversifying",
    
    
    HOG == "N5.HOG0046344" & 
      species_gene_name %in% 
      c("Paramormyrops_kingsleyae_rna_XM_023822671_1") ~ "Diversifying",
    
    HOG == "N5.HOG0010622" & 
      species_gene_name %in% 
      c("Paramormyrops_kingsleyae_rna_XM_023831520_1") ~ "Relaxation",

  
    HOG == "N5.HOG0030200" & 
      species_gene_name %in% 
      c("Node14",
        "Node15",
        "Pangasius_djambal_PDJAM_T00103720") ~ "Diversifying",
    
    
    HOG == "N5.HOG0034670" & 
      species_gene_name %in% 
      c("Pangasianodon_gigas_PGIGA_T00102290",
        "Pangasianodon_gigas_PGIGA_T00260930") ~ "Diversifying",
    

    HOG == "N5.HOG0019111" & 
      species_gene_name %in% 
      c("Clupea_harengus_rna_XM_031577667_2") ~ "Diversifying"
    
  ))

ALL_omega_df <- 
  ALL_omega_df %>%
  mutate(Selection_type = ifelse(is.na(Selection_type), "Other", Selection_type))



#Let's plot :)

pdf("~/Horizontal_transfer_project/Raw_plots/Selection_violinPlots.pdf",width = 10,  height = 4.61)


ALL_omega_df %>%
  ggplot(., aes(x=HOG, y=MLE)) +
  geom_violin() +
  geom_jitter(position=position_jitter(0.2), 
              aes(color= Branch_type, size=Branch_type, shape=Selection_type, alpha=Branch_type)) +
  scale_color_manual(values = c("Other" = "#CEC7C7", "Donors" = "#CEC7C7", 
                                "Recipient_1" = "#0072B2", "Recipient_2" = "#FFB000")) +
  scale_alpha_manual(values = c("Other" = 0.4, "Donors" = 0.4, 
                                "Recipient_1" = 1, "Recipient_2" = 1)) +
  scale_size_manual(values = c("Other" = 1, "Donors" = 1, 
                               "Recipient_1" = 3, "Recipient_2" = 3)) +
  scale_shape_manual(values = c("Other" = 16, "Diversifying" = 17, "Relaxation" = 18)) +
  ylab("dN/dS") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none")



ALL_omega_df %>%
  ggplot(., aes(x=HOG, y=MLE)) +
  geom_violin() +
  geom_jitter(position=position_jitter(0.2), 
              aes(color= Branch_type, size=Branch_type, shape=Selection_type, alpha=Branch_type)) +
  scale_color_manual(values = c("Other" = "#CEC7C7", "Donors" = "#CEC7C7", 
                                "Recipient_1" = "#0072B2", "Recipient_2" = "#FFB000")) +
  scale_alpha_manual(values = c("Other" = 0.4, "Donors" = 0.4, 
                                "Recipient_1" = 1, "Recipient_2" = 1)) +
  scale_size_manual(values = c("Other" = 1, "Donors" = 1, 
                               "Recipient_1" = 3, "Recipient_2" = 3)) +
  scale_shape_manual(values = c("Other" = 16, "Diversifying" = 17, "Relaxation" = 18)) +
  ylab("dN/dS") +
  theme_classic() +
  theme(axis.text=element_text(size=5),
        axis.title=element_text(size=5),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none")

ALL_omega_df %>%
  filter(HOG %in% c(
    "N5.HOG0047987",
    "N5.HOG0001647",
    "N5.HOG0013514",
    "N5.HOG0018901",
    "N5.HOG0008885",
    "N5.HOG0004480",
    "N5.HOG0028899",
    "N5.HOG0004450",
    "N5.HOG0029633")) %>%
  ggplot(., aes(x=HOG, y=MLE)) +
  geom_violin() +
  geom_jitter(position=position_jitter(0.2), 
              aes(color= Branch_type, size=Branch_type, shape=Selection_type, alpha=Branch_type)) +
  scale_color_manual(values = c("Other" = "#CEC7C7", "Donors" = "#CEC7C7", 
                                "Recipient_1" = "#0072B2", "Recipient_2" = "#FFB000")) +
  scale_alpha_manual(values = c("Other" = 0.4, "Donors" = 0.4, 
                                "Recipient_1" = 1, "Recipient_2" = 1)) +
  scale_size_manual(values = c("Other" = 1, "Donors" = 1, 
                               "Recipient_1" = 3, "Recipient_2" = 3)) +
  scale_shape_manual(values = c("Other" = 16, "Diversifying" = 17, "Relaxation" = 18)) +
  ylab("dN/dS") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none")


ALL_omega_df %>%
  filter(! HOG %in% c(
    "N5.HOG0047987",
    "N5.HOG0001647",
    "N5.HOG0013514",
    "N5.HOG0018901",
    "N5.HOG0008885",
    "N5.HOG0004480",
    "N5.HOG0028899",
    "N5.HOG0004450",
    "N5.HOG0029633")) %>%
  ggplot(., aes(x=HOG, y=MLE)) +
  geom_violin() +
  geom_jitter(position=position_jitter(0.2), 
              aes(color= Branch_type, size=Branch_type, shape=Selection_type, alpha=Branch_type)) +
  scale_color_manual(values = c("Other" = "#CEC7C7", "Donors" = "#CEC7C7", 
                                "Recipient_1" = "#0072B2", "Recipient_2" = "#FFB000")) +
  scale_alpha_manual(values = c("Other" = 0.4, "Donors" = 0.4, 
                                "Recipient_1" = 1, "Recipient_2" = 1)) +
  scale_size_manual(values = c("Other" = 1, "Donors" = 1, 
                               "Recipient_1" = 3, "Recipient_2" = 3)) +
  scale_shape_manual(values = c("Other" = 16, "Diversifying" = 17, "Relaxation" = 18)) +
  ylab("dN/dS") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none")

dev.off()

ALL_omega_df %>%
  filter(Branch_type %in% c("Recipient_1", "Recipient_2")) %>%
  group_by(Selection_type) %>%
  summarise(count = n())



