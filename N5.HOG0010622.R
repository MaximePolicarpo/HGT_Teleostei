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
library(phylolm)
library(ggtreeExtra)
library(ggtree)
library(treeio)
library(tidytree)
library(RColorBrewer)
library(ggnewscale)
library(gggenomes)

##### Import HGT data  ---------------------------------

args = commandArgs(trailingOnly=TRUE)

HGT_candidate_df <- 
  fread("HGT_candidates_dataframe.filtered.csv",
        sep=",",
        header=TRUE)

HGT_candidate_df <- 
  HGT_candidate_df %>%
  filter(ks < quantile_value_005_ks_OGG) %>%
  filter(same_scaff_nbgene_seq1 > 1 & same_scaff_nbgene_seq2 > 1) %>%
  filter(micro_score == 0)

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



putative_TE_genes <- 
  unique(c(putative_genes_TE_diamond_repbase,putative_genes_TE_blast_repbase,
           putative_genes_TE_diamond_replib, putative_genes_TE_blast_replib,
           putative_genes_TE_pfam))

putative_OGG_TE_pfam_blast <- 
  scan("N5_OGG_Transposable_elements_candidates_PFAM_BLAST_20perc.diamond_blast_pfam.txt", what="character") 


##### Import species tree ---------------------------------

#import teleost species list
list_teleost_sp <- scan("teleost_annotated_species.txt",
                        what="character")

#Import the teleost species tree
species_tree <- read.tree("AMAS_concatenated_alignment_BUSCO.fa.timetree.nwk")
teleost_tree <- keep.tip(species_tree, list_teleost_sp)



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

species_order <-  species_taxonomy %>% dplyr::select(species, order)


#### Define order colors  ---------------------------------

list_order <- species_order %>% pull(order) %>% unique()

#Define a color palette
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector <- unique(col_vector)
col_vector <- col_vector[ !col_vector == '#FFFF99']
col_vector <- col_vector[ !col_vector == '#FFF2AE']
col_vector <- col_vector[ !col_vector == '#FFFFCC']

orders_colors <- c()
i=1
for(curr_orders in list_order){
  
  orders_colors <- c(orders_colors, col_vector[i])
  i=i+1
  
}

list_order <- c(list_order, "noHGT")
orders_colors <- c(orders_colors, "black")
names(orders_colors) <- list_order



#### Draw  phylogenies with alignment -- OGG   ---------------------------------

curr_OGG <- "N5.HOG0010622"


current_tree <- 
  read.tree("N5.HOG0010622.prot.aln.treefile")

#Root tree at midpoint 
current_tree_rooted <- midpoint.root(current_tree)

#Associate each gene name with its species and with its order

tips_df <- as.data.frame(current_tree_rooted$tip.label)
colnames(tips_df) <- c("label")


list_teleost_sp[list_teleost_sp == "Salvelinus_sp_IW2-2015"] <- "Salvelinus_sp_IW2_2015"
tips_df_species <- as.data.frame(NULL)
for(curr_sp in list_teleost_sp){
  
  
  curr_df <- 
    as.data.frame(
      tips_df %>%
        rowwise() %>%
        filter(grepl(curr_sp, label))
    ) %>%
    mutate(species = curr_sp)
  
  tips_df_species <-
    rbind(tips_df_species, curr_df)
}

#add the order to the dataframe
tips_df_species_order <- left_join(tips_df_species, species_order, by="species")


#Define wich gene are part of candidate HGT events, and with which ks filter ... 

curr_OGG_table <- 
  HGT_candidate_df %>%
  filter(OGG == curr_OGG) %>%
  filter(ks < quantile_value_005_ks_OGG) %>%
  mutate(ks_pass = case_when(
    ks < quantile_value_005_ks_OGG & ks > quantile_value_004_ks_OGG ~ "ks_05",
    ks < quantile_value_004_ks_OGG & ks > quantile_value_003_ks_OGG ~ "ks_04",
    ks < quantile_value_003_ks_OGG & ks > quantile_value_002_ks_OGG ~ "ks_03",
    ks < quantile_value_002_ks_OGG & ks > quantile_value_001_ks_OGG ~ "ks_02",
    ks < quantile_value_001_ks_OGG ~ "ks_01"
  ))



tips_HGT_ks5 <- c(
  curr_OGG_table %>% filter(ks_pass == "ks_05") %>% pull(seq1),
  curr_OGG_table %>% filter(ks_pass == "ks_05") %>% pull(seq2))
tips_HGT_ks4 <- c(
  curr_OGG_table %>% filter(ks_pass == "ks_04") %>% pull(seq1),
  curr_OGG_table %>% filter(ks_pass == "ks_04") %>% pull(seq2))
tips_HGT_ks3 <- c(
  curr_OGG_table %>% filter(ks_pass == "ks_03") %>% pull(seq1),
  curr_OGG_table %>% filter(ks_pass == "ks_03") %>% pull(seq2))
tips_HGT_ks2 <- c(
  curr_OGG_table %>% filter(ks_pass == "ks_02") %>% pull(seq1),
  curr_OGG_table %>% filter(ks_pass == "ks_02") %>% pull(seq2))
tips_HGT_ks1 <- c(
  curr_OGG_table %>% filter(ks_pass == "ks_01") %>% pull(seq1),
  curr_OGG_table %>% filter(ks_pass == "ks_01") %>% pull(seq2))
all_tips_HTG <- c(tips_HGT_ks5, tips_HGT_ks4, tips_HGT_ks3, tips_HGT_ks2, tips_HGT_ks1)


tips_df_species_order <- 
  tips_df_species_order %>%
  mutate(ks_color_cat = case_when(
    label %in% tips_HGT_ks5 ~ "ks_05",
    label %in% tips_HGT_ks4 ~ "ks_04",
    label %in% tips_HGT_ks3 ~ "ks_03",
    label %in% tips_HGT_ks2 ~ "ks_02",
    label %in% tips_HGT_ks1 ~ "ks_01",
    ! label %in% all_tips_HTG ~ "noHGT"
  )) %>%
  mutate(HGT_or_not = case_when(
    label %in% all_tips_HTG ~ "Yes",
    ! label %in% all_tips_HTG ~ "No"
  ))



# Plot the rooted gene trees

ks_colors <- 
  c("ks_01" = "#009E73",
    "ks_02" = "#56B4E9",
    "ks_03" = "#E69F00",
    "ks_04" = "#CC79A7",
    "ks_05" = "#000000",
    "noHGT" = "white")


#Add the alignment
prot_alignment = "N5.HOG0010622.prot.aln"
prot_alignment_trim = "N5.HOG0010622.prot.trimmed.aln"


p0 <- 
  ggtree(current_tree_rooted) %<+% tips_df_species_order +
  geom_tiplab(size = 3) +
  aes(color=as.character(order), alpha=as.character(HGT_or_not)) +
  scale_color_manual(values = orders_colors) + 
  scale_alpha_manual(values = c("No" = 0.5, "Yes" = 1)) + 
  new_scale_color() + 
  geom_tippoint(aes(color=ks_color_cat), size=0.5) +
  scale_color_manual(values = ks_colors) + 
  new_scale_color() + 
  theme(legend.position = "none") 
p1 <- msaplot(p0, prot_alignment, offset=1, width=3)
p2 <- msaplot(p0, prot_alignment_trim, offset=1, width=3)




#Save trees as PDF


ggsave(paste("N5.HOG0010622/",
             curr_OGG,
             ".gene.prot_alignment.pdf",
             sep=""), 
       p1, width=15.34, height=24.34, 
       units="in", scale=2)

ggsave(paste("N5.HOG0010622/",
             curr_OGG,
             ".gene.prot_alignment_trim.pdf",
             sep=""), 
       p2, width=15.34, height=24.34, 
       units="in", scale=2)


### Now draw a species phylogeny 

#Define which species have at-least one gene in the curr OGG

species_OGG <- tips_df_species_order %>% pull(species) %>% unique()

presence_species_df <- as.data.frame(list_teleost_sp)
colnames(presence_species_df) <- c("species")

presence_species_df <- 
  presence_species_df %>%
  rowwise() %>%
  mutate(presence_absence = if_else(
    species %in% species_OGG,
    "present",
    "absent"
  ))


species_order <-  species_taxonomy %>% dplyr::select(species, order)
species_taxonomy[(species_taxonomy$species == "Salvelinus_sp_IW2-2015"),"species"] <- "Salvelinus_sp_IW2_2015"
species_order[(species_order$species == "Salvelinus_sp_IW2-2015"),"species"] <- "Salvelinus_sp_IW2_2015"


species_order <- left_join(presence_species_df, species_order, by="species")

#Define candidate transfers species, and at which ks .. 

curr_OGG_table_sp <-
  curr_OGG_table %>%
  dplyr::select(species_1, species_2, ks_pass) %>%
  distinct()

#Lets draw a tree with link between candidate HGT species

teleost_tree$tip.label[teleost_tree$tip.label == "Salvelinus_sp_IW2-2015"] <- "Salvelinus_sp_IW2_2015"

p2 <-
  ggtree(teleost_tree, layout="inward_circular", xlim=0.5) %<+% species_order +
  aes(color=as.character(order), alpha=as.character(presence_absence)) +
  scale_color_manual(values = orders_colors, guide = "none") +
  scale_alpha_manual(values = c("absent" = 0.1, "present" = 1), guide = "none") +
  new_scale_color() +
  geom_tippoint(aes(color=presence_absence), size=1) +
  scale_color_manual(values = c("presence_absence" = "white", "present" = "black")) +
  new_scale_color() +
  geom_taxalink(data=curr_OGG_table_sp, 
                mapping=aes(taxa1=species_1, 
                            taxa2=species_2, 
                            color=ks_pass), 
                size=1,
                ncp=10,
                offset=0.02) +
  scale_color_manual(values = ks_colors) + 
  ggtitle(curr_OGG) 

#Save trees as PDF


ggsave(paste("N5.HOG0010622/",
             curr_OGG,
             ".species.pdf",
             sep=""), 
       p2, width=8, height=6, 
       units="in", scale=2)

#### Draw  phylogenies with alignment -- Blastp   ---------------------------------


current_tree <- 
  read.tree("Uniprot_plus_closefish.aln.treefile")

#Root tree at midpoint 
current_tree_rooted <- midpoint.root(current_tree)

#Associate each gene name with its species and with its order

tips_df <- as.data.frame(current_tree_rooted$tip.label)
colnames(tips_df) <- c("label")


list_teleost_sp[list_teleost_sp == "Salvelinus_sp_IW2-2015"] <- "Salvelinus_sp_IW2_2015"
tips_df_species <- as.data.frame(NULL)
for(curr_sp in list_teleost_sp){
  
  
  curr_df <- 
    as.data.frame(
      tips_df %>%
        rowwise() %>%
        filter(grepl(curr_sp, label))
    ) %>%
    mutate(species = curr_sp)
  
  tips_df_species <-
    rbind(tips_df_species, curr_df)
}

#add the order to the dataframe
tips_df_species_order <- left_join(tips_df_species, species_order, by="species")



#Add the alignment
prot_alignment = "Uniprot_plus_closefish.aln"
prot_alignment_trim= "Uniprot_plus_closefish.aln.trimmed"

p0 <- 
  ggtree(current_tree_rooted) %<+% tips_df_species_order +
  geom_tiplab(size = 2) +
  aes(color=as.character(order)) +
  scale_color_manual(values = orders_colors) + 
  new_scale_color() + 
  new_scale_color() + 
  theme(legend.position = "none")

p2 <- msaplot(p0, prot_alignment, offset=2, width=3)
p3 <- msaplot(p0, prot_alignment_trim, offset=2, width=3)




#Save trees as PDF

ggsave(paste("N5.HOG0010622/",
             curr_OGG,
             ".uniprot.prot_alignment.pdf",
             sep=""), 
       p2, width=9, height=12, 
       units="in", scale=2)

ggsave(paste("N5.HOG0010622/",
             curr_OGG,
             ".uniprot.prot_alignment_trim.pdf",
             sep=""), 
       p3, width=9, height=12, 
       units="in", scale=2)



#### DNA subtree  ---------------------------------


current_tree <- 
  read.tree("N5.HOG0010622/HGT_clade.cds.aln.treefile")

current_tree_rooted <- midpoint.root(current_tree)

tips_df <- as.data.frame(current_tree_rooted$tip.label)
colnames(tips_df) <- c("label")


list_teleost_sp[list_teleost_sp == "Salvelinus_sp_IW2-2015"] <- "Salvelinus_sp_IW2_2015"
tips_df_species <- as.data.frame(NULL)
for(curr_sp in list_teleost_sp){
  
  
  curr_df <- 
    as.data.frame(
      tips_df %>%
        rowwise() %>%
        filter(grepl(curr_sp, label))
    ) %>%
    mutate(species = curr_sp)
  
  tips_df_species <-
    rbind(tips_df_species, curr_df)
}

#add the order to the dataframe
tips_df_species_order <- left_join(tips_df_species, species_order, by="species")


DNA_tree <- 
  ggtree(current_tree_rooted, size=1) %<+% tips_df_species_order +
  geom_tiplab(size = 5) +
  aes(color=as.character(order)) +
  scale_color_manual(values = orders_colors) +
  new_scale_color() + 
  theme(legend.position = "none")  +
  geom_text2(aes(subset = !isTip, label = label), hjust = 1.2, vjust = 1.6,size=2) +
  xlim(0, 0.6) +
  geom_treescale(x=0.02, y=12, width=0.05, color='black') 


DNA_tree





#### Draw  phylogenies with alignment -- Alternative assembly + other paramormyrops   ---------------------------------


current_tree <- 
  read.tree("OGG_plus_ManuallyFound.aln.treefile")

#Root tree at midpoint 
current_tree_rooted <- midpoint.root(current_tree)

#Associate each gene name with its species and with its order

tips_df <- as.data.frame(current_tree_rooted$tip.label)
colnames(tips_df) <- c("label")


list_teleost_sp[list_teleost_sp == "Salvelinus_sp_IW2-2015"] <- "Salvelinus_sp_IW2_2015"
tips_df_species <- as.data.frame(NULL)
for(curr_sp in list_teleost_sp){
  
  
  curr_df <- 
    as.data.frame(
      tips_df %>%
        rowwise() %>%
        filter(grepl(curr_sp, label))
    ) %>%
    mutate(species = curr_sp)
  
  tips_df_species <-
    rbind(tips_df_species, curr_df)
}

#add the order to the dataframe
tips_df_species_order <- left_join(tips_df_species, species_order, by="species")

toadd <- as.data.frame(cbind("Paramormyrops_hopkinsi", "Paramormyrops_kingsleyae", "present", "Osteoglossiformes"))
colnames(toadd) <- colnames(tips_df_species_order)
tips_df_species_order <- rbind(tips_df_species_order, toadd)

#Add the alignment
prot_alignment = "OGG_plus_ManuallyFound.aln"

p0 <- 
  ggtree(current_tree_rooted) %<+% tips_df_species_order +
  geom_tiplab(size = 2) +
  aes(color=as.character(order)) +
  scale_color_manual(values = orders_colors) + 
  new_scale_color() + 
  new_scale_color() + 
  theme(legend.position = "none")

p1 <- msaplot(p0, prot_alignment, offset=1, width=3)


pdf("N5.HOG0010622/Manually_annoated_genes.pdf",width = 8.34,  height = 4.61)

p1

dev.off()



#### Species Ks distribution -- Paramormyrops_kingsleyae vs Ictalurus_punctatus   ---------------------------------


current_ks_df <- 
  fread("HGT_stats_Results_per_sp/Paramormyrops_kingsleyae.cds.stats",
        header=FALSE,
        sep=",")

curr_stats_df <-
  fread("KS_Stats_per_species/Paramormyrops_kingsleyae.csv",
        header=TRUE,
        sep=",")



species1 <- "Paramormyrops_kingsleyae"
species2 <- "Ictalurus_punctatus"

colnames(current_ks_df) <- 
  c("OGG","seq1", "seq2", "ks", "ka", "vks", "vka", "nt_used", "codon_used", "cds_identity")

current_ks_df <- 
  current_ks_df %>%
  filter(ks != 9.99999900) %>% 
  filter(ks > 0) %>%
  filter(codon_used > 100) 

putative_TE_genes <- gsub("---", "_", putative_TE_genes)
putative_TE_genes <- gsub("\\.", "_", putative_TE_genes)
putative_TE_genes <- gsub("-", "_", putative_TE_genes)

current_ks_df <- 
  current_ks_df %>%
  filter(! seq1 %in% putative_TE_genes) %>%
  filter(! seq2 %in% putative_TE_genes) %>%
  filter(! OGG %in% putative_OGG_TE_pfam_blast)


current_ks_df <- 
  current_ks_df %>% 
  filter(grepl(species2, seq2))

curr_stats_df_long <- 
  curr_stats_df %>% 
  filter(species_2 == species2) %>%
  dplyr::select(species_1, species_2, quantile_OGG_001, quantile_OGG_002, quantile_OGG_003,
                quantile_OGG_004, quantile_OGG_005) %>%
  pivot_longer(!c(species_1, species_2), names_to = "quantile", values_to = "ks")


mean_ks <- 
  curr_stats_df %>% 
  filter(species_2 == species2) %>%
  pull(mean_OGG)

current_ks_HGT <-
  current_ks_df %>%
  filter(seq1 %in% c("Paramormyrops_kingsleyae_rna_XM_023831520_1")) 

pdf("N5.HOG0010622/ks_distribution-Paramormyrops_kingsleya-Ictalurus_punctatus.pdf",width = 8.34,  height = 4.61)

current_ks_df %>%
  ggplot(., aes(x=ks)) +
  geom_histogram(bins=50) +
  xlab("Ks") +
  geom_vline(data = curr_stats_df_long, aes(xintercept = ks), linetype = "dashed", color = "black") +
  geom_vline(xintercept = mean_ks, color = "black") +
  geom_vline(data = current_ks_HGT, aes(xintercept = ks), color = "red", alpha=0.5)  +
  theme_minimal() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none")

dev.off()

#### Species Ks distribution -- Paramormyrops_kingsleyae vs Pygocentrus_nattereri   ---------------------------------


current_ks_df <- 
  fread("HGT_stats_Results_per_sp/Paramormyrops_kingsleyae.cds.stats",
        header=FALSE,
        sep=",")

curr_stats_df <-
  fread("KS_Stats_per_species/Paramormyrops_kingsleyae.csv",
        header=TRUE,
        sep=",")



species1 <- "Paramormyrops_kingsleyae"
species2 <- "Pygocentrus_nattereri"

colnames(current_ks_df) <- 
  c("OGG","seq1", "seq2", "ks", "ka", "vks", "vka", "nt_used", "codon_used", "cds_identity")

current_ks_df <- 
  current_ks_df %>%
  filter(ks != 9.99999900) %>% 
  filter(ks > 0) %>%
  filter(codon_used > 100) 

putative_TE_genes <- gsub("---", "_", putative_TE_genes)
putative_TE_genes <- gsub("\\.", "_", putative_TE_genes)
putative_TE_genes <- gsub("-", "_", putative_TE_genes)

current_ks_df <- 
  current_ks_df %>%
  filter(! seq1 %in% putative_TE_genes) %>%
  filter(! seq2 %in% putative_TE_genes) %>%
  filter(! OGG %in% putative_OGG_TE_pfam_blast)


current_ks_df <- 
  current_ks_df %>% 
  filter(grepl(species2, seq2))

curr_stats_df_long <- 
  curr_stats_df %>% 
  filter(species_2 == species2) %>%
  dplyr::select(species_1, species_2, quantile_OGG_001, quantile_OGG_002, quantile_OGG_003,
                quantile_OGG_004, quantile_OGG_005) %>%
  pivot_longer(!c(species_1, species_2), names_to = "quantile", values_to = "ks")


mean_ks <- 
  curr_stats_df %>% 
  filter(species_2 == species2) %>%
  pull(mean_OGG)

current_ks_HGT <-
  current_ks_df %>%
  filter(seq1 %in% c("Paramormyrops_kingsleyae_rna_XM_023831520_1")) 

pdf("N5.HOG0010622/ks_distribution-Paramormyrops_kingsleya-Pygocentrus_nattereri.pdf",width = 8.34,  height = 4.61)

current_ks_df %>%
  ggplot(., aes(x=ks)) +
  geom_histogram(bins=50) +
  xlab("Ks") +
  geom_vline(data = curr_stats_df_long, aes(xintercept = ks), linetype = "dashed", color = "black") +
  geom_vline(xintercept = mean_ks, color = "black") +
  geom_vline(data = current_ks_HGT, aes(xintercept = ks), color = "red", alpha=0.5)  +
  theme_minimal() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none")

dev.off()

#### Species Ks distribution -- Paramormyrops_kingsleyae vs Electrophorus_electricus   ---------------------------------


current_ks_df <- 
  fread("HGT_stats_Results_per_sp/Paramormyrops_kingsleyae.cds.stats",
        header=FALSE,
        sep=",")

curr_stats_df <-
  fread("KS_Stats_per_species/Paramormyrops_kingsleyae.csv",
        header=TRUE,
        sep=",")



species1 <- "Paramormyrops_kingsleyae"
species2 <- "Electrophorus_electricus"

colnames(current_ks_df) <- 
  c("OGG","seq1", "seq2", "ks", "ka", "vks", "vka", "nt_used", "codon_used", "cds_identity")

current_ks_df <- 
  current_ks_df %>%
  filter(ks != 9.99999900) %>% 
  filter(ks > 0) %>%
  filter(codon_used > 100) 

putative_TE_genes <- gsub("---", "_", putative_TE_genes)
putative_TE_genes <- gsub("\\.", "_", putative_TE_genes)
putative_TE_genes <- gsub("-", "_", putative_TE_genes)

current_ks_df <- 
  current_ks_df %>%
  filter(! seq1 %in% putative_TE_genes) %>%
  filter(! seq2 %in% putative_TE_genes) %>%
  filter(! OGG %in% putative_OGG_TE_pfam_blast)


current_ks_df <- 
  current_ks_df %>% 
  filter(grepl(species2, seq2))

curr_stats_df_long <- 
  curr_stats_df %>% 
  filter(species_2 == species2) %>%
  dplyr::select(species_1, species_2, quantile_OGG_001, quantile_OGG_002, quantile_OGG_003,
                quantile_OGG_004, quantile_OGG_005) %>%
  pivot_longer(!c(species_1, species_2), names_to = "quantile", values_to = "ks")


mean_ks <- 
  curr_stats_df %>% 
  filter(species_2 == species2) %>%
  pull(mean_OGG)

current_ks_HGT <-
  current_ks_df %>%
  filter(seq1 %in% c("Paramormyrops_kingsleyae_rna_XM_023831520_1")) 

pdf("N5.HOG0010622/ks_distribution-Paramormyrops_kingsleya-Electrophorus_electricus.pdf",width = 8.34,  height = 4.61)

current_ks_df %>%
  ggplot(., aes(x=ks)) +
  geom_histogram(bins=50) +
  xlab("Ks") +
  geom_vline(data = curr_stats_df_long, aes(xintercept = ks), linetype = "dashed", color = "black") +
  geom_vline(xintercept = mean_ks, color = "black") +
  geom_vline(data = current_ks_HGT, aes(xintercept = ks), color = "red", alpha=0.5)  +
  theme_minimal() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none")

dev.off()

#### Orthogroup Ks distribution   ---------------------------------


GLS_function <- function(gls_rslt) {
  GLS_cc <- coef(gls_rslt)
  GLS_f <- function(x) GLS_cc[1] + GLS_cc[2]*x
  return(GLS_f)
}

all_comb <- apply(combn(list_teleost_sp,2),2,paste,collapse=',')
pairs_df <- as.data.frame(t(sapply(all_comb, function(x) unlist(strsplit(x, ",")))))
row.names(pairs_df) <- NULL
colnames(pairs_df) <- c("species_1", "species_2")



HOG_stats_df <- 
  fread("Coding_sequences_alignments/N5.HOG0010622.cds.stats",
        sep=",",
        header=TRUE)


HOG_stats_df_sp <- as.data.frame(NULL)
for(curr_sp in list_teleost_sp){
  
  curr_df <- 
    HOG_stats_df %>%
    filter(grepl(curr_sp, seq1)) %>%
    mutate(species_1 = curr_sp)
  
  HOG_stats_df_sp <- rbind(HOG_stats_df_sp, curr_df)
  
}

HOG_stats_df_species <- as.data.frame(NULL)
for(curr_sp in list_teleost_sp){
  
  curr_df <- 
    HOG_stats_df_sp %>%
    filter(grepl(curr_sp, seq2)) %>%
    mutate(species_2 = curr_sp)
  
  HOG_stats_df_species <- rbind(HOG_stats_df_species, curr_df)
  
}


#Remove redundant comparisons

HOG_stats_df_species <- 
  left_join(pairs_df, HOG_stats_df_species, by=c("species_1", "species_2")) %>%
  filter(! is.na(ks))


#Add divergence times

HOG_stats_df_species$species_1 <- 
  gsub("Salvelinus_sp_IW2-2015", "Salvelinus_sp_IW2_2015", HOG_stats_df_species$species_1)
HOG_stats_df_species$species_2 <- 
  gsub("Salvelinus_sp_IW2-2015", "Salvelinus_sp_IW2_2015", HOG_stats_df_species$species_2)

HOG_stats_df_species <- 
  as.data.frame(HOG_stats_df_species %>%
                  rowwise() %>%
                  mutate(div_time = 
                           get_pairwise_distances(teleost_tree, 
                                                  species_1, 
                                                  species_2, 
                                                  as_edge_counts=FALSE, 
                                                  check_input=TRUE)/2))


#Filter for aberrant ks values

HOG_stats_df_species <- 
  HOG_stats_df_species %>%
  filter(ks != 9.99999900) %>% 
  filter(ks > 0) %>%
  filter(codon_used > 100) 


#Find the slope of the linear regression  and plot

HOG_stats_df_species$divergence_time <- HOG_stats_df_species$div_time * 1000
ks_vs_time_lm <- lm(ks ~ divergence_time, 
                    data = HOG_stats_df_species)
summary_lm <- summary(ks_vs_time_lm)
slope_lm_ks <- summary_lm$coefficients[2]
inter_lm_ks <- summary_lm$coefficients[1]
curr_function_lm <- GLS_function(summary_lm)

my_seq1 <-  c("Paramormyrops_kingsleyae_rna_XM_023831520_1")

my_seq2 <- c("Neoarius_graeffei_rna_XM_060902382_1",
             "Neoarius_graeffei_rna_XM_060901346_1",
             "Silurus_meridionalis_rna_XM_046865888_1",
             "Silurus_asotus_mrna_C0J50_21101",
             "Ictalurus_furcatus_rna_XM_053653479_1",
             "Ictalurus_punctatus_rna_XM_017491464_3",
             "Ameiurus_melas_AMELA_T00023050",
             "Pangasianodon_hypophthalmus_rna_XM_026916402_3",
             "Hemibagrus_wyckioides_rna_XM_058413206_1",
             "Hemibagrus_guttatus_male_06655",
             "Tachysurus_vachellii_rna_XM_060892018_1",
             "Tachysurus_fulvidraco_rna_XM_027141518_2",
             "Bagarius_yarrelli_mrna_Baya_4298",
             "Clarias_magur_DAT39_007587_RA_mrna",
             "Clarias_gariepinus_rna_XM_053514871_1",
             "Trichomycterus_rosablanca_rna_XM_062984713_1",
             "Hemibagrus_wyckioides_rna_XM_058394735_1",
             "Hemibagrus_guttatus_male_16671",
             "Tachysurus_vachellii_rna_XM_060873936_1",
             "Silurus_asotus_mrna_C0J50_8883",
             "Ameiurus_melas_AMELA_T00177240",
             "Ictalurus_punctatus_rna_XM_017463853_3",
             "Pangasianodon_hypophthalmus_rna_XM_053232595_1",
             "Trichomycterus_rosablanca_rna_XM_062995414_1")


my_seq3 <- c("Pygocentrus_nattereri_rna_XM_017714086_2",
             "Colossoma_macropomum_rna_XM_036589671_1",
             "Astyanax_mexicanus_rna_XM_007241775_4")

my_seq4 <- c("Electrophorus_electricus_rna_XM_027006769_2")

HOG_stats_df_species <- 
  HOG_stats_df_species %>%
  mutate(group_hgt = case_when(
    (seq1 %in% my_seq1 & seq2 %in% my_seq2) | (seq1 %in% my_seq2 & seq2 %in% my_seq1) ~ "HGT",
    (seq1 %in% my_seq1 & seq2 %in% my_seq3) | (seq1 %in% my_seq3 & seq2 %in% my_seq1) ~ "HGT_2",
    (seq1 %in% my_seq1 & seq2 %in% my_seq4) | (seq1 %in% my_seq4 & seq2 %in% my_seq1) ~ "HGT_3",
  ))

HOG_stats_df_species <- HOG_stats_df_species %>%
  mutate(group_hgt = ifelse(is.na(group_hgt), "noHGT", group_hgt))




pdf("N5.HOG0010622/HOG_ks_distribution.pdf",width = 8.34,  height = 4.61)

HOG_stats_df_species %>%
  ggplot(., aes(x=(divergence_time), y=ks, color=group_hgt)) +
  geom_point() + 
  scale_color_manual(values = c("no_HGT" = "#8A8686", "HGT" = "#D55E00", "HGT_2" = "#1E88E5", "HGT_3" = "#004D40")) + 
  stat_function(fun = curr_function_lm, color="black") +
  theme_classic() +
  ylim(0, 5) +
  xlab("Divergence time (Mya)") +
  ylab("Ks") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none")

dev.off()


#### Draw synteny  ---------------------------------


my_seqs.clade1 <- read.table("cluster_id_arranged.clade1.txt", sep=",")
colnames(my_seqs.clade1) <- c("bin_id", "seq_id", "length")
my_seqs.clade1 <- tibble(my_seqs.clade1)

my_genes.clade1 <- read.table("seq_clustered_infos_ogg_num.clade1.csv", sep=",")
colnames(my_genes.clade1) <-  c("seq_id", "start", "end","gene_name","OGG", "strand", "species")
my_genes.clade1 <- tibble(my_genes.clade1)


my_links.clade1 <- read.table("link_table.clade1.txt", sep=",")
colnames(my_links.clade1) <-  c("species1","seq_id", "start", "end", "species2","seq_id2", "start2","end2", "OGG")
my_links.clade1 <- tibble(my_links.clade1)


pair_to_fix_onlyscaff <- 
  c("Ictalurus_punctatus,Pangasianodon_hypophthalmus")

for(species_pair in pair_to_fix_onlyscaff){
  
  my_sp1 <- gsub(",.*","",species_pair)
  my_sp2 <- gsub(".*,","",species_pair)
  
  
  test <- 
    my_links.clade1 %>%
    filter(species1 == my_sp1) %>%
    filter(species2 == my_sp2) %>%
    filter(seq_id == "NC_030417_2") %>%
    mutate(start3 = end2) %>%
    mutate(end3 = start2) %>%
    dplyr::select(species1, seq_id, start, end, species2, seq_id2, start3, end3, OGG)
  colnames(test) <- colnames(my_links.clade1)
  
  my_links.clade1_new <- 
    my_links.clade1 %>% 
    filter(! (species1 == my_sp1 & species2 == my_sp2 & seq_id == "NC_030417_2"))
  
  my_links.clade1 <- rbind(my_links.clade1_new, test)
  
  my_links.clade1 <- tibble(my_links.clade1)
  
  
}

pair_to_fix_onlyscaff <- 
  c("Electrophorus_electricus,Pygocentrus_nattereri")

for(species_pair in pair_to_fix_onlyscaff){
  
  my_sp1 <- gsub(",.*","",species_pair)
  my_sp2 <- gsub(".*,","",species_pair)
  
  
  test <- 
    my_links.clade1 %>%
    filter(species1 == my_sp1) %>%
    filter(species2 == my_sp2) %>%
    filter(seq_id == "NC_049535_1") %>%
    mutate(start3 = end2) %>%
    mutate(end3 = start2) %>%
    dplyr::select(species1, seq_id, start, end, species2, seq_id2, start3, end3, OGG)
  colnames(test) <- colnames(my_links.clade1)
  
  my_links.clade1_new <- 
    my_links.clade1 %>% 
    filter(! (species1 == my_sp1 & species2 == my_sp2 & seq_id == "NC_049535_1"))
  
  my_links.clade1 <- rbind(my_links.clade1_new, test)
  
  my_links.clade1 <- tibble(my_links.clade1)
  
  
}


pair_to_fix_onlyscaff <- 
  c("Brienomyrus_brachyistius,Paramormyrops_kingsleyae")

for(species_pair in pair_to_fix_onlyscaff){
  
  my_sp1 <- gsub(",.*","",species_pair)
  my_sp2 <- gsub(".*,","",species_pair)
  
  
  test <- 
    my_links.clade1 %>%
    filter(species1 == my_sp1) %>%
    filter(species2 == my_sp2) %>%
    filter(seq_id == "NC_064539_1") %>%
    mutate(start3 = end2) %>%
    mutate(end3 = start2) %>%
    dplyr::select(species1, seq_id, start, end, species2, seq_id2, start3, end3, OGG)
  colnames(test) <- colnames(my_links.clade1)
  
  my_links.clade1_new <- 
    my_links.clade1 %>% 
    filter(! (species1 == my_sp1 & species2 == my_sp2 & seq_id == "NC_064539_1"))
  
  my_links.clade1 <- rbind(my_links.clade1_new, test)
  
  my_links.clade1 <- tibble(my_links.clade1)
  
  
}

p.clade1 <-  
  gggenomes(seqs=my_seqs.clade1, genes=my_genes.clade1, links=my_links.clade1) %>%
  flip_seqs("NC_030417_2") %>%
  flip_seqs("NC_051224_1") %>%
  flip_seqs("NW_019713398_1")

#Lets make a  color plaette, where everything is gray expect the HGT OGG

all_OGG <- my_genes.clade1$OGG %>% unique()
all_OGG_gray_vec <- setNames(rep("gray", length(all_OGG)), all_OGG)
all_OGG_gray_vec["N5_HOG0010622"] <- "red"


p1 <- 
  p.clade1 + geom_seq() + 
  geom_gene(aes(fill = OGG)) + 
  geom_link(aes(fill=OGG, color=OGG)) +
  scale_fill_manual(values = all_OGG_gray_vec) +
  scale_color_manual(values = all_OGG_gray_vec) +
  theme(legend.position = "none") 

#Save as PDF
ggsave("Micro_synteny.pdf", p1, width=8.34, height=6.34, units="in", scale=2)





#### Draw TEs locations  ---------------------------------

my_seqs.clade1 <- read.table("clusters_ID_TE.txt", sep=",")
colnames(my_seqs.clade1) <- c("bin_id", "seq_id", "length")
my_seqs.clade1 <- tibble(my_seqs.clade1)

my_genes.clade1 <- read.table("seq_clustered_infos_ogg.TE.txt", sep=",")
colnames(my_genes.clade1) <-  c("seq_id", "start", "end","exon_name","OGG", "strand", "species")
my_genes.clade1 <- tibble(my_genes.clade1)


all_OGG <- my_genes.clade1$OGG %>% unique()
all_OGG_gray_vec <- setNames(rep("gray", length(all_OGG)), all_OGG)
all_OGG_gray_vec["N5_HOG0010622"] <- "red"
#all_OGG_gray_vec["L1-19_SSa_2p:ClassI:LINE:L1:L1"] <- "#FFC107"

p.clade1 <-  
  gggenomes(seqs=my_seqs.clade1, genes=my_genes.clade1)


synt_TE_p <- 
  p.clade1 + geom_seq(arrow = NULL) + 
  geom_gene(aes(fill = OGG, color=OGG)) +
  theme(legend.position = "none")  +
  scale_fill_manual(values = all_OGG_gray_vec) +
  scale_color_manual(values = all_OGG_gray_vec) 



#Save as PDF
ggsave("Micro_synteny.TE.pdf", synt_TE_p, width=8, height=3, units="in", scale=2)


#### Region alignment -- Paramormyrops versus Ictalurus  ---------------------------------

library(dotplot)

paramormyrops_region_1 <- seqinr::read.fasta("N5.HOG0010622.Paramormyrops_kingsleyae.extended.1.fa", seqtype = "DNA", as.string = TRUE)[[1]]
ictalurus_region_1 <- seqinr::read.fasta("N5.HOG0010622.Ictalurus_punctatus.extended.1.fa", seqtype = "DNA", as.string = TRUE)[[1]]




dotplot.p1 <- 
  dotPlotg(paramormyrops_region_1, ictalurus_region_1, wsize = 60, wstep = 1, nmatch = 38)  +
  geom_rect(aes(xmin=5221, xmax=9782, ymin=0, ymax=Inf), fill=NA, color="indianred1", linewidth=0.2, alpha=0.4) + 
  geom_rect(aes(xmin=9890, xmax=9893, ymin=0, ymax=Inf), fill=NA, color="indianred1", linewidth=0.2, alpha=0.4) + 
  
  geom_rect(aes(ymin=2872, ymax=7472, xmin=0, xmax=Inf), fill=NA, color="indianred1", linewidth=0.2, alpha=0.4) + 
  geom_rect(aes(ymin=7561, ymax=7564, xmin=0, xmax=Inf), fill=NA, color="indianred1", linewidth=0.2, alpha=0.4) + 
  
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") +
  xlab("P. kingsleyae (NW_019713398.1:637791-652729)")+
  ylab("I. punctatus (NC_030417.2:23089172-23102712)")


ictalurus_region_2 <- seqinr::read.fasta("N5.HOG0010622.Ictalurus_punctatus.extended.2.fa", seqtype = "DNA", as.string = TRUE)[[1]]

dotplot.p2 <- 
  dotPlotg(paramormyrops_region_1, ictalurus_region_2, wsize = 60, wstep = 1, nmatch = 38)  +
  geom_rect(aes(xmin=5221, xmax=9782, ymin=0, ymax=Inf), fill=NA, color="indianred1", linewidth=0.2, alpha=0.4) + 
  geom_rect(aes(xmin=9890, xmax=9893, ymin=0, ymax=Inf), fill=NA, color="indianred1", linewidth=0.2, alpha=0.4) + 
  
  geom_rect(aes(ymin=5389, ymax=9869, xmin=0, xmax=Inf), fill=NA, color="indianred1", linewidth=0.2, alpha=0.4) + 
  geom_rect(aes(ymin=11354, ymax=11378, xmin=0, xmax=Inf), fill=NA, color="indianred1", linewidth=0.2, alpha=0.4) + 
  
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") +
  xlab("P. kingsleyae (NW_019713398.1:637791-652729)")+
  ylab("I. punctatus (NW_026521115.1:118979-133731)")



ggsave("dotplot-Paramormyrops_kingsleyae-Ictalurus_punctatus.1.pdf", 
       dotplot.p1, width=8.34, height=6.34, units="in", scale=2)

ggsave("dotplot-Paramormyrops_kingsleyae-Ictalurus_punctatus.2.pdf", 
       dotplot.p2, width=8.34, height=6.34, units="in", scale=2)


