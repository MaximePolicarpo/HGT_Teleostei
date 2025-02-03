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

curr_OGG <- "N5.HOG0013514"


current_tree <- 
  read.tree("N5.HOG0013514.prot.aln.treefile")

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


#Prune tree for sequences present in prot_alignment_trim
setwd("N5.HOG0013514")
seq_retained <-
  scan("seq_retained.id", what="character")

current_tree_rooted_retained <- keep.tip(current_tree_rooted, seq_retained)

#Add the alignment
prot_alignment = "N5.HOG0013514.prot.aln.retained"
prot_alignment_trim = "N5.HOG0013514.prot.trimmed.aln"




p0 <- 
  ggtree(current_tree_rooted_retained) %<+% tips_df_species_order +
  geom_tiplab(size = 1) +
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


ggsave(paste("N5.HOG0013514/",
             curr_OGG,
             ".gene.prot_alignment.pdf",
             sep=""), 
       p1, width=15.34, height=24.34, 
       units="in", scale=2)

ggsave(paste("N5.HOG0013514/",
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


ggsave(paste("N5.HOG0013514/",
             curr_OGG,
             ".species.pdf",
             sep=""), 
       p2, width=8, height=6, 
       units="in", scale=2)

#### Draw  phylogenies with alignment -- Blastp   ---------------------------------


current_tree <- 
  read.tree("Uniprot_plus_closefish.aln.trimmed.treefile")

#Root tree at midpoint 
outgroup_id <- scan("closest_nonfish_seq.id", what="character")
outgroup_id_true <- 
  outgroup_id[! outgroup_id %in% c("sp_D3YWJ0_SLIP_MOUSE_Nuclear_GTPase_SLIP-GC_OS_Mus_musculus_OX_10090_GN_Nuggc_PE_3_SV_1",
                                 "sp_Q68CJ6_SLIP_HUMAN_Nuclear_GTPase_SLIP-GC_OS_Homo_sapiens_OX_9606_GN_NUGGC_PE_2_SV_3")]


MRCA_outgroup <- findMRCA(current_tree, tips=outgroup_id_true, type="node")
current_tree_rooted <-  root(current_tree, node=MRCA_outgroup, resolve.root= TRUE)


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

ggsave(paste("N5.HOG0013514/",
             curr_OGG,
             ".uniprot.prot_alignment.pdf",
             sep=""), 
       p2, width=9, height=12, 
       units="in", scale=2)

ggsave(paste("N5.HOG0013514/",
             curr_OGG,
             ".uniprot.prot_alignment_trim.pdf",
             sep=""), 
       p3, width=9, height=12, 
       units="in", scale=2)



#### DNA subtree   ---------------------------------


current_tree <- 
  read.tree("N5.HOG0013514/only_Clupea_Osmeriformes.cds.aln.treefile")

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
  xlim(0, 2) +
  geom_treescale(x=0.1, y=15, width=0.2, color='black') 


DNA_tree







#### Species Ks distribution -- Hypomesus vs Clupea   ---------------------------------


current_ks_df <- 
  fread("HGT_stats_Results_per_sp/Hypomesus_transpacificus.cds.stats",
        header=FALSE,
        sep=",")

curr_stats_df <-
  fread("KS_Stats_per_species/Hypomesus_transpacificus.csv",
        header=TRUE,
        sep=",")



species1 <- "Hypomesus_transpacificus"
species2 <- "Clupea_harengus"

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
  filter(seq1 %in% c("Hypomesus_transpacificus_rna_XM_047045767_1", "Hypomesus_transpacificus_rna_XM_047045786_1")) 


pdf("N5.HOG0013514/ks_distribution.hypo_clupea.pdf",width = 8.34,  height = 4.61)

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

#### Species Ks distribution -- Osmerus vs Clupea   ---------------------------------


current_ks_df <- 
  fread("HGT_stats_Results_per_sp/Osmerus_eperlanus.cds.stats",
        header=FALSE,
        sep=",")

curr_stats_df <-
  fread("KS_Stats_per_species/Osmerus_eperlanus.csv",
        header=TRUE,
        sep=",")



species1 <- "Osmerus_eperlanus"
species2 <- "Clupea_harengus"

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
  filter(seq1 %in% c("Osmerus_eperlanus_rna_XM_062468715_1", "Osmerus_eperlanus_rna_XM_062468701_1", "Osmerus_eperlanus_rna_XM_062468699_1")) 


pdf("N5.HOG0013514/ks_distribution.osmerus_clupea.pdf",width = 8.34,  height = 4.61)

current_ks_df %>%
  ggplot(., aes(x=ks)) +
  geom_histogram(bins=50) +
  xlab("Ks") +
  geom_vline(data = curr_stats_df_long, aes(xintercept = ks), linetype = "dashed", color = "black") +
  geom_vline(xintercept = mean_ks, color = "black") +
  geom_vline(data = current_ks_HGT, aes(xintercept = ks), color = "red")  +
  theme_minimal() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none")

dev.off()


#### Species Ks distribution -- Danio_rerio vs Clupea   ---------------------------------


current_ks_df <- 
  fread("HGT_stats_Results_per_sp/Danio_rerio.cds.stats",
        header=FALSE,
        sep=",")

curr_stats_df <-
  fread("KS_Stats_per_species/Danio_rerio.csv",
        header=TRUE,
        sep=",")



species1 <- "Osmerus_eperlanus"
species2 <- "Clupea_harengus"

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
  filter(OGG == "N5.HOG0013514") 


pdf("N5.HOG0013514/ks_distribution.danio_clupea.pdf",width = 8.34,  height = 4.61)

current_ks_df %>%
  ggplot(., aes(x=ks)) +
  geom_histogram(bins=50) +
  xlab("Ks") +
  geom_vline(data = curr_stats_df_long, aes(xintercept = ks), linetype = "dashed", color = "black") +
  geom_vline(xintercept = mean_ks, color = "black") +
  geom_vline(data = current_ks_HGT, aes(xintercept = ks), color = "red")  +
  theme_minimal() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none")

dev.off()



#### Species Ks distribution -- Astyanax vs Clupea   ---------------------------------


current_ks_df <- 
  fread("HGT_stats_Results_per_sp/Clupea_harengus.cds.stats",
        header=FALSE,
        sep=",")

curr_stats_df <-
  fread("KS_Stats_per_species/Clupea_harengus.csv",
        header=TRUE,
        sep=",")



species1 <- "Clupea_harengus"
species2 <- "Astyanax_mexicanus"

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
  filter(OGG == "N5.HOG0013514") 


pdf("N5.HOG0013514/ks_distribution.astyanax_clupea.pdf",width = 8.34,  height = 4.61)

current_ks_df %>%
  ggplot(., aes(x=ks)) +
  geom_histogram(bins=50) +
  xlab("Ks") +
  geom_vline(data = curr_stats_df_long, aes(xintercept = ks), linetype = "dashed", color = "black") +
  geom_vline(xintercept = mean_ks, color = "black") +
  geom_vline(data = current_ks_HGT, aes(xintercept = ks), color = "red")  +
  theme_minimal() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none")

dev.off()




#### Species Ks distribution -- Electrophorus vs Clupea   ---------------------------------


current_ks_df <- 
  fread("HGT_stats_Results_per_sp/Clupea_harengus.cds.stats",
        header=FALSE,
        sep=",")

curr_stats_df <-
  fread("KS_Stats_per_species/Clupea_harengus.csv",
        header=TRUE,
        sep=",")



species1 <- "Clupea_harengus"
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
  filter(OGG == "N5.HOG0013514") 


pdf("N5.HOG0013514/ks_distribution.electrophorus_clupea.pdf",width = 8.34,  height = 4.61)

current_ks_df %>%
  ggplot(., aes(x=ks)) +
  geom_histogram(bins=50) +
  xlab("Ks") +
  geom_vline(data = curr_stats_df_long, aes(xintercept = ks), linetype = "dashed", color = "black") +
  geom_vline(xintercept = mean_ks, color = "black") +
  geom_vline(data = current_ks_HGT, aes(xintercept = ks), color = "red")  +
  theme_minimal() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none")

dev.off()
#### Species Ks distribution -- Engraulis vs Clupea_harengus   ---------------------------------


current_ks_df <- 
  fread("HGT_stats_Results_per_sp/Engraulis_encrasicolus.cds.stats",
        header=FALSE,
        sep=",")

curr_stats_df <-
  fread("KS_Stats_per_species/Engraulis_encrasicolus.csv",
        header=TRUE,
        sep=",")



species1 <- "Engraulis_encrasicolus"
species2 <- "Clupea_harengus"

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

ks_quantile_1 <- quantile(current_ks_df$ks, 0.001) 
ks_quantile_2 <- quantile(current_ks_df$ks, 0.002) 
ks_quantile_3 <- quantile(current_ks_df$ks, 0.003) 
ks_quantile_4 <- quantile(current_ks_df$ks, 0.004) 
ks_quantile_5 <- quantile(current_ks_df$ks, 0.005) 
mean_ks <- mean(current_ks_df$ks)

curr_stats_df_long <- 
  as.data.frame(
    c(ks_quantile_1, ks_quantile_2, ks_quantile_3, ks_quantile_4, ks_quantile_5))
colnames(curr_stats_df_long) <- c("ks")


current_ks_HGT <-
  current_ks_df %>%
  filter(OGG == "N5.HOG0013514") 

pdf("N5.HOG0013514/ks_distribution.engraulis_clupea.pdf",width = 8.34,  height = 4.61)

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
  fread("Coding_sequences_alignments/N5.HOG0013514.cds.stats",
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

species1 <- "Hypomesus_transpacificus"
species1_sec <- "Osmerus_eperlanus"
species2 <- "Clupea_harengus"

HOG_stats_df_species <- 
  HOG_stats_df_species %>%
  mutate(group_hgt = case_when(
    (species_1 == species1 & species_2 == species2) | (species_1 == species2 & species_2 == species1) ~ "HGT_1",
    (species_1 == species1_sec & species_2 == species2) | (species_1 == species2 & species_2 == species1_sec) ~ "HGT_2"  ))

HOG_stats_df_species <- HOG_stats_df_species %>%
  mutate(group_hgt = ifelse(is.na(group_hgt), "noHGT", group_hgt))


pdf("N5.HOG0013514/HOG_ks_distribution.pdf",width = 8.34,  height = 4.61)

HOG_stats_df_species %>%
  ggplot(., aes(x=(divergence_time), y=ks, color=group_hgt)) +
  geom_point() + 
  scale_color_manual(values = c("no_HGT" = "#8A8686", "HGT_1" = "#D55E00", "HGT_2" = "#009E73")) + 
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

setwd("N5.HOG0013514")
library(gggenomes)

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
  c("Alosa_alosa,Sardina_pilchardus",
    "Sardina_pilchardus,Clupea_harengus")

for(species_pair in pair_to_fix_onlyscaff){
  
  my_sp1 <- gsub(",.*","",species_pair)
  my_sp2 <- gsub(".*,","",species_pair)
  
  
  test <- 
    my_links.clade1 %>%
    filter(species1 == my_sp1) %>%
    filter(species2 == my_sp2) %>%
    mutate(start3 = end2) %>%
    mutate(end3 = start2) %>%
    dplyr::select(species1, seq_id, start, end, species2, seq_id2, start3, end3, OGG)
  colnames(test) <- colnames(my_links.clade1)
  
  my_links.clade1_new <- 
    my_links.clade1 %>% 
    filter(! (species1 == my_sp1 & species2 == my_sp2))
  
  my_links.clade1 <- rbind(my_links.clade1_new, test)
  
  my_links.clade1 <- tibble(my_links.clade1)
  
  
}


p.clade1 <-  
  gggenomes(seqs=my_seqs.clade1, genes=my_genes.clade1, links=my_links.clade1) %>%
  flip_seqs("NC_085000_1")
  





#Lets make a  color plaette, where everything is gray expect the HGT OGG

all_OGG <- my_genes.clade1$OGG %>% unique()
all_OGG_gray_vec <- setNames(rep("gray", length(all_OGG)), all_OGG)
all_OGG_gray_vec["N5_HOG0013514"] <- "red"


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
all_OGG_gray_vec["N5_HOG0013514"] <- "red"
all_OGG_gray_vec["HERO-1_AFC_polLINE/R2-Hero"] <- "#FFC107"

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





#### Draw Full length alignment for Clupea   ---------------------------------


current_tree <- 
  read.tree("Full_length_clupea.aln.treefile")

current_tree_rooted <- midpoint.root(current_tree)


prot_alignment <- "Full_length_clupea.aln"

p0 <- 
  ggtree(current_tree_rooted) + 
  geom_tiplab(size = 3) +
  theme(legend.position = "none")

p1 <- msaplot(p0, prot_alignment, offset=0.5, width=4) + theme(legend.position = "none")



#Save trees as PDF


ggsave("Full_length_alignment.pdf", p1, width=8.34, height=4, units="in", scale=2)

#### Species Ks distribution FULL LENGTH -- Hypomesus vs Clupea   ---------------------------------


current_ks_df <- 
  fread("HGT_stats_Results_per_sp/Hypomesus_transpacificus.cds.stats",
        header=FALSE,
        sep=",")

curr_stats_df <-
  fread("KS_Stats_per_species/Hypomesus_transpacificus.csv",
        header=TRUE,
        sep=",")



species1 <- "Hypomesus_transpacificus"
species2 <- "Clupea_harengus"

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
  filter(seq1 %in% c("Hypomesus_transpacificus_rna_XM_047045767_1", "Hypomesus_transpacificus_rna_XM_047045786_1")) 


current_ks_FULL_length <- 
  fread("Full_length.stats.filt", header=TRUE, sep=",")
current_ks_FULL_length <- 
  current_ks_FULL_length %>% 
  filter(! grepl("NC_045155.1", seq1) & grepl("NC_045155.1", seq2))




pdf("N5.HOG0013514/ks_distribution.hypo_clupea.TRUE.pdf",width = 8.34,  height = 4.61)

current_ks_df %>%
  ggplot(., aes(x=ks)) +
  geom_histogram(bins=50) +
  xlab("Ks") +
  geom_vline(data = curr_stats_df_long, aes(xintercept = ks), linetype = "dashed", color = "black") +
  geom_vline(xintercept = mean_ks, color = "black") +
  geom_vline(data = current_ks_HGT, aes(xintercept = ks), color = "red", alpha=0.5)  +
  geom_vline(data = current_ks_FULL_length, aes(xintercept = ks), color = "#004D40", alpha=0.5)  +
  theme_minimal() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none")

dev.off()


#### Region alignment -- Hypomesus versus Clupea -- first region  ---------------------------------


hypomesus_region <- seqinr::read.fasta("N5.HOG0013514.Hypomesus_transpacificus.extended.1.fa", seqtype = "DNA", as.string = TRUE)[[1]]
clupea_region_1 <- seqinr::read.fasta("N5.HOG0013514.Clupea_harengus.extended.1.fa.rev", seqtype = "DNA", as.string = TRUE)[[1]]


dotplot.p1 <- dotPlotg(hypomesus_region, clupea_region_1, wsize = 40, wstep = 1, nmatch = 28) +
  
  geom_rect(aes(xmin=3709, xmax=3709, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 osm gene
  geom_rect(aes(xmin=6103, xmax=6780, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=7247, xmax=7367, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 osm gene
  geom_rect(aes(xmin=7513, xmax=7735, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=7837, xmax=8049, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 osm gene
  geom_rect(aes(xmin=8168, xmax=8362, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=8987, xmax=9174, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 osm gene
  geom_rect(aes(xmin=9331, xmax=9367, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=9465, xmax=9551, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 osm gene
  geom_rect(aes(xmin=9635, xmax=9722, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=9919, xmax=10052, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 osm gene
  geom_rect(aes(xmin=10148, xmax=10327, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=10450, xmax=10607, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 osm gene
  geom_rect(aes(xmin=10688, xmax=10817, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=11516, xmax=11624, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=11727, xmax=11863, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=12020, xmax=12104, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=12297, xmax=12352, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  
  geom_rect(aes(xmin=31307, xmax=31307, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 osm gene
  geom_rect(aes(xmin=33780, xmax=34463, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=34921, xmax=35041, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 osm gene
  geom_rect(aes(xmin=35191, xmax=35413, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=35515, xmax=35727, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 osm gene
  geom_rect(aes(xmin=35847, xmax=36041, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=36719, xmax=36906, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 osm gene
  geom_rect(aes(xmin=37062, xmax=37098, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=37197, xmax=37301, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 osm gene
  geom_rect(aes(xmin=37386, xmax=37473, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=37674, xmax=37807, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 osm gene
  geom_rect(aes(xmin=37903, xmax=38082, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=38205, xmax=38362, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 osm gene
  geom_rect(aes(xmin=38445, xmax=38574, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=38995, xmax=39103, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=39201, xmax=39334, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=39491, xmax=39575, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=39768, xmax=39823, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  
  geom_rect(aes(xmin=30607, xmax=30945, ymin=0, ymax=Inf, color="#FFC107"), fill=NA, linewidth=0.2, alpha=0.4) +
  
  geom_rect(aes(ymin=59681, ymax=60359, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=60820, ymax=60940, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=61104, ymax=61326, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=61428, ymax=61640, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=61759, ymax=61953, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=62617, ymax=62804, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=62968, ymax=63004, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=63103, ymax=63195, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=63279, ymax=63366, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=63574, ymax=63707, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=63803, ymax=63982, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=64105, ymax=64262, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=64343, ymax=64472, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=65157, ymax=65265, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=65370, ymax=65503, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=65660, ymax=65744, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=65937, ymax=65989, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  
  geom_rect(aes(ymin=45134, ymax=45809, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=46188, ymax=46308, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=46473, ymax=46695, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=46795, ymax=47007, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=47115, ymax=47309, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=47842, ymax=48029, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=48190, ymax=48226, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=48324, ymax=48416, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=48501, ymax=48588, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=48778, ymax=48911, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=49007, ymax=49186, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=49309, ymax=49466, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=49549, ymax=49684, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=50709, ymax=50817, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=50913, ymax=51046, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=51203, ymax=51287, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=51488, ymax=51540, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  
  geom_rect(aes(ymin=6173,  ymax=6857, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=7368,  ymax=7497, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=7664,  ymax=7877, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=7979,  ymax=8191, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=9442,  ymax=9636, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=10235, ymax=10422, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=10591, ymax=10627, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=10726, ymax=10818, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=10904, ymax=10991, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=11194, ymax=11327, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=11423, ymax=11599, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=11726, ymax=11880, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=11964, ymax=12099, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=12848, ymax=12956, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=13052, ymax=13185, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=13298, ymax=13382, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=13584, ymax=13636, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  
  geom_rect(aes(ymin=30344, ymax=30617, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=31274, ymax=31403, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=31578, ymax=31791, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=31893, ymax=32105, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=33896, ymax=34090, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=34687, ymax=34874, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=35044, ymax=35080, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=35179, ymax=35271, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=35357, ymax=35444, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=35647, ymax=35780, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=35876, ymax=36051, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=36178, ymax=36332, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=36415, ymax=36550, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=39164, ymax=39272, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=39368, ymax=39501, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=39614, ymax=39698, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=39905, ymax=39960, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  
  
  geom_rect(aes(ymin=57915, ymax=58477, xmin=0, xmax=Inf, color="#FFC107"), fill=NA, linewidth=0.2, alpha=0.4) +

  
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") +
  xlab("H. transpacificus (NC_061062.1:6295541-6339487)")+
  ylab("C. harengus (NC_045155.1:5717415-5788992)")


ggsave("dotplot_hypomesus_clupea.1.pdf", dotplot.p1, width=8.34, height=6.34, units="in", scale=2)





dotPlotg(hypomesus_region, clupea_region_1, wsize = 40, wstep = 5, nmatch = 28) +
  
  geom_rect(aes(xmin=3709, xmax=3709, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 osm gene
  geom_rect(aes(xmin=6103, xmax=6780, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=7247, xmax=7367, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 osm gene
  geom_rect(aes(xmin=7513, xmax=7735, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=7837, xmax=8049, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 osm gene
  geom_rect(aes(xmin=8168, xmax=8362, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=8987, xmax=9174, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 osm gene
  geom_rect(aes(xmin=9331, xmax=9367, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=9465, xmax=9551, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 osm gene
  geom_rect(aes(xmin=9635, xmax=9722, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=9919, xmax=10052, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 osm gene
  geom_rect(aes(xmin=10148, xmax=10327, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=10450, xmax=10607, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 osm gene
  geom_rect(aes(xmin=10688, xmax=10817, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=11516, xmax=11624, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=11727, xmax=11863, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=12020, xmax=12104, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=12297, xmax=12352, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  
  geom_rect(aes(xmin=31307, xmax=31307, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 osm gene
  geom_rect(aes(xmin=33780, xmax=34463, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=34921, xmax=35041, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 osm gene
  geom_rect(aes(xmin=35191, xmax=35413, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=35515, xmax=35727, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 osm gene
  geom_rect(aes(xmin=35847, xmax=36041, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=36719, xmax=36906, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 osm gene
  geom_rect(aes(xmin=37062, xmax=37098, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=37197, xmax=37301, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 osm gene
  geom_rect(aes(xmin=37386, xmax=37473, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=37674, xmax=37807, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 osm gene
  geom_rect(aes(xmin=37903, xmax=38082, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=38205, xmax=38362, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 osm gene
  geom_rect(aes(xmin=38445, xmax=38574, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=38995, xmax=39103, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=39201, xmax=39334, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=39491, xmax=39575, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=39768, xmax=39823, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  
  geom_rect(aes(xmin=30607, xmax=30945, ymin=0, ymax=Inf, color="#FFC107"), fill=NA, linewidth=0.2, alpha=0.4) +
  
  geom_rect(aes(ymin=59681, ymax=60359, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=60820, ymax=60940, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=61104, ymax=61326, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=61428, ymax=61640, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=61759, ymax=61953, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=62617, ymax=62804, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=62968, ymax=63004, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=63103, ymax=63195, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=63279, ymax=63366, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=63574, ymax=63707, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=63803, ymax=63982, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=64105, ymax=64262, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=64343, ymax=64472, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=65157, ymax=65265, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=65370, ymax=65503, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=65660, ymax=65744, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=65937, ymax=65989, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  
  geom_rect(aes(ymin=45134, ymax=45809, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=46188, ymax=46308, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=46473, ymax=46695, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=46795, ymax=47007, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=47115, ymax=47309, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=47842, ymax=48029, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=48190, ymax=48226, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=48324, ymax=48416, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=48501, ymax=48588, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=48778, ymax=48911, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=49007, ymax=49186, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=49309, ymax=49466, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=49549, ymax=49684, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=50709, ymax=50817, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=50913, ymax=51046, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=51203, ymax=51287, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=51488, ymax=51540, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  
  geom_rect(aes(ymin=6173,  ymax=6857, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=7368,  ymax=7497, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=7664,  ymax=7877, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=7979,  ymax=8191, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=9442,  ymax=9636, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=10235, ymax=10422, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=10591, ymax=10627, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=10726, ymax=10818, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=10904, ymax=10991, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=11194, ymax=11327, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=11423, ymax=11599, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=11726, ymax=11880, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=11964, ymax=12099, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=12848, ymax=12956, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=13052, ymax=13185, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=13298, ymax=13382, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=13584, ymax=13636, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  
  geom_rect(aes(ymin=30344, ymax=30617, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=31274, ymax=31403, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=31578, ymax=31791, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=31893, ymax=32105, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=33896, ymax=34090, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=34687, ymax=34874, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=35044, ymax=35080, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=35179, ymax=35271, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=35357, ymax=35444, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=35647, ymax=35780, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=35876, ymax=36051, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=36178, ymax=36332, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=36415, ymax=36550, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=39164, ymax=39272, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=39368, ymax=39501, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=39614, ymax=39698, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=39905, ymax=39960, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene

  geom_rect(aes(ymin=57915, ymax=58477, xmin=0, xmax=Inf, color="#FFC107"), fill=NA, linewidth=0.2, alpha=0.4) +
  
  
  geom_rect(aes(xmin=12352, xmax=13352, ymin=0, ymax=Inf), color="black", fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 osm gene
  geom_rect(aes(xmin=39823, xmax=40823, ymin=0, ymax=Inf), color="black", fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 osm gene
  
  geom_rect(aes(ymin=13636, ymax=14636, xmin=0, xmax=Inf), color="black", fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=39960, ymax=40960, xmin=0, xmax=Inf), color="black", fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=51540, ymax=52540, xmin=0, xmax=Inf), color="black", fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=65989, ymax=66989, xmin=0, xmax=Inf), color="black", fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  
  
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") +
  xlab("H. transpacificus (NC_061062.1:6295541-6339487)")+
  ylab("C. harengus (NC_045155.1:5717415-5788992)")




#### Region alignment -- Hypomesus versus Clupea -- second region  ---------------------------------

library(dotplot)
setwd("N5.HOG0013514")


hypomesus_region <- seqinr::read.fasta("N5.HOG0013514.Hypomesus_transpacificus.extended.1.fa", seqtype = "DNA", as.string = TRUE)[[1]]
clupea_region_2 <- seqinr::read.fasta("N5.HOG0013514.Clupea_harengus.extended.2.fa.rev", seqtype = "DNA", as.string = TRUE)[[1]]



dotplot.p2 <- dotPlotg(hypomesus_region, clupea_region_2, wsize = 40, wstep = 1, nmatch = 28) +
  
  geom_rect(aes(xmin=3709, xmax=3709, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 osm gene
  geom_rect(aes(xmin=6103, xmax=6780, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=7247, xmax=7367, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 osm gene
  geom_rect(aes(xmin=7513, xmax=7735, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=7837, xmax=8049, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 osm gene
  geom_rect(aes(xmin=8168, xmax=8362, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=8987, xmax=9174, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 osm gene
  geom_rect(aes(xmin=9331, xmax=9367, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=9465, xmax=9551, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 osm gene
  geom_rect(aes(xmin=9635, xmax=9722, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=9919, xmax=10052, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 osm gene
  geom_rect(aes(xmin=10148, xmax=10327, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=10450, xmax=10607, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 osm gene
  geom_rect(aes(xmin=10688, xmax=10817, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=11516, xmax=11624, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=11727, xmax=11863, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=12020, xmax=12104, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=12297, xmax=12352, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  
  geom_rect(aes(xmin=31307, xmax=31307, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 osm gene
  geom_rect(aes(xmin=33780, xmax=34463, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=34921, xmax=35041, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 osm gene
  geom_rect(aes(xmin=35191, xmax=35413, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=35515, xmax=35727, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 osm gene
  geom_rect(aes(xmin=35847, xmax=36041, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=36719, xmax=36906, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 osm gene
  geom_rect(aes(xmin=37062, xmax=37098, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=37197, xmax=37301, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 osm gene
  geom_rect(aes(xmin=37386, xmax=37473, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=37674, xmax=37807, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 osm gene
  geom_rect(aes(xmin=37903, xmax=38082, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=38205, xmax=38362, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 osm gene
  geom_rect(aes(xmin=38445, xmax=38574, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=38995, xmax=39103, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=39201, xmax=39334, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=39491, xmax=39575, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  geom_rect(aes(xmin=39768, xmax=39823, ymin=0, ymax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 2 hypo gene
  
  geom_rect(aes(xmin=30607, xmax=30945, ymin=0, ymax=Inf, color="#FFC107"), fill=NA, linewidth=0.2, alpha=0.4) +
  
  geom_rect(aes(ymin=4072,  ymax=4739, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=5419,  ymax=5532, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=8657,  ymax=8879, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=8961,  ymax=9173, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=9294,  ymax=9488, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=10061, ymax=10248, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=10390, ymax=10432, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=10525, ymax=10611, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=10698, ymax=10785, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=11017, ymax=11150, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=11248, ymax=11424, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=11664, ymax=11821, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=11904, ymax=12033, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=12727, ymax=12835, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=12945, ymax=13078, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=13229, ymax=13313, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  geom_rect(aes(ymin=13515, ymax=13564, xmin=0, xmax=Inf, color="#D6B554"), fill=NA, linewidth=0.1, alpha=0.1) + #exon 1 clupea gene
  
  
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") +
  xlab("H. transpacificus (NC_061062.1:6295541-6339487)")+
  ylab("C. harengus (NC_045155.1:1092148-1106954)")


ggsave("dotplot_hypomesus_clupea.2.pdf", dotplot.p2, width=8.34, height=6.34, units="in", scale=2)







