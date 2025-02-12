livrary(dplyr)
library("UniprotR")

#Import HOG involved in HGT
HGT_HOG <- 
  c("N5.HOG0047987", "N5.HOG0001647", "N5.HOG0013514", "N5.HOG0018901", "N5.HOG0008885", "N5.HOG0004480",
    "N5.HOG0028899", "N5.HOG0004450", "N5.HOG0029633", "N5.HOG0010622", "N5.HOG0013500", "N5.HOG0019032",
    "N5.HOG0004451", "N5.HOG0030200", "N5.HOG0034670", "N5.HOG0019111", "N5.HOG0046344")


#Import diamond blastp table
diamond_df <- 
  read.table("HOG_vs_Uniprot.diamond.tsv",
             sep="\t", 
             header=FALSE)

colnames(diamond_df) <- 
  c("HOG","sseqid" ,"pident", "length" ,"mismatch", "gapopen", "qstart" ,"qend" ,
    "sstart", "send" ,"evalue", "bitscore")
  
diamond_df <- 
  diamond_df %>% 
  filter(evalue < 1e-05)


#Combine diamond and HOG names.
HOG_uniprot <- 
  diamond_df %>%
  mutate(uniprot_acc = 
           gsub("sp\\|", "", sseqid)) %>%
  mutate(uniprot_gene = 
           gsub("\\|.*", "", uniprot_acc)) %>%
  dplyr::select(HOG, uniprot_gene) %>%
  distinct(HOG, .keep_all = TRUE)



#Extract the list of uniprot accessions
list_uniprot_gene <- HOG_uniprot %>% pull(uniprot_gene) %>% unique()


#Find the ontologies for each gene

process_chunk <- function(gene_chunk) {
  attempt <- 1
  max_attempts <- 8  # Set max retries for connection issues
  
  while (attempt <= max_attempts) {
    tryCatch({
      result <- GetProteinGOInfo(gene_chunk, directorypath = NULL)
      return(result)  # Success, return the result
    }, error = function(e) {
      message(sprintf("Error in attempt %d for chunk. Retrying...", attempt))
      Sys.sleep(5)  # Wait before retrying
      attempt <<- attempt + 1
    })
  }
  
  message("Failed after multiple attempts. Skipping this chunk.")
  return(NULL)  # Return NULL if all attempts fail
}


chunk_size <- 200
uniprot_ontology_df <- as.data.frame(NULL)
for (i in seq(1, length(list_uniprot_gene), by = chunk_size)){
  curr_chunk <- list_uniprot_gene[i:min(i + chunk_size - 1, length(list_uniprot_gene))] 
  
  
  uniprot_ontology_df <- rbind(uniprot_ontology_df, process_chunk(curr_chunk))
}



uniprot_ontology_df_GO <- 
  uniprot_ontology_df %>%
  rownames_to_column(var = "uniprot_gene") %>%
  dplyr::select(uniprot_gene, Gene.Ontology.IDs) %>%
  distinct()



colnames(uniprot_ontology_df_GO) <- c("uniprot_gene", "GO")



#Load parent go terms

go_parent_df <- read.table("go_parent_child.tsv", sep="\t", header=FALSE)
colnames(go_parent_df) <- c("par_GO", "desc", "GO")


#Combine the HOG names with ontologies

HOG_uniprot_GO <- left_join(uniprot_ontology_df_GO, HOG_uniprot, by="uniprot_gene")


HOG_uniprot_GO %>% 
  filter(HOG %in% HGT_HOG) %>%
  pull(HOG) %>%
  unique()

#Go from wide to long

HOG_uniprot_GO <- 
  as.data.frame(HOG_uniprot_GO %>%
  separate_rows(GO, sep = ";") %>%
  mutate(GO_nospace = gsub(" ", "", GO)) %>%
  dplyr::select(-GO)) %>%
  dplyr::select(-uniprot_gene)

colnames(HOG_uniprot_GO) <- c("HOG", "GO")


  
#Add the parent go terms

temp_parent_GO <- 
  left_join(HOG_uniprot_GO, go_parent_df, by="GO") %>%
  dplyr::select(-desc) %>%
  filter(! is.na(par_GO)) %>%
  dplyr::select(HOG, par_GO) %>%
  distinct()

colnames(temp_parent_GO) <- c("HOG", "GO")

HOG_uniprot_GO <- rbind(HOG_uniprot_GO , temp_parent_GO)


#Import background list

OGG_good_ks_bad_synt <- scan("OGG_good_ks_bad_synt.txt", what="character")

#Keep only biological processes goterms

GO_BP <- 
  uniprot_ontology_df %>%
  dplyr::select(Gene.Ontology..biological.process.) %>%
  separate_rows(Gene.Ontology..biological.process., sep = ";") %>%
  mutate(GO = gsub(".*\\[", "",Gene.Ontology..biological.process. )) %>%
  mutate(GO_sec = gsub("\\]", "", GO)) %>%
  dplyr::pull(GO_sec) %>%
  unique()

GO_BP <- c(GO_BP, go_parent_df %>% pull(par_GO) %>% unique())
GO_BP <- unique(GO_BP)


HOG_uniprot_GO_BP <- 
  HOG_uniprot_GO %>%
  filter(GO %in% GO_BP) %>%
  filter(! is.na(GO))

#Make a Fisher test to test the enrichment of goterms among HGT genes


HOG_uniprot_GO_BP %>%
  filter(HOG %in% HGT_HOG) %>%
  group_by(GO) %>%
  summarise(count = n()) %>%
  arrange(desc(count))


length(HGT_HOG[HGT_HOG %in% HOG_uniprot_GO_BP$HOG])
length(OGG_good_ks_bad_synt[OGG_good_ks_bad_synt %in% HOG_uniprot_GO_BP$HOG])


GO_term_prez_subset <- 
  HOG_uniprot_GO_BP %>%
  filter(HOG %in% HGT_HOG) %>%
  pull(GO)



HGT_HOG_subset.badsynt <- 
  HOG_uniprot_GO_BP %>%
  filter(HOG %in% c(HGT_HOG, OGG_good_ks_bad_synt)) %>%
  filter(GO %in% GO_term_prez_subset)
  



go_counts.badsynt <- 
  HGT_HOG_subset.badsynt %>%
  mutate(InTestSet = HOG %in% HGT_HOG) %>%
  group_by(GO) %>%
  summarise(
    count_InTest = sum(InTestSet),
    count_InBackground = n() - sum(InTestSet),
    count_TotalTest = 11,
    count_TotalBackground = 4641)



go_counts.badsynt <- 
  go_counts.badsynt %>%
  rowwise() %>%
  mutate(p_value =
           fisher.test(matrix(c(count_InTest, count_InBackground, 
                                count_TotalTest - count_InTest,
                                count_TotalBackground - count_InBackground), 
                              nrow = 2))$p.value) %>%
  ungroup() 






go_counts.badsynt <-
  as.data.frame(
    go_counts.badsynt %>%
  mutate(p_adjusted = p.adjust(p_value, method = "fdr")) %>%
  arrange(p_value)
  )


#Add the HOG associated to goterms

HOG_per_GO <- 
  as.data.frame(
    HOG_uniprot_GO_BP %>%
  filter(HOG %in% HGT_HOG) %>%
  group_by(GO) %>%
  summarise(HOG = paste(unique(HOG), collapse = ";"), .groups = "drop")
  )


go_counts.badsynt <- left_join(go_counts.badsynt, HOG_per_GO, by="GO")


#Add the description of goterms to the fisher table


desc_col <- 
  uniprot_ontology_df %>%
  dplyr::select(Gene.Ontology..biological.process.) %>%
  separate_rows(Gene.Ontology..biological.process., sep = ";") %>%
  mutate(desc = gsub("^ *", "", Gene.Ontology..biological.process.)) %>%
  mutate(desc_good = gsub(" $", "", desc)) %>%
  dplyr::select(desc_good)
  
GO_col <- 
  uniprot_ontology_df %>%
  dplyr::select(Gene.Ontology..biological.process.) %>%
  separate_rows(Gene.Ontology..biological.process., sep = ";") %>%
  mutate(GO = gsub(".*\\[", "",Gene.Ontology..biological.process. )) %>%
  mutate(GO_sec = gsub("\\]", "", GO)) %>%
  dplyr::select(GO_sec)
  
desc_GO_df <- 
  cbind(desc_col, GO_col) %>%
  distinct() 

colnames(desc_GO_df) <- c("description", "GO")


parent_desc <- go_parent_df %>% dplyr::select(desc, par_GO) %>% distinct()
colnames(parent_desc) <- c("description", "GO")

desc_GO_df <- rbind(desc_GO_df, parent_desc)


go_counts.badsynt <- left_join(go_counts.badsynt, desc_GO_df, by="GO")



#Extract significant pvalues

go_counts.badsynt <- 
  go_counts.badsynt %>%
  filter(p_adjusted < 0.05) %>%
  arrange(desc(count_InTest))


Final_GO_df <- 
  go_counts.badsynt %>%
  filter(count_InTest >= 2)

