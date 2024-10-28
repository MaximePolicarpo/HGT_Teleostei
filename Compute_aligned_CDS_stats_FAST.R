library("bio3d")
library("tidyr")
library("reshape2")
library(dplyr)
library("seqinr")
library(MSA2dist)
library(Biostrings)

args = commandArgs(trailingOnly=TRUE)

OG_align <- seqinr::read.alignment(args[1], "fasta", forceToLower = TRUE)
ka_ks_df <- kaks(OG_align, verbose = FALSE, debug = FALSE, forceUpperCase = TRUE, rmgap = FALSE)
ks_OG <- ka_ks_df$ks
ka_OG <- ka_ks_df$ka
vks_OG <- ka_ks_df$vks
vka_OG <- ka_ks_df$vka
ks_OG_df <- melt(as.matrix(ks_OG), varnames = c("row", "col"))
ka_OG_df <- melt(as.matrix(ka_OG), varnames = c("row", "col"))
vks_OG_df <- melt(as.matrix(vks_OG), varnames = c("row", "col"))
vka_OG_df <- melt(as.matrix(vka_OG), varnames = c("row", "col"))
colnames(ks_OG_df) <- c("seq1", "seq2", "ks")
colnames(ka_OG_df) <- c("seq1", "seq2", "ka")
colnames(vks_OG_df) <- c("seq1", "seq2", "vks")
colnames(vka_OG_df) <- c("seq1", "seq2", "vka")
ks_ka_df <- left_join(ks_OG_df, ka_OG_df, by=c("seq1", "seq2"))
ks_ka_vks_df <- left_join(ks_ka_df, vks_OG_df, by=c("seq1", "seq2"))
ks_ka_vks_vka_df <- left_join(ks_ka_vks_df, vka_OG_df, by=c("seq1", "seq2"))


#compute the number of nucleotides compared
OG_align_string <- aln2dnastring(OG_align)
dna_dist <- OG_align_string |> dnastring2dist(model="raw", threads = 10)
dna_dist_sitesUsed <- dna_dist$sitesUsed
dna_dist_sitesUsed_df <- melt(as.matrix(dna_dist_sitesUsed), varnames = c("row", "col"))
colnames(dna_dist_sitesUsed_df) <- c("seq1", "seq2", "nt_used")
dna_dist_sitesUsed_df <- dna_dist_sitesUsed_df %>% mutate(codon_used = nt_used/3)


#compute sequence identity
OG_align <- bio3d::read.fasta(args[1], rm.dup = TRUE, to.upper = FALSE, to.dash=TRUE)
OG_align_ident <- as.data.frame(seqidentity(OG_align, normalize=TRUE, similarity=FALSE, ncore=10, nseg.scale=1))
OG_align_ident <- tibble::rownames_to_column(OG_align_ident, "Seq1")
OG_align_ident_reshape <- as.data.frame(OG_align_ident %>% pivot_longer(!Seq1,names_to = "Seq2", values_to = "Identity"))
colnames(OG_align_ident_reshape) <- c("seq1", "seq2", "identity")


#merge tables
ks_ka_vks_vka_used <- left_join(ks_ka_vks_vka_df, dna_dist_sitesUsed_df, by=c("seq1", "seq2"))
ks_ka_vks_vka_used_ident <- left_join(ks_ka_vks_vka_used, OG_align_ident_reshape, by=c("seq1", "seq2"))


write.table(ks_ka_vks_vka_used_ident, file = args[2], sep=",", row.names=FALSE, col.names=TRUE, quote=FALSE)
