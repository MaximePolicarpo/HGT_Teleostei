##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##========================================== Now generate an alignment and tree with best blastp matches =============================
##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================

samtools faidx Proteomes_BUSCO80/Hypomesus_transpacificus.fa Hypomesus_transpacificus---rna-XM_047049897.1 > Hypomesus_transpacificus---rna-XM_047049897.1.prot

#blastp against the fish proteome (50 best hits) + Uniprot database without any fish species (10 best hits). Here I will only take 20 catffishes otherwise there are way too much of them

blastp -query Hypomesus_transpacificus---rna-XM_047049897.1.prot -db concatenated_proteomes.fa -outfmt 6 -out Gene_vs_FishProteome.blastp -num_threads 8 -max_target_seqs 51
blastp -query Hypomesus_transpacificus---rna-XM_047049897.1.prot -db non_actino_uniprot.fa -outfmt 6 -out Gene_vs_Uniprot.blastp -num_threads 8  -max_target_seqs 10 -evalue 1e-5

cut -f2 Gene_vs_FishProteome.blastp | sort | uniq  > closest_fish_seq.id
cut -f2 Gene_vs_Uniprot.blastp | sort | uniq > closest_nonfish_seq.id 

xargs samtools faidx concatenated_proteomes.fa < closest_fish_seq.id > closest_fish_seq.fa
xargs samtools faidx non_actino_uniprot.fa  < closest_nonfish_seq.id > closest_nonfish_seq.fa

#Align with muscle and trim with trimal

cat closest_fish_seq.fa closest_nonfish_seq.fa > Uniprot_plus_closefish.fa

muscle5.1 -align closest_fish_seq.fa -output closest_fish_seq.aln
mafft --add closest_nonfish_seq.fa --keeplength closest_fish_seq.aln > Uniprot_plus_closefish.aln
trimal -in Uniprot_plus_closefish.aln -gt 0.7 -out Uniprot_plus_closefish.aln.trimmed

#Make a ML tree with IQ-TREE
iqtree -s Uniprot_plus_closefish.aln -st AA -nt 8 -m TEST -mrate G4 --redo
iqtree -s Uniprot_plus_closefish.aln.trimmed -st AA -nt 8 -m TEST -mrate G4 --redo






##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##========================================== Search the gene in an alternative assembly of Hypomesus  ==============================
##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================

exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs"  Hypomesus_transpacificus---rna-XM_047049897.1.prot GCA_021870715.1_mHypTra1_genomic.fna > Hypomesus_transpacificus.alternative.exo
samtools faidx Hypomesus_transpacificus.cds rna-XM_047049897.1 > Hypomesus_transpacificus---rna-XM_047049897.1.fa

transeq Hypomesus.alternative.fa Hypomesus.alternative.prot ; sed -i 's/_1$//g' Hypomesus.alternative.prot


needle -asequence Hypomesus_transpacificus---rna-XM_047049897.1.fa -bsequence Hypomesus.alternative.fa -outfile Hypomesus.combined.aln -gapopen 10.0 -gapextend 0.5

samtools faidx Coding_sequences_alignments/N5.HOG0019111.prot Thunnus_maccoyii_rna_XM_042435142_1 > N5.HOG0019111.samples.prot
samtools faidx Coding_sequences_alignments/N5.HOG0019111.prot Thunnus_maccoyii_rna_XM_042435120_1 >> N5.HOG0019111.samples.prot
samtools faidx Coding_sequences_alignments/N5.HOG0019111.prot Thunnus_maccoyii_rna_XM_042435094_1 >> N5.HOG0019111.samples.prot
samtools faidx Coding_sequences_alignments/N5.HOG0019111.prot Clupea_harengus_rna_XM_031577666_1 >> N5.HOG0019111.samples.prot
samtools faidx Coding_sequences_alignments/N5.HOG0019111.prot Clupea_harengus_rna_XM_031578285_2 >> N5.HOG0019111.samples.prot
samtools faidx Coding_sequences_alignments/N5.HOG0019111.prot Clupea_harengus_rna_XM_031577667_2 >> N5.HOG0019111.samples.prot
samtools faidx Coding_sequences_alignments/N5.HOG0019111.prot Alosa_sapidissima_rna_XM_042075309_1 >> N5.HOG0019111.samples.prot
samtools faidx Coding_sequences_alignments/N5.HOG0019111.prot Alosa_alosa_rna_XM_048231371_1 >> N5.HOG0019111.samples.prot
samtools faidx Coding_sequences_alignments/N5.HOG0019111.prot Sardina_pilchardus_rna_XM_062535392_1 >> N5.HOG0019111.samples.prot
samtools faidx Coding_sequences_alignments/N5.HOG0019111.prot Thalassophryne_amazonica_rna_XM_034187378_1 >> N5.HOG0019111.samples.prot

cat N5.HOG0019111.samples.prot Hypomesus.alternative.prot Hypomesus_transpacificus---rna-XM_047049897.1.prot > N5.HOG0019111.combined.prot

muscle5.1 -align N5.HOG0019111.combined.prot -output N5.HOG0019111.combined.prot.aln

trimal -in N5.HOG0019111.combined.prot.aln -gt 0.7 -out N5.HOG0019111.combined.prot.aln.trimmed


iqtree -s N5.HOG0019111.combined.prot.aln.trimmed -st AA -nt 8 -m TEST -mrate G4 --redo




#Hypomesus.alternative.fa
>Hypomesus_transpacificus.alternative.GCA_021870715.1.CM038951.1-14039283-14040745
GTGTCTACATACGCCCAGGAGGCCAAGGATGATCTACCATACATAAACAACTACAAGGACAAGATTATCA
GGGTTTCGAACAATTGTAATGTGCAGCCTTCTGTCGTGGCTGGCATCATCTCCAGAGAGACACGCGGAGG
CAGAGGGGCAGGTCTTGATTCCCATGGCTATGGAGACAATGGGAACGGCTATGGACTAATGCAGGTTGAC
AAACGCTACCACACGCTTCAGGGTGCGTGGGACAGCCAGACCAATATTGAGCAAGGAATTGGGATTCTAC
AATACTTCCTCTCAGACTTTGATAACAAGTGTCCTAGTTGGAGCCGTGTGCAGAAGCTGAAAGGAGCGCT
GGCTGCCTACAACAAGGGTACAGAAAAAATGGATTACAACAACGTGGACGGTCACACTACTTGGGGGAAC
TACTCCACTGATGTCCTTGCCAGGGCCGAGTTCTTTAAGAGGAATGGGTA





##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##========================================== Now make a synteny plot  #1  ============================================================
##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================


#N5.HOG0019111.clade_hgt.txt=> clade where the HGT is extracted from iTOL (from N5.HOG0019111.prot.aln.treefile)

cut -f1,2 -d "_" N5.HOG0019111.clade_hgt.txt | sort | uniq  > species_to_draw.clade1
#Choose close species without the gene .. 



echo "Engraulis_encrasicolus" > non_transfer_species_1
echo "Alosa_alosa" >> non_transfer_species_1
echo "Sardina_pilchardus" >> non_transfer_species_1
echo "Osmerus_eperlanus" > non_transfer_species_2
echo "Borostomias_antarcticus" >> non_transfer_species_2
echo "Gadus_morhua" >> non_transfer_species_2



rm -r Syn_tables_dir ; mkdir Syn_tables_dir


#First order the species names 

rm species_to_draw.clade1.ordered

echo "Clupea_harengus" >> species_to_draw.clade1.ordered
echo "Hypomesus_transpacificus" >> species_to_draw.clade1.ordered


#If there are more than 4 species from the same order, keep only 4 


#First extract genes around 

curr_OGG=N5.HOG0019111

for curr_sp in `cat species_to_draw.clade1.ordered` ; do

  grep "$curr_sp" $curr_OGG.clade_hgt.txt > Syn_tables_dir/gene_names.$curr_sp.txt
  
  echo "" > Syn_tables_dir/$curr_sp.synt.df
  
  for gene in `cat Syn_tables_dir/gene_names.$curr_sp.txt` ; do 

    gene_name=`echo "$gene" | sed "s/$curr_sp\_//g"`
  
    if grep -q "$gene_name" Syn_tables_dir/$curr_sp.synt.df ; then 
      echo "gene already in the table"
    else
      
      grep -A10 -B10 "$gene_name" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df
    
    fi
  
  done

  if [ $curr_sp == "Clupea_harengus" ] ; then
  	grep -A7 -B7 "N5_HOG0046495" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df
  fi


  sed -i '/^$/d' Syn_tables_dir/$curr_sp.synt.df

  awk '{ if (NF < 5 || $5 == "") $5 = "no_OGG"; print $0;}' Syn_tables_dir/$curr_sp.synt.df > Syn_tables_dir/$curr_sp.synt.final.df
  

  sed -i 's/	/,/g' Syn_tables_dir/$curr_sp.synt.final.df
  sed -i 's/ no_OGG/,no_OGG/g' Syn_tables_dir/$curr_sp.synt.final.df
  sed -i 's/ /,/g' Syn_tables_dir/$curr_sp.synt.final.df

  #Add the strand of each gene ! 

  rm Syn_tables_dir/$curr_sp.synt.final.df.strand
  for syn_line in `cat Syn_tables_dir/$curr_sp.synt.final.df` ; do 

    curr_name=`echo "$syn_line" | cut -f4 -d "," | sed 's/rna.//g' | sed 's/..$//g' | sed 's/_$//g'`

    if grep -q "$curr_name" GFF3_files_per_species/$curr_sp.gff ; then
      curr_strand=`grep "$curr_name" GFF3_files_per_species/$curr_sp.gff | cut -f7 | sort | uniq | head -1`
    else 
      curr_strand="+"
    fi

    echo "$syn_line,$curr_strand,$curr_sp" >> Syn_tables_dir/$curr_sp.synt.final.df.strand
  
  done
  mv Syn_tables_dir/$curr_sp.synt.final.df.strand Syn_tables_dir/$curr_sp.synt.final.df

  rm Syn_tables_dir/gene_names.$curr_sp.txt

done



#First add Sardina_pilchardus

curr_OGG=N5.HOG0019111
curr_sp=Sardina_pilchardus
ref_sp=Clupea_harengus

grep -A10 -B10 "XM_062535392" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0015783" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df


awk '{ if (NF < 5 || $5 == "") $5 = "no_OGG"; print $0;}' Syn_tables_dir/$curr_sp.synt.df > Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/	/,/g' Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/ no_OGG/,no_OGG/g' Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/ /,/g' Syn_tables_dir/$curr_sp.synt.final.df
grep -v ",,,," Syn_tables_dir/$curr_sp.synt.final.df > temp ; mv temp Syn_tables_dir/$curr_sp.synt.final.df


rm Syn_tables_dir/$curr_sp.synt.final.df.strand
for syn_line in `cat Syn_tables_dir/$curr_sp.synt.final.df` ; do 

  curr_name=`echo "$syn_line" | cut -f4 -d "," | sed 's/rna.//g' | sed 's/..$//g' | sed 's/_$//g'`
  
  if grep -q "$curr_name" GFF3_files_per_species/$curr_sp.gff ; then
    curr_strand=`grep "$curr_name" GFF3_files_per_species/$curr_sp.gff | cut -f7 | sort | uniq | head -1`
  else 
    curr_strand="+"
  fi

  echo "$syn_line,$curr_strand,$curr_sp" >> Syn_tables_dir/$curr_sp.synt.final.df.strand

done
mv Syn_tables_dir/$curr_sp.synt.final.df.strand Syn_tables_dir/$curr_sp.synt.final.df




#First add Alosa_alosa

curr_OGG=N5.HOG0019111
curr_sp=Alosa_alosa
ref_sp=Clupea_harengus

grep -A10 -B10 "XM_048231371" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0040207" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0023643" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df


awk '{ if (NF < 5 || $5 == "") $5 = "no_OGG"; print $0;}' Syn_tables_dir/$curr_sp.synt.df > Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/	/,/g' Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/ no_OGG/,no_OGG/g' Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/ /,/g' Syn_tables_dir/$curr_sp.synt.final.df
grep -v ",,,," Syn_tables_dir/$curr_sp.synt.final.df > temp ; mv temp Syn_tables_dir/$curr_sp.synt.final.df


rm Syn_tables_dir/$curr_sp.synt.final.df.strand
for syn_line in `cat Syn_tables_dir/$curr_sp.synt.final.df` ; do 

  curr_name=`echo "$syn_line" | cut -f4 -d "," | sed 's/rna.//g' | sed 's/..$//g' | sed 's/_$//g'`
  
  if grep -q "$curr_name" GFF3_files_per_species/$curr_sp.gff ; then
    curr_strand=`grep "$curr_name" GFF3_files_per_species/$curr_sp.gff | cut -f7 | sort | uniq | head -1`
  else 
    curr_strand="+"
  fi

  echo "$syn_line,$curr_strand,$curr_sp" >> Syn_tables_dir/$curr_sp.synt.final.df.strand

done
mv Syn_tables_dir/$curr_sp.synt.final.df.strand Syn_tables_dir/$curr_sp.synt.final.df




#Now add Engraulis_encrasicolus

curr_OGG=N5.HOG0019111
curr_sp=Engraulis_encrasicolus
ref_sp=Clupea_harengus

grep -A10 -B10 "N5_HOG0032487" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0032111" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df
grep -A5 -B5 "N5_HOG0006299" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df
grep -A5 -B5 "N5_HOG0036590" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df



awk '{ if (NF < 5 || $5 == "") $5 = "no_OGG"; print $0;}' Syn_tables_dir/$curr_sp.synt.df > Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/	/,/g' Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/ no_OGG/,no_OGG/g' Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/ /,/g' Syn_tables_dir/$curr_sp.synt.final.df
grep -v ",,,," Syn_tables_dir/$curr_sp.synt.final.df > temp ; mv temp Syn_tables_dir/$curr_sp.synt.final.df


rm Syn_tables_dir/$curr_sp.synt.final.df.strand
for syn_line in `cat Syn_tables_dir/$curr_sp.synt.final.df` ; do 

  curr_name=`echo "$syn_line" | cut -f4 -d "," | sed 's/rna.//g' | sed 's/..$//g' | sed 's/_$//g'`
  
  if grep -q "$curr_name" GFF3_files_per_species/$curr_sp.gff ; then
    curr_strand=`grep "$curr_name" GFF3_files_per_species/$curr_sp.gff | cut -f7 | sort | uniq | head -1`
  else 
    curr_strand="+"
  fi

  echo "$syn_line,$curr_strand,$curr_sp" >> Syn_tables_dir/$curr_sp.synt.final.df.strand

done
mv Syn_tables_dir/$curr_sp.synt.final.df.strand Syn_tables_dir/$curr_sp.synt.final.df




#Now add Osmerus_eperlanus

curr_OGG=N5.HOG0019111
curr_sp=Osmerus_eperlanus
ref_sp=Hypomesus_transpacificus

grep -A25 -B10 "N5_HOG0015817" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df


awk '{ if (NF < 5 || $5 == "") $5 = "no_OGG"; print $0;}' Syn_tables_dir/$curr_sp.synt.df > Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/	/,/g' Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/ no_OGG/,no_OGG/g' Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/ /,/g' Syn_tables_dir/$curr_sp.synt.final.df
grep -v ",,,," Syn_tables_dir/$curr_sp.synt.final.df > temp ; mv temp Syn_tables_dir/$curr_sp.synt.final.df


rm Syn_tables_dir/$curr_sp.synt.final.df.strand
for syn_line in `cat Syn_tables_dir/$curr_sp.synt.final.df` ; do 

  curr_name=`echo "$syn_line" | cut -f4 -d "," | sed 's/rna.//g' | sed 's/..$//g' | sed 's/_$//g'`
  
  if grep -q "$curr_name" GFF3_files_per_species/$curr_sp.gff ; then
    curr_strand=`grep "$curr_name" GFF3_files_per_species/$curr_sp.gff | cut -f7 | sort | uniq | head -1`
  else 
    curr_strand="+"
  fi

  echo "$syn_line,$curr_strand,$curr_sp" >> Syn_tables_dir/$curr_sp.synt.final.df.strand

done
mv Syn_tables_dir/$curr_sp.synt.final.df.strand Syn_tables_dir/$curr_sp.synt.final.df



#Now add Borostomias_antarcticus

curr_OGG=N5.HOG0019111
curr_sp=Borostomias_antarcticus
ref_sp=Hypomesus_transpacificus

grep -A10 -B15 "N5_HOG0015817" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0024350" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df


awk '{ if (NF < 5 || $5 == "") $5 = "no_OGG"; print $0;}' Syn_tables_dir/$curr_sp.synt.df > Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/	/,/g' Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/ no_OGG/,no_OGG/g' Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/ /,/g' Syn_tables_dir/$curr_sp.synt.final.df
grep -v ",,,," Syn_tables_dir/$curr_sp.synt.final.df > temp ; mv temp Syn_tables_dir/$curr_sp.synt.final.df


rm Syn_tables_dir/$curr_sp.synt.final.df.strand
for syn_line in `cat Syn_tables_dir/$curr_sp.synt.final.df` ; do 

  curr_name=`echo "$syn_line" | cut -f4 -d "," | sed 's/rna.//g' | sed 's/..$//g' | sed 's/_$//g'`
  
  if grep -q "$curr_name" GFF3_files_per_species/$curr_sp.gff ; then
    curr_strand=`grep "$curr_name" GFF3_files_per_species/$curr_sp.gff | cut -f7 | sort | uniq | head -1`
  else 
    curr_strand="+"
  fi

  echo "$syn_line,$curr_strand,$curr_sp" >> Syn_tables_dir/$curr_sp.synt.final.df.strand

done
mv Syn_tables_dir/$curr_sp.synt.final.df.strand Syn_tables_dir/$curr_sp.synt.final.df



#Finally add Gadus_morhua

curr_OGG=N5.HOG0019111
curr_sp=Gadus_morhua
ref_sp=Hypomesus_transpacificus

grep -A7 -B7 "N5_HOG0015817" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df
grep -A7 -B7 "N5_HOG0023474" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df
grep -A7 -B7 "N5_HOG0031439" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df
grep -A9 -B9 "N5_HOG0024350" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df


awk '{ if (NF < 5 || $5 == "") $5 = "no_OGG"; print $0;}' Syn_tables_dir/$curr_sp.synt.df > Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/	/,/g' Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/ no_OGG/,no_OGG/g' Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/ /,/g' Syn_tables_dir/$curr_sp.synt.final.df
grep -v ",,,," Syn_tables_dir/$curr_sp.synt.final.df > temp ; mv temp Syn_tables_dir/$curr_sp.synt.final.df


rm Syn_tables_dir/$curr_sp.synt.final.df.strand
for syn_line in `cat Syn_tables_dir/$curr_sp.synt.final.df` ; do 

  curr_name=`echo "$syn_line" | cut -f4 -d "," | sed 's/rna.//g' | sed 's/..$//g' | sed 's/_$//g'`
  
  if grep -q "$curr_name" GFF3_files_per_species/$curr_sp.gff ; then
    curr_strand=`grep "$curr_name" GFF3_files_per_species/$curr_sp.gff | cut -f7 | sort | uniq | head -1`
  else 
    curr_strand="+"
  fi

  echo "$syn_line,$curr_strand,$curr_sp" >> Syn_tables_dir/$curr_sp.synt.final.df.strand

done
mv Syn_tables_dir/$curr_sp.synt.final.df.strand Syn_tables_dir/$curr_sp.synt.final.df




#Lets create a seq info ogg file for both clades


cat non_transfer_species_1 species_to_draw.clade1.ordered non_transfer_species_2 > temp ; mv temp species_to_draw.clade1.ordered


rm seq_clustered_infos_ogg.clade1.csv
for curr_sp in `cat species_to_draw.clade1.ordered` ; do 
  cat Syn_tables_dir/$curr_sp.synt.final.df >> seq_clustered_infos_ogg.clade1.csv
done



sed -i 's/ /,/g' seq_clustered_infos_ogg.clade1.csv


### Now lets crate a seqID file that list all the clusters and their length per species, and arrange it for the graph representation
### We also create a good geneID file 

cut -f1,7 -d "," seq_clustered_infos_ogg.clade1.csv  | sort | uniq > uniq_clusters.clade1


for species in `cat species_to_draw.clade1.ordered` ; do 
  grep "$species" uniq_clusters.clade1 >> temp.txt
done ; mv temp.txt uniq_clusters.clade1




rm cluster_id_arranged.clade1.txt
rm seq_clustered_infos_ogg_num.clade1.csv
for cluster in `cat uniq_clusters.clade1` ; do 

  scaffold=`echo "$cluster" | cut -f1 -d ","`
  species=`echo "$cluster" | cut -f2 -d ","`

  start=`grep "$scaffold.*$species" seq_clustered_infos_ogg.clade1.csv | cut -f2 -d "," | sort -n | head -1`
  end=`grep "$scaffold.*$species" seq_clustered_infos_ogg.clade1.csv | cut -f3 -d "," | sort -n | tail -1`
  length=$((end - start))

  echo "$species,$scaffold,$length" >> cluster_id_arranged.clade1.txt

  start_min_1=$((start - 1))
  grep "$scaffold.*$species" seq_clustered_infos_ogg.clade1.csv | awk -F, -v start_min_1="$start_min_1" '{$2 = $2 - start_min_1; print}' | awk -F" " -v start_min_1="$start_min_1" '{$3 = $3 - start_min_1; print}' | sed 's/ /,/g'  >> seq_clustered_infos_ogg_num.clade1.csv

done 
 



##Now create a link table for clade 1




i=1
j=`wc -l  < species_to_draw.clade1.ordered`
rm Species_Comparison_table.clade1.txt
while [ $i -lt $j ] ; do 

  n=$(( i + 1 )) 


  sp1=`head -$i species_to_draw.clade1.ordered | tail -1`
  sp2=`head -$n species_to_draw.clade1.ordered | tail -1`

  echo "$sp1,$sp2" >> Species_Comparison_table.clade1.txt

  i=$(( i + 1 )) 


done





rm link_table.clade1.txt
for comparisons in `cat Species_Comparison_table.clade1.txt` ; do  #Species_Comparison_table.txt = table with two columns with each pair of species being compared

  species1=`echo "$comparisons" | cut -f1 -d ","`
  species2=`echo "$comparisons" | cut -f2 -d ","`

  grep "$species1" seq_clustered_infos_ogg_num.clade1.csv  | cut -f5 -d "," | sort | uniq > uniq_OGG_sp1
  grep "$species2" seq_clustered_infos_ogg_num.clade1.csv  | cut -f5 -d "," | sort | uniq > uniq_OGG_sp2
  cat uniq_OGG_sp1 uniq_OGG_sp2 | sort | uniq -c | sed 's/^ *//g' | grep -v "^1 " | grep -v "no_OGG" | cut -f2 -d " " | sort | uniq > OGG_list.txt 

  for OGG in `cat OGG_list.txt` ; do

    grep "$OGG" seq_clustered_infos_ogg_num.clade1.csv > current_OGG.txt 
    grep "$species1" current_OGG.txt | cut -f7,1,2,3 -d "," | awk -v OFS="," -F"," '{print($4,$1,$2,$3)}' > species1.ogg.txt
    grep "$species2" current_OGG.txt | cut -f7,1,2,3 -d "," | awk -v OFS="," -F"," '{print($4,$1,$2,$3)}' > species2.ogg.txt

    for line_species1 in `cat species1.ogg.txt` ; do

      for line_species2 in `cat species2.ogg.txt` ; do 

        echo "$line_species1,$line_species2,$OGG" >> link_table.clade1.txt


      done
    done

    rm species1.ogg.txt ; rm species2.ogg.txt

  done
done


mv cluster_id_arranged.clade1.txt cluster_id_arranged.clade1.ClupOsm.txt
mv seq_clustered_infos_ogg_num.clade1.csv seq_clustered_infos_ogg_num.clade1.Cluposm.csv
mv link_table.clade1.txt link_table.clade1.ClupOsm.txt



##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##========================================== Now make a synteny plot  #2  ============================================================
##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================


#N5.HOG0019111.clade_hgt.txt=> clade where the HGT is extracted from iTOL (from N5.HOG0019111.prot.aln.treefile)

cut -f1,2 -d "_" N5.HOG0019111.clade_hgt.txt | sort | uniq  > species_to_draw.clade1
#Choose close species without the gene .. 



echo "Engraulis_encrasicolus" > non_transfer_species_1
echo "Alosa_alosa" >> non_transfer_species_1
echo "Sardina_pilchardus" >> non_transfer_species_1

echo "Scomber_scombrus" > non_transfer_species_2
echo "Thalassophryne_amazonica" >> non_transfer_species_2
echo "Centroberyx_gerrardi" >> non_transfer_species_2



rm -r Syn_tables_dir ; mkdir Syn_tables_dir


#First order the species names 

rm species_to_draw.clade1.ordered

echo "Clupea_harengus" >> species_to_draw.clade1.ordered
echo "Thunnus_maccoyii" >> species_to_draw.clade1.ordered


#If there are more than 4 species from the same order, keep only 4 


#First extract genes around 

curr_OGG=N5.HOG0019111

for curr_sp in `cat species_to_draw.clade1.ordered` ; do

  grep "$curr_sp" $curr_OGG.clade_hgt.txt > Syn_tables_dir/gene_names.$curr_sp.txt
  
  echo "" > Syn_tables_dir/$curr_sp.synt.df
  
  for gene in `cat Syn_tables_dir/gene_names.$curr_sp.txt` ; do 

    gene_name=`echo "$gene" | sed "s/$curr_sp\_//g"`
  
    if grep -q "$gene_name" Syn_tables_dir/$curr_sp.synt.df ; then 
      echo "gene already in the table"
    else
      
      if [ $curr_sp == "Thunnus_maccoyii" ] ; then
      	grep -A15 -B15 "$gene_name" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df
      else

      	grep -A10 -B10 "$gene_name" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df
   	  fi
  
    fi
  
  done



  sed -i '/^$/d' Syn_tables_dir/$curr_sp.synt.df

  awk '{ if (NF < 5 || $5 == "") $5 = "no_OGG"; print $0;}' Syn_tables_dir/$curr_sp.synt.df > Syn_tables_dir/$curr_sp.synt.final.df
  

  sed -i 's/	/,/g' Syn_tables_dir/$curr_sp.synt.final.df
  sed -i 's/ no_OGG/,no_OGG/g' Syn_tables_dir/$curr_sp.synt.final.df
  sed -i 's/ /,/g' Syn_tables_dir/$curr_sp.synt.final.df

  #Add the strand of each gene ! 

  rm Syn_tables_dir/$curr_sp.synt.final.df.strand
  for syn_line in `cat Syn_tables_dir/$curr_sp.synt.final.df` ; do 

    curr_name=`echo "$syn_line" | cut -f4 -d "," | sed 's/rna.//g' | sed 's/..$//g' | sed 's/_$//g'`

    if grep -q "$curr_name" GFF3_files_per_species/$curr_sp.gff ; then
      curr_strand=`grep "$curr_name" GFF3_files_per_species/$curr_sp.gff | cut -f7 | sort | uniq | head -1`
    else 
      curr_strand="+"
    fi

    echo "$syn_line,$curr_strand,$curr_sp" >> Syn_tables_dir/$curr_sp.synt.final.df.strand
  
  done
  mv Syn_tables_dir/$curr_sp.synt.final.df.strand Syn_tables_dir/$curr_sp.synt.final.df

  rm Syn_tables_dir/gene_names.$curr_sp.txt

done



#First add Sardina_pilchardus

curr_OGG=N5.HOG0019111
curr_sp=Sardina_pilchardus
ref_sp=Clupea_harengus

grep -A10 -B10 "XM_062535392" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0015783" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df


awk '{ if (NF < 5 || $5 == "") $5 = "no_OGG"; print $0;}' Syn_tables_dir/$curr_sp.synt.df > Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/	/,/g' Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/ no_OGG/,no_OGG/g' Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/ /,/g' Syn_tables_dir/$curr_sp.synt.final.df
grep -v ",,,," Syn_tables_dir/$curr_sp.synt.final.df > temp ; mv temp Syn_tables_dir/$curr_sp.synt.final.df


rm Syn_tables_dir/$curr_sp.synt.final.df.strand
for syn_line in `cat Syn_tables_dir/$curr_sp.synt.final.df` ; do 

  curr_name=`echo "$syn_line" | cut -f4 -d "," | sed 's/rna.//g' | sed 's/..$//g' | sed 's/_$//g'`
  
  if grep -q "$curr_name" GFF3_files_per_species/$curr_sp.gff ; then
    curr_strand=`grep "$curr_name" GFF3_files_per_species/$curr_sp.gff | cut -f7 | sort | uniq | head -1`
  else 
    curr_strand="+"
  fi

  echo "$syn_line,$curr_strand,$curr_sp" >> Syn_tables_dir/$curr_sp.synt.final.df.strand

done
mv Syn_tables_dir/$curr_sp.synt.final.df.strand Syn_tables_dir/$curr_sp.synt.final.df




#First add Alosa_alosa

curr_OGG=N5.HOG0019111
curr_sp=Alosa_alosa
ref_sp=Clupea_harengus

grep -A10 -B10 "XM_048231371" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0040207" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0023643" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df


awk '{ if (NF < 5 || $5 == "") $5 = "no_OGG"; print $0;}' Syn_tables_dir/$curr_sp.synt.df > Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/	/,/g' Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/ no_OGG/,no_OGG/g' Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/ /,/g' Syn_tables_dir/$curr_sp.synt.final.df
grep -v ",,,," Syn_tables_dir/$curr_sp.synt.final.df > temp ; mv temp Syn_tables_dir/$curr_sp.synt.final.df


rm Syn_tables_dir/$curr_sp.synt.final.df.strand
for syn_line in `cat Syn_tables_dir/$curr_sp.synt.final.df` ; do 

  curr_name=`echo "$syn_line" | cut -f4 -d "," | sed 's/rna.//g' | sed 's/..$//g' | sed 's/_$//g'`
  
  if grep -q "$curr_name" GFF3_files_per_species/$curr_sp.gff ; then
    curr_strand=`grep "$curr_name" GFF3_files_per_species/$curr_sp.gff | cut -f7 | sort | uniq | head -1`
  else 
    curr_strand="+"
  fi

  echo "$syn_line,$curr_strand,$curr_sp" >> Syn_tables_dir/$curr_sp.synt.final.df.strand

done
mv Syn_tables_dir/$curr_sp.synt.final.df.strand Syn_tables_dir/$curr_sp.synt.final.df




#Now add Engraulis_encrasicolus

curr_OGG=N5.HOG0019111
curr_sp=Engraulis_encrasicolus
ref_sp=Clupea_harengus

grep -A10 -B10 "N5_HOG0032487" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0032111" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df
grep -A5 -B5 "N5_HOG0006299" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df
grep -A5 -B5 "N5_HOG0036590" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df



awk '{ if (NF < 5 || $5 == "") $5 = "no_OGG"; print $0;}' Syn_tables_dir/$curr_sp.synt.df > Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/	/,/g' Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/ no_OGG/,no_OGG/g' Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/ /,/g' Syn_tables_dir/$curr_sp.synt.final.df
grep -v ",,,," Syn_tables_dir/$curr_sp.synt.final.df > temp ; mv temp Syn_tables_dir/$curr_sp.synt.final.df


rm Syn_tables_dir/$curr_sp.synt.final.df.strand
for syn_line in `cat Syn_tables_dir/$curr_sp.synt.final.df` ; do 

  curr_name=`echo "$syn_line" | cut -f4 -d "," | sed 's/rna.//g' | sed 's/..$//g' | sed 's/_$//g'`
  
  if grep -q "$curr_name" GFF3_files_per_species/$curr_sp.gff ; then
    curr_strand=`grep "$curr_name" GFF3_files_per_species/$curr_sp.gff | cut -f7 | sort | uniq | head -1`
  else 
    curr_strand="+"
  fi

  echo "$syn_line,$curr_strand,$curr_sp" >> Syn_tables_dir/$curr_sp.synt.final.df.strand

done
mv Syn_tables_dir/$curr_sp.synt.final.df.strand Syn_tables_dir/$curr_sp.synt.final.df




#Now add Scomber_scombrus

curr_OGG=N5.HOG0019111
curr_sp=Scomber_scombrus
ref_sp=Thunnus_maccoyii

grep -A10 -B10 "HOG0019111" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "NC_084971" > Syn_tables_dir/$curr_sp.synt.df


awk '{ if (NF < 5 || $5 == "") $5 = "no_OGG"; print $0;}' Syn_tables_dir/$curr_sp.synt.df > Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/	/,/g' Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/ no_OGG/,no_OGG/g' Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/ /,/g' Syn_tables_dir/$curr_sp.synt.final.df
grep -v ",,,," Syn_tables_dir/$curr_sp.synt.final.df > temp ; mv temp Syn_tables_dir/$curr_sp.synt.final.df


rm Syn_tables_dir/$curr_sp.synt.final.df.strand
for syn_line in `cat Syn_tables_dir/$curr_sp.synt.final.df` ; do 

  curr_name=`echo "$syn_line" | cut -f4 -d "," | sed 's/rna.//g' | sed 's/..$//g' | sed 's/_$//g'`
  
  if grep -q "$curr_name" GFF3_files_per_species/$curr_sp.gff ; then
    curr_strand=`grep "$curr_name" GFF3_files_per_species/$curr_sp.gff | cut -f7 | sort | uniq | head -1`
  else 
    curr_strand="+"
  fi

  echo "$syn_line,$curr_strand,$curr_sp" >> Syn_tables_dir/$curr_sp.synt.final.df.strand

done
mv Syn_tables_dir/$curr_sp.synt.final.df.strand Syn_tables_dir/$curr_sp.synt.final.df





#Now add Thalassophryne_amazonica

curr_OGG=N5.HOG0019111
curr_sp=Thalassophryne_amazonica
ref_sp=Thunnus_maccoyii

grep -A15 -B15 "N5_HOG0019111" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "NC_047117" > Syn_tables_dir/$curr_sp.synt.df

awk '{ if (NF < 5 || $5 == "") $5 = "no_OGG"; print $0;}' Syn_tables_dir/$curr_sp.synt.df > Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/	/,/g' Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/ no_OGG/,no_OGG/g' Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/ /,/g' Syn_tables_dir/$curr_sp.synt.final.df
grep -v ",,,," Syn_tables_dir/$curr_sp.synt.final.df > temp ; mv temp Syn_tables_dir/$curr_sp.synt.final.df


rm Syn_tables_dir/$curr_sp.synt.final.df.strand
for syn_line in `cat Syn_tables_dir/$curr_sp.synt.final.df` ; do 

  curr_name=`echo "$syn_line" | cut -f4 -d "," | sed 's/rna.//g' | sed 's/..$//g' | sed 's/_$//g'`
  
  if grep -q "$curr_name" GFF3_files_per_species/$curr_sp.gff ; then
    curr_strand=`grep "$curr_name" GFF3_files_per_species/$curr_sp.gff | cut -f7 | sort | uniq | head -1`
  else 
    curr_strand="+"
  fi

  echo "$syn_line,$curr_strand,$curr_sp" >> Syn_tables_dir/$curr_sp.synt.final.df.strand

done
mv Syn_tables_dir/$curr_sp.synt.final.df.strand Syn_tables_dir/$curr_sp.synt.final.df






#Finally add  Centroberyx_gerrardi

curr_OGG=N5.HOG0019111
curr_sp=Centroberyx_gerrardi
ref_sp=Thunnus_maccoyii

grep -A8 -B8 "N5_HOG0019111" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "JAPMTB010005717" > Syn_tables_dir/$curr_sp.synt.df
grep -A8 -B9 "N5_HOG0016604" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df
grep -A8 -B9 "N5_HOG0040864" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df



awk '{ if (NF < 5 || $5 == "") $5 = "no_OGG"; print $0;}' Syn_tables_dir/$curr_sp.synt.df > Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/	/,/g' Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/ no_OGG/,no_OGG/g' Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/ /,/g' Syn_tables_dir/$curr_sp.synt.final.df
grep -v ",,,," Syn_tables_dir/$curr_sp.synt.final.df > temp ; mv temp Syn_tables_dir/$curr_sp.synt.final.df


rm Syn_tables_dir/$curr_sp.synt.final.df.strand
for syn_line in `cat Syn_tables_dir/$curr_sp.synt.final.df` ; do 

  curr_name=`echo "$syn_line" | cut -f4 -d "," | sed 's/rna.//g' | sed 's/..$//g' | sed 's/_$//g'`
  
  if grep -q "$curr_name" GFF3_files_per_species/$curr_sp.gff ; then
    curr_strand=`grep "$curr_name" GFF3_files_per_species/$curr_sp.gff | cut -f7 | sort | uniq | head -1`
  else 
    curr_strand="+"
  fi

  echo "$syn_line,$curr_strand,$curr_sp" >> Syn_tables_dir/$curr_sp.synt.final.df.strand

done
mv Syn_tables_dir/$curr_sp.synt.final.df.strand Syn_tables_dir/$curr_sp.synt.final.df




#Lets create a seq info ogg file for both clades


cat non_transfer_species_1 species_to_draw.clade1.ordered non_transfer_species_2 > temp ; mv temp species_to_draw.clade1.ordered


rm seq_clustered_infos_ogg.clade1.csv
for curr_sp in `cat species_to_draw.clade1.ordered` ; do 
  cat Syn_tables_dir/$curr_sp.synt.final.df >> seq_clustered_infos_ogg.clade1.csv
done



sed -i 's/ /,/g' seq_clustered_infos_ogg.clade1.csv


### Now lets crate a seqID file that list all the clusters and their length per species, and arrange it for the graph representation
### We also create a good geneID file 

cut -f1,7 -d "," seq_clustered_infos_ogg.clade1.csv  | sort | uniq > uniq_clusters.clade1


for species in `cat species_to_draw.clade1.ordered` ; do 
  grep "$species" uniq_clusters.clade1 >> temp.txt
done ; mv temp.txt uniq_clusters.clade1




rm cluster_id_arranged.clade1.txt
rm seq_clustered_infos_ogg_num.clade1.csv
for cluster in `cat uniq_clusters.clade1` ; do 

  scaffold=`echo "$cluster" | cut -f1 -d ","`
  species=`echo "$cluster" | cut -f2 -d ","`

  start=`grep "$scaffold.*$species" seq_clustered_infos_ogg.clade1.csv | cut -f2 -d "," | sort -n | head -1`
  end=`grep "$scaffold.*$species" seq_clustered_infos_ogg.clade1.csv | cut -f3 -d "," | sort -n | tail -1`
  length=$((end - start))

  echo "$species,$scaffold,$length" >> cluster_id_arranged.clade1.txt

  start_min_1=$((start - 1))
  grep "$scaffold.*$species" seq_clustered_infos_ogg.clade1.csv | awk -F, -v start_min_1="$start_min_1" '{$2 = $2 - start_min_1; print}' | awk -F" " -v start_min_1="$start_min_1" '{$3 = $3 - start_min_1; print}' | sed 's/ /,/g'  >> seq_clustered_infos_ogg_num.clade1.csv

done 
 



##Now create a link table for clade 1




i=1
j=`wc -l  < species_to_draw.clade1.ordered`
rm Species_Comparison_table.clade1.txt
while [ $i -lt $j ] ; do 

  n=$(( i + 1 )) 


  sp1=`head -$i species_to_draw.clade1.ordered | tail -1`
  sp2=`head -$n species_to_draw.clade1.ordered | tail -1`

  echo "$sp1,$sp2" >> Species_Comparison_table.clade1.txt

  i=$(( i + 1 )) 


done





rm link_table.clade1.txt
for comparisons in `cat Species_Comparison_table.clade1.txt` ; do  #Species_Comparison_table.txt = table with two columns with each pair of species being compared

  species1=`echo "$comparisons" | cut -f1 -d ","`
  species2=`echo "$comparisons" | cut -f2 -d ","`

  grep "$species1" seq_clustered_infos_ogg_num.clade1.csv  | cut -f5 -d "," | sort | uniq > uniq_OGG_sp1
  grep "$species2" seq_clustered_infos_ogg_num.clade1.csv  | cut -f5 -d "," | sort | uniq > uniq_OGG_sp2
  cat uniq_OGG_sp1 uniq_OGG_sp2 | sort | uniq -c | sed 's/^ *//g' | grep -v "^1 " | grep -v "no_OGG" | cut -f2 -d " " | sort | uniq > OGG_list.txt 

  for OGG in `cat OGG_list.txt` ; do

    grep "$OGG" seq_clustered_infos_ogg_num.clade1.csv > current_OGG.txt 
    grep "$species1" current_OGG.txt | cut -f7,1,2,3 -d "," | awk -v OFS="," -F"," '{print($4,$1,$2,$3)}' > species1.ogg.txt
    grep "$species2" current_OGG.txt | cut -f7,1,2,3 -d "," | awk -v OFS="," -F"," '{print($4,$1,$2,$3)}' > species2.ogg.txt

    for line_species1 in `cat species1.ogg.txt` ; do

      for line_species2 in `cat species2.ogg.txt` ; do 

        echo "$line_species1,$line_species2,$OGG" >> link_table.clade1.txt


      done
    done

    rm species1.ogg.txt ; rm species2.ogg.txt

  done
done



mv cluster_id_arranged.clade1.txt cluster_id_arranged.clade1.ClupSco.txt
mv seq_clustered_infos_ogg_num.clade1.csv seq_clustered_infos_ogg_num.clade1.ClupSco.csv
mv link_table.clade1.txt link_table.clade1.ClupSco.txt





##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##========================================== Now make a synteny plot  #3  ============================================================
##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================


#N5.HOG0019111.clade_hgt.txt=> clade where the HGT is extracted from iTOL (from N5.HOG0019111.prot.aln.treefile)

cut -f1,2 -d "_" N5.HOG0019111.clade_hgt.txt | sort | uniq  > species_to_draw.clade1
#Choose close species without the gene .. 



echo "Gadus_morhua" > non_transfer_species_1
echo "Borostomias_antarcticus" >> non_transfer_species_1
echo "Osmerus_eperlanus" >> non_transfer_species_1

echo "Scomber_scombrus" > non_transfer_species_2
echo "Thalassophryne_amazonica" >> non_transfer_species_2
echo "Centroberyx_gerrardi" >> non_transfer_species_2



rm -r Syn_tables_dir ; mkdir Syn_tables_dir


#First order the species names 

rm species_to_draw.clade1.ordered

echo "Hypomesus_transpacificus" >> species_to_draw.clade1.ordered
echo "Thunnus_maccoyii" >> species_to_draw.clade1.ordered


#If there are more than 4 species from the same order, keep only 4 


#First extract genes around 

curr_OGG=N5.HOG0019111

for curr_sp in `cat species_to_draw.clade1.ordered` ; do

  grep "$curr_sp" $curr_OGG.clade_hgt.txt > Syn_tables_dir/gene_names.$curr_sp.txt
  
  echo "" > Syn_tables_dir/$curr_sp.synt.df
  
  for gene in `cat Syn_tables_dir/gene_names.$curr_sp.txt` ; do 

    gene_name=`echo "$gene" | sed "s/$curr_sp\_//g"`
  
    if grep -q "$gene_name" Syn_tables_dir/$curr_sp.synt.df ; then 
      echo "gene already in the table"
    else
      
      if [ $curr_sp == "Thunnus_maccoyii" ] ; then
      	grep -A15 -B15 "$gene_name" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df
      else

      	grep -A10 -B10 "$gene_name" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df
   	  fi
  
    fi
  
  done


  if [ $curr_sp == "Clupea_harengus" ] ; then
  	grep -A7 -B7 "N5_HOG0046495" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df
  fi



  sed -i '/^$/d' Syn_tables_dir/$curr_sp.synt.df

  awk '{ if (NF < 5 || $5 == "") $5 = "no_OGG"; print $0;}' Syn_tables_dir/$curr_sp.synt.df > Syn_tables_dir/$curr_sp.synt.final.df
  

  sed -i 's/	/,/g' Syn_tables_dir/$curr_sp.synt.final.df
  sed -i 's/ no_OGG/,no_OGG/g' Syn_tables_dir/$curr_sp.synt.final.df
  sed -i 's/ /,/g' Syn_tables_dir/$curr_sp.synt.final.df

  #Add the strand of each gene ! 

  rm Syn_tables_dir/$curr_sp.synt.final.df.strand
  for syn_line in `cat Syn_tables_dir/$curr_sp.synt.final.df` ; do 

    curr_name=`echo "$syn_line" | cut -f4 -d "," | sed 's/rna.//g' | sed 's/..$//g' | sed 's/_$//g'`

    if grep -q "$curr_name" GFF3_files_per_species/$curr_sp.gff ; then
      curr_strand=`grep "$curr_name" GFF3_files_per_species/$curr_sp.gff | cut -f7 | sort | uniq | head -1`
    else 
      curr_strand="+"
    fi

    echo "$syn_line,$curr_strand,$curr_sp" >> Syn_tables_dir/$curr_sp.synt.final.df.strand
  
  done
  mv Syn_tables_dir/$curr_sp.synt.final.df.strand Syn_tables_dir/$curr_sp.synt.final.df

  rm Syn_tables_dir/gene_names.$curr_sp.txt

done


#Now add Osmerus_eperlanus

curr_OGG=N5.HOG0019111
curr_sp=Osmerus_eperlanus
ref_sp=Hypomesus_transpacificus

grep -A25 -B10 "N5_HOG0015817" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df


awk '{ if (NF < 5 || $5 == "") $5 = "no_OGG"; print $0;}' Syn_tables_dir/$curr_sp.synt.df > Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/	/,/g' Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/ no_OGG/,no_OGG/g' Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/ /,/g' Syn_tables_dir/$curr_sp.synt.final.df
grep -v ",,,," Syn_tables_dir/$curr_sp.synt.final.df > temp ; mv temp Syn_tables_dir/$curr_sp.synt.final.df


rm Syn_tables_dir/$curr_sp.synt.final.df.strand
for syn_line in `cat Syn_tables_dir/$curr_sp.synt.final.df` ; do 

  curr_name=`echo "$syn_line" | cut -f4 -d "," | sed 's/rna.//g' | sed 's/..$//g' | sed 's/_$//g'`
  
  if grep -q "$curr_name" GFF3_files_per_species/$curr_sp.gff ; then
    curr_strand=`grep "$curr_name" GFF3_files_per_species/$curr_sp.gff | cut -f7 | sort | uniq | head -1`
  else 
    curr_strand="+"
  fi

  echo "$syn_line,$curr_strand,$curr_sp" >> Syn_tables_dir/$curr_sp.synt.final.df.strand

done
mv Syn_tables_dir/$curr_sp.synt.final.df.strand Syn_tables_dir/$curr_sp.synt.final.df



#Now add Borostomias_antarcticus

curr_OGG=N5.HOG0019111
curr_sp=Borostomias_antarcticus
ref_sp=Hypomesus_transpacificus

grep -A10 -B15 "N5_HOG0015817" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0024350" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df


awk '{ if (NF < 5 || $5 == "") $5 = "no_OGG"; print $0;}' Syn_tables_dir/$curr_sp.synt.df > Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/	/,/g' Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/ no_OGG/,no_OGG/g' Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/ /,/g' Syn_tables_dir/$curr_sp.synt.final.df
grep -v ",,,," Syn_tables_dir/$curr_sp.synt.final.df > temp ; mv temp Syn_tables_dir/$curr_sp.synt.final.df


rm Syn_tables_dir/$curr_sp.synt.final.df.strand
for syn_line in `cat Syn_tables_dir/$curr_sp.synt.final.df` ; do 

  curr_name=`echo "$syn_line" | cut -f4 -d "," | sed 's/rna.//g' | sed 's/..$//g' | sed 's/_$//g'`
  
  if grep -q "$curr_name" GFF3_files_per_species/$curr_sp.gff ; then
    curr_strand=`grep "$curr_name" GFF3_files_per_species/$curr_sp.gff | cut -f7 | sort | uniq | head -1`
  else 
    curr_strand="+"
  fi

  echo "$syn_line,$curr_strand,$curr_sp" >> Syn_tables_dir/$curr_sp.synt.final.df.strand

done
mv Syn_tables_dir/$curr_sp.synt.final.df.strand Syn_tables_dir/$curr_sp.synt.final.df



#Finally add Gadus_morhua

curr_OGG=N5.HOG0019111
curr_sp=Gadus_morhua
ref_sp=Hypomesus_transpacificus

grep -A7 -B7 "N5_HOG0015817" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df
grep -A7 -B7 "N5_HOG0023474" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df
grep -A7 -B7 "N5_HOG0031439" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df
grep -A9 -B9 "N5_HOG0024350" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df


awk '{ if (NF < 5 || $5 == "") $5 = "no_OGG"; print $0;}' Syn_tables_dir/$curr_sp.synt.df > Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/	/,/g' Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/ no_OGG/,no_OGG/g' Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/ /,/g' Syn_tables_dir/$curr_sp.synt.final.df
grep -v ",,,," Syn_tables_dir/$curr_sp.synt.final.df > temp ; mv temp Syn_tables_dir/$curr_sp.synt.final.df


rm Syn_tables_dir/$curr_sp.synt.final.df.strand
for syn_line in `cat Syn_tables_dir/$curr_sp.synt.final.df` ; do 

  curr_name=`echo "$syn_line" | cut -f4 -d "," | sed 's/rna.//g' | sed 's/..$//g' | sed 's/_$//g'`
  
  if grep -q "$curr_name" GFF3_files_per_species/$curr_sp.gff ; then
    curr_strand=`grep "$curr_name" GFF3_files_per_species/$curr_sp.gff | cut -f7 | sort | uniq | head -1`
  else 
    curr_strand="+"
  fi

  echo "$syn_line,$curr_strand,$curr_sp" >> Syn_tables_dir/$curr_sp.synt.final.df.strand

done
mv Syn_tables_dir/$curr_sp.synt.final.df.strand Syn_tables_dir/$curr_sp.synt.final.df



#Now add Scomber_scombrus

curr_OGG=N5.HOG0019111
curr_sp=Scomber_scombrus
ref_sp=Thunnus_maccoyii

grep -A10 -B10 "HOG0019111" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "NC_084971" > Syn_tables_dir/$curr_sp.synt.df


awk '{ if (NF < 5 || $5 == "") $5 = "no_OGG"; print $0;}' Syn_tables_dir/$curr_sp.synt.df > Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/	/,/g' Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/ no_OGG/,no_OGG/g' Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/ /,/g' Syn_tables_dir/$curr_sp.synt.final.df
grep -v ",,,," Syn_tables_dir/$curr_sp.synt.final.df > temp ; mv temp Syn_tables_dir/$curr_sp.synt.final.df


rm Syn_tables_dir/$curr_sp.synt.final.df.strand
for syn_line in `cat Syn_tables_dir/$curr_sp.synt.final.df` ; do 

  curr_name=`echo "$syn_line" | cut -f4 -d "," | sed 's/rna.//g' | sed 's/..$//g' | sed 's/_$//g'`
  
  if grep -q "$curr_name" GFF3_files_per_species/$curr_sp.gff ; then
    curr_strand=`grep "$curr_name" GFF3_files_per_species/$curr_sp.gff | cut -f7 | sort | uniq | head -1`
  else 
    curr_strand="+"
  fi

  echo "$syn_line,$curr_strand,$curr_sp" >> Syn_tables_dir/$curr_sp.synt.final.df.strand

done
mv Syn_tables_dir/$curr_sp.synt.final.df.strand Syn_tables_dir/$curr_sp.synt.final.df





#Now add Thalassophryne_amazonica

curr_OGG=N5.HOG0019111
curr_sp=Thalassophryne_amazonica
ref_sp=Thunnus_maccoyii

grep -A15 -B15 "N5_HOG0019111" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "NC_047117" > Syn_tables_dir/$curr_sp.synt.df

awk '{ if (NF < 5 || $5 == "") $5 = "no_OGG"; print $0;}' Syn_tables_dir/$curr_sp.synt.df > Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/	/,/g' Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/ no_OGG/,no_OGG/g' Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/ /,/g' Syn_tables_dir/$curr_sp.synt.final.df
grep -v ",,,," Syn_tables_dir/$curr_sp.synt.final.df > temp ; mv temp Syn_tables_dir/$curr_sp.synt.final.df


rm Syn_tables_dir/$curr_sp.synt.final.df.strand
for syn_line in `cat Syn_tables_dir/$curr_sp.synt.final.df` ; do 

  curr_name=`echo "$syn_line" | cut -f4 -d "," | sed 's/rna.//g' | sed 's/..$//g' | sed 's/_$//g'`
  
  if grep -q "$curr_name" GFF3_files_per_species/$curr_sp.gff ; then
    curr_strand=`grep "$curr_name" GFF3_files_per_species/$curr_sp.gff | cut -f7 | sort | uniq | head -1`
  else 
    curr_strand="+"
  fi

  echo "$syn_line,$curr_strand,$curr_sp" >> Syn_tables_dir/$curr_sp.synt.final.df.strand

done
mv Syn_tables_dir/$curr_sp.synt.final.df.strand Syn_tables_dir/$curr_sp.synt.final.df






#Finally add  Centroberyx_gerrardi

curr_OGG=N5.HOG0019111
curr_sp=Centroberyx_gerrardi
ref_sp=Thunnus_maccoyii

grep -A8 -B8 "N5_HOG0019111" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "JAPMTB010005717" > Syn_tables_dir/$curr_sp.synt.df
grep -A8 -B9 "N5_HOG0016604" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df
grep -A8 -B9 "N5_HOG0040864" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df



awk '{ if (NF < 5 || $5 == "") $5 = "no_OGG"; print $0;}' Syn_tables_dir/$curr_sp.synt.df > Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/	/,/g' Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/ no_OGG/,no_OGG/g' Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/ /,/g' Syn_tables_dir/$curr_sp.synt.final.df
grep -v ",,,," Syn_tables_dir/$curr_sp.synt.final.df > temp ; mv temp Syn_tables_dir/$curr_sp.synt.final.df


rm Syn_tables_dir/$curr_sp.synt.final.df.strand
for syn_line in `cat Syn_tables_dir/$curr_sp.synt.final.df` ; do 

  curr_name=`echo "$syn_line" | cut -f4 -d "," | sed 's/rna.//g' | sed 's/..$//g' | sed 's/_$//g'`
  
  if grep -q "$curr_name" GFF3_files_per_species/$curr_sp.gff ; then
    curr_strand=`grep "$curr_name" GFF3_files_per_species/$curr_sp.gff | cut -f7 | sort | uniq | head -1`
  else 
    curr_strand="+"
  fi

  echo "$syn_line,$curr_strand,$curr_sp" >> Syn_tables_dir/$curr_sp.synt.final.df.strand

done
mv Syn_tables_dir/$curr_sp.synt.final.df.strand Syn_tables_dir/$curr_sp.synt.final.df




#Lets create a seq info ogg file for both clades


cat non_transfer_species_1 species_to_draw.clade1.ordered non_transfer_species_2 > temp ; mv temp species_to_draw.clade1.ordered


rm seq_clustered_infos_ogg.clade1.csv
for curr_sp in `cat species_to_draw.clade1.ordered` ; do 
  cat Syn_tables_dir/$curr_sp.synt.final.df >> seq_clustered_infos_ogg.clade1.csv
done



sed -i 's/ /,/g' seq_clustered_infos_ogg.clade1.csv


### Now lets crate a seqID file that list all the clusters and their length per species, and arrange it for the graph representation
### We also create a good geneID file 

cut -f1,7 -d "," seq_clustered_infos_ogg.clade1.csv  | sort | uniq > uniq_clusters.clade1


for species in `cat species_to_draw.clade1.ordered` ; do 
  grep "$species" uniq_clusters.clade1 >> temp.txt
done ; mv temp.txt uniq_clusters.clade1




rm cluster_id_arranged.clade1.txt
rm seq_clustered_infos_ogg_num.clade1.csv
for cluster in `cat uniq_clusters.clade1` ; do 

  scaffold=`echo "$cluster" | cut -f1 -d ","`
  species=`echo "$cluster" | cut -f2 -d ","`

  start=`grep "$scaffold.*$species" seq_clustered_infos_ogg.clade1.csv | cut -f2 -d "," | sort -n | head -1`
  end=`grep "$scaffold.*$species" seq_clustered_infos_ogg.clade1.csv | cut -f3 -d "," | sort -n | tail -1`
  length=$((end - start))

  echo "$species,$scaffold,$length" >> cluster_id_arranged.clade1.txt

  start_min_1=$((start - 1))
  grep "$scaffold.*$species" seq_clustered_infos_ogg.clade1.csv | awk -F, -v start_min_1="$start_min_1" '{$2 = $2 - start_min_1; print}' | awk -F" " -v start_min_1="$start_min_1" '{$3 = $3 - start_min_1; print}' | sed 's/ /,/g'  >> seq_clustered_infos_ogg_num.clade1.csv

done 
 



##Now create a link table for clade 1




i=1
j=`wc -l  < species_to_draw.clade1.ordered`
rm Species_Comparison_table.clade1.txt
while [ $i -lt $j ] ; do 

  n=$(( i + 1 )) 


  sp1=`head -$i species_to_draw.clade1.ordered | tail -1`
  sp2=`head -$n species_to_draw.clade1.ordered | tail -1`

  echo "$sp1,$sp2" >> Species_Comparison_table.clade1.txt

  i=$(( i + 1 )) 


done





rm link_table.clade1.txt
for comparisons in `cat Species_Comparison_table.clade1.txt` ; do  #Species_Comparison_table.txt = table with two columns with each pair of species being compared

  species1=`echo "$comparisons" | cut -f1 -d ","`
  species2=`echo "$comparisons" | cut -f2 -d ","`

  grep "$species1" seq_clustered_infos_ogg_num.clade1.csv  | cut -f5 -d "," | sort | uniq > uniq_OGG_sp1
  grep "$species2" seq_clustered_infos_ogg_num.clade1.csv  | cut -f5 -d "," | sort | uniq > uniq_OGG_sp2
  cat uniq_OGG_sp1 uniq_OGG_sp2 | sort | uniq -c | sed 's/^ *//g' | grep -v "^1 " | grep -v "no_OGG" | cut -f2 -d " " | sort | uniq > OGG_list.txt 

  for OGG in `cat OGG_list.txt` ; do

    grep "$OGG" seq_clustered_infos_ogg_num.clade1.csv > current_OGG.txt 
    grep "$species1" current_OGG.txt | cut -f7,1,2,3 -d "," | awk -v OFS="," -F"," '{print($4,$1,$2,$3)}' > species1.ogg.txt
    grep "$species2" current_OGG.txt | cut -f7,1,2,3 -d "," | awk -v OFS="," -F"," '{print($4,$1,$2,$3)}' > species2.ogg.txt

    for line_species1 in `cat species1.ogg.txt` ; do

      for line_species2 in `cat species2.ogg.txt` ; do 

        echo "$line_species1,$line_species2,$OGG" >> link_table.clade1.txt


      done
    done

    rm species1.ogg.txt ; rm species2.ogg.txt

  done
done



#Reduce size: Gadus_morhua,NC_044068_1,7245838



grep "Gadus_morhua" seq_clustered_infos_ogg_num.clade1.csv | grep "NC_044068_1" | cut -f2 -d ","  | sort -n
grep "Gadus_morhua" seq_clustered_infos_ogg_num.clade1.csv | grep "NC_044068_1" | cut -f3 -d ","  | sort -n

grep "Gadus_morhua" seq_clustered_infos_ogg_num.clade1.csv | grep "NC_044068_1" | cut -f2 -d ","  | awk '$1 > 2000000' > coord_to_reduce
grep "Gadus_morhua" seq_clustered_infos_ogg_num.clade1.csv | grep "NC_044068_1" | cut -f3 -d ","  | awk '$1 > 2000000' >> coord_to_reduce


for curr_coord in `cat coord_to_reduce` ; do
  new_coord=$(( curr_coord -  5487822))

  sed -i "s/$curr_coord/$new_coord/g" seq_clustered_infos_ogg_num.clade1.csv
  sed -i "s/$curr_coord/$new_coord/g" link_table.clade1.txt

done

new_max=`grep "Gadus_morhua" seq_clustered_infos_ogg_num.clade1.csv | cut -f2,3 -d "," | tr "," "\n" | sort -n | tail -1`
new_max=$(( new_max + 2000 ))
sed -i "s/7245838/$new_max/g" cluster_id_arranged.clade1.txt 



mv cluster_id_arranged.clade1.txt cluster_id_arranged.clade1.OsmSco.txt
mv seq_clustered_infos_ogg_num.clade1.csv seq_clustered_infos_ogg_num.clade1.OsmSco.csv
mv link_table.clade1.txt link_table.clade1.OsmSco.txt




##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##========================================== Verify that the gene is absent in close species ==========================================
##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================




samtools faidx Proteomes_BUSCO80/Hypomesus_transpacificus.fa Hypomesus_transpacificus---rna-XM_047049897.1 > Hypomesus_transpacificus---rna-XM_047049897.1.prot

exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs" Hypomesus_transpacificus---rna-XM_047049897.1.prot GCF_963692335.1_fOsmEpe2.1_genomic.fna > Osmerus_eperlanus.whole.exo
exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs" Hypomesus_transpacificus---rna-XM_047049897.1.prot GCA_949987555.1_fBorAnt1.1_genomic.fna > Borostomias_antarcticus.whole.exo
exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs" Hypomesus_transpacificus---rna-XM_047049897.1.prot GCF_902167405.1_gadMor3.0_genomic.fna > Gadus_morhua.whole.exo


samtools faidx Proteomes_BUSCO80/Alosa_alosa.fa Alosa_alosa---rna-XM_048231371.1 > Alosa_alosa---rna-XM_048231371.1.prot

exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs" Alosa_alosa---rna-XM_048231371.1.prot GCF_034702125.1_IST_EnEncr_1.0_genomic.fna > Engraulis_encrasicolus.whole.exo

#Check in the other hypomesus

exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs" Hypomesus_transpacificus---rna-XM_047049897.1.prot GCA_018346875.1_ASM1834687v1_genomic.fna  > Hypomesus_nipponensis.whole.exo



##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##========================================== Now let's find TE in the region  ===================================================
##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================


#Extract the region 

#Clupea_harengus
grep -A3 -B3 "N5_HOG0019111" GFF3_N5_OGGs/Clupea_harengus.gff.simplified.sorted.OGG.tiret 
samtools faidx GCF_900700415.2_Ch_v2.0.2_genomic.fna NC_045163.1:6841927-6871444 > N5.HOG0019111.Clupea_harengus.extended.1.fa
sed -i 's/:/-/g' N5.HOG0019111.Clupea_harengus.extended.1.fa
makeblastdb -in N5.HOG0019111.Clupea_harengus.extended.1.fa -dbtype nucl


#Thunnus_maccoyii
grep -A3 -B3 "N5_HOG0019111" GFF3_N5_OGGs/Thunnus_maccoyii.gff.simplified.sorted.OGG.tiret 
samtools faidx GCF_910596095.1_fThuMac1.1_genomic.fna NC_056534.1:35550502-35592373 > N5.HOG0019111.Thunnus_maccoyii.extended.1.fa
sed -i 's/:/-/g' N5.HOG0019111.Thunnus_maccoyii.extended.1.fa
makeblastdb -in N5.HOG0019111.Thunnus_maccoyii.extended.1.fa -dbtype nucl



#Hypomesus_transpacificus
grep -A3 -B3 "N5_HOG0019111"  GFF3_N5_OGGs/Hypomesus_transpacificus.gff.simplified.sorted.OGG.tiret
samtools faidx GCF_021917145.1_fHypTra1_genomic.fna NC_061085.1:4537408-4557518 > N5.HOG0019111.Hypomesus_transpacificus.extended.1.fa
sed -i 's/:/-/g' N5.HOG0019111.Hypomesus_transpacificus.extended.1.fa
makeblastdb -in N5.HOG0019111.Hypomesus_transpacificus.extended.1.fa -dbtype nucl




#Now lets find TEs using the repbase and dfam databases

tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0019111.Clupea_harengus.extended.1.fa -evalue 1e-5 -outfmt 6 -out TE.Clupea_harengus.1.tblastn -num_threads 8
sed -i 's/#//g' TE.Clupea_harengus.1.tblastn

tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0019111.Thunnus_maccoyii.extended.1.fa -evalue 1e-5 -outfmt 6 -out TE.Thunnus_maccoyii.1.tblastn -num_threads 8
sed -i 's/#//g' TE.Thunnus_maccoyii.1.tblastn

tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0019111.Hypomesus_transpacificus.extended.1.fa -evalue 1e-5 -outfmt 6 -out TE.Hypomesus_transpacificus.1.tblastn -num_threads 8
sed -i 's/#//g' TE.Hypomesus_transpacificus.1.tblastn

#merge tblastn hits and find the best TE match by doing a blastx

Rscript Rscript_merge_blast_hits.R TE.Clupea_harengus.1.tblastn TE.Clupea_harengus.1.tblastn.merged
Rscript Rscript_merge_blast_hits.R TE.Thunnus_maccoyii.1.tblastn TE.Thunnus_maccoyii.1.tblastn.merged
Rscript Rscript_merge_blast_hits.R TE.Hypomesus_transpacificus.1.tblastn TE.Hypomesus_transpacificus.1.tblastn.merged


xargs samtools faidx N5.HOG0019111.Clupea_harengus.extended.1.fa < TE.Clupea_harengus.1.tblastn.merged > TE.Clupea_harengus.1.BEST.fa
blastx -query TE.Clupea_harengus.1.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Clupea_harengus.1.BEST.blastx -max_target_seqs 1
cut -f1 TE.Clupea_harengus.1.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Clupea_harengus.1.BEST.blastx >> temp  ; done ; mv temp TE.Clupea_harengus.1.BEST.blastx

xargs samtools faidx N5.HOG0019111.Thunnus_maccoyii.extended.1.fa < TE.Thunnus_maccoyii.1.tblastn.merged > TE.Thunnus_maccoyii.1.BEST.fa
blastx -query TE.Thunnus_maccoyii.1.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Thunnus_maccoyii.1.BEST.blastx -max_target_seqs 1
cut -f1 TE.Thunnus_maccoyii.1.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Thunnus_maccoyii.1.BEST.blastx >> temp  ; done ; mv temp TE.Thunnus_maccoyii.1.BEST.blastx

xargs samtools faidx N5.HOG0019111.Hypomesus_transpacificus.extended.1.fa < TE.Hypomesus_transpacificus.1.tblastn.merged > TE.Hypomesus_transpacificus.1.BEST.fa
blastx -query TE.Hypomesus_transpacificus.1.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Hypomesus_transpacificus.1.BEST.blastx -max_target_seqs 1
cut -f1  TE.Hypomesus_transpacificus.1.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Hypomesus_transpacificus.1.BEST.blastx >> temp  ; done ; mv temp TE.Hypomesus_transpacificus.1.BEST.blastx


#Now find shared elements

cut -f2 TE.Clupea_harengus.1.BEST.blastx | sort | uniq > TE.Clupea_harengus.1.uniqTE
cut -f2 TE.Thunnus_maccoyii.1.BEST.blastx | sort | uniq > TE.Thunnus_maccoyii.1.uniqTE
cut -f2 TE.Hypomesus_transpacificus.1.BEST.blastx | sort | uniq > TE.Hypomesus_transpacificus.1.uniqTE


comm -12 TE.Clupea_harengus.1.uniqTE TE.Hypomesus_transpacificus.1.uniqTE
comm -12 TE.Clupea_harengus.1.uniqTE TE.Thunnus_maccoyii.1.uniqTE
comm -12 TE.Thunnus_maccoyii.1.uniqTE TE.Hypomesus_transpacificus.1.uniqTE

## No shared TEs


### In a FINAL step we will make another synteny plot more zoomed




## Clupea_harengus

rm seq_clustered_infos_ogg.TE.txt ; rm clusters_ID_TE.txt


scaffold=`grep ">" N5.HOG0019111.Clupea_harengus.extended.1.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0019111.Clupea_harengus.extended.1.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0019111.Clupea_harengus.extended.1.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Clupea_harengus,$scaffold,$length" > clusters_ID_TE.txt

grep "N5_HOG0019111"  GFF3_N5_OGGs/Clupea_harengus.gff.simplified.sorted.OGG.tiret

grep "XM_031577666" GFF3_files_per_species/Clupea_harengus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons
grep "XM_031577667" GFF3_files_per_species/Clupea_harengus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_2.exons
grep "XM_031578285" GFF3_files_per_species/Clupea_harengus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_3.exons


nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene1.$nbexon,N5_HOG0019111,+,Clupea_harengus" >> seq_clustered_infos_ogg.TE.txt
done
nbexon=0
for line in `cat gene_2.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene2.$nbexon,N5_HOG0019111,+,Clupea_harengus" >> seq_clustered_infos_ogg.TE.txt
done
nbexon=0
for line in `cat gene_3.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene3.$nbexon,N5_HOG0019111,+,Clupea_harengus" >> seq_clustered_infos_ogg.TE.txt
done

cut -f1 TE.Clupea_harengus.1.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Clupea_harengus.1.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0019111.Clupea_harengus.extended.1.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0019111.Clupea_harengus.extended.1.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold,/g" | sed "s/$/,$strand,Clupea_harengus/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done



## Hypomesus_transpacificus


scaffold=`grep ">" N5.HOG0019111.Hypomesus_transpacificus.extended.1.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0019111.Hypomesus_transpacificus.extended.1.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0019111.Hypomesus_transpacificus.extended.1.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Hypomesus_transpacificus,$scaffold,$length" >> clusters_ID_TE.txt

grep "N5_HOG0019111"  GFF3_N5_OGGs/Hypomesus_transpacificus.gff.simplified.sorted.OGG.tiret

grep "XM_047049897" GFF3_files_per_species/Hypomesus_transpacificus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons


nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene1.$nbexon,N5_HOG0019111,+,Hypomesus_transpacificus" >> seq_clustered_infos_ogg.TE.txt
done


cut -f1 TE.Hypomesus_transpacificus.1.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Hypomesus_transpacificus.1.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0019111.Hypomesus_transpacificus.extended.1.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0019111.Hypomesus_transpacificus.extended.1.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold,/g" | sed "s/$/,$strand,Hypomesus_transpacificus/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done



## Thunnus_maccoyii


scaffold=`grep ">" N5.HOG0019111.Thunnus_maccoyii.extended.1.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0019111.Thunnus_maccoyii.extended.1.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0019111.Thunnus_maccoyii.extended.1.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Thunnus_maccoyii,$scaffold,$length" >> clusters_ID_TE.txt

grep "N5_HOG0019111"  GFF3_N5_OGGs/Thunnus_maccoyii.gff.simplified.sorted.OGG.tiret

grep "XM_042435094" GFF3_files_per_species/Thunnus_maccoyii.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons
grep "XM_042435120" GFF3_files_per_species/Thunnus_maccoyii.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_2.exons
grep "XM_042435142" GFF3_files_per_species/Thunnus_maccoyii.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_3.exons


nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene1.$nbexon,N5_HOG0019111,-,Thunnus_maccoyii" >> seq_clustered_infos_ogg.TE.txt
done
nbexon=0
for line in `cat gene_2.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene2.$nbexon,N5_HOG0019111,-,Thunnus_maccoyii" >> seq_clustered_infos_ogg.TE.txt
done
nbexon=0
for line in `cat gene_3.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene3.$nbexon,N5_HOG0019111,-,Thunnus_maccoyii" >> seq_clustered_infos_ogg.TE.txt
done

cut -f1 TE.Thunnus_maccoyii.1.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Thunnus_maccoyii.1.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0019111.Thunnus_maccoyii.extended.1.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0019111.Thunnus_maccoyii.extended.1.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold,/g" | sed "s/$/,$strand,Thunnus_maccoyii/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done





### Check the observed non-coding conserved region

samtools faidx N5.HOG0019111.Hypomesus_transpacificus.extended.1.fa NC_061085.1-4537408-4557518:2700-4200 > region1.Hypomesus_transpacificus.fa
samtools faidx N5.HOG0019111.Clupea_harengus.extended.1.fa NC_045163.1-6841927-6871444:3300-4800 > region1.Clupea_harengus.fa

sed -i 's/>.*/>H.transpacificus/g' region1.Hypomesus_transpacificus.fa
sed -i 's/>.*/>C.harengus/g' region1.Clupea_harengus.fa

needle -asequence region1.Hypomesus_transpacificus.fa -bsequence region1.Clupea_harengus.fa -outfile region1.aln -gapopen 10.0 -gapextend 0.5


cp ../launch_blast_allgenomes.sh ./

sbatch --qos=6hours -c 10 --mem=20G launch_blast_allgenomes.sh region1.Hypomesus_transpacificus.fa 1e-20 region1.blastn.tsv



##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##========================================== Selective pressure analysis ============================================================
##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================



#Compute the dN/dS on every branches
sbatch --qos=1day -c 4 --mem=10G -e error.fitMG4.out -o slurm.fitMG4.out launch_fitMG4.sh N5.HOG0019111


Recipient branches : 

Hypomesus_transpacificus_rna_XM_047049897_1 #> Osmeriformes recipient
Clupea_harengus_rna_XM_031577667_2 #> Clupeiformes recipient
Node1109 #> Clupeiformes recipient
Clupea_harengus_rna_XM_031577666_1 #> Clupeiformes recipient
Node1108 #> Clupeiformes recipient
Clupea_harengus_rna_XM_031578285_2 #> Clupeiformes recipient
Node1107 #> Clupeiformes recipient
Alosa_alosa_rna_XM_048231371_1 #> Clupeiformes recipient
Alosa_sapidissima_rna_XM_042075309_1 #> Clupeiformes recipient
Node1118 #> Clupeiformes recipient
Sardina_pilchardus_rna_XM_062535392_1 #> Clupeiformes recipient
Node1117 #> Clupeiformes recipient

#Test positive selection and relaxed selection on receiver branch

sbatch --qos=1week -c 4 --mem=10G -e error.absrel.cand.out -o slurm.absrel.cand.out --job-name=HOG0019111 launch_absrel_cand.sh N5.HOG0019111


#Launch RELAX

### Adaptive branch site random effects likelihood test 
#Likelihood ratio test for episodic diversifying positive selection at Holm-Bonferroni corrected _p =   0.0500_ found **1** branches under selection among **12** tested.
#* Clupea_harengus_rna_XM_031577667_2, p-value =  0.00000


===> $HOG.prot.aln.treefile.SelecMarked.relax => select branches detected with accelerated evolution
sbatch --qos=1week -c 6 --mem=10G -e error.relax.out -o slurm.relax.out --job-name=HOG0019111 launch_RELAX.sh N5.HOG0019111







