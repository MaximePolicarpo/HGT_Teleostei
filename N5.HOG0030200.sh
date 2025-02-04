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

samtools faidx Proteomes_BUSCO80/Pangasianodon_hypophthalmus.fa Pangasianodon_hypophthalmus---rna-XM_026945431.3 > Pangasianodon_hypophthalmus---rna-XM_026945431.3.prot

#blastp against the fish proteome (50 best hits) + Uniprot database without any fish species (10 best hits). Here I will only take 20 catffishes otherwise there are way too much of them

blastp -query Pangasianodon_hypophthalmus---rna-XM_026945431.3.prot -db concatenated_proteomes.fa -outfmt 6 -out Gene_vs_FishProteome.blastp -num_threads 8 -max_target_seqs 52
blastp -query Pangasianodon_hypophthalmus---rna-XM_026945431.3.prot -db non_actino_uniprot.fa -outfmt 6 -out Gene_vs_Uniprot.blastp -num_threads 8  -max_target_seqs 10 -evalue 1e-5

cut -f2 Gene_vs_FishProteome.blastp | sort | uniq  > closest_fish_seq.id
cut -f2 Gene_vs_Uniprot.blastp | sort | uniq > closest_nonfish_seq.id 

xargs samtools faidx concatenated_proteomes.fa < closest_fish_seq.id > closest_fish_seq.fa
xargs samtools faidx non_actino_uniprot.fa  < closest_nonfish_seq.id > closest_nonfish_seq.fa

#Align with muscle and trim with trimal

sed -i 's/sp_Q00973_B4GN1_HUMAN_Beta-1,4_N-acetylgalactosaminyltransferase_1_OS_Homo_sapiens_OX_9606_GN_B4GALNT1_PE_1_SV_2/sp_Q00973_B4GN1_HUMAN_Beta-1/g' closest_nonfish_seq.fa
sed -i 's/sp_Q09199_B4GN2_MOUSE_Beta-1,4_N-acetylgalactosaminyltransferase_2_OS_Mus_musculus_OX_10090_GN_B4galnt2_PE_2_SV_1/sp_Q09199_B4GN2_MOUSE_Beta-1/g' closest_nonfish_seq.fa
sed -i 's/sp_Q09200_B4GN1_MOUSE_Beta-1,4_N-acetylgalactosaminyltransferase_1_OS_Mus_musculus_OX_10090_GN_B4galnt1_PE_1_SV_1/sp_Q09200_B4GN1_MOUSE_Beta-1/g' closest_nonfish_seq.fa
sed -i 's/sp_Q10468_B4GN1_RAT_Beta-1,4_N-acetylgalactosaminyltransferase_1_OS_Rattus_norvegicus_OX_10116_GN_B4galnt1_PE_1_SV_1/sp_Q10468_B4GN1_RAT_Beta-1/g' closest_nonfish_seq.fa
sed -i 's/sp_Q8NHY0_B4GN2_HUMAN_Beta-1,4_N-acetylgalactosaminyltransferase_2_OS_Homo_sapiens_OX_9606_GN_B4GALNT2_PE_1_SV_2/sp_Q8NHY0_B4GN2_HUMAN_Beta-1/g' closest_nonfish_seq.fa


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
##========================================== Now make a synteny plot  ================================================================
##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================


#N5.HOG0030200.clade_hgt.txt=> clade where the HGT is extracted from iTOL (from N5.HOG0030200.prot.aln.treefile)

cut -f1,2 -d "_" N5.HOG0030200.clade_hgt.txt | sort | uniq  > species_to_draw.clade1
#Choose close species without the gene .. 


echo "Silurus_meridionalis" > non_transfer_species_1
echo "Ictalurus_punctatus" >> non_transfer_species_1
echo "Pangasius_djambal" >> non_transfer_species_1

echo "Cyprinus_carpio" > non_transfer_species_2
echo "Onychostoma_macrolepis" >> non_transfer_species_2
echo "Danionella_translucida" >> non_transfer_species_2




rm -r Syn_tables_dir ; mkdir Syn_tables_dir


#First order the species names 

rm species_to_draw.clade1.ordered

echo "Pangasianodon_hypophthalmus" >> species_to_draw.clade1.ordered
echo "Carassius_auratus" >> species_to_draw.clade1.ordered


#If there are more than 4 species from the same order, keep only 4 


#First extract genes around 

curr_OGG=N5.HOG0030200

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



#First add Pangasius_djambal

curr_OGG=N5.HOG0030200
curr_sp=Pangasius_djambal
ref_sp=Pangasianodon_hypophthalmus

grep -A10 -B10 "N5_HOG0030200" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "CM040976\|CM040975" > Syn_tables_dir/$curr_sp.synt.df


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




#Now add Ictalurus_punctatus

curr_OGG=N5.HOG0030200
curr_sp=Ictalurus_punctatus
ref_sp=Pangasianodon_hypophthalmus

grep -A10 -B10 "N5_HOG0029507" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df
grep -A5 -B5 "N5_HOG0001046" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "NW_026521123" >> Syn_tables_dir/$curr_sp.synt.df
grep -A5 -B5 "N5_HOG0020794" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df

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



#Now add Silurus_meridionalis

curr_OGG=N5.HOG0030200
curr_sp=Silurus_meridionalis
ref_sp=Pangasianodon_hypophthalmus

grep -A10 -B10 "N5_HOG0029507" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df
grep -A5 -B5 "N5_HOG0001046" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "NC_060900" >> Syn_tables_dir/$curr_sp.synt.df
grep -A5 -B5 "N5_HOG0020794" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df

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




#Now add Cyprinus_carpio

curr_OGG=N5.HOG0030200
curr_sp=Cyprinus_carpio
ref_sp=Carassius_auratus

grep -A10 -B10 "HOG0030200" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "NC_056574" > Syn_tables_dir/$curr_sp.synt.df

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




#Now add Onychostoma_macrolepis

curr_OGG=N5.HOG0030200
curr_sp=Onychostoma_macrolepis
ref_sp=Carassius_auratus

grep -A10 -B10 "HOG0030200" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "NC_081157_1	4" > Syn_tables_dir/$curr_sp.synt.df

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




#Now add Danionella_translucida

curr_OGG=N5.HOG0030200
curr_sp=Danionella_translucida
ref_sp=Carassius_auratus

grep -A10 -B10 "HOG0030200" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "SRMA01026995" > Syn_tables_dir/$curr_sp.synt.df
grep -A5 -B5 "N5_HOG0031656" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df

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

exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs" Pangasianodon_hypophthalmus---rna-XM_026945431.3.prot GCF_014805685.1_ASM1480568v1_genomic.fna > Silurus_meridionalis.whole.exo
exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs" Pangasianodon_hypophthalmus---rna-XM_026945431.3.prot GCF_001660625.3_Coco_2.0_genomic.fna > Ictalurus_punctatus.whole.exo



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

#Carassius_auratus
grep -A3 -B3 "N5_HOG0030200" GFF3_N5_OGGs/Carassius_auratus.gff.simplified.sorted.OGG.tiret | grep "NC_039245"

samtools faidx GCF_003368295.1_ASM336829v1_genomic.fna NC_039245.1:4797615-4817324 > N5.HOG0030200.Carassius_auratus.extended.1.fa

sed -i 's/:/-/g' N5.HOG0030200.Carassius_auratus.extended.1.fa
makeblastdb -in N5.HOG0030200.Carassius_auratus.extended.1.fa -dbtype nucl



#Pangasianodon_hypophthalmus
grep -A3 -B3 "N5_HOG0030200"  GFF3_N5_OGGs/Pangasianodon_hypophthalmus.gff.simplified.sorted.OGG.tiret
samtools faidx GCF_027358585.1_fPanHyp1.pri_genomic.fna NC_069711.1:3821864-3837403 > N5.HOG0030200.Pangasianodon_hypophthalmus.extended.1.fa

sed -i 's/:/-/g' N5.HOG0030200.Pangasianodon_hypophthalmus.extended.1.fa
makeblastdb -in N5.HOG0030200.Pangasianodon_hypophthalmus.extended.1.fa -dbtype nucl




#Now lets find TEs using the repbase and dfam databases

tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0030200.Carassius_auratus.extended.1.fa -evalue 1e-5 -outfmt 6 -out TE.Carassius_auratus.1.tblastn -num_threads 8
sed -i 's/#//g' TE.Carassius_auratus.1.tblastn

tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0030200.Pangasianodon_hypophthalmus.extended.1.fa -evalue 1e-5 -outfmt 6 -out TE.Pangasianodon_hypophthalmus.1.tblastn -num_threads 8
sed -i 's/#//g' TE.Pangasianodon_hypophthalmus.1.tblastn

#merge tblastn hits and find the best TE match by doing a blastx

Rscript Rscript_merge_blast_hits.R TE.Carassius_auratus.1.tblastn TE.Carassius_auratus.1.tblastn.merged
Rscript Rscript_merge_blast_hits.R TE.Pangasianodon_hypophthalmus.1.tblastn TE.Pangasianodon_hypophthalmus.1.tblastn.merged


xargs samtools faidx N5.HOG0030200.Carassius_auratus.extended.1.fa < TE.Carassius_auratus.1.tblastn.merged > TE.Carassius_auratus.1.BEST.fa
blastx -query TE.Carassius_auratus.1.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Carassius_auratus.1.BEST.blastx -max_target_seqs 1
cut -f1 TE.Carassius_auratus.1.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Carassius_auratus.1.BEST.blastx >> temp  ; done ; mv temp TE.Carassius_auratus.1.BEST.blastx

xargs samtools faidx N5.HOG0030200.Pangasianodon_hypophthalmus.extended.1.fa < TE.Pangasianodon_hypophthalmus.1.tblastn.merged > TE.Pangasianodon_hypophthalmus.1.BEST.fa
blastx -query TE.Pangasianodon_hypophthalmus.1.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Pangasianodon_hypophthalmus.1.BEST.blastx -max_target_seqs 1
cut -f1  TE.Pangasianodon_hypophthalmus.1.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Pangasianodon_hypophthalmus.1.BEST.blastx >> temp  ; done ; mv temp TE.Pangasianodon_hypophthalmus.1.BEST.blastx


#Now find shared elements

cut -f2 TE.Carassius_auratus.1.BEST.blastx | sort | uniq > TE.Carassius_auratus.1.uniqTE
cut -f2 TE.Pangasianodon_hypophthalmus.1.BEST.blastx | sort | uniq > TE.Pangasianodon_hypophthalmus.1.uniqTE
comm -12 TE.Carassius_auratus.1.uniqTE TE.Pangasianodon_hypophthalmus.1.uniqTE

## No shared TEs



### In a FINAL step we will make another synteny plot more zoomed


## Carassius_auratus

rm seq_clustered_infos_ogg.TE.txt ; rm clusters_ID_TE.txt


scaffold=`grep ">" N5.HOG0030200.Carassius_auratus.extended.1.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0030200.Carassius_auratus.extended.1.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0030200.Carassius_auratus.extended.1.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Carassius_auratus,$scaffold,$length" > clusters_ID_TE.txt

grep "N5_HOG0030200"  GFF3_N5_OGGs/Carassius_auratus.gff.simplified.sorted.OGG.tiret

grep "XM_026200108" GFF3_files_per_species/Carassius_auratus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons


nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene1.$nbexon,N5_HOG0030200,+,Carassius_auratus" >> seq_clustered_infos_ogg.TE.txt
done

cut -f1 TE.Carassius_auratus.1.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Carassius_auratus.1.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0030200.Carassius_auratus.extended.1.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0030200.Carassius_auratus.extended.1.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold,/g" | sed "s/$/,$strand,Carassius_auratus/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done




## Pangasianodon_hypophthalmus

scaffold=`grep ">" N5.HOG0030200.Pangasianodon_hypophthalmus.extended.1.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0030200.Pangasianodon_hypophthalmus.extended.1.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0030200.Pangasianodon_hypophthalmus.extended.1.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Pangasianodon_hypophthalmus,$scaffold,$length" >> clusters_ID_TE.txt

grep "N5_HOG0030200"  GFF3_N5_OGGs/Pangasianodon_hypophthalmus.gff.simplified.sorted.OGG.tiret

grep "XM_026945431" GFF3_files_per_species/Pangasianodon_hypophthalmus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons


nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene1.$nbexon,N5_HOG0030200,-,Pangasianodon_hypophthalmus" >> seq_clustered_infos_ogg.TE.txt
done

cut -f1 TE.Pangasianodon_hypophthalmus.1.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Pangasianodon_hypophthalmus.1.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0030200.Pangasianodon_hypophthalmus.extended.1.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0030200.Pangasianodon_hypophthalmus.extended.1.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold,/g" | sed "s/$/,$strand,Pangasianodon_hypophthalmus/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done



### Check the conserved non-coding region observed

samtools faidx N5.HOG0030200.Pangasianodon_hypophthalmus.extended.1.fa.rev NC_069711.1-3821864-3837403:5100-5900 > region1.Pangasianodon_hypophthalmus.fa
samtools faidx N5.HOG0030200.Carassius_auratus.extended.1.fa NC_039245.1-4797615-4817324:4900-6100 > region1.Carassius_auratus.fa

sed -i 's/>.*/>P.hypophthalmus/g' region1.Pangasianodon_hypophthalmus.fa
sed -i 's/>.*/>C.auratus/g' region1.Carassius_auratus.fa

needle -asequence region1.Pangasianodon_hypophthalmus.fa -bsequence region1.Carassius_auratus.fa -outfile region1.aln -gapopen 10.0 -gapextend 0.5


sbatch --qos=6hours -c 10 --mem=20G launch_blast_allgenomes.sh region1.Pangasianodon_hypophthalmus.fa 1e-2 region1.blastn.tsv


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
sbatch --qos=1day -c 4 --mem=10G -e error.fitMG4.out -o slurm.fitMG4.out launch_fitMG4.sh N5.HOG0030200


Recipient branches : 
Node14
Node15
Pangasianodon hypophthalmus rna XM 026945431 3
Pangasius djambal PDJAM T00103720
Pangasius djambal PDJAM T00012360


#Test positive selection and relaxed selection on receiver branch

sbatch --qos=1week -c 4 --mem=10G -e error.absrel.cand.out -o slurm.absrel.cand.out --job-name=HOG0030200 launch_absrel_cand.sh N5.HOG0030200

### Adaptive branch site random effects likelihood test 
#Likelihood ratio test for episodic diversifying positive selection at Holm-Bonferroni corrected _p =   0.0500_ found **3** branches under selection among **5** tested.
#* Node15, p-value =  0.00000
#* Pangasius_djambal_PDJAM_T00103720, p-value =  0.00053
#* Node14, p-value =  0.01181


#Now , lets launch RELAX

===> $HOG.prot.aln.treefile.SelecMarked.relax => select branches detected with accelerated evolution
sbatch --qos=1week -c 6 --mem=10G -e error.relax.out -o slurm.relax.out --job-name=HOG0030200 launch_RELAX.sh N5.HOG0030200



##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##========================================== DNA phylogeny of the HGT clade ==========================================================
##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================




grep ">" N5.HOG0030200.prot | grep "Myxocyprinus_asiaticus_rna_XM_051690293_1\|Myxocyprinus_asiaticus_rna_XM_051690292_1\|Xyrauchen_texanus_rna_XM_052119494_1\|Carassius_gibelio_rna_XM_052538799_1\|Carassius_auratus_rna_XM_026200108_1\|Carassius_carassius_rna_XM_059557239_1\|Cyprinus_carpio_rna_XM_042728143_1\|Sinocyclocheilus_grahami_rna_XM_016268950_1\|Sinocyclocheilus_anshuiensis_rna_XM_016459463_1\|Pangasius_djambal_PDJAM_T00103720\|Pangasius_djambal_PDJAM_T00012360\|Pangasianodon_hypophthalmus_rna_XM_026945431_3\|Onychostoma_macrolepis_rna_XM_058770295_1\|Onychostoma_macrolepis_rna_XM_058768245_1\|Puntigrus_tetrazona_rna_XM_043229772_1\|Cirrhinus_molitorella_evm_model_Chr01_785\|Labeo_rohita_rna_XM_051103582_1\|Ctenopharyngodon_idella_rna_XM_051888784_1\|Anabarilius_grahami_mrna_DPX16_16220\|Megalobrama_amblycephala_rna_XM_048197940_1\|Rhinichthys_klamathensis_goyatoka_rna_XM_056234689_1\|Pimephales_promelas_rna_XM_039687464_1\|Danionella_translucida_DT02135_RA\|Xyrauchen_texanus_rna_XM_052138467_1\|Myxocyprinus_asiaticus_rna_XM_051717578_1\|Xyrauchen_texanus_rna_XM_052131238_1\|Myxocyprinus_asiaticus_rna_XM_051715317_1\|Misgurnus_anguillicaudatus_rna_XM_055194518_1\|Misgurnus_anguillicaudatus_rna_XM_055194515_1"  | sed 's/>//g' > curr_ID.txt
xargs samtools faidx N5.HOG0030200.prot  < curr_ID.txt > HGT_clade.prot
xargs samtools faidx N5.HOG0030200.cds  < curr_ID.txt > HGT_clade.cds

muscle5.1.linux_intel64 -align HGT_clade.prot -output HGT_clade.aln

trimal -in HGT_clade.aln -gt 0.6 -cons 60 -backtrans HGT_clade.cds -out HGT_clade.cds.aln

iqtree -s HGT_clade.cds.aln -st DNA -nt 8 -bb 1000 --redo




##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##=================================== Extract micro-synteny scores between recipient species ==========================================
##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================


grep -B10 -A10 "XM_026945431" GFF3_N5_OGGs/Pangasianodon_hypophthalmus.gff.simplified.sorted.OGG.prez.f145.clustered.tiret | sed 's/.*N5_HOG/N5.HOG/g' | sed 's/,.*//g' | grep "[A-Z]" | grep -v "N5.HOG0030200" > Pangasianodon_hypophthalmus.synt.1.txt

grep -B10 -A10 "T00012360" GFF3_N5_OGGs/Pangasius_djambal.gff.simplified.sorted.OGG.prez.f145.clustered.tiret | sed 's/.*N5_HOG/N5.HOG/g' | sed 's/,.*//g' | grep "[A-Z]" | grep -v "N5.HOG0030200"  > Pangasius_djambal.synt.1.txt
grep -B10 -A10 "T00103720" GFF3_N5_OGGs/Pangasius_djambal.gff.simplified.sorted.OGG.prez.f145.clustered.tiret | sed 's/.*N5_HOG/N5.HOG/g' | sed 's/,.*//g' | grep "[A-Z]" | grep -v "N5.HOG0030200"  > Pangasius_djambal.synt.2.txt


comm -12 <(sort Pangasianodon_hypophthalmus.synt.1.txt) <(sort Pangasius_djambal.synt.1.txt) | wc -l
comm -12 <(sort Pangasianodon_hypophthalmus.synt.1.txt) <(sort Pangasius_djambal.synt.2.txt) | wc -l



