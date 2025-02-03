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


samtools faidx Proteomes_BUSCO80/Lates_calcarifer.fa Lates_calcarifer---rna-XM_018688510.2 > Lates_calcarifer---rna-XM_018688510.2.prot

#blastp against the fish proteome (50 best hits) + Uniprot database without any fish species (10 best hits). Here I will only take 20 catffishes otherwise there are way too much of them

blastp -query Lates_calcarifer---rna-XM_018688510.2.prot -db concatenated_proteomes.fa -outfmt 6 -out Gene_vs_FishProteome.blastp -num_threads 8 -max_target_seqs 51
blastp -query Lates_calcarifer---rna-XM_018688510.2.prot -db non_actino_uniprot.fa -outfmt 6 -out Gene_vs_Uniprot.blastp -num_threads 8  -max_target_seqs 10 -evalue 1e-5

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
##========================================== Now make a synteny plot  ================================================================
##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================


#N5.HOG0004450.clade_hgt.txt=> clade where the HGT is extracted from iTOL (from N5.HOG0004450.prot.aln.treefile)

cut -f1,2 -d "_" N5.HOG0004450.clade_hgt.txt | sort | uniq  > species_to_draw.clade1
#Choose close species without the gene .. 


echo "Engraulis_encrasicolus" > non_transfer_species_1
echo "Alosa_alosa" >> non_transfer_species_1
echo "Sardina_pilchardus" >> non_transfer_species_1

echo "Caranx_melampygus" > non_transfer_species_2
echo "Solea_solea" >> non_transfer_species_2
echo "Betta_splendens" >> non_transfer_species_2




rm -r Syn_tables_dir ; mkdir Syn_tables_dir


#First order the species names 

rm species_to_draw.clade1.ordered

echo "Clupea_harengus" >> species_to_draw.clade1.ordered
echo "Lates_calcarifer" >> species_to_draw.clade1.ordered


#If there are more than 4 species from the same order, keep only 4 


#First extract genes around 

curr_OGG=N5.HOG0004450

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


grep "NW_026115662_1\|NW_026115863" Syn_tables_dir/Lates_calcarifer.synt.final.df > temp ; mv temp Syn_tables_dir/Lates_calcarifer.synt.final.df

#First add Sardina_pilchardus

curr_OGG=N5.HOG0004450
curr_sp=Sardina_pilchardus
ref_sp=Clupea_harengus

grep -A10 -B10 "N5_HOG0004450" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0034945" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df


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



#Now add Alosa_alosa

curr_OGG=N5.HOG0004450
curr_sp=Alosa_alosa
ref_sp=Clupea_harengus

grep -A10 -B10 "N5_HOG0004450" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "NC_063189\|NC_063191" > Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0046016" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df
grep -A5 -B5 "N5_HOG0054012" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df


sort Syn_tables_dir/$curr_sp.synt.df | uniq > temp ; mv temp Syn_tables_dir/$curr_sp.synt.df

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

curr_OGG=N5.HOG0004450
curr_sp=Engraulis_encrasicolus
ref_sp=Clupea_harengus

grep -A10 -B10 "N5_HOG0004450" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "NW_026945808\|NW_026946109" > Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0046016" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df


sort Syn_tables_dir/$curr_sp.synt.df | uniq > temp ; mv temp Syn_tables_dir/$curr_sp.synt.df

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





#Now add Caranx_melampygus

curr_OGG=N5.HOG0004450
curr_sp=Caranx_melampygus
ref_sp=Lates_calcarifer

grep -A5 -B5 "N5_HOG0004450" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret  > Syn_tables_dir/$curr_sp.synt.df
grep -A5 -B5 "N5_HOG0033160" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "JAFELL010001300"  >> Syn_tables_dir/$curr_sp.synt.df
grep -A5 -B5 "N5_HOG0004482" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "JAFELL010003188"  >> Syn_tables_dir/$curr_sp.synt.df


sort Syn_tables_dir/$curr_sp.synt.df | uniq > temp ; mv temp Syn_tables_dir/$curr_sp.synt.df

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




#Now add Solea_solea

curr_OGG=N5.HOG0004450
curr_sp=Solea_solea
ref_sp=Lates_calcarifer

grep -A5 -B5 "N5_HOG0005227" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df
grep -A5 -B5 "N5_HOG0046157" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df
grep -A5 -B5 "N5_HOG0026425" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "NC_081153" >> Syn_tables_dir/$curr_sp.synt.df
grep -A5 -B5 "N5_HOG0003377" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df
grep -A5 -B5 "N5_HOG0031761" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df




sort Syn_tables_dir/$curr_sp.synt.df | uniq > temp ; mv temp Syn_tables_dir/$curr_sp.synt.df

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





#Finally, add Betta_splendens

curr_OGG=N5.HOG0004450
curr_sp=Betta_splendens
ref_sp=Lates_calcarifer

grep -A5 -B5 "N5_HOG0005227" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df
grep -A5 -B5 "N5_HOG0046157" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df
grep -A5 -B5 "N5_HOG0026425" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "NC_040891" >> Syn_tables_dir/$curr_sp.synt.df
grep -A5 -B5 "N5_HOG0003377" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df
grep -A5 -B5 "N5_HOG0031761" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df




sort Syn_tables_dir/$curr_sp.synt.df | uniq > temp ; mv temp Syn_tables_dir/$curr_sp.synt.df

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

exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs" Lates_calcarifer---rna-XM_018688510.2.prot GCF_958295425.1_fSolSol10.1_genomic.fna > Solea_solea.whole.exo
exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs" Lates_calcarifer---rna-XM_018688510.2.prot GCF_900634795.4_fBetSpl5.4_genomic.fna > Betta_splendens.whole.exo
exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs" Lates_calcarifer---rna-XM_018688510.2.prot GCF_016859285.1_ASM1685928v1_genomic.fna > Xiphias_gladius.whole.exo
exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs" Lates_calcarifer---rna-XM_018688510.2.prot GCF_002814215.2_Sedor1_genomic.fna > test.pos.whole.exo
exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs" Lates_calcarifer---rna-XM_018688510.2.prot GCA_022829145.1_AGOR_1.0_genomic.fna > Albula.whole.exo


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

#Caranx_melampygus
grep -A3 -B3 "N5_HOG0004450"  GFF3_N5_OGGs/Caranx_melampygus.gff.simplified.sorted.OGG.tiret

samtools faidx GCA_019059645.1_BYU_Cmel_1.0_genomic.fna JAFELL010000420.1:7875-23713 > N5.HOG0004450.Caranx_melampygus.extended.1.fa
samtools faidx GCA_019059645.1_BYU_Cmel_1.0_genomic.fna JAFELL010002116.1:27006-41038 > N5.HOG0004450.Caranx_melampygus.extended.2.fa
samtools faidx GCA_019059645.1_BYU_Cmel_1.0_genomic.fna JAFELL010002462.1:248780-261895 > N5.HOG0004450.Caranx_melampygus.extended.3.fa

sed -i 's/:/-/g' N5.HOG0004450.Caranx_melampygus.extended.1.fa
makeblastdb -in N5.HOG0004450.Caranx_melampygus.extended.1.fa -dbtype nucl

sed -i 's/:/-/g' N5.HOG0004450.Caranx_melampygus.extended.2.fa
makeblastdb -in N5.HOG0004450.Caranx_melampygus.extended.2.fa -dbtype nucl

sed -i 's/:/-/g' N5.HOG0004450.Caranx_melampygus.extended.3.fa
makeblastdb -in N5.HOG0004450.Caranx_melampygus.extended.3.fa -dbtype nucl


#Sardina_pilchardus
grep -A3 -B3 "N5_HOG0004450"  GFF3_N5_OGGs/Sardina_pilchardus.gff.simplified.sorted.OGG.tiret
samtools faidx GCF_963854185.1_fSarPil1.1_genomic.fna NC_085013.1:30433942-30446476 > N5.HOG0004450.Sardina_pilchardus.extended.1.fa
samtools faidx GCF_963854185.1_fSarPil1.1_genomic.fna NC_085013.1:30501317-30527634 > N5.HOG0004450.Sardina_pilchardus.extended.2.fa

sed -i 's/:/-/g' N5.HOG0004450.Sardina_pilchardus.extended.1.fa
makeblastdb -in N5.HOG0004450.Sardina_pilchardus.extended.1.fa -dbtype nucl


sed -i 's/:/-/g' N5.HOG0004450.Sardina_pilchardus.extended.2.fa
makeblastdb -in N5.HOG0004450.Sardina_pilchardus.extended.2.fa -dbtype nucl




#Now lets find TEs using the repbase and dfam databases

tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0004450.Caranx_melampygus.extended.1.fa -evalue 1e-5 -outfmt 6 -out TE.Caranx_melampygus.1.tblastn -num_threads 8
sed -i 's/#//g' TE.Caranx_melampygus.1.tblastn
tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0004450.Caranx_melampygus.extended.2.fa -evalue 1e-5 -outfmt 6 -out TE.Caranx_melampygus.2.tblastn -num_threads 8
sed -i 's/#//g' TE.Caranx_melampygus.2.tblastn
tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0004450.Caranx_melampygus.extended.3.fa -evalue 1e-5 -outfmt 6 -out TE.Caranx_melampygus.3.tblastn -num_threads 8
sed -i 's/#//g' TE.Caranx_melampygus.3.tblastn


tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0004450.Sardina_pilchardus.extended.1.fa -evalue 1e-5 -outfmt 6 -out TE.Sardina_pilchardus.1.tblastn -num_threads 8
sed -i 's/#//g' TE.Sardina_pilchardus.1.tblastn
tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0004450.Sardina_pilchardus.extended.2.fa -evalue 1e-5 -outfmt 6 -out TE.Sardina_pilchardus.2.tblastn -num_threads 8
sed -i 's/#//g' TE.Sardina_pilchardus.2.tblastn

#merge tblastn hits and find the best TE match by doing a blastx

Rscript Rscript_merge_blast_hits.R TE.Caranx_melampygus.1.tblastn TE.Caranx_melampygus.1.tblastn.merged
Rscript Rscript_merge_blast_hits.R TE.Caranx_melampygus.2.tblastn TE.Caranx_melampygus.2.tblastn.merged
Rscript Rscript_merge_blast_hits.R TE.Caranx_melampygus.3.tblastn TE.Caranx_melampygus.3.tblastn.merged

Rscript Rscript_merge_blast_hits.R TE.Sardina_pilchardus.1.tblastn TE.Sardina_pilchardus.1.tblastn.merged
Rscript Rscript_merge_blast_hits.R TE.Sardina_pilchardus.2.tblastn TE.Sardina_pilchardus.2.tblastn.merged


xargs samtools faidx N5.HOG0004450.Caranx_melampygus.extended.1.fa < TE.Caranx_melampygus.1.tblastn.merged > TE.Caranx_melampygus.1.BEST.fa
blastx -query TE.Caranx_melampygus.1.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Caranx_melampygus.1.BEST.blastx -max_target_seqs 1
cut -f1 TE.Caranx_melampygus.1.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Caranx_melampygus.1.BEST.blastx >> temp  ; done ; mv temp TE.Caranx_melampygus.1.BEST.blastx

xargs samtools faidx N5.HOG0004450.Caranx_melampygus.extended.2.fa < TE.Caranx_melampygus.2.tblastn.merged > TE.Caranx_melampygus.2.BEST.fa
blastx -query TE.Caranx_melampygus.2.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Caranx_melampygus.2.BEST.blastx -max_target_seqs 1
cut -f1 TE.Caranx_melampygus.2.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Caranx_melampygus.2.BEST.blastx >> temp  ; done ; mv temp TE.Caranx_melampygus.2.BEST.blastx

xargs samtools faidx N5.HOG0004450.Caranx_melampygus.extended.3.fa < TE.Caranx_melampygus.3.tblastn.merged > TE.Caranx_melampygus.3.BEST.fa
blastx -query TE.Caranx_melampygus.3.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Caranx_melampygus.3.BEST.blastx -max_target_seqs 1
cut -f1 TE.Caranx_melampygus.3.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Caranx_melampygus.3.BEST.blastx >> temp  ; done ; mv temp TE.Caranx_melampygus.3.BEST.blastx



xargs samtools faidx N5.HOG0004450.Sardina_pilchardus.extended.1.fa < TE.Sardina_pilchardus.1.tblastn.merged > TE.Sardina_pilchardus.1.BEST.fa
blastx -query TE.Sardina_pilchardus.1.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Sardina_pilchardus.1.BEST.blastx -max_target_seqs 1
cut -f1  TE.Sardina_pilchardus.1.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Sardina_pilchardus.1.BEST.blastx >> temp  ; done ; mv temp TE.Sardina_pilchardus.1.BEST.blastx

xargs samtools faidx N5.HOG0004450.Sardina_pilchardus.extended.2.fa < TE.Sardina_pilchardus.2.tblastn.merged > TE.Sardina_pilchardus.2.BEST.fa
blastx -query TE.Sardina_pilchardus.2.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Sardina_pilchardus.2.BEST.blastx -max_target_seqs 1
cut -f1  TE.Sardina_pilchardus.2.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Sardina_pilchardus.2.BEST.blastx >> temp  ; done ; mv temp TE.Sardina_pilchardus.2.BEST.blastx


#Now find shared elements

cut -f2 TE.Caranx_melampygus.1.BEST.blastx | sort | uniq > TE.Caranx_melampygus.1.uniqTE
cut -f2 TE.Caranx_melampygus.2.BEST.blastx | sort | uniq > TE.Caranx_melampygus.2.uniqTE
cut -f2 TE.Caranx_melampygus.3.BEST.blastx | sort | uniq > TE.Caranx_melampygus.3.uniqTE

cut -f2 TE.Sardina_pilchardus.1.BEST.blastx | sort | uniq > TE.Sardina_pilchardus.1.uniqTE
cut -f2 TE.Sardina_pilchardus.2.BEST.blastx | sort | uniq > TE.Sardina_pilchardus.2.uniqTE


comm -12 TE.Caranx_melampygus.1.uniqTE TE.Sardina_pilchardus.1.uniqTE
comm -12 TE.Caranx_melampygus.2.uniqTE TE.Sardina_pilchardus.1.uniqTE
comm -12 TE.Caranx_melampygus.3.uniqTE TE.Sardina_pilchardus.1.uniqTE

comm -12 TE.Caranx_melampygus.1.uniqTE TE.Sardina_pilchardus.2.uniqTE
comm -12 TE.Caranx_melampygus.2.uniqTE TE.Sardina_pilchardus.2.uniqTE
comm -12 TE.Caranx_melampygus.3.uniqTE TE.Sardina_pilchardus.2.uniqTE

## No shared TEs




### In a FINAL step we will make another synteny plot more zoomed


## Sardina_pilchardus

rm seq_clustered_infos_ogg.TE.txt ; rm clusters_ID_TE.txt


scaffold=`grep ">" N5.HOG0004450.Sardina_pilchardus.extended.1.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0004450.Sardina_pilchardus.extended.1.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0004450.Sardina_pilchardus.extended.1.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Sardina_pilchardus,$scaffold,$length" > clusters_ID_TE.txt

grep "N5_HOG0004450"  GFF3_N5_OGGs/Sardina_pilchardus.gff.simplified.sorted.OGG.tiret

grep "XM_062523630" GFF3_files_per_species/Sardina_pilchardus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons


nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene1.$nbexon,N5_HOG0004450,+,Sardina_pilchardus" >> seq_clustered_infos_ogg.TE.txt
done



## Sardina_pilchardus ### second scaffold


scaffold=`grep ">" N5.HOG0004450.Sardina_pilchardus.extended.2.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0004450.Sardina_pilchardus.extended.2.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0004450.Sardina_pilchardus.extended.2.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Sardina_pilchardus,$scaffold.sec,$length" >> clusters_ID_TE.txt

grep "N5_HOG0004450"  GFF3_N5_OGGs/Sardina_pilchardus.gff.simplified.sorted.OGG.tiret

grep "XM_062523290" GFF3_files_per_species/Sardina_pilchardus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons


nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold.sec,$real_exon_start,$real_exon_stop,gene2.$nbexon,N5_HOG0004450,+,Sardina_pilchardus" >> seq_clustered_infos_ogg.TE.txt
done

cut -f1 TE.Sardina_pilchardus.2.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Sardina_pilchardus.2.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0004450.Sardina_pilchardus.extended.2.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0004450.Sardina_pilchardus.extended.2.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold.sec,/g" | sed "s/$/,$strand,Sardina_pilchardus/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done



## Caranx_melampygus #1



scaffold=`grep ">" N5.HOG0004450.Caranx_melampygus.extended.1.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0004450.Caranx_melampygus.extended.1.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0004450.Caranx_melampygus.extended.1.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Caranx_melampygus,$scaffold,$length" >> clusters_ID_TE.txt

grep "N5_HOG0004450"  GFF3_N5_OGGs/Caranx_melampygus.gff.simplified.sorted.OGG.tiret

grep "mRNA25461" GFF3_files_per_species/Caranx_melampygus.gff | grep "CDS	" | cut -f4,5 | tr "\t" "," > gene_1.exons


nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene1.$nbexon,N5_HOG0004450,-,Caranx_melampygus" >> seq_clustered_infos_ogg.TE.txt
done


cut -f1 TE.Caranx_melampygus.1.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Caranx_melampygus.1.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0004450.Caranx_melampygus.extended.1.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0004450.Caranx_melampygus.extended.1.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold,/g" | sed "s/$/,$strand,Caranx_melampygus/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done



## Caranx_melampygus #2



scaffold=`grep ">" N5.HOG0004450.Caranx_melampygus.extended.2.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0004450.Caranx_melampygus.extended.2.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0004450.Caranx_melampygus.extended.2.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Caranx_melampygus,$scaffold,$length" >> clusters_ID_TE.txt

grep "N5_HOG0004450"  GFF3_N5_OGGs/Caranx_melampygus.gff.simplified.sorted.OGG.tiret

grep "mRNA11584" GFF3_files_per_species/Caranx_melampygus.gff | grep "CDS	" | cut -f4,5 | tr "\t" "," > gene_1.exons


nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene2.$nbexon,N5_HOG0004450,-,Caranx_melampygus" >> seq_clustered_infos_ogg.TE.txt
done





## Caranx_melampygus #3



scaffold=`grep ">" N5.HOG0004450.Caranx_melampygus.extended.3.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0004450.Caranx_melampygus.extended.3.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0004450.Caranx_melampygus.extended.3.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Caranx_melampygus,$scaffold,$length" >> clusters_ID_TE.txt

grep "N5_HOG0004450"  GFF3_N5_OGGs/Caranx_melampygus.gff.simplified.sorted.OGG.tiret

grep "mRNA27263" GFF3_files_per_species/Caranx_melampygus.gff | grep "CDS	" | cut -f4,5 | tr "\t" "," > gene_1.exons


nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene3.$nbexon,N5_HOG0004450,+,Caranx_melampygus" >> seq_clustered_infos_ogg.TE.txt
done


cut -f1 TE.Caranx_melampygus.3.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Caranx_melampygus.3.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0004450.Caranx_melampygus.extended.3.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0004450.Caranx_melampygus.extended.3.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold,/g" | sed "s/$/,$strand,Caranx_melampygus/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done





### Check the conserved non-coding region

samtools faidx N5.HOG0004450.Caranx_melampygus.extended.3.fa JAFELL010002462.1-248780-261895:6641-7000 > region1.Caranx_melampygus.fa
samtools faidx N5.HOG0004450.Sardina_pilchardus.extended.1.fa NC_085013.1-30433942-30446476:12243-12643 > region1.Sardina_pilchardus.fa

sed -i 's/>.*/>C.gerrardi/g' region1.Caranx_melampygus.fa
sed -i 's/>.*/>M.cyprinoides/g' region1.Sardina_pilchardus.fa

needle -asequence region1.Caranx_melampygus.fa -bsequence region1.Sardina_pilchardus.fa -outfile region1.aln -gapopen 10.0 -gapextend 0.5

sbatch --qos=6hours -c 10 --mem=20G launch_blast_allgenomes.sh region1.Caranx_melampygus.fa 1e-20 region1.blastn.tsv





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
sbatch --qos=1day -c 4 --mem=10G -e error.fitMG4.out -o slurm.fitMG4.out launch_fitMG4.sh N5.HOG0004450

Recipient branches : 
Node381
Node377
Lates_calcarifer_rna_XM_018688510_2
Lates_calcarifer_rna_XM_018686387_2
Lates_japonicus_gene_AKAME5_002034700
Seriola_lalandi_dorsalis_rna_XM_023399112_1
Caranx_melampygus_mRNA27263
Caranx_melampygus_mRNA11584
Caranx_melampygus_mRNA25461




#Test positive selection and relaxed selection on receiver branch

sbatch --qos=1week -c 4 --mem=10G -e error.absrel.cand.out -o slurm.absrel.cand.out --job-name=HOG0004450 launch_absrel_cand.sh N5.HOG0004450



#Launch RELAX

### Adaptive branch site random effects likelihood test 
### Likelihood ratio test for episodic diversifying positive selection at Holm-Bonferroni corrected _p =   0.0500_ found **1** branches under selection among **10** tested.
### * Node377, p-value =  0.00000


===> $HOG.prot.aln.treefile.SelecMarked.relax => select branches detected with accelerated evolution
sbatch --qos=1week -c 6 --mem=10G -e error.relax.out -o slurm.relax.out --job-name=HOG0004450 launch_RELAX.sh N5.HOG0004450




#Extract dN/dS to table

grep "LB\":" N5.HOG0004450.cds.aln.FITTER.json | sed 's/.*://g' | sed 's/,.*//g' > curr_LB_values.txt
grep "MLE\":" N5.HOG0004450.cds.aln.FITTER.json | sed 's/.*://g' | sed 's/,.*//g' > curr_MLE_values.txt
grep "UB\":" N5.HOG0004450.cds.aln.FITTER.json | sed 's/.*://g' | sed 's/,.*//g' > curr_UB_values.txt
grep "\"dN\"" N5.HOG0004450.cds.aln.FITTER.json  | grep -v "{" | sed 's/.*://g' | sed 's/,//g' > curr_dN_values.txt
grep "\"dS\"" N5.HOG0004450.cds.aln.FITTER.json  | grep -v "{" | sed 's/.*://g' | sed 's/,//g' > curr_dS_values.txt
grep -B2 "LB\":" N5.HOG0004450.cds.aln.FITTER.json | grep -v "\-\-" | grep -v "Confidence" | grep -v "LB\":"  | sed 's/\"//g' | sed 's/:.*//g' | sed 's/^ *//g' > curr_labels

paste -d "," curr_labels curr_LB_values.txt curr_MLE_values.txt curr_UB_values.txt curr_dN_values.txt curr_dS_values.txt > N5.HOG0004450.dN_dS.csv





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


grep ">" N5.HOG0004450.prot | grep "Alosa_sapidissima_rna_XM_042093290_1\|Alosa_alosa_rna_XM_048256059_1\|Sardina_pilchardus_rna_XM_062523630_1\|Sardina_pilchardus_rna_XM_062523290_1\|Alosa_sapidissima_rna_XM_042103484_1\|Alosa_alosa_rna_XM_048250590_1\|Alosa_alosa_rna_XM_048239489_1\|Clupea_harengus_rna_XM_042704194_1\|Lates_calcarifer_rna_XM_018688510_2\|Lates_calcarifer_rna_XM_018686387_2\|Lates_japonicus_gene_AKAME5_002034700\|Seriola_lalandi_dorsalis_rna_XM_023399112_1\|Caranx_melampygus_mRNA27263\|Caranx_melampygus_mRNA11584\|Caranx_melampygus_mRNA25461\|Engraulis_encrasicolus_rna_XM_063193640_1\|Engraulis_encrasicolus_rna_XM_063193637_1\|Engraulis_encrasicolus_rna_XM_063194025_1\|Engraulis_encrasicolus_rna_XM_063193635_1"  | sed 's/>//g' > curr_ID.txt
xargs samtools faidx N5.HOG0004450.prot  < curr_ID.txt > HGT_clade.prot
xargs samtools faidx N5.HOG0004450.cds  < curr_ID.txt > HGT_clade.cds

muscle5.1.linux_intel64 -align HGT_clade.prot -output HGT_clade.aln

trimal -in HGT_clade.aln -gt 0.6 -cons 60 -backtrans HGT_clade.cds -out HGT_clade.cds.aln

iqtree -s HGT_clade.cds.aln -st DNA -nt 8 -bb 1000 --redo











