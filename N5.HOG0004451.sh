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


samtools faidx Proteomes_BUSCO80/Centroberyx_gerrardi.fa Centroberyx_gerrardi---g23848.t1_1 > Centroberyx_gerrardi---g23848.t1_1.prot

#blastp against the fish proteome (50 best hits) + Uniprot database without any fish species (10 best hits). Here I will only take 20 catffishes otherwise there are way too much of them

blastp -query Centroberyx_gerrardi---g23848.t1_1.prot -db concatenated_proteomes.fa -outfmt 6 -out Gene_vs_FishProteome.blastp -num_threads 8 -max_target_seqs 53
blastp -query Centroberyx_gerrardi---g23848.t1_1.prot -db non_actino_uniprot.fa -outfmt 6 -out Gene_vs_Uniprot.blastp -num_threads 8  -max_target_seqs 10 -evalue 1e-5


cut -f2 Gene_vs_FishProteome.blastp | sort | uniq  > closest_fish_seq.id
cut -f2 Gene_vs_Uniprot.blastp | sort | uniq > closest_nonfish_seq.id 

xargs samtools faidx concatenated_proteomes.fa < closest_fish_seq.id > closest_fish_seq.fa
xargs samtools faidx non_actino_uniprot.fa  < closest_nonfish_seq.id > closest_nonfish_seq.fa

#Align with muscle and trim with trimal

cat closest_fish_seq.fa closest_nonfish_seq.fa > Uniprot_plus_closefish.fa

muscle5.1 -align closest_fish_seq.fa -output closest_fish_seq.aln
cp closest_fish_seq.aln Uniprot_plus_closefish.aln
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


#N5.HOG0004451.clade_hgt.txt=> clade where the HGT is extracted from iTOL (from N5.HOG0004451.prot.aln.treefile)

cut -f1,2 -d "_" N5.HOG0004451.clade_hgt.txt | sort | uniq  > species_to_draw.clade1
#Choose close species without the gene .. 


echo "Anguilla_anguilla" > non_transfer_species_1
echo "Albula_glossodonta" >> non_transfer_species_1
echo "Aldrovandia_affinis" >> non_transfer_species_1

echo "Myripristis_murdjan" > non_transfer_species_2
echo "Hoplostethus_atlanticus" >> non_transfer_species_2
echo "Lampris_incognitus" >> non_transfer_species_2




rm -r Syn_tables_dir ; mkdir Syn_tables_dir


#First order the species names 

rm species_to_draw.clade1.ordered

echo "Megalops_cyprinoides" >> species_to_draw.clade1.ordered
echo "Centroberyx_gerrardi" >> species_to_draw.clade1.ordered



#If there are more than 4 species from the same order, keep only 4 


#First extract genes around 

curr_OGG=N5.HOG0004451

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

grep -v "NC_050585_1" Syn_tables_dir/Megalops_cyprinoides.synt.final.df > temp ; mv temp Syn_tables_dir/Megalops_cyprinoides.synt.final.df 


#First add Aldrovandia_affinis

curr_OGG=N5.HOG0004451
curr_sp=Aldrovandia_affinis
ref_sp=Megalops_cyprinoides

grep -A10 -B10 "N5_HOG0041299" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df
grep -A5 -B5 "N5_HOG0004451" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret  >> Syn_tables_dir/$curr_sp.synt.df
grep -A5 -B5 "N5_HOG0001304" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df




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



#Now add Albula_glossodonta

curr_OGG=N5.HOG0004451
curr_sp=Albula_glossodonta
ref_sp=Megalops_cyprinoides

grep -A10 -B10 "N5_HOG0004451" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "JAFBMS010000030" > Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0004451" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "JAFBMS010000145" >> Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0041299" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df




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




#Now add Anguilla_anguilla

curr_OGG=N5.HOG0004451
curr_sp=Anguilla_anguilla
ref_sp=Megalops_cyprinoides

grep -A10 -B10 "N5_HOG0004451" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0041299" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df
grep -A5 -B5 "N5_HOG0051296" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df


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



#Now add Myripristis_murdjan

curr_OGG=N5.HOG0004451
curr_sp=Myripristis_murdjan
ref_sp=Centroberyx_gerrardi

grep -A10 -B10 "N5_HOG0014892" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0035294" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df


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




#Now add Hoplostethus_atlanticus

curr_OGG=N5.HOG0004451
curr_sp=Hoplostethus_atlanticus
ref_sp=Centroberyx_gerrardi

grep -A10 -B10 "N5_HOG0042201" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0035294" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df


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



#Finally add Lampris_incognitus

curr_OGG=N5.HOG0004451
curr_sp=Lampris_incognitus
ref_sp=Centroberyx_gerrardi

grep -A10 -B10 "N5_HOG0042201" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0035294" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df

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

exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs" Centroberyx_gerrardi---g23848.t1_1.prot GCF_902150065.1_fMyrMur1.1_genomic.fna > Myripristis_murdjan.whole.exo
exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs" Centroberyx_gerrardi---g23848.t1_1.prot GCA_034670405.1_ASM3467040v1_genomic.fna > Hoplostethus_atlanticus.whole.exo
exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs" Centroberyx_gerrardi---g23848.t1_1.prot GCF_029633865.1_fLamInc1.hap2_genomic.fna > Lampris_incognitus.whole.exo



## Search the gene in the related species Beryx_splendens

exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs" Centroberyx_gerrardi---g23848.t1_1.prot GCA_900312565.1_ASM90031256v1_genomic.fna > Beryx_splendens.whole.exo

#Beryx_splendens.HOG0004451.fa
>Beryx_splendens.OMOV01022145.1-1931-5459
ATGGACTTGAAGGACATTGACATGCCGGCCCTCGGACGCCCCGTGCACATCGGCATGCTGTACAACATTT
TCACTGACAACTTCCTCTCAGGAGTCTTACTGTGGGATGAAGTGACCATTAAGAAAAGCACACATGTGGA
CCAAAGTCCCAGCACAAGGTTTACGATAACTGCGGAGGACTCCCTCACTGAGAAGACCAACATGCTGAAC
GTAAGCGCCTCCCTGAAGGCCAGTTTCCTGGGGGGgctggtggaggtgggggggtcGGCCAGTTACCTGA
AGAACAGAGTGTCCTCCGTAAAGCAATGCAGGGTCACAATGAAATACGAACAGGAAACTGAGCACAAGTC
TATACACTTTGATGAGCTGGGACAGGTGACGCACCCCGAGGTGTTTGAGCAGCAGAACGCCACCCACGTT
GTCATAGGAGTGACGTACGGAGCCAAGGCCTTCATGGTCTTCGATCGGACGACTTCAGATGTAGAAAATA
AGCAAGATATTCAGGGGCAACTGGGTGTGATGGTCAAGAAGATCCCTTCAACGGAAATCAGCGGTGAAGG
AAAAGTGGCCCTGACAGACGAGGACaaggagaaggtggagaggTTCACGTGCACATTCTTCGGAGACTTC
AGACTGGATCACAACCCCTCCACTTATGAGGAAGCTGTACTGGTGTACAAAAGCCTCCCAAAGCTGctgg
gggagaagggagagaaagcagtGCCTGTGAAAATGCTTCTGCATCCTCTGAAGAAGCTGGACCCCAAAGC
TGCCACGCTGGTGAGAGAAATCAGACCGGAGCTGGTGTCAGATGTGGAGGATGTTATGGAGGAGCTTCAT
AAGGCCAAAATTAGAGCCAATGACCTGATTAGCCAGTGCAAGCTCATCAAGGCTGCTGCCATAGAGAAAA
AACTGACAACATTCGAGGACGAGCTGGACGATTACACCGTCCCTCTCCAGAAGAATCTGAGCAAGGTCCT
GCCAGCCATCCGGGCGGGTACGGAGGATGAGCAGAGGCTGGAAAACATCCTGAAGTTCCACAATGAATCA
TCCTTCAGCTTCAAAAAGATGAATAAATGGCTGGATGAGAAAGAGACGGAAATACGGGTGCTGAGAACCC
ACATCGGTCCTCTTTGCAGCAACCCTGCAGTTTCCATTGCCCCTCCAGGCCCAGAGCTCGACACTTTCAT
CAGCAATTCCAGAACTGATTGGCTGTATGTGCTGAACTTCACCTCACTGGACTGCAAGGAGGCCTACCTG
TCCACCCTGTCTGAGTGCAAGGCCATGGACGAGTTCAAGAAAATGGAGGACATCTCTGTAGCTCATCACT
GCAGCGCCACAGAGACCCGGCCGTGGTACAAAGACCCCCAGTTCAAAGACAGGCTGAATGGCGTTTTGCA
TAATTTTGAGTATTATGTAGGGAGTTACAAAACCATCATCAGCTACATGCCAGACGCTGTGAAACCTGGA
GTCAGTCTGTACCTGTACAATTATGGGAAACTGTCCAGT


transeq Beryx_splendens.HOG0004451.fa Beryx_splendens.HOG0004451.prot ; sed -i 's/_1$//g' Beryx_splendens.HOG0004451.prot
grep "Albula_glossodonta_mRNA15989\|Albula_glossodonta_mRNA15988\|Albula_glossodonta_mRNA15990\|Albula_goreensis_AGOR_T00122240\|Albula_glossodonta_mRNA7121\|Aldrovandia_affinis_mrna_AAFF_T00314980\|Aldrovandia_affinis_mrna_AAFF_T00314970\|Aldrovandia_affinis_mrna_AAFF_T00102350\|Centroberyx_gerrardi_g23848_t1\|Megalops_cyprinoides_rna_XM_036520634_1\|Conger_conger_rna_XM_061230612_1\|Conger_conger_rna_XM_061230611_1\|Conger_conger_rna_XM_061230614_1\|Conger_conger_rna_XM_061230613_1\|Conger_conger_rna_XM_061232904_1\|Conger_conger_rna_XM_061233141_1\|Conger_conger_rna_XM_061232361_1\|Conger_conger_rna_XM_061230610_1\|Synaphobranchus_kaupii_mrna_SKAU_T00116890\|Synaphobranchus_kaupii_mrna_SKAU_T00116910\|Synaphobranchus_kaupii_mrna_SKAU_T00116870\|Synaphobranchus_kaupii_mrna_SKAU_T00116840\|Synaphobranchus_kaupii_mrna_SKAU_T00116810\|Synaphobranchus_kaupii_mrna_SKAU_T00116900\|Synaphobranchus_kaupii_mrna_SKAU_T00116880\|Synaphobranchus_kaupii_mrna_SKAU_T00116920\|Synaphobranchus_kaupii_mrna_SKAU_T00116860\|Synaphobranchus_kaupii_mrna_SKAU_T00116790\|Synaphobranchus_kaupii_mrna_SKAU_T00116770\|Synaphobranchus_kaupii_mrna_SKAU_T00116800\|Anguilla_anguilla_rna_XM_035407033_1\|Anguilla_anguilla_rna_XM_035403412_1\|Anguilla_anguilla_rna_XM_035407034_1" Coding_sequences_alignments/N5.HOG0004451.prot | sed 's/>//g' > curr_seq.id

xargs samtools faidx Coding_sequences_alignments/N5.HOG0004451.prot < curr_seq.id > curr_seq.prot 


cat Beryx_splendens.HOG0004451.prot curr_seq.prot > HOG0004451.clade.Beryx.prot


muscle5.1 -align HOG0004451.clade.Beryx.prot -output HOG0004451.clade.Beryx.aln
trimal -in HOG0004451.clade.Beryx.aln -gt 0.7 -out HOG0004451.clade.Beryx.aln.trimmed

#Make a ML tree with IQ-TREE
iqtree -s HOG0004451.clade.Beryx.aln.trimmed -st AA -nt 8 -m TEST -mrate G4 --redo


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

#Centroberyx_gerrardi
grep -A3 -B3 "N5_HOG0004451"  GFF3_N5_OGGs/Centroberyx_gerrardi.gff.simplified.sorted.OGG.tiret

samtools faidx GCA_026898505.1_AGI_CSIRO_Cger_v1_genomic.fna JAPMTB010005118.1:498411-512515 > N5.HOG0004451.Centroberyx_gerrardi.extended.1.fa
sed -i 's/:/-/g' N5.HOG0004451.Centroberyx_gerrardi.extended.1.fa
makeblastdb -in N5.HOG0004451.Centroberyx_gerrardi.extended.1.fa -dbtype nucl


#Megalops_cyprinoides
grep -A3 -B3 "N5_HOG0004451"  GFF3_N5_OGGs/Megalops_cyprinoides.gff.simplified.sorted.OGG.tiret
samtools faidx GCF_013368585.1_fMegCyp1.pri_genomic.fna NC_050584.1:69815062-69826893 > N5.HOG0004451.Megalops_cyprinoides.extended.1.fa

sed -i 's/:/-/g' N5.HOG0004451.Megalops_cyprinoides.extended.1.fa
makeblastdb -in N5.HOG0004451.Megalops_cyprinoides.extended.1.fa -dbtype nucl





#Now lets find TEs using the repbase and dfam databases

tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0004451.Centroberyx_gerrardi.extended.1.fa -evalue 1e-5 -outfmt 6 -out TE.Centroberyx_gerrardi.1.tblastn -num_threads 8
sed -i 's/#//g' TE.Centroberyx_gerrardi.1.tblastn

tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0004451.Megalops_cyprinoides.extended.1.fa -evalue 1e-5 -outfmt 6 -out TE.Megalops_cyprinoides.1.tblastn -num_threads 8
sed -i 's/#//g' TE.Megalops_cyprinoides.1.tblastn

#merge tblastn hits and find the best TE match by doing a blastx

Rscript Rscript_merge_blast_hits.R TE.Centroberyx_gerrardi.1.tblastn TE.Centroberyx_gerrardi.1.tblastn.merged
Rscript Rscript_merge_blast_hits.R TE.Megalops_cyprinoides.1.tblastn TE.Megalops_cyprinoides.1.tblastn.merged


xargs samtools faidx N5.HOG0004451.Centroberyx_gerrardi.extended.1.fa < TE.Centroberyx_gerrardi.1.tblastn.merged > TE.Centroberyx_gerrardi.1.BEST.fa
blastx -query TE.Centroberyx_gerrardi.1.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Centroberyx_gerrardi.1.BEST.blastx -max_target_seqs 1
cut -f1 TE.Centroberyx_gerrardi.1.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Centroberyx_gerrardi.1.BEST.blastx >> temp  ; done ; mv temp TE.Centroberyx_gerrardi.1.BEST.blastx


xargs samtools faidx N5.HOG0004451.Megalops_cyprinoides.extended.1.fa < TE.Megalops_cyprinoides.1.tblastn.merged > TE.Megalops_cyprinoides.1.BEST.fa
blastx -query TE.Megalops_cyprinoides.1.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Megalops_cyprinoides.1.BEST.blastx -max_target_seqs 1
cut -f1  TE.Megalops_cyprinoides.1.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Megalops_cyprinoides.1.BEST.blastx >> temp  ; done ; mv temp TE.Megalops_cyprinoides.1.BEST.blastx


#Now find shared elements

cut -f2 TE.Centroberyx_gerrardi.1.BEST.blastx | sort | uniq > TE.Centroberyx_gerrardi.1.uniqTE
cut -f2 TE.Megalops_cyprinoides.1.BEST.blastx | sort | uniq > TE.Megalops_cyprinoides.1.uniqTE

comm -12 TE.Centroberyx_gerrardi.1.uniqTE TE.Megalops_cyprinoides.1.uniqTE

## No shared TEs

### In a FINAL step we will make another synteny plot more zoomed


## Megalops_cyprinoides

rm seq_clustered_infos_ogg.TE.txt ; rm clusters_ID_TE.txt


scaffold=`grep ">" N5.HOG0004451.Megalops_cyprinoides.extended.1.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0004451.Megalops_cyprinoides.extended.1.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0004451.Megalops_cyprinoides.extended.1.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Megalops_cyprinoides,$scaffold,$length" > clusters_ID_TE.txt

grep "N5_HOG0004451"  GFF3_N5_OGGs/Megalops_cyprinoides.gff.simplified.sorted.OGG.tiret

grep "XM_036520634" GFF3_files_per_species/Megalops_cyprinoides.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons


nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene1.$nbexon,N5_HOG0004451,-,Megalops_cyprinoides" >> seq_clustered_infos_ogg.TE.txt
done

cut -f1 TE.Megalops_cyprinoides.1.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Megalops_cyprinoides.1.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0004451.Megalops_cyprinoides.extended.1.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0004451.Megalops_cyprinoides.extended.1.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold,/g" | sed "s/$/,$strand,Megalops_cyprinoides/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done




## Centroberyx_gerrardi

scaffold=`grep ">" N5.HOG0004451.Centroberyx_gerrardi.extended.1.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0004451.Centroberyx_gerrardi.extended.1.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0004451.Centroberyx_gerrardi.extended.1.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Centroberyx_gerrardi,$scaffold,$length" >> clusters_ID_TE.txt

grep "N5_HOG0004451"  GFF3_N5_OGGs/Centroberyx_gerrardi.gff.simplified.sorted.OGG.tiret

grep "g23848" GFF3_files_per_species/Centroberyx_gerrardi.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons


nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene1.$nbexon,N5_HOG0004451,+,Centroberyx_gerrardi" >> seq_clustered_infos_ogg.TE.txt
done

cut -f1 TE.Centroberyx_gerrardi.1.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Centroberyx_gerrardi.1.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0004451.Centroberyx_gerrardi.extended.1.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0004451.Centroberyx_gerrardi.extended.1.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold,/g" | sed "s/$/,$strand,Centroberyx_gerrardi/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done






### Check the conserved non-coding region

samtools faidx N5.HOG0004451.Centroberyx_gerrardi.extended.1.fa.rev JAPMTB010005118.1-498411-512515:3800-5000 > region1.Centroberyx_gerrardi.fa
samtools faidx N5.HOG0004451.Megalops_cyprinoides.extended.1.fa NC_050584.1-69815062-69826893:4000-5000 > region1.Megalops_cyprinoides.fa

sed -i 's/>.*/>C.gerrardi/g' region1.Centroberyx_gerrardi.fa
sed -i 's/>.*/>M.cyprinoides/g' region1.Megalops_cyprinoides.fa

needle -asequence region1.Centroberyx_gerrardi.fa -bsequence region1.Megalops_cyprinoides.fa -outfile region1.aln -gapopen 10.0 -gapextend 0.5

sbatch --qos=6hours -c 10 --mem=20G launch_blast_allgenomes.sh region1.Centroberyx_gerrardi.fa 1e-2 region1.blastn.tsv 
sbatch --qos=6hours -c 10 --mem=20G launch_blast_allgenomes.sh region1.Megalops_cyprinoides.fa 1e-2 region1.blastn.2.tsv 



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
sbatch --qos=1day -c 4 --mem=10G -e error.fitMG4.out -o slurm.fitMG4.out launch_fitMG4.sh N5.HOG0004451


Recipient branches : 
Centroberyx_gerrardi_g23848_t1

#Test positive selection and relaxed selection on receiver branch

sbatch --qos=1week -c 4 --mem=10G -e error.absrel.cand.out -o slurm.absrel.cand.out --job-name=HOG0004451 launch_absrel_cand.sh N5.HOG0004451


## Not need to launch RELAX - No positive selection detected



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




grep ">" N5.HOG0004451.prot | grep "Synaphobranchus_kaupii_mrna_SKAU_T00116920\|Synaphobranchus_kaupii_mrna_SKAU_T00116880\|Synaphobranchus_kaupii_mrna_SKAU_T00116900\|Synaphobranchus_kaupii_mrna_SKAU_T00116860\|Synaphobranchus_kaupii_mrna_SKAU_T00116840\|Synaphobranchus_kaupii_mrna_SKAU_T00116910\|Synaphobranchus_kaupii_mrna_SKAU_T00116890\|Synaphobranchus_kaupii_mrna_SKAU_T00116870\|Synaphobranchus_kaupii_mrna_SKAU_T00116810\|Synaphobranchus_kaupii_mrna_SKAU_T00116790\|Synaphobranchus_kaupii_mrna_SKAU_T00116770\|Synaphobranchus_kaupii_mrna_SKAU_T00116800\|Anguilla_anguilla_rna_XM_035407033_1\|Anguilla_anguilla_rna_XM_035403412_1\|Anguilla_anguilla_rna_XM_035407034_1\|Conger_conger_rna_XM_061230614_1\|Conger_conger_rna_XM_061230611_1\|Conger_conger_rna_XM_061230613_1\|Conger_conger_rna_XM_061232904_1\|Conger_conger_rna_XM_061233141_1\|Conger_conger_rna_XM_061232361_1\|Conger_conger_rna_XM_061230610_1\|Conger_conger_rna_XM_061230612_1\|Albula_glossodonta_mRNA15989\|Albula_glossodonta_mRNA15988\|Albula_goreensis_AGOR_T00122240\|Albula_glossodonta_mRNA15990\|Albula_glossodonta_mRNA7121\|Aldrovandia_affinis_mrna_AAFF_T00314980\|Aldrovandia_affinis_mrna_AAFF_T00314970\|Aldrovandia_affinis_mrna_AAFF_T00102350\|Megalops_cyprinoides_rna_XM_036520634_1\|Centroberyx_gerrardi_g23848_t1"  | sed 's/>//g' > curr_ID.txt
xargs samtools faidx N5.HOG0004451.prot  < curr_ID.txt > HGT_clade.prot
xargs samtools faidx N5.HOG0004451.cds  < curr_ID.txt > HGT_clade.cds

muscle5.1.linux_intel64 -align HGT_clade.prot -output HGT_clade.aln

trimal -in HGT_clade.aln -gt 0.6 -cons 60 -backtrans HGT_clade.cds -out HGT_clade.cds.aln

iqtree -s HGT_clade.cds.aln -st DNA -nt 8 -bb 1000 --redo





