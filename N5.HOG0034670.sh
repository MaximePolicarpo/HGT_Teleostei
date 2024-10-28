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

samtools faidx Proteomes_BUSCO80/Pangasianodon_hypophthalmus.fa Pangasianodon_hypophthalmus---rna-XM_053236680.1 > Pangasianodon_hypophthalmus---rna-XM_053236680.1.prot

#blastp against the fish proteome (50 best hits) + Uniprot database without any fish species (10 best hits). Here I will only take 20 catffishes otherwise there are way too much of them

blastp -query Pangasianodon_hypophthalmus---rna-XM_053236680.1.prot -db concatenated_proteomes.fa -outfmt 6 -out Gene_vs_FishProteome.blastp -num_threads 8 -max_target_seqs 51
blastp -query Pangasianodon_hypophthalmus---rna-XM_053236680.1.prot -db non_actino_uniprot.fa -outfmt 6 -out Gene_vs_Uniprot.blastp -num_threads 8  -max_target_seqs 10 -evalue 1e-5

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
iqtree -s Uniprot_plus_closefish.aln.trimmed -st AA -nt 8 -m TEST -mrate G4 --redo
iqtree -s Uniprot_plus_closefish.aln -st AA -nt 8 -m TEST -mrate G4 --redo



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


#N5.HOG0034670.clade_hgt.txt=> clade where the HGT is extracted from iTOL (from N5.HOG0034670.prot.aln.treefile)

cut -f1,2 -d "_" N5.HOG0034670.clade_hgt.txt | sort | uniq  > species_to_draw.clade1
#Choose close species without the gene .. 


echo "Silurus_meridionalis" > non_transfer_species_1
echo "Ictalurus_punctatus" >> non_transfer_species_1
echo "Pangasius_djambal" >> non_transfer_species_1

echo "Cyprinus_carpio" > non_transfer_species_2
echo "Onychostoma_macrolepis" >> non_transfer_species_2
echo "Danio_rerio" >> non_transfer_species_2




rm -r Syn_tables_dir ; mkdir Syn_tables_dir


#First order the species names 

rm species_to_draw.clade1.ordered

echo "Pangasianodon_hypophthalmus" >> species_to_draw.clade1.ordered
echo "Carassius_auratus" >> species_to_draw.clade1.ordered


#If there are more than 4 species from the same order, keep only 4 


#First extract genes around 

curr_OGG=N5.HOG0034670

for curr_sp in `cat species_to_draw.clade1.ordered` ; do

  grep "$curr_sp" $curr_OGG.clade_hgt.txt > Syn_tables_dir/gene_names.$curr_sp.txt
  
  echo "" > Syn_tables_dir/$curr_sp.synt.df
  
  for gene in `cat Syn_tables_dir/gene_names.$curr_sp.txt` ; do 

    gene_name=`echo "$gene" | sed "s/$curr_sp\_//g"`
  
    if grep -q "$gene_name" Syn_tables_dir/$curr_sp.synt.df ; then 
      echo "gene already in the table"
    else
      
      grep -A15 -B15 "$gene_name" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df
    
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


grep "NC_039266_1\|NW_020523543_1\|NW_020525806_1\|NW_020528938_1" Syn_tables_dir/Carassius_auratus.synt.final.df > temp ; mv temp Syn_tables_dir/Carassius_auratus.synt.final.df


#First add Pangasius_djambal

curr_OGG=N5.HOG0034670
curr_sp=Pangasius_djambal
ref_sp=Pangasianodon_hypophthalmus

grep -A15 -B15 "N5_HOG0034670" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "CM040983" > Syn_tables_dir/$curr_sp.synt.df
grep -A5 -B5 "N5_HOG0004869" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "CM040975" >> Syn_tables_dir/$curr_sp.synt.df
grep -A5 -B5 "N5_HOG0001048" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "CM040983" >> Syn_tables_dir/$curr_sp.synt.df


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

curr_OGG=N5.HOG0034670
curr_sp=Ictalurus_punctatus
ref_sp=Pangasianodon_hypophthalmus

grep -A7 -B7 "N5_HOG0001048" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret  > Syn_tables_dir/$curr_sp.synt.df
grep -A7 -B7 "N5_HOG0011349" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "NC_030441"  >> Syn_tables_dir/$curr_sp.synt.df


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

curr_OGG=N5.HOG0034670
curr_sp=Silurus_meridionalis
ref_sp=Pangasianodon_hypophthalmus

grep -A7 -B7 "N5_HOG0001048" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret  > Syn_tables_dir/$curr_sp.synt.df
grep -A7 -B7 "N5_HOG0011349" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "NC_060907"  >> Syn_tables_dir/$curr_sp.synt.df

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

curr_OGG=N5.HOG0034670
curr_sp=Cyprinus_carpio
ref_sp=Carassius_auratus

grep -A6 -B6 "N5_HOG0005063" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "NC_056599" > Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B8 "N5_HOG0002942" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "NC_056603"  >> Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B8 "N5_HOG0005798" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "NC_056620" >> Syn_tables_dir/$curr_sp.synt.df
grep -A5 -B5 "N5_HOG0034670" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "NC_056620" >> Syn_tables_dir/$curr_sp.synt.df

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

curr_OGG=N5.HOG0034670
curr_sp=Onychostoma_macrolepis
ref_sp=Carassius_auratus



grep -A6 -B6 "N5_HOG0005063" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "NC_081157" > Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B8 "N5_HOG0002942" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "NC_081161"  >> Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B8 "N5_HOG0005798" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df
grep -A5 -B5 "N5_HOG0034670" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "NC_081178" >> Syn_tables_dir/$curr_sp.synt.df





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




#Now add Danio_rerio

curr_OGG=N5.HOG0034670
curr_sp=Danio_rerio
ref_sp=Carassius_auratus


grep -A10 -B10 "N5_HOG0034670" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "NC_007135" > Syn_tables_dir/$curr_sp.synt.df
grep -A6 -B6 "N5_HOG0000493" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "NC_007115" | awk -F'\t' '$2 > 40000000' >> Syn_tables_dir/$curr_sp.synt.df
grep -A8 -B8 "N5_HOG0021582" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret  >> Syn_tables_dir/$curr_sp.synt.df


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


grep "NC_007115" Syn_tables_dir/$curr_sp.synt.final.df > temp




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

exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs" Pangasianodon_hypophthalmus---rna-XM_053236680.1.prot GCF_014805685.1_ASM1480568v1_genomic.fna > Silurus_meridionalis.whole.exo
exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs" Pangasianodon_hypophthalmus---rna-XM_053236680.1.prot GCF_001660625.3_Coco_2.0_genomic.fna > Ictalurus_punctatus.whole.exo

## No gene !

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
grep -A3 -B3 "N5_HOG0034670" GFF3_N5_OGGs/Carassius_auratus.gff.simplified.sorted.OGG.tiret 

grep -A3 -B3 "XM_026265006\|XM_026243545\|XM_026243633" GFF3_N5_OGGs/Carassius_auratus.gff.simplified.sorted.OGG.tiret 

samtools faidx GCF_003368295.1_ASM336829v1_genomic.fna NW_020523543.1:1402666-1433535 > N5.HOG0034670.Carassius_auratus.extended.1.fa
samtools faidx GCF_003368295.1_ASM336829v1_genomic.fna NW_020523543.1:1497263-1511665 > N5.HOG0034670.Carassius_auratus.extended.2.fa
samtools faidx GCF_003368295.1_ASM336829v1_genomic.fna NW_020528938.1:159523-177872 > N5.HOG0034670.Carassius_auratus.extended.3.fa

sed -i 's/:/-/g' N5.HOG0034670.Carassius_auratus.extended.1.fa
makeblastdb -in N5.HOG0034670.Carassius_auratus.extended.1.fa -dbtype nucl
sed -i 's/:/-/g' N5.HOG0034670.Carassius_auratus.extended.2.fa
makeblastdb -in N5.HOG0034670.Carassius_auratus.extended.2.fa -dbtype nucl
sed -i 's/:/-/g' N5.HOG0034670.Carassius_auratus.extended.3.fa
makeblastdb -in N5.HOG0034670.Carassius_auratus.extended.3.fa -dbtype nucl



#Pangasianodon_hypophthalmus
grep -A3 -B3 "N5_HOG0034670"  GFF3_N5_OGGs/Pangasianodon_hypophthalmus.gff.simplified.sorted.OGG.tiret
samtools faidx GCF_027358585.1_fPanHyp1.pri_genomic.fna NC_069718.1:29758409-29791397 > N5.HOG0034670.Pangasianodon_hypophthalmus.extended.1.fa
samtools faidx GCF_027358585.1_fPanHyp1.pri_genomic.fna NC_069718.1:29954005-29986143 > N5.HOG0034670.Pangasianodon_hypophthalmus.extended.2.fa
samtools faidx GCF_027358585.1_fPanHyp1.pri_genomic.fna NC_069718.1:30023632-30057517 > N5.HOG0034670.Pangasianodon_hypophthalmus.extended.3.fa
samtools faidx GCF_027358585.1_fPanHyp1.pri_genomic.fna NC_069718.1:30096928-30128711 > N5.HOG0034670.Pangasianodon_hypophthalmus.extended.4.fa

sed -i 's/:/-/g' N5.HOG0034670.Pangasianodon_hypophthalmus.extended.1.fa
makeblastdb -in N5.HOG0034670.Pangasianodon_hypophthalmus.extended.1.fa -dbtype nucl
sed -i 's/:/-/g' N5.HOG0034670.Pangasianodon_hypophthalmus.extended.2.fa
makeblastdb -in N5.HOG0034670.Pangasianodon_hypophthalmus.extended.2.fa -dbtype nucl
sed -i 's/:/-/g' N5.HOG0034670.Pangasianodon_hypophthalmus.extended.3.fa
makeblastdb -in N5.HOG0034670.Pangasianodon_hypophthalmus.extended.3.fa -dbtype nucl
sed -i 's/:/-/g' N5.HOG0034670.Pangasianodon_hypophthalmus.extended.4.fa
makeblastdb -in N5.HOG0034670.Pangasianodon_hypophthalmus.extended.4.fa -dbtype nucl


#Now lets find TEs using the repbase and dfam databases

tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0034670.Carassius_auratus.extended.1.fa -evalue 1e-5 -outfmt 6 -out TE.Carassius_auratus.1.tblastn -num_threads 8
sed -i 's/#//g' TE.Carassius_auratus.1.tblastn
tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0034670.Carassius_auratus.extended.2.fa -evalue 1e-5 -outfmt 6 -out TE.Carassius_auratus.2.tblastn -num_threads 8
sed -i 's/#//g' TE.Carassius_auratus.2.tblastn
tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0034670.Carassius_auratus.extended.3.fa -evalue 1e-5 -outfmt 6 -out TE.Carassius_auratus.3.tblastn -num_threads 8
sed -i 's/#//g' TE.Carassius_auratus.3.tblastn

tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0034670.Pangasianodon_hypophthalmus.extended.1.fa -evalue 1e-5 -outfmt 6 -out TE.Pangasianodon_hypophthalmus.1.tblastn -num_threads 8
sed -i 's/#//g' TE.Pangasianodon_hypophthalmus.1.tblastn
tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0034670.Pangasianodon_hypophthalmus.extended.2.fa -evalue 1e-5 -outfmt 6 -out TE.Pangasianodon_hypophthalmus.2.tblastn -num_threads 8
sed -i 's/#//g' TE.Pangasianodon_hypophthalmus.2.tblastn
tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0034670.Pangasianodon_hypophthalmus.extended.3.fa -evalue 1e-5 -outfmt 6 -out TE.Pangasianodon_hypophthalmus.3.tblastn -num_threads 8
sed -i 's/#//g' TE.Pangasianodon_hypophthalmus.3.tblastn
tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0034670.Pangasianodon_hypophthalmus.extended.4.fa -evalue 1e-5 -outfmt 6 -out TE.Pangasianodon_hypophthalmus.4.tblastn -num_threads 8
sed -i 's/#//g' TE.Pangasianodon_hypophthalmus.4.tblastn


#merge tblastn hits and find the best TE match by doing a blastx

Rscript Rscript_merge_blast_hits.R TE.Carassius_auratus.1.tblastn TE.Carassius_auratus.1.tblastn.merged
Rscript Rscript_merge_blast_hits.R TE.Carassius_auratus.2.tblastn TE.Carassius_auratus.2.tblastn.merged
Rscript Rscript_merge_blast_hits.R TE.Carassius_auratus.3.tblastn TE.Carassius_auratus.3.tblastn.merged

Rscript Rscript_merge_blast_hits.R TE.Pangasianodon_hypophthalmus.1.tblastn TE.Pangasianodon_hypophthalmus.1.tblastn.merged
Rscript Rscript_merge_blast_hits.R TE.Pangasianodon_hypophthalmus.2.tblastn TE.Pangasianodon_hypophthalmus.2.tblastn.merged
Rscript Rscript_merge_blast_hits.R TE.Pangasianodon_hypophthalmus.3.tblastn TE.Pangasianodon_hypophthalmus.3.tblastn.merged
Rscript Rscript_merge_blast_hits.R TE.Pangasianodon_hypophthalmus.4.tblastn TE.Pangasianodon_hypophthalmus.4.tblastn.merged


xargs samtools faidx N5.HOG0034670.Carassius_auratus.extended.1.fa < TE.Carassius_auratus.1.tblastn.merged > TE.Carassius_auratus.1.BEST.fa
blastx -query TE.Carassius_auratus.1.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Carassius_auratus.1.BEST.blastx -max_target_seqs 1
cut -f1 TE.Carassius_auratus.1.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Carassius_auratus.1.BEST.blastx >> temp  ; done ; mv temp TE.Carassius_auratus.1.BEST.blastx

xargs samtools faidx N5.HOG0034670.Carassius_auratus.extended.2.fa < TE.Carassius_auratus.2.tblastn.merged > TE.Carassius_auratus.2.BEST.fa
blastx -query TE.Carassius_auratus.2.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Carassius_auratus.2.BEST.blastx -max_target_seqs 1
cut -f1 TE.Carassius_auratus.2.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Carassius_auratus.2.BEST.blastx >> temp  ; done ; mv temp TE.Carassius_auratus.2.BEST.blastx

xargs samtools faidx N5.HOG0034670.Carassius_auratus.extended.3.fa < TE.Carassius_auratus.3.tblastn.merged > TE.Carassius_auratus.3.BEST.fa
blastx -query TE.Carassius_auratus.3.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Carassius_auratus.3.BEST.blastx -max_target_seqs 1
cut -f1 TE.Carassius_auratus.3.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Carassius_auratus.3.BEST.blastx >> temp  ; done ; mv temp TE.Carassius_auratus.3.BEST.blastx


xargs samtools faidx N5.HOG0034670.Pangasianodon_hypophthalmus.extended.1.fa < TE.Pangasianodon_hypophthalmus.1.tblastn.merged > TE.Pangasianodon_hypophthalmus.1.BEST.fa
blastx -query TE.Pangasianodon_hypophthalmus.1.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Pangasianodon_hypophthalmus.1.BEST.blastx -max_target_seqs 1
cut -f1  TE.Pangasianodon_hypophthalmus.1.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Pangasianodon_hypophthalmus.1.BEST.blastx >> temp  ; done ; mv temp TE.Pangasianodon_hypophthalmus.1.BEST.blastx

xargs samtools faidx N5.HOG0034670.Pangasianodon_hypophthalmus.extended.2.fa < TE.Pangasianodon_hypophthalmus.2.tblastn.merged > TE.Pangasianodon_hypophthalmus.2.BEST.fa
blastx -query TE.Pangasianodon_hypophthalmus.2.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Pangasianodon_hypophthalmus.2.BEST.blastx -max_target_seqs 1
cut -f1  TE.Pangasianodon_hypophthalmus.2.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Pangasianodon_hypophthalmus.2.BEST.blastx >> temp  ; done ; mv temp TE.Pangasianodon_hypophthalmus.2.BEST.blastx

xargs samtools faidx N5.HOG0034670.Pangasianodon_hypophthalmus.extended.3.fa < TE.Pangasianodon_hypophthalmus.3.tblastn.merged > TE.Pangasianodon_hypophthalmus.3.BEST.fa
blastx -query TE.Pangasianodon_hypophthalmus.3.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Pangasianodon_hypophthalmus.3.BEST.blastx -max_target_seqs 1
cut -f1  TE.Pangasianodon_hypophthalmus.3.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Pangasianodon_hypophthalmus.3.BEST.blastx >> temp  ; done ; mv temp TE.Pangasianodon_hypophthalmus.3.BEST.blastx

xargs samtools faidx N5.HOG0034670.Pangasianodon_hypophthalmus.extended.4.fa < TE.Pangasianodon_hypophthalmus.4.tblastn.merged > TE.Pangasianodon_hypophthalmus.4.BEST.fa
blastx -query TE.Pangasianodon_hypophthalmus.4.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Pangasianodon_hypophthalmus.4.BEST.blastx -max_target_seqs 1
cut -f1  TE.Pangasianodon_hypophthalmus.4.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Pangasianodon_hypophthalmus.4.BEST.blastx >> temp  ; done ; mv temp TE.Pangasianodon_hypophthalmus.4.BEST.blastx

#Now find shared elements

cut -f2 TE.Carassius_auratus.1.BEST.blastx | sort | uniq > TE.Carassius_auratus.1.uniqTE
cut -f2 TE.Carassius_auratus.2.BEST.blastx | sort | uniq > TE.Carassius_auratus.2.uniqTE
cut -f2 TE.Carassius_auratus.3.BEST.blastx | sort | uniq > TE.Carassius_auratus.3.uniqTE

cut -f2 TE.Pangasianodon_hypophthalmus.1.BEST.blastx | sort | uniq > TE.Pangasianodon_hypophthalmus.1.uniqTE
cut -f2 TE.Pangasianodon_hypophthalmus.2.BEST.blastx | sort | uniq > TE.Pangasianodon_hypophthalmus.2.uniqTE
cut -f2 TE.Pangasianodon_hypophthalmus.3.BEST.blastx | sort | uniq > TE.Pangasianodon_hypophthalmus.3.uniqTE
cut -f2 TE.Pangasianodon_hypophthalmus.4.BEST.blastx | sort | uniq > TE.Pangasianodon_hypophthalmus.4.uniqTE


comm -12 TE.Carassius_auratus.1.uniqTE TE.Pangasianodon_hypophthalmus.1.uniqTE
comm -12 TE.Carassius_auratus.2.uniqTE TE.Pangasianodon_hypophthalmus.1.uniqTE
comm -12 TE.Carassius_auratus.3.uniqTE TE.Pangasianodon_hypophthalmus.1.uniqTE

comm -12 TE.Carassius_auratus.1.uniqTE TE.Pangasianodon_hypophthalmus.2.uniqTE
comm -12 TE.Carassius_auratus.2.uniqTE TE.Pangasianodon_hypophthalmus.2.uniqTE
comm -12 TE.Carassius_auratus.3.uniqTE TE.Pangasianodon_hypophthalmus.2.uniqTE

comm -12 TE.Carassius_auratus.1.uniqTE TE.Pangasianodon_hypophthalmus.3.uniqTE
comm -12 TE.Carassius_auratus.2.uniqTE TE.Pangasianodon_hypophthalmus.3.uniqTE
comm -12 TE.Carassius_auratus.3.uniqTE TE.Pangasianodon_hypophthalmus.3.uniqTE

comm -12 TE.Carassius_auratus.1.uniqTE TE.Pangasianodon_hypophthalmus.4.uniqTE
comm -12 TE.Carassius_auratus.2.uniqTE TE.Pangasianodon_hypophthalmus.4.uniqTE
comm -12 TE.Carassius_auratus.3.uniqTE TE.Pangasianodon_hypophthalmus.4.uniqTE

## Two shared TEs => L2-5_DRep:ClassI:LINE:Jockey:L2 and REX1-7_XT_pol#LINE/Rex-Babar

samtools faidx Dfam_plus_Repbase.cdhit80.prot L2-5_DRep:ClassI:LINE:Jockey:L2 > L2-5.prot
samtools faidx Dfam_plus_Repbase.cdhit80.prot REX1-7_XT_pol#LINE/Rex-Babar > REX1-7.prot



###### FIRST L2-5

#Silurus_meridionalis
tblastn -query L2-5.prot -db GCF_014805685.1_ASM1480568v1_genomic.fna -evalue 1e-1 -outfmt 6 -out L2-5.Silurus_meridionalis.tblastn -num_threads 8
sed -i 's/#//g' L2-5.Silurus_meridionalis.tblastn
Rscript Rscript_merge_blast_hits.R L2-5.Silurus_meridionalis.tblastn L2-5.Silurus_meridionalis.tblastn.merged
xargs samtools faidx GCF_014805685.1_ASM1480568v1_genomic.fna < L2-5.Silurus_meridionalis.tblastn.merged > L2-5.Silurus_meridionalis.BEST.fa
diamond blastx --query L2-5.Silurus_meridionalis.BEST.fa --db ../../BetweenActino_HGT/Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out L2-5.Silurus_meridionalis.BEST.blastx --max-target-seqs 1
grep "L2-5" L2-5.Silurus_meridionalis.BEST.blastx > L2-5.Silurus_meridionalis.list

#Ictalurus_punctatus
tblastn -query L2-5.prot -db GCF_001660625.3_Coco_2.0_genomic.fna -evalue 1e-1 -outfmt 6 -out L2-5.Ictalurus_punctatus.tblastn -num_threads 8
sed -i 's/#//g' L2-5.Ictalurus_punctatus.tblastn
Rscript Rscript_merge_blast_hits.R L2-5.Ictalurus_punctatus.tblastn L2-5.Ictalurus_punctatus.tblastn.merged
xargs samtools faidx GCF_001660625.3_Coco_2.0_genomic.fna  < L2-5.Ictalurus_punctatus.tblastn.merged > L2-5.Ictalurus_punctatus.BEST.fa
diamond blastx --query L2-5.Ictalurus_punctatus.BEST.fa --db ../../BetweenActino_HGT/Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out L2-5.Ictalurus_punctatus.BEST.blastx --max-target-seqs 1
grep "L2-5" L2-5.Ictalurus_punctatus.BEST.blastx > L2-5.Ictalurus_punctatus.list


#Pangasius_djambal
tblastn -query L2-5.prot -db GCA_022985145.1_PDJA_genomic.fna -evalue 1e-1 -outfmt 6 -out L2-5.Pangasius_djambal.tblastn -num_threads 8
sed -i 's/#//g' L2-5.Pangasius_djambal.tblastn
Rscript Rscript_merge_blast_hits.R L2-5.Pangasius_djambal.tblastn L2-5.Pangasius_djambal.tblastn.merged
xargs samtools faidx GCA_022985145.1_PDJA_genomic.fna  < L2-5.Pangasius_djambal.tblastn.merged > L2-5.Pangasius_djambal.BEST.fa
diamond blastx --query L2-5.Pangasius_djambal.BEST.fa --db ../../BetweenActino_HGT/Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out L2-5.Pangasius_djambal.BEST.blastx --max-target-seqs 1
grep "L2-5" L2-5.Pangasius_djambal.BEST.blastx > L2-5.Pangasius_djambal.list

#Pangasianodon_hypophthalmus
tblastn -query L2-5.prot -db GCF_027358585.1_fPanHyp1.pri_genomic.fna -evalue 1e-1 -outfmt 6 -out L2-5.Pangasianodon_hypophthalmus.tblastn -num_threads 8
sed -i 's/#//g' L2-5.Pangasianodon_hypophthalmus.tblastn
Rscript Rscript_merge_blast_hits.R L2-5.Pangasianodon_hypophthalmus.tblastn L2-5.Pangasianodon_hypophthalmus.tblastn.merged
xargs samtools faidx GCF_027358585.1_fPanHyp1.pri_genomic.fna < L2-5.Pangasianodon_hypophthalmus.tblastn.merged > L2-5.Pangasianodon_hypophthalmus.BEST.fa
diamond blastx --query L2-5.Pangasianodon_hypophthalmus.BEST.fa --db ../../BetweenActino_HGT/Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out L2-5.Pangasianodon_hypophthalmus.BEST.blastx --max-target-seqs 1
grep "L2-5" L2-5.Pangasianodon_hypophthalmus.BEST.blastx > L2-5.Pangasianodon_hypophthalmus.list


#Carassius_auratus
tblastn -query L2-5.prot -db GCF_003368295.1_ASM336829v1_genomic.fna -evalue 1e-1 -outfmt 6 -out L2-5.Carassius_auratus.tblastn -num_threads 8
sed -i 's/#//g' L2-5.Carassius_auratus.tblastn
Rscript Rscript_merge_blast_hits.R L2-5.Carassius_auratus.tblastn L2-5.Carassius_auratus.tblastn.merged
xargs samtools faidx GCF_003368295.1_ASM336829v1_genomic.fna  < L2-5.Carassius_auratus.tblastn.merged > L2-5.Carassius_auratus.BEST.fa
diamond blastx --query L2-5.Carassius_auratus.BEST.fa --db ../../BetweenActino_HGT/Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out L2-5.Carassius_auratus.BEST.blastx --max-target-seqs 1
grep "L2-5" L2-5.Carassius_auratus.BEST.blastx > L2-5.Carassius_auratus.list


#Cyprinus_carpio
tblastn -query L2-5.prot -db GCF_018340385.1_ASM1834038v1_genomic.fna -evalue 1e-1 -outfmt 6 -out L2-5.Cyprinus_carpio.tblastn -num_threads 8
sed -i 's/#//g' L2-5.Cyprinus_carpio.tblastn
Rscript Rscript_merge_blast_hits.R L2-5.Cyprinus_carpio.tblastn L2-5.Cyprinus_carpio.tblastn.merged
xargs samtools faidx GCF_018340385.1_ASM1834038v1_genomic.fna < L2-5.Cyprinus_carpio.tblastn.merged > L2-5.Cyprinus_carpio.BEST.fa
diamond blastx --query L2-5.Cyprinus_carpio.BEST.fa --db ../../BetweenActino_HGT/Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out L2-5.Cyprinus_carpio.BEST.blastx --max-target-seqs 1
grep "L2-5" L2-5.Cyprinus_carpio.BEST.blastx > L2-5.Cyprinus_carpio.list

#Onychostoma_macrolepis
tblastn -query L2-5.prot -db GCF_012432095.1_ASM1243209v1_genomic.fna -evalue 1e-1 -outfmt 6 -out L2-5.Onychostoma_macrolepis.tblastn -num_threads 8
sed -i 's/#//g' L2-5.Onychostoma_macrolepis.tblastn
Rscript Rscript_merge_blast_hits.R L2-5.Onychostoma_macrolepis.tblastn L2-5.Onychostoma_macrolepis.tblastn.merged
xargs samtools faidx GCF_012432095.1_ASM1243209v1_genomic.fna  < L2-5.Onychostoma_macrolepis.tblastn.merged > L2-5.Onychostoma_macrolepis.BEST.fa
diamond blastx --query L2-5.Onychostoma_macrolepis.BEST.fa --db ../../BetweenActino_HGT/Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out L2-5.Onychostoma_macrolepis.BEST.blastx --max-target-seqs 1
grep "L2-5" L2-5.Onychostoma_macrolepis.BEST.blastx > L2-5.Onychostoma_macrolepis.list

#Danio_rerio
tblastn -query L2-5.prot -db GCF_000002035.6_GRCz11_genomic.fna.fasfas -evalue 1e-1 -outfmt 6 -out L2-5.Danio_rerio.tblastn -num_threads 8
sed -i 's/#//g' L2-5.Danio_rerio.tblastn
Rscript Rscript_merge_blast_hits.R L2-5.Danio_rerio.tblastn L2-5.Danio_rerio.tblastn.merged
xargs samtools faidx GCF_000002035.6_GRCz11_genomic.fna.fasfas < L2-5.Danio_rerio.tblastn.merged > L2-5.Danio_rerio.BEST.fa
diamond blastx --query L2-5.Danio_rerio.BEST.fa --db ../../BetweenActino_HGT/Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out L2-5.Danio_rerio.BEST.blastx --max-target-seqs 1
grep "L2-5" L2-5.Danio_rerio.BEST.blastx > L2-5.Danio_rerio.list


wc -l L2-5.Silurus_meridionalis.list
wc -l L2-5.Ictalurus_punctatus.list
wc -l L2-5.Pangasius_djambal.list
wc -l L2-5.Pangasianodon_hypophthalmus.list
wc -l L2-5.Carassius_auratus.list
wc -l L2-5.Cyprinus_carpio.list
wc -l L2-5.Onychostoma_macrolepis.list
wc -l L2-5.Danio_rerio.list




###### THEN REX1-7


tblastn -query REX1-7.prot -db GCF_014805685.1_ASM1480568v1_genomic.fna -evalue 1e-1 -outfmt 6 -out REX1-7.Silurus_meridionalis.tblastn -num_threads 8
sed -i 's/#//g' REX1-7.Silurus_meridionalis.tblastn
Rscript Rscript_merge_blast_hits.R REX1-7.Silurus_meridionalis.tblastn REX1-7.Silurus_meridionalis.tblastn.merged
xargs samtools faidx GCF_014805685.1_ASM1480568v1_genomic.fna < REX1-7.Silurus_meridionalis.tblastn.merged > REX1-7.Silurus_meridionalis.BEST.fa
diamond blastx --query REX1-7.Silurus_meridionalis.BEST.fa --db ../../BetweenActino_HGT/Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out REX1-7.Silurus_meridionalis.BEST.blastx --max-target-seqs 1
grep "REX1-7" REX1-7.Silurus_meridionalis.BEST.blastx > REX1-7.Silurus_meridionalis.list

#Ictalurus_punctatus
tblastn -query REX1-7.prot -db GCF_001660625.3_Coco_2.0_genomic.fna -evalue 1e-1 -outfmt 6 -out REX1-7.Ictalurus_punctatus.tblastn -num_threads 8
sed -i 's/#//g' REX1-7.Ictalurus_punctatus.tblastn
Rscript Rscript_merge_blast_hits.R REX1-7.Ictalurus_punctatus.tblastn REX1-7.Ictalurus_punctatus.tblastn.merged
xargs samtools faidx GCF_001660625.3_Coco_2.0_genomic.fna  < REX1-7.Ictalurus_punctatus.tblastn.merged > REX1-7.Ictalurus_punctatus.BEST.fa
diamond blastx --query REX1-7.Ictalurus_punctatus.BEST.fa --db ../../BetweenActino_HGT/Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out REX1-7.Ictalurus_punctatus.BEST.blastx --max-target-seqs 1
grep "REX1-7" REX1-7.Ictalurus_punctatus.BEST.blastx > REX1-7.Ictalurus_punctatus.list


#Pangasius_djambal
tblastn -query REX1-7.prot -db GCA_022985145.1_PDJA_genomic.fna -evalue 1e-1 -outfmt 6 -out REX1-7.Pangasius_djambal.tblastn -num_threads 8
sed -i 's/#//g' REX1-7.Pangasius_djambal.tblastn
Rscript Rscript_merge_blast_hits.R REX1-7.Pangasius_djambal.tblastn REX1-7.Pangasius_djambal.tblastn.merged
xargs samtools faidx GCA_022985145.1_PDJA_genomic.fna  < REX1-7.Pangasius_djambal.tblastn.merged > REX1-7.Pangasius_djambal.BEST.fa
diamond blastx --query REX1-7.Pangasius_djambal.BEST.fa --db ../../BetweenActino_HGT/Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out REX1-7.Pangasius_djambal.BEST.blastx --max-target-seqs 1
grep "REX1-7" REX1-7.Pangasius_djambal.BEST.blastx > REX1-7.Pangasius_djambal.list

#Pangasianodon_hypophthalmus
tblastn -query REX1-7.prot -db GCF_027358585.1_fPanHyp1.pri_genomic.fna -evalue 1e-1 -outfmt 6 -out REX1-7.Pangasianodon_hypophthalmus.tblastn -num_threads 8
sed -i 's/#//g' REX1-7.Pangasianodon_hypophthalmus.tblastn
Rscript Rscript_merge_blast_hits.R REX1-7.Pangasianodon_hypophthalmus.tblastn REX1-7.Pangasianodon_hypophthalmus.tblastn.merged
xargs samtools faidx GCF_027358585.1_fPanHyp1.pri_genomic.fna < REX1-7.Pangasianodon_hypophthalmus.tblastn.merged > REX1-7.Pangasianodon_hypophthalmus.BEST.fa
diamond blastx --query REX1-7.Pangasianodon_hypophthalmus.BEST.fa --db ../../BetweenActino_HGT/Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out REX1-7.Pangasianodon_hypophthalmus.BEST.blastx --max-target-seqs 1
grep "REX1-7" REX1-7.Pangasianodon_hypophthalmus.BEST.blastx > REX1-7.Pangasianodon_hypophthalmus.list


#Carassius_auratus
tblastn -query REX1-7.prot -db GCF_003368295.1_ASM336829v1_genomic.fna -evalue 1e-1 -outfmt 6 -out REX1-7.Carassius_auratus.tblastn -num_threads 8
sed -i 's/#//g' REX1-7.Carassius_auratus.tblastn
Rscript Rscript_merge_blast_hits.R REX1-7.Carassius_auratus.tblastn REX1-7.Carassius_auratus.tblastn.merged
xargs samtools faidx GCF_003368295.1_ASM336829v1_genomic.fna  < REX1-7.Carassius_auratus.tblastn.merged > REX1-7.Carassius_auratus.BEST.fa
diamond blastx --query REX1-7.Carassius_auratus.BEST.fa --db ../../BetweenActino_HGT/Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out REX1-7.Carassius_auratus.BEST.blastx --max-target-seqs 1
grep "REX1-7" REX1-7.Carassius_auratus.BEST.blastx > REX1-7.Carassius_auratus.list


#Cyprinus_carpio
tblastn -query REX1-7.prot -db GCF_018340385.1_ASM1834038v1_genomic.fna -evalue 1e-1 -outfmt 6 -out REX1-7.Cyprinus_carpio.tblastn -num_threads 8
sed -i 's/#//g' REX1-7.Cyprinus_carpio.tblastn
Rscript Rscript_merge_blast_hits.R REX1-7.Cyprinus_carpio.tblastn REX1-7.Cyprinus_carpio.tblastn.merged
xargs samtools faidx GCF_018340385.1_ASM1834038v1_genomic.fna < REX1-7.Cyprinus_carpio.tblastn.merged > REX1-7.Cyprinus_carpio.BEST.fa
diamond blastx --query REX1-7.Cyprinus_carpio.BEST.fa --db ../../BetweenActino_HGT/Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out REX1-7.Cyprinus_carpio.BEST.blastx --max-target-seqs 1
grep "REX1-7" REX1-7.Cyprinus_carpio.BEST.blastx > REX1-7.Cyprinus_carpio.list

#Onychostoma_macrolepis
tblastn -query REX1-7.prot -db GCF_012432095.1_ASM1243209v1_genomic.fna -evalue 1e-1 -outfmt 6 -out REX1-7.Onychostoma_macrolepis.tblastn -num_threads 8
sed -i 's/#//g' REX1-7.Onychostoma_macrolepis.tblastn
Rscript Rscript_merge_blast_hits.R REX1-7.Onychostoma_macrolepis.tblastn REX1-7.Onychostoma_macrolepis.tblastn.merged
xargs samtools faidx GCF_012432095.1_ASM1243209v1_genomic.fna  < REX1-7.Onychostoma_macrolepis.tblastn.merged > REX1-7.Onychostoma_macrolepis.BEST.fa
diamond blastx --query REX1-7.Onychostoma_macrolepis.BEST.fa --db ../../BetweenActino_HGT/Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out REX1-7.Onychostoma_macrolepis.BEST.blastx --max-target-seqs 1
grep "REX1-7" REX1-7.Onychostoma_macrolepis.BEST.blastx > REX1-7.Onychostoma_macrolepis.list

#Danio_rerio
tblastn -query REX1-7.prot -db GCF_000002035.6_GRCz11_genomic.fna.fasfas -evalue 1e-1 -outfmt 6 -out REX1-7.Danio_rerio.tblastn -num_threads 8
sed -i 's/#//g' REX1-7.Danio_rerio.tblastn
Rscript Rscript_merge_blast_hits.R REX1-7.Danio_rerio.tblastn REX1-7.Danio_rerio.tblastn.merged
xargs samtools faidx GCF_000002035.6_GRCz11_genomic.fna.fasfas < REX1-7.Danio_rerio.tblastn.merged > REX1-7.Danio_rerio.BEST.fa
diamond blastx --query REX1-7.Danio_rerio.BEST.fa --db ../../BetweenActino_HGT/Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out REX1-7.Danio_rerio.BEST.blastx --max-target-seqs 1
grep "REX1-7" REX1-7.Danio_rerio.BEST.blastx > REX1-7.Danio_rerio.list


wc -l REX1-7.Silurus_meridionalis.list
wc -l REX1-7.Ictalurus_punctatus.list
wc -l REX1-7.Pangasius_djambal.list
wc -l REX1-7.Pangasianodon_hypophthalmus.list
wc -l REX1-7.Carassius_auratus.list
wc -l REX1-7.Cyprinus_carpio.list
wc -l REX1-7.Onychostoma_macrolepis.list
wc -l REX1-7.Danio_rerio.list


### In a FINAL step we will make another synteny plot more zoomed


## Carassius_auratus #1

rm seq_clustered_infos_ogg.TE.txt ; rm clusters_ID_TE.txt


scaffold=`grep ">" N5.HOG0034670.Carassius_auratus.extended.1.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0034670.Carassius_auratus.extended.1.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0034670.Carassius_auratus.extended.1.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Carassius_auratus,$scaffold,$length" > clusters_ID_TE.txt

grep "XM_026243633" GFF3_files_per_species/Carassius_auratus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons


nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene1.$nbexon,N5_HOG0034670,-,Carassius_auratus" >> seq_clustered_infos_ogg.TE.txt
done

cut -f1 TE.Carassius_auratus.1.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Carassius_auratus.1.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0034670.Carassius_auratus.extended.1.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0034670.Carassius_auratus.extended.1.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold,/g" | sed "s/$/,$strand,Carassius_auratus/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done


## Carassius_auratus #2


scaffold=`grep ">" N5.HOG0034670.Carassius_auratus.extended.2.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0034670.Carassius_auratus.extended.2.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0034670.Carassius_auratus.extended.2.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Carassius_auratus,$scaffold.sec,$length" >> clusters_ID_TE.txt

grep "XM_026243545" GFF3_files_per_species/Carassius_auratus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons


nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold.sec,$real_exon_start,$real_exon_stop,gene2.$nbexon,N5_HOG0034670,-,Carassius_auratus" >> seq_clustered_infos_ogg.TE.txt
done

cut -f1 TE.Carassius_auratus.2.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Carassius_auratus.2.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0034670.Carassius_auratus.extended.2.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0034670.Carassius_auratus.extended.2.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold.sec,/g" | sed "s/$/,$strand,Carassius_auratus/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done


## Carassius_auratus #3


scaffold=`grep ">" N5.HOG0034670.Carassius_auratus.extended.3.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0034670.Carassius_auratus.extended.3.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0034670.Carassius_auratus.extended.3.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Carassius_auratus,$scaffold,$length" >> clusters_ID_TE.txt

grep "XM_026265006" GFF3_files_per_species/Carassius_auratus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons


nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene3.$nbexon,N5_HOG0034670,+,Carassius_auratus" >> seq_clustered_infos_ogg.TE.txt
done

cut -f1 TE.Carassius_auratus.3.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Carassius_auratus.3.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0034670.Carassius_auratus.extended.3.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0034670.Carassius_auratus.extended.3.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold,/g" | sed "s/$/,$strand,Carassius_auratus/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done




## Pangasianodon_hypophthalmus #1

scaffold=`grep ">" N5.HOG0034670.Pangasianodon_hypophthalmus.extended.1.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0034670.Pangasianodon_hypophthalmus.extended.1.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0034670.Pangasianodon_hypophthalmus.extended.1.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Pangasianodon_hypophthalmus,$scaffold,$length" >> clusters_ID_TE.txt

grep -A3 -B3 "N5_HOG0034670"  GFF3_N5_OGGs/Pangasianodon_hypophthalmus.gff.simplified.sorted.OGG.tiret

grep "XM_053237209" GFF3_files_per_species/Pangasianodon_hypophthalmus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons


nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene1.$nbexon,N5_HOG0034670,-,Pangasianodon_hypophthalmus" >> seq_clustered_infos_ogg.TE.txt
done

cut -f1 TE.Pangasianodon_hypophthalmus.1.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Pangasianodon_hypophthalmus.1.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0034670.Pangasianodon_hypophthalmus.extended.1.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0034670.Pangasianodon_hypophthalmus.extended.1.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold,/g" | sed "s/$/,$strand,Pangasianodon_hypophthalmus/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done




## Pangasianodon_hypophthalmus #2

scaffold=`grep ">" N5.HOG0034670.Pangasianodon_hypophthalmus.extended.2.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0034670.Pangasianodon_hypophthalmus.extended.2.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0034670.Pangasianodon_hypophthalmus.extended.2.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Pangasianodon_hypophthalmus,$scaffold.sec,$length" >> clusters_ID_TE.txt

grep -A3 -B3 "N5_HOG0034670"  GFF3_N5_OGGs/Pangasianodon_hypophthalmus.gff.simplified.sorted.OGG.tiret

grep "XM_053237136" GFF3_files_per_species/Pangasianodon_hypophthalmus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons


nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold.sec,$real_exon_start,$real_exon_stop,gene2.$nbexon,N5_HOG0034670,+,Pangasianodon_hypophthalmus" >> seq_clustered_infos_ogg.TE.txt
done

cut -f1 TE.Pangasianodon_hypophthalmus.2.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Pangasianodon_hypophthalmus.2.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0034670.Pangasianodon_hypophthalmus.extended.2.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0034670.Pangasianodon_hypophthalmus.extended.2.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold.sec,/g" | sed "s/$/,$strand,Pangasianodon_hypophthalmus/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done



## Pangasianodon_hypophthalmus #3

scaffold=`grep ">" N5.HOG0034670.Pangasianodon_hypophthalmus.extended.3.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0034670.Pangasianodon_hypophthalmus.extended.3.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0034670.Pangasianodon_hypophthalmus.extended.3.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Pangasianodon_hypophthalmus,$scaffold.third,$length" >> clusters_ID_TE.txt

grep -A3 -B3 "N5_HOG0034670"  GFF3_N5_OGGs/Pangasianodon_hypophthalmus.gff.simplified.sorted.OGG.tiret

grep "XM_053236680" GFF3_files_per_species/Pangasianodon_hypophthalmus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons


nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold.third,$real_exon_start,$real_exon_stop,gene3.$nbexon,N5_HOG0034670,+,Pangasianodon_hypophthalmus" >> seq_clustered_infos_ogg.TE.txt
done

cut -f1 TE.Pangasianodon_hypophthalmus.3.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Pangasianodon_hypophthalmus.3.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0034670.Pangasianodon_hypophthalmus.extended.3.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0034670.Pangasianodon_hypophthalmus.extended.3.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold.third,/g" | sed "s/$/,$strand,Pangasianodon_hypophthalmus/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done

## Pangasianodon_hypophthalmus #4

scaffold=`grep ">" N5.HOG0034670.Pangasianodon_hypophthalmus.extended.4.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0034670.Pangasianodon_hypophthalmus.extended.4.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0034670.Pangasianodon_hypophthalmus.extended.4.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Pangasianodon_hypophthalmus,$scaffold.fourth,$length" >> clusters_ID_TE.txt

grep -A3 -B3 "N5_HOG0034670"  GFF3_N5_OGGs/Pangasianodon_hypophthalmus.gff.simplified.sorted.OGG.tiret

grep "XM_053236679" GFF3_files_per_species/Pangasianodon_hypophthalmus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons


nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold.fourth,$real_exon_start,$real_exon_stop,gene4.$nbexon,N5_HOG0034670,+,Pangasianodon_hypophthalmus" >> seq_clustered_infos_ogg.TE.txt
done

cut -f1 TE.Pangasianodon_hypophthalmus.4.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Pangasianodon_hypophthalmus.4.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0034670.Pangasianodon_hypophthalmus.extended.4.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0034670.Pangasianodon_hypophthalmus.extended.4.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold.fourth,/g" | sed "s/$/,$strand,Pangasianodon_hypophthalmus/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done



### Look at in more details the observed conserved region


samtools faidx N5.HOG0034670.Pangasianodon_hypophthalmus.extended.1.fa NC_069718.1-29758409-29791397:5000-5294 > region1.Pangasianodon_hypophthalmus.fa
samtools faidx N5.HOG0034670.Pangasianodon_hypophthalmus.extended.1.fa NC_069718.1-29758409-29791397:25558-25858 > region2.Pangasianodon_hypophthalmus.fa

samtools faidx N5.HOG0034670.Carassius_auratus.extended.1.fa NW_020523543.1-1402666-1433535:5000-5382 > region1.Carassius_auratus.fa
samtools faidx N5.HOG0034670.Carassius_auratus.extended.1.fa NW_020523543.1-1402666-1433535:24907-25207 > region2.Carassius_auratus.fa


sed -i 's/>.*/>P.hypophthalmus/g' region1.Pangasianodon_hypophthalmus.fa
sed -i 's/>.*/>P.hypophthalmus/g' region2.Pangasianodon_hypophthalmus.fa

sed -i 's/>.*/>C.auratus/g' region1.Carassius_auratus.fa
sed -i 's/>.*/>C.auratus/g' region2.Carassius_auratus.fa

needle -asequence region1.Pangasianodon_hypophthalmus.fa -bsequence region1.Carassius_auratus.fa -outfile region1.aln -gapopen 10.0 -gapextend 0.5
needle -asequence region2.Pangasianodon_hypophthalmus.fa -bsequence region2.Carassius_auratus.fa -outfile region2.aln -gapopen 10.0 -gapextend 0.5


sbatch --qos=6hours -c 10 --mem=20G launch_blast_allgenomes.sh region1.Pangasianodon_hypophthalmus.fa 1e-5 region1.blastn.tsv
sbatch --qos=6hours -c 10 --mem=20G launch_blast_allgenomes.sh region2.Pangasianodon_hypophthalmus.fa 1e-5 region2.blastn.tsv



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
sbatch --qos=1day -c 4 --mem=10G -e error.fitMG4.out -o slurm.fitMG4.out launch_fitMG4.sh N5.HOG0034670




Recipient branches : 
Pangasianodon_gigas_PGIGA_T00260940
Pangasianodon_gigas_PGIGA_T00102290
Node73
Pangasianodon_hypophthalmus_rna_XM_053237209_1
Pangasianodon_hypophthalmus_rna_XM_053236680_1
Node80
Pangasianodon_hypophthalmus_rna_XM_053236679_1
Node79
Pangasianodon_hypophthalmus_rna_XM_053237136_1
Node78
Pangasianodon_gigas_PGIGA_T00260840
Pangasius_djambal_PDJAM_T00258610
Node85
Node77
Pangasianodon_gigas_PGIGA_T00260930
Node76
Node72

#Test positive selection and relaxed selection on receiver branch

sbatch --qos=1week -c 4 --mem=10G -e error.absrel.cand.out -o slurm.absrel.cand.out --job-name=HOG0034670 launch_absrel_cand.sh N5.HOG0034670


### Adaptive branch site random effects likelihood test 
#Likelihood ratio test for episodic diversifying positive selection at Holm-Bonferroni corrected _p =   0.0500_ found **2** branches under selection among **15** tested.
#* Pangasianodon_gigas_PGIGA_T00102290, p-value =  0.00000
#* Pangasianodon_gigas_PGIGA_T00260930, p-value =  0.00028

#Now , lets launch RELAX


===> $HOG.prot.aln.treefile.SelecMarked.relax => select branches detected with accelerated evolution
sbatch --qos=1week -c 6 --mem=10G -e error.relax.out -o slurm.relax.out --job-name=HOG0034670 launch_RELAX.sh N5.HOG0034670


