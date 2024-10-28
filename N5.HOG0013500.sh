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


samtools faidx Proteomes_BUSCO80/Thalassophryne_amazonica.fa Thalassophryne_amazonica---rna-XM_034166331.1 > Thalassophryne_amazonica---rna-XM_034166331.1.prot

#blastp against the fish proteome (50 best hits) + Uniprot database without any fish species (10 best hits). Here I will only take 20 catffishes otherwise there are way too much of them

blastp -query Thalassophryne_amazonica---rna-XM_034166331.1.prot -db concatenated_proteomes.fa -outfmt 6 -out Gene_vs_FishProteome.blastp -num_threads 8 -max_target_seqs 51
blastp -query Thalassophryne_amazonica---rna-XM_034166331.1.prot -db non_actino_uniprot.fa -outfmt 6 -out Gene_vs_Uniprot.blastp -num_threads 8  -max_target_seqs 10 -evalue 1e-5


cut -f2 Gene_vs_FishProteome.blastp | sort | uniq  > closest_fish_seq.id
cut -f2 Gene_vs_Uniprot.blastp | sort | uniq > closest_nonfish_seq.id 

xargs samtools faidx concatenated_proteomes.fa < closest_fish_seq.id > closest_fish_seq.fa
xargs samtools faidx non_actino_uniprot.fa  < closest_nonfish_seq.id > closest_nonfish_seq.fa

#Align with muscle and trim with trimal

cat closest_fish_seq.fa closest_nonfish_seq.fa > Uniprot_plus_closefish.fa

muscle5.1 -align closest_fish_seq.fa -output closest_fish_seq.aln
cp closest_fish_seq.aln Uniprot_plus_closefish.aln
trimal -in Uniprot_plus_closefish.aln -automated1 -out Uniprot_plus_closefish.aln.trimmed

#Make a ML tree with IQ-TREE
iqtree -s Uniprot_plus_closefish.aln -st AA -nt 8 -m TEST -mrate G4 --redo
iqtree -s Uniprot_plus_closefish.aln.trimmed -st AA -nt 8 -m TEST -mrate G4 --redo



##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##========================================== Re-annotation of the T. amazonica truncated gene  =======================================
##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================


grep "XM_034167571" GFF3_N5_OGGs/Thalassophryne_amazonica.gff.simplified.sorted.OGG.tiret

samtools faidx GCF_902500255.1_fThaAma1.1_genomic.fna NC_047105.1:54418817-54442016 > Thalassophryne_amazonica.region.XM_034167571.fa

exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs"  Thalassophryne_amazonica---rna-XM_034166331.1.prot  Thalassophryne_amazonica.region.XM_034167571.fa > Thalassophryne_amazonica.region.XM_034167571.exo

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


#N5.HOG0013500.clade_hgt.txt=> clade where the HGT is extracted from iTOL (from N5.HOG0013500.prot.aln.treefile)

cut -f1,2 -d "_" N5.HOG0013500.clade_hgt.txt | sort | uniq  > species_to_draw.clade1
#Choose close species without the gene .. 


echo "Engraulis_encrasicolus" > non_transfer_species_1
echo "Alosa_alosa" >> non_transfer_species_1
echo "Sardina_pilchardus" >> non_transfer_species_1

echo "Sphaeramia_orbicularis" > non_transfer_species_2
echo "Thunnus_maccoyii" >> non_transfer_species_2
echo "Lucifuga_dentata" >> non_transfer_species_2



rm -r Syn_tables_dir ; mkdir Syn_tables_dir


#First order the species names 

rm species_to_draw.clade1.ordered

echo "Clupea_harengus" >> species_to_draw.clade1.ordered
echo "Thalassophryne_amazonica" >> species_to_draw.clade1.ordered



#If there are more than 4 species from the same order, keep only 4 


#First extract genes around 

curr_OGG=N5.HOG0013500




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


	grep -A10 -B10  "N5_HOG0046052" GFF3_N5_OGGs/Clupea_harengus.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/Clupea_harengus.synt.df

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


grep "NC_045157_1\|NW_024880005_1" Syn_tables_dir/Clupea_harengus.synt.final.df  > temp ; mv temp Syn_tables_dir/Clupea_harengus.synt.final.df
sort Syn_tables_dir/Clupea_harengus.synt.final.df | uniq > temp ; mv temp Syn_tables_dir/Clupea_harengus.synt.final.df

#First add Sardina_pilchardus

curr_OGG=N5.HOG0013500
curr_sp=Sardina_pilchardus
ref_sp=Clupea_harengus


grep -A10 -B10 "N5_HOG0013500" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "NC_085004"  > Syn_tables_dir/$curr_sp.synt.df
grep -A25 -B25 "N5_HOG0022738" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "NC_085004" >> Syn_tables_dir/$curr_sp.synt.df

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

curr_OGG=N5.HOG0013500
curr_sp=Alosa_alosa
ref_sp=Clupea_harengus


grep -A25 -B25 "N5_HOG0022738" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "NC_063202" > Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10  "N5_HOG0046052" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "NC_063202" >> Syn_tables_dir/$curr_sp.synt.df




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

curr_OGG=N5.HOG0013500
curr_sp=Engraulis_encrasicolus
ref_sp=Clupea_harengus

grep -A25 -B25 "N5_HOG0022738" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "NC_085878" > Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10  "N5_HOG0046052" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df


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


#Now add Sphaeramia_orbicularis

curr_OGG=N5.HOG0013500
curr_sp=Sphaeramia_orbicularis
ref_sp=Thalassophryne_amazonica

grep -A10 -B10 "N5_HOG0040950" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret  > Syn_tables_dir/$curr_sp.synt.df


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



#Now add Thunnus_maccoyii

curr_OGG=N5.HOG0013500
curr_sp=Thunnus_maccoyii
ref_sp=Thalassophryne_amazonica

grep -A10 -B10 "N5_HOG0040950" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret  > Syn_tables_dir/$curr_sp.synt.df


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




#Now add Lucifuga_dentata

curr_OGG=N5.HOG0013500
curr_sp=Lucifuga_dentata
ref_sp=Thalassophryne_amazonica

grep -A10 -B10 "N5_HOG0040950" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret  > Syn_tables_dir/$curr_sp.synt.df


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


exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs" Thalassophryne_amazonica---rna-XM_034166331.1.prot GCF_902148855.1_fSphaOr1.1_genomic.fna > Sphaeramia_orbicularis.whole.fa
exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs" Thalassophryne_amazonica---rna-XM_034166331.1.prot GCF_910596095.1_fThuMac1.1_genomic.fna > Thunnus_maccoyii.whole.fa
exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs" Thalassophryne_amazonica---rna-XM_034166331.1.prot GCA_014773175.1_Ldentata1.0_genomic.fna > Lucifuga_dentata.whole.fa


best_score=`grep "Raw score:" Sphaeramia_orbicularis.whole.fa | sed 's/.*: //g' | sort -n | tail -1`
grep -B5 -A100 "Raw score: $best_score" Sphaeramia_orbicularis.whole.fa

best_score=`grep "Raw score:" Thunnus_maccoyii.whole.fa | sed 's/.*: //g' | sort -n | tail -1`
grep -B5 -A100 "Raw score: $best_score" Thunnus_maccoyii.whole.fa

best_score=`grep "Raw score:" Lucifuga_dentata.whole.fa | sed 's/.*: //g' | sort -n | tail -1`
grep -B5 -A100 "Raw score: $best_score" Lucifuga_dentata.whole.fa

==> no gene we look for





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

#Thalassophryne_amazonica
grep -A3 -B3 "N5_HOG0013500"  GFF3_N5_OGGs/Thalassophryne_amazonica.gff.simplified.sorted.OGG.tiret
samtools faidx GCF_902500255.1_fThaAma1.1_genomic.fna NC_047105.1:54401744-54437016 > N5.HOG0013500.Thalassophryne_amazonica.extended.1.fa
sed -i 's/:/-/g' N5.HOG0013500.Thalassophryne_amazonica.extended.1.fa
makeblastdb -in N5.HOG0013500.Thalassophryne_amazonica.extended.1.fa -dbtype nucl

#Clupea_harengus
grep -A3 -B3 "N5_HOG0013500"  GFF3_N5_OGGs/Clupea_harengus.gff.simplified.sorted.OGG.tiret
samtools faidx GCF_900700415.2_Ch_v2.0.2_genomic.fna NC_045157.1:8227306-8322129 > N5.HOG0013500.Clupea_harengus.extended.1.fa
samtools faidx GCF_900700415.2_Ch_v2.0.2_genomic.fna NC_045157.1:9199529-9212069 > N5.HOG0013500.Clupea_harengus.extended.2.fa
samtools faidx GCF_900700415.2_Ch_v2.0.2_genomic.fna NW_024880005.1:25144-42711 > N5.HOG0013500.Clupea_harengus.extended.3.fa

sed -i 's/:/-/g' N5.HOG0013500.Clupea_harengus.extended.1.fa
makeblastdb -in N5.HOG0013500.Clupea_harengus.extended.1.fa -dbtype nucl
sed -i 's/:/-/g' N5.HOG0013500.Clupea_harengus.extended.2.fa
makeblastdb -in N5.HOG0013500.Clupea_harengus.extended.2.fa -dbtype nucl
sed -i 's/:/-/g' N5.HOG0013500.Clupea_harengus.extended.3.fa
makeblastdb -in N5.HOG0013500.Clupea_harengus.extended.3.fa -dbtype nucl


#Now lets find TEs using the repbase and dfam databases

tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0013500.Thalassophryne_amazonica.extended.1.fa -evalue 1e-5 -outfmt 6 -out TE.Thalassophryne_amazonica.1.tblastn -num_threads 8
sed -i 's/#//g' TE.Thalassophryne_amazonica.1.tblastn

tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0013500.Clupea_harengus.extended.1.fa -evalue 1e-5 -outfmt 6 -out TE.Clupea_harengus.1.tblastn -num_threads 8
sed -i 's/#//g' TE.Clupea_harengus.1.tblastn

tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0013500.Clupea_harengus.extended.2.fa -evalue 1e-5 -outfmt 6 -out TE.Clupea_harengus.2.tblastn -num_threads 8
sed -i 's/#//g' TE.Clupea_harengus.2.tblastn

tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0013500.Clupea_harengus.extended.3.fa -evalue 1e-5 -outfmt 6 -out TE.Clupea_harengus.3.tblastn -num_threads 8
sed -i 's/#//g' TE.Clupea_harengus.3.tblastn

#merge tblastn hits and find the best TE match by doing a blastx

Rscript Rscript_merge_blast_hits.R TE.Thalassophryne_amazonica.1.tblastn TE.Thalassophryne_amazonica.1.tblastn.merged
Rscript Rscript_merge_blast_hits.R TE.Clupea_harengus.1.tblastn TE.Clupea_harengus.1.tblastn.merged
Rscript Rscript_merge_blast_hits.R TE.Clupea_harengus.2.tblastn TE.Clupea_harengus.2.tblastn.merged
Rscript Rscript_merge_blast_hits.R TE.Clupea_harengus.3.tblastn TE.Clupea_harengus.3.tblastn.merged


xargs samtools faidx N5.HOG0013500.Thalassophryne_amazonica.extended.1.fa < TE.Thalassophryne_amazonica.1.tblastn.merged > TE.Thalassophryne_amazonica.1.BEST.fa
blastx -query TE.Thalassophryne_amazonica.1.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Thalassophryne_amazonica.1.BEST.blastx -max_target_seqs 1
cut -f1 TE.Thalassophryne_amazonica.1.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Thalassophryne_amazonica.1.BEST.blastx >> temp  ; done ; mv temp TE.Thalassophryne_amazonica.1.BEST.blastx


xargs samtools faidx N5.HOG0013500.Clupea_harengus.extended.1.fa < TE.Clupea_harengus.1.tblastn.merged > TE.Clupea_harengus.1.BEST.fa
blastx -query TE.Clupea_harengus.1.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Clupea_harengus.1.BEST.blastx -max_target_seqs 1
cut -f1  TE.Clupea_harengus.1.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Clupea_harengus.1.BEST.blastx >> temp  ; done ; mv temp TE.Clupea_harengus.1.BEST.blastx


xargs samtools faidx  N5.HOG0013500.Clupea_harengus.extended.2.fa < TE.Clupea_harengus.2.tblastn.merged > TE.Clupea_harengus.2.BEST.fa
blastx -query TE.Clupea_harengus.2.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Clupea_harengus.2.BEST.blastx -max_target_seqs 1
cut -f1  TE.Clupea_harengus.2.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Clupea_harengus.2.BEST.blastx >> temp  ; done ; mv temp TE.Clupea_harengus.2.BEST.blastx

xargs samtools faidx  N5.HOG0013500.Clupea_harengus.extended.3.fa < TE.Clupea_harengus.3.tblastn.merged > TE.Clupea_harengus.3.BEST.fa
blastx -query TE.Clupea_harengus.3.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Clupea_harengus.3.BEST.blastx -max_target_seqs 1
cut -f1  TE.Clupea_harengus.3.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Clupea_harengus.3.BEST.blastx >> temp  ; done ; mv temp TE.Clupea_harengus.3.BEST.blastx


#Now find shared elements

cut -f2 TE.Thalassophryne_amazonica.1.BEST.blastx | sort | uniq > TE.Thalassophryne_amazonica.1.uniqTE
cut -f2 TE.Clupea_harengus.1.BEST.blastx | sort | uniq > TE.Clupea_harengus.1.uniqTE
cut -f2 TE.Clupea_harengus.2.BEST.blastx | sort | uniq > TE.Clupea_harengus.2.uniqTE
cut -f2 TE.Clupea_harengus.3.BEST.blastx | sort | uniq > TE.Clupea_harengus.3.uniqTE

comm -12 TE.Thalassophryne_amazonica.1.uniqTE TE.Clupea_harengus.1.uniqTE
comm -12 TE.Thalassophryne_amazonica.1.uniqTE TE.Clupea_harengus.2.uniqTE
comm -12 TE.Thalassophryne_amazonica.1.uniqTE TE.Clupea_harengus.3.uniqTE




#NO SHARED TEs


### In a FINAL step we will make another synteny plot focused on the two species


## Paramormyrops

rm seq_clustered_infos_ogg.TE.txt ; rm clusters_ID_TE.txt


scaffold=`grep ">" N5.HOG0013500.Thalassophryne_amazonica.extended.1.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0013500.Thalassophryne_amazonica.extended.1.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0013500.Thalassophryne_amazonica.extended.1.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Thalassophryne_amazonica,$scaffold,$length" > clusters_ID_TE.txt

grep "N5_HOG0013500"  GFF3_N5_OGGs/Thalassophryne_amazonica.gff.simplified.sorted.OGG.tiret

grep "XM_034166331" GFF3_files_per_species/Thalassophryne_amazonica.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons
grep "XM_034167571" GFF3_files_per_species/Thalassophryne_amazonica.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_2.exons


nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene1.$nbexon,N5_HOG0013500,-,Thalassophryne_amazonica" >> seq_clustered_infos_ogg.TE.txt
done
nbexon=0
for line in `cat gene_2.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene2.$nbexon,N5_HOG0013500,-,Thalassophryne_amazonica" >> seq_clustered_infos_ogg.TE.txt
done

cut -f1 TE.Thalassophryne_amazonica.1.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Thalassophryne_amazonica.1.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0013500.Thalassophryne_amazonica.extended.1.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0013500.Thalassophryne_amazonica.extended.1.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold,/g" | sed "s/$/,$strand,Thalassophryne_amazonica/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done




## Clupea_harengus ## SCAFF 1


scaffold=`grep ">" N5.HOG0013500.Clupea_harengus.extended.1.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0013500.Clupea_harengus.extended.1.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0013500.Clupea_harengus.extended.1.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Clupea_harengus,$scaffold,$length" >> clusters_ID_TE.txt

grep -A3 -B3 "N5_HOG0013500"  GFF3_N5_OGGs/Clupea_harengus.gff.simplified.sorted.OGG.tiret


grep "XM_042708011" GFF3_files_per_species/Clupea_harengus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons
grep "XM_031568919" GFF3_files_per_species/Clupea_harengus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_2.exons
grep "XM_031568920" GFF3_files_per_species/Clupea_harengus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_3.exons

nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene1.$nbexon,N5_HOG0013500,+,Clupea_harengus" >> seq_clustered_infos_ogg.TE.txt
done
nbexon=0
for line in `cat gene_2.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene2.$nbexon,N5_HOG0013500,+,Clupea_harengus" >> seq_clustered_infos_ogg.TE.txt
done
nbexon=0
for line in `cat gene_3.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene3.$nbexon,N5_HOG0013500,+,Clupea_harengus" >> seq_clustered_infos_ogg.TE.txt
done


cut -f1 TE.Clupea_harengus.1.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Clupea_harengus.1.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0013500.Clupea_harengus.extended.1.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0013500.Clupea_harengus.extended.1.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold,/g" | sed "s/$/,$strand,Clupea_harengus/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done



#Clupea_harengus ## SCAFF 2


scaffold=`grep ">" N5.HOG0013500.Clupea_harengus.extended.2.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0013500.Clupea_harengus.extended.2.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0013500.Clupea_harengus.extended.2.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Clupea_harengus,$scaffold.sec,$length" >> clusters_ID_TE.txt

grep "N5_HOG0013500"  GFF3_N5_OGGs/Clupea_harengus.gff.simplified.sorted.OGG.tiret

grep "XM_042707968" GFF3_files_per_species/Clupea_harengus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons

nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold.sec,$real_exon_start,$real_exon_stop,gene4.$nbexon,N5_HOG0013500,-,Clupea_harengus" >> seq_clustered_infos_ogg.TE.txt
done


cut -f1 TE.Clupea_harengus.2.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Clupea_harengus.2.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0013500.Clupea_harengus.extended.2.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0013500.Clupea_harengus.extended.2.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold.sec,/g" | sed "s/$/,$strand,Clupea_harengus/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done


#Clupea_harengus ## SCAFF 3


scaffold=`grep ">" N5.HOG0013500.Clupea_harengus.extended.3.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0013500.Clupea_harengus.extended.3.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0013500.Clupea_harengus.extended.3.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Clupea_harengus,$scaffold,$length" >> clusters_ID_TE.txt

grep "N5_HOG0013500"  GFF3_N5_OGGs/Clupea_harengus.gff.simplified.sorted.OGG.tiret

grep "XM_042705999" GFF3_files_per_species/Clupea_harengus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons
grep "XM_042706001" GFF3_files_per_species/Clupea_harengus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_2.exons

nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene5.$nbexon,N5_HOG0013500,+,Clupea_harengus" >> seq_clustered_infos_ogg.TE.txt
done
nbexon=0
for line in `cat gene_2.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene6.$nbexon,N5_HOG0013500,+,Clupea_harengus" >> seq_clustered_infos_ogg.TE.txt
done

cut -f1 TE.Clupea_harengus.3.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Clupea_harengus.3.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0013500.Clupea_harengus.extended.3.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0013500.Clupea_harengus.extended.3.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold,/g" | sed "s/$/,$strand,Clupea_harengus/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done





### Check conserved non-coding regions


samtools faidx N5.HOG0013500.Thalassophryne_amazonica.extended.1.fa NC_047105.1-54401744-54437016:8600-9000 > region1.Thalassophryne_amazonica.fa
samtools faidx N5.HOG0013500.Thalassophryne_amazonica.extended.1.fa NC_047105.1-54401744-54437016:30000-30500 > region2.Thalassophryne_amazonica.fa

samtools faidx N5.HOG0013500.Clupea_harengus.extended.2.fa NC_045157.1-9199529-9212069:7200-7450 > region1.Clupea_harengus.fa

sed -i 's/>.*/>T.amazonica/g' region1.Thalassophryne_amazonica.fa
sed -i 's/>.*/>T.amazonica/g' region2.Thalassophryne_amazonica.fa

sed -i 's/>.*/>C.harengus/g' region1.Clupea_harengus.fa

needle -asequence region1.Thalassophryne_amazonica.fa -bsequence region1.Clupea_harengus.fa -outfile region1.aln -gapopen 10.0 -gapextend 0.5
needle -asequence region2.Thalassophryne_amazonica.fa -bsequence region1.Clupea_harengus.fa -outfile region2.aln -gapopen 10.0 -gapextend 0.5


blastn -query region1.Thalassophryne_amazonica.fa -db Concatenated_assemblies.fa -outfmt 6 -out regionA.blastn.tsv -num_threads 10 
blastn -query region1.Clupea_harengus.fa -dbConcatenated_assemblies.fa -outfmt 6 -out regionA.blastn.2.tsv -num_threads 10





##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##======================Search for the gene in other Batrachoididae non-annotated genomes  ============================================
##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================

exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" Thalassophryne_amazonica---rna-XM_034166331.1.prot GCA_900302635.1_ASM90030263v1_genomic.fna > Chatrabus_melanurus.whole.fa
exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" Thalassophryne_amazonica---rna-XM_034166331.1.prot GCA_900660325.1_Opsanus_beta_assembly_genomic.fna > Opsanus_beta.whole.fa


samtools faidx GCA_900660325.1_Opsanus_beta_assembly_genomic.fna CABFOU010245896.1:0-10000 > CABFOU010245896.fa
exonerate -E True --showtargetgff TRUE --model protein2genome --ryo "%tcs" Thalassophryne_amazonica---rna-XM_034166331.1.prot CABFOU010245896.fa > Opsanus_beta.focused.fa


#Batrachoididae_N5.HOG0013500.fa

>Chatrabus_melanurus.GCA_900302635.1.OMMV01052601.1-6252-5536
GTATTTATTGATGGCAGTGACTTGAAAAGTCATTCTGACCTCAATCCATACACCTTCGATGACAACCTCA
ACCTTGAGTGGCAGTTGTTTACAGGGTCTGTTCCCGATGGAGCAGTGGGCTTCTGGAACAGCTACGCAAA
TCGCTACGATTACATGTGCAGGCTATCTGATGGCTGTGAGAGTGGATTCTACTACAGTGGCTTTGGTAAA
TGCCCTTTCCAACATGGTCACCCTCGTTTGTACAGCACTCAAAGTTTTTTCATAAACAAAGATGAGTTTG
AAGTTCTGGAGTGGAAATGGGAATCCCATGGGTCTGTTCCTCAAAACTTGGTCAGAACATGTTTCGGATC
AACTCAGGATtatgttggaaaaaaactgtatGGTTTCAGATTGCTTTATTCGGGATTTTTATTCTATTTG
GCTTGGGTAGAAAAATATTCTAGTGGTGAAATCTTTGGGTAGAGTACATGGTCCAAAAGGTCCTACGAGG
CTTTGACCATAAATACAGACAAACGTAAGCAGGAAATCAAAGATGTGGAGTATTTCATTGAACAAGGTGA
GATCATAGAGGCTTCCCCATACGGCATTATTGAATACACTTCAAACAATAGAGAATGCAATCCTGTAACC
GAGAAAACTGCACGATGTACATCAAAAACTTGGCAAGTCTCTTACTCCATGACCCTTGCAATCAGCACCA
CTGTTACTGCT
>Opsanus_beta.GCA_900660325.CABFOU010245896.1-2795-3950
ACCTTTGCTGAAGAGCCTGCTCCATCAGCTCAGATTGCAGCCCACAAGGCACTGATTGATGACA
GCGCCTCGAGAAGTCATTCTCACGTCCCTCCATACACCTTTGATGACAACCCCAACCTTGAGTGGCAGTT
GTTTAATGGGTATCTTCCCAACGGAGCAGTGGGCATCTGGAACAGTTACACAAATCGCTACGACTATGTG
TGCAGAACATATGATGGTTGTGAGAGCGGATTCTACAACAGTGACTTTGGTACATGCTTATTCCCACGTT
ACCCTGTCTTGTGTAGCACTCAGAGTTTTTACATCCTCGTAAACAAAGAtgagtttgtatttttggagtG
GAAATCGGGATCCGGTGGGTCTGTTCCCCAAAACTCAGTCAGAACATGTTTTAAGTCAGATAAGGTTTAT
GTTGGAAAAAATGAGTATGGTTTAGGGTTGCTTTCCTCTGGAGATTCATTCTATTTGCCATGGGTTCAAG
AATATTCCGGTGGTACAATCTATGGGTACAGTACATGGTATAGACGCACTTACCAGGCTTTGACAGTAAA
TACAGACAAATATAAGCAGGAAATCAAAGACATAAAGTATTTCACTGATCAAGCTAAGATCATAGGCACT
CCCCCATTCGGTATTGTCAAAACCAGTTTAAACAACAGAGAATGCCATCCTGTAACTTTGACGACCACAC
TATCTACatcccaaacaaaaacaaataattggGAATTCTCTTACTCTATGAGCCTTTCCATCAGTACCAC
TGTATCTGCTGGGATTCCTGACATAGTTGACTTTAGTGTTACCATTGGTGTGGAGCAGACTTTTACAGTA
ACAAAGGAAATATCCCTCTCTCAAACTCAGACTACTAGCCTCGAGGTCCAAGTTACTGTGCCCCCAAACA
TGACGTGCACTATTGAGATGGAAGGCAGGAAATTCACATCAAATATACCTTTTCAGGCTCGCCTTGGCCG
TACTTACAGCAATGGTGACACAAAGTGGACCACCATAACCGGTATATATGATGGACTTCAGGTGTCTGAT
TTTCAGGCTACCACGAAGAAGTGTGAGCCTGTGGATGACCTTGTGACTTGTGGGAATGATGTCTATGTGT
CTGATATATCCTCCTCACATCAACAACCCCCTGTT
>Thalassophryne_amazonica_rna_XM_034167571_1-Manually_Corrected
ATGAGGTTCTACCTGTCTGTTGTTTTCACCATTCTGGTCATGGCTGCGGACCACACCCGTGATTGTGAAC
CACTTACCGAAACAGTGAGACTACCCTCAGTCAAGGTTTCCAAGGAAATCAAAGAAACTCAGATCTCGAC
CGTTGCTGAAGAGCCTGCTCCATCAgctcagattgcatttcacaagGCATTGATTGATGGCAGCACCTCG
AGAAGTGATTCTCACGTCCCTCCATACACCTTTGATGACAACCCCAACCTTGAGTGGCAGTTGTTTAATG
GGTCTCTTCCCAACGGAGCAGTGGGCATCTGGAACAGTTACACAAATCGCTATGACTATGTGTGCAGAGC
ATATGATGGTTGTGAGAGCGGATTCTACAGCAGTGACTTTGGTAAATGCTTGTTCCCAAGTTACCCTTTG
TTGTTGATTACTCAGAGTTTTTACATCCTCGTAAACAAAGACGACTTTGTATTTCTGGAGTGGAAATGGG
GATCAGGTGGGTCTGTTCCCCAAAACTCAATCAGGACATGTTTTAAGCCACAGTTTTATGTTGGAAAAAA
TAAGTATGGTTTAGGGTTGCTTTCCTCTGGAGATTCATTCTATTTGCCATGGGTTGAAGAATATTCCAAT
GGTACAATCTATGGGTATAGTACATGGTATAGAAGCTCTTACCAGGCTTTGACAAATACAGACAAATATA
AGCAGGAAATCAAAGATGTAAAGTATTTCACGGATCAAGCTAAGATCATAGGCACTCCCCCATTCGGTAT
TGTCAAAAACACTTTAAGCAACAGAGAATGCCATCCTGTAACTTTGACGACCACATTGTCTACATctcaa
acaaaaacaagcaattGGCAATTCTCTTACTCCATGACCCTTTCCATCAGCACCACTGTATCTGCTAGCA
TCCCTGACATCGTTGACTTTAGTGTCACCATTGGTGTAGAGCAGTTTACAGTAACAAAGGAAATATCCCT
GTCTCAAACTCAAACTACTAGCCTTGTGGTTGAAGTCACTGTGCCCCCAAACATGAGGTGCACTATTGAG
ATGGAAGGCAAGAAATTCACAAATATACCTTTTCAGGCTCGCCTTGGCCGTACCTACAGCGGTGACACAA
AGTGGACCACCATAACCGGCATATATGATGGGCTTCAGGTGTCTGATTTTCAGGCTACTGCAAAGACTTG
TGAACCTGTGGATGAACTTGTGTCATGTGGGAATGATATCTATGTTTCTGATCTATCCTCCTTACATCAA
CAATCCCCTGTT


samtools faidx Thalassophryne_amazonica.cds rna-XM_034166331.1 > Thalassophryne_amazonica---rna-XM_034166331.1.fa
samtools faidx Batrachoididae_N5.HOG0013500.fa Chatrabus_melanurus.GCA_900302635.1.OMMV01052601.1-6252-5536 > Chatrabus_melanurus.HOG0013500.fa
needle -asequence Thalassophryne_amazonica---rna-XM_034166331.1.fa -bsequence Chatrabus_melanurus.HOG0013500.fa -outfile Thalassophryne_amazonica.vs.Chatrabus_melanurus.aln -gapopen 10.0 -gapextend 0.5
needle -asequence test.fa -bsequence Chatrabus_melanurus.HOG0013500.fa -outfile Thalassophryne_amazonica.vs.Chatrabus_melanurus.aln -gapopen 10.0 -gapextend 0.5 

samtools faidx Batrachoididae_N5.HOG0013500.fa Opsanus_beta.GCA_900660325.CABFOU010245896.1:2795-3950 > Opsanus_beta.HOG0013500.fa
needle -asequence Thalassophryne_amazonica---rna-XM_034166331.1.fa -bsequence Opsanus_beta.HOG0013500.fa -outfile Thalassophryne_amazonica.vs.Opsanus_beta.aln -gapopen 10.0 -gapextend 0.5
needle -asequence prout.fa -bsequence Opsanus_beta.HOG0013500.fa -outfile Thalassophryne_amazonica.vs.Opsanus_beta.aln -gapopen 10.0 -gapextend 0.5



transeq Batrachoididae_N5.HOG0013500.fa Batrachoididae_N5.HOG0013500.prot ; sed -i 's/_1$//g' Batrachoididae_N5.HOG0013500.prot

samtools faidx Coding_sequences_alignments/N5.HOG0013500.prot Clupea_harengus_rna_XM_031568919_1 > Clupeiformes_HOG0013500.prot
samtools faidx Coding_sequences_alignments/N5.HOG0013500.prot Clupea_harengus_rna_XM_042705999_1 >> Clupeiformes_HOG0013500.prot
samtools faidx Coding_sequences_alignments/N5.HOG0013500.prot Clupea_harengus_rna_XM_031568920_2 >> Clupeiformes_HOG0013500.prot
samtools faidx Coding_sequences_alignments/N5.HOG0013500.prot Clupea_harengus_rna_XM_042708011_1 >> Clupeiformes_HOG0013500.prot
samtools faidx Coding_sequences_alignments/N5.HOG0013500.prot Clupea_harengus_rna_XM_042706001_1 >> Clupeiformes_HOG0013500.prot
samtools faidx Coding_sequences_alignments/N5.HOG0013500.prot Clupea_harengus_rna_XM_042707968_1 >> Clupeiformes_HOG0013500.prot
samtools faidx Coding_sequences_alignments/N5.HOG0013500.prot Thalassophryne_amazonica_rna_XM_034167571_1 >> Clupeiformes_HOG0013500.prot
samtools faidx Coding_sequences_alignments/N5.HOG0013500.prot Thalassophryne_amazonica_rna_XM_034166331_1 >> Clupeiformes_HOG0013500.prot
samtools faidx Coding_sequences_alignments/N5.HOG0013500.prot Sardina_pilchardus_rna_XM_062549915_1 >> Clupeiformes_HOG0013500.prot
samtools faidx Coding_sequences_alignments/N5.HOG0013500.prot Sardina_pilchardus_rna_XM_062549788_1 >> Clupeiformes_HOG0013500.prot


cat Batrachoididae_N5.HOG0013500.prot Thalassophryne_amazonica---rna-XM_034166331.1.prot Clupeiformes_HOG0013500.prot  > HOG0013500.additionalTree.prot


muscle5.1 -align HOG0013500.additionalTree.prot -output HOG0013500.additionalTree.prot.aln
trimal -in HOG0013500.additionalTree.prot.aln -gt 0.9 -out HOG0013500.additionalTree.prot.aln.trimmed


samtools faidx Coding_sequences_alignments/N5.HOG0013500.prot Scleropages_formosus_rna_XM_018734537_1 > outgroup.prot
samtools faidx Coding_sequences_alignments/N5.HOG0013500.prot Scleropages_formosus_rna_XM_029259819_1 >> outgroup.prot
samtools faidx concatenated_proteomes.fa Gasterosteus_aculeatus_aculeatus---rna-XM_040194669.1 >> outgroup.prot
samtools faidx concatenated_proteomes.fa Austrofundulus_limnaeus---rna-XM_014014583.1 >> outgroup.prot
samtools faidx concatenated_proteomes.fa Periophthalmus_magnuspinnatus---rna-XM_055224356.1 >> outgroup.prot
samtools faidx concatenated_proteomes.fa Anarrhichthys_ocellatus---rna-XM_031840404.1 >> outgroup.prot


mafft --add outgroup.prot --keeplength HOG0013500.additionalTree.prot.aln.trimmed > temp ; mv temp HOG0013500.additionalTree.prot.aln.trimmed


iqtree -s HOG0013500.additionalTree.prot.aln.trimmed -st AA -nt 8 -bb 1000 --redo -m LG+F+G4



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
sbatch --qos=1day -c 4 --mem=10G -e error.fitMG4.out -o slurm.fitMG4.out launch_fitMG4.sh N5.HOG0013500


Recipient branches : 
Node21
Thalassophryne_amazonica_rna_XM_034167571_1
Thalassophryne_amazonica_rna_XM_034166331_1


#Test positive selection and relaxed selection on receiver branch

sbatch --qos=1week -c 4 --mem=10G -e error.absrel.cand.out -o slurm.absrel.cand.out --job-name=HOG0013500 launch_absrel_cand.sh N5.HOG0013500

##Not needed to run RELAX, no branch detected under accelerated evolution
