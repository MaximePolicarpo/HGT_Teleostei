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


samtools faidx Proteomes_BUSCO80/Neoarius_graeffei.fa Neoarius_graeffei---rna-XM_060939486.1 > Neoarius_graeffei---rna-XM_060939486.1.prot

#blastp against the fish proteome (50 best hits) + Uniprot database without any fish species (10 best hits). Here I will only take 20 catffishes otherwise there are way too much of them

blastp -query Neoarius_graeffei---rna-XM_060939486.1.prot -db concatenated_proteomes.fa -outfmt 6 -out Gene_vs_FishProteome.blastp -num_threads 8 -max_target_seqs 53
blastp -query Neoarius_graeffei---rna-XM_060939486.1.prot -db non_actino_uniprot.fa -outfmt 6 -out Gene_vs_Uniprot.blastp -num_threads 8  -max_target_seqs 10 -evalue 1e-5


cut -f2 Gene_vs_FishProteome.blastp | sort | uniq  > closest_fish_seq.id
cut -f2 Gene_vs_Uniprot.blastp | sort | uniq > closest_nonfish_seq.id 

xargs samtools faidx concatenated_proteomes.fa < closest_fish_seq.id > closest_fish_seq.fa
samtools faidx Proteomes_BUSCO80/Neoarius_graeffei.fa Neoarius_graeffei---rna-XM_060938829.1 >> closest_fish_seq.fa

xargs samtools faidx non_actino_uniprot.fa  < closest_nonfish_seq.id > closest_nonfish_seq.fa

#Align with muscle and trim with trimal

cat closest_fish_seq.fa closest_nonfish_seq.fa > Uniprot_plus_closefish.fa

muscle5.1 -align closest_fish_seq.fa -output closest_fish_seq.aln
mafft --add closest_nonfish_seq.fa --keeplength closest_fish_seq.aln > Uniprot_plus_closefish.aln
trimal -in Uniprot_plus_closefish.aln -gt 0.4 -out Uniprot_plus_closefish.aln.trimmed

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


#N5.HOG0019032.clade_hgt.txt=> clade where the HGT is extracted from iTOL (from N5.HOG0019032.prot.aln.treefile)

cut -f1,2 -d "_" N5.HOG0019032.clade_hgt.txt | sort | uniq  > species_to_draw.clade1
#Choose close species without the gene .. 


echo "Megalops_atlanticus" > non_transfer_species_1
echo "Anguilla_anguilla" >> non_transfer_species_1
echo "Clupea_harengus" >> non_transfer_species_1

echo "Ictalurus_punctatus" > non_transfer_species_2
echo "Silurus_meridionalis" >> non_transfer_species_2
echo "Trichomycterus_rosablanca" >> non_transfer_species_2




rm -r Syn_tables_dir ; mkdir Syn_tables_dir


#First order the species names 

rm species_to_draw.clade1.ordered

echo "Chanos_chanos" >> species_to_draw.clade1.ordered
echo "Neoarius_graeffei" >> species_to_draw.clade1.ordered



#If there are more than 4 species from the same order, keep only 4 


#First extract genes around 

curr_OGG=N5.HOG0019032

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


grep "NC_045157_1\|NW_024880005_1" Syn_tables_dir/Clupea_harengus.synt.final.df  > temp ; mv temp Syn_tables_dir/Clupea_harengus.synt.final.df
sort Syn_tables_dir/Clupea_harengus.synt.final.df | uniq > temp ; mv temp Syn_tables_dir/Clupea_harengus.synt.final.df



#First add Clupea_harengus

curr_OGG=N5.HOG0019032
curr_sp=Clupea_harengus
ref_sp=Chanos_chanos


grep -A10 -B10 "N5_HOG0026677" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0027283" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0018366" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "NC_045159" >> Syn_tables_dir/$curr_sp.synt.df

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

curr_OGG=N5.HOG0019032
curr_sp=Anguilla_anguilla
ref_sp=Chanos_chanos


grep -A15 -B15 "N5_HOG0026677" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0018366" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "NC_049209" >> Syn_tables_dir/$curr_sp.synt.df



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




#First add Megalops_atlanticus

curr_OGG=N5.HOG0019032
curr_sp=Megalops_atlanticus
ref_sp=Chanos_chanos


grep -A15 -B15 "N5_HOG0026677" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0018366" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "CM032883" >> Syn_tables_dir/$curr_sp.synt.df


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

curr_OGG=N5.HOG0019032
curr_sp=Ictalurus_punctatus
ref_sp=Neoarius_graeffei


grep -A10 -B10 "N5_HOG0027865" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df

grep -A10 -B10 "N5_HOG0049184" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret





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

curr_OGG=N5.HOG0019032
curr_sp=Silurus_meridionalis
ref_sp=Neoarius_graeffei


grep -A10 -B10 "N5_HOG0027865" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df

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




#Now add Trichomycterus_rosablanca

curr_OGG=N5.HOG0019032
curr_sp=Trichomycterus_rosablanca
ref_sp=Neoarius_graeffei


grep -A10 -B10 "N5_HOG0025334" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df

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


samtools faidx Proteomes_BUSCO80/Chanos_chanos.fa Chanos_chanos---rna-XM_030788612.1 > Chanos_chanos---rna-XM_030788612.1.prot


exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs" Neoarius_graeffei---rna-XM_060939486.1.prot GCF_001660625.3_Coco_2.0_genomic.fna > Ictalurus_punctatus.whole.fa
exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs" Neoarius_graeffei---rna-XM_060939486.1.prot GCF_014805685.1_ASM1480568v1_genomic.fna > Silurus_meridionalis.whole.fa
exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs" Neoarius_graeffei---rna-XM_060939486.1.prot GCF_030014385.1_fTriRos1.hap1_genomic.fna > Trichomycterus_rosablanca.whole.fa


exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs" Chanos_chanos---rna-XM_030788612.1.prot GCA_019176425.1_MATL_1.0_genomic.fna > Megalops_atlanticus.whole.fa
exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs" Chanos_chanos---rna-XM_030788612.1.prot GCF_013347855.1_fAngAng1.pri_genomic.fna > Anguilla_anguilla.whole.fa
exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs" Chanos_chanos---rna-XM_030788612.1.prot GCF_900700415.2_Ch_v2.0.2_genomic.fna > Clupea_harengus.whole.fa



=> no related gene.. 

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

#Neoarius_graeffei
grep -A3 -B3 "N5_HOG0019032"  GFF3_N5_OGGs/Neoarius_graeffei.gff.simplified.sorted.OGG.tiret

samtools faidx GCF_027579695.1_fNeoGra1.pri_genomic.fna NC_083582.1:27011828-27049295 > N5.HOG0019032.Neoarius_graeffei.extended.1.fa
sed -i 's/:/-/g' N5.HOG0019032.Neoarius_graeffei.extended.1.fa
makeblastdb -in N5.HOG0019032.Neoarius_graeffei.extended.1.fa -dbtype nucl

samtools faidx GCF_027579695.1_fNeoGra1.pri_genomic.fna NC_083582.1:27282822-27351468 > N5.HOG0019032.Neoarius_graeffei.extended.2.fa
sed -i 's/:/-/g' N5.HOG0019032.Neoarius_graeffei.extended.2.fa
makeblastdb -in N5.HOG0019032.Neoarius_graeffei.extended.2.fa -dbtype nucl

samtools faidx GCF_027579695.1_fNeoGra1.pri_genomic.fna NC_083582.1:27848466-27905574 > N5.HOG0019032.Neoarius_graeffei.extended.3.fa
sed -i 's/:/-/g' N5.HOG0019032.Neoarius_graeffei.extended.3.fa
makeblastdb -in N5.HOG0019032.Neoarius_graeffei.extended.3.fa -dbtype nucl


#Chanos_chanos
grep -A3 -B3 "N5_HOG0019032"  GFF3_N5_OGGs/Chanos_chanos.gff.simplified.sorted.OGG.tiret
samtools faidx GCF_902362185.1_fChaCha1.1_genomic.fna NC_044506.1:22802197-22929311 > N5.HOG0019032.Chanos_chanos.extended.1.fa
samtools faidx GCF_902362185.1_fChaCha1.1_genomic.fna NC_044509.1:5148601-5167851 > N5.HOG0019032.Chanos_chanos.extended.2.fa

sed -i 's/:/-/g' N5.HOG0019032.Chanos_chanos.extended.1.fa
makeblastdb -in N5.HOG0019032.Chanos_chanos.extended.1.fa -dbtype nucl
sed -i 's/:/-/g' N5.HOG0019032.Chanos_chanos.extended.2.fa
makeblastdb -in N5.HOG0019032.Chanos_chanos.extended.2.fa -dbtype nucl


#Now lets find TEs using the repbase and dfam databases

tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0019032.Neoarius_graeffei.extended.1.fa -evalue 1e-5 -outfmt 6 -out TE.Neoarius_graeffei.1.tblastn -num_threads 8
sed -i 's/#//g' TE.Neoarius_graeffei.1.tblastn

tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0019032.Neoarius_graeffei.extended.2.fa -evalue 1e-5 -outfmt 6 -out TE.Neoarius_graeffei.2.tblastn -num_threads 8
sed -i 's/#//g' TE.Neoarius_graeffei.2.tblastn

tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0019032.Neoarius_graeffei.extended.3.fa -evalue 1e-5 -outfmt 6 -out TE.Neoarius_graeffei.3.tblastn -num_threads 8
sed -i 's/#//g' TE.Neoarius_graeffei.3.tblastn

tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0019032.Chanos_chanos.extended.1.fa -evalue 1e-5 -outfmt 6 -out TE.Chanos_chanos.1.tblastn -num_threads 8
sed -i 's/#//g' TE.Chanos_chanos.1.tblastn

tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0019032.Chanos_chanos.extended.2.fa -evalue 1e-5 -outfmt 6 -out TE.Chanos_chanos.2.tblastn -num_threads 8
sed -i 's/#//g' TE.Chanos_chanos.2.tblastn

#merge tblastn hits and find the best TE match by doing a blastx


Rscript Rscript_merge_blast_hits.R TE.Neoarius_graeffei.1.tblastn TE.Neoarius_graeffei.1.tblastn.merged
Rscript Rscript_merge_blast_hits.R TE.Neoarius_graeffei.1.tblastn TE.Neoarius_graeffei.2.tblastn.merged
Rscript Rscript_merge_blast_hits.R TE.Neoarius_graeffei.1.tblastn TE.Neoarius_graeffei.3.tblastn.merged

Rscript Rscript_merge_blast_hits.R TE.Chanos_chanos.1.tblastn TE.Chanos_chanos.1.tblastn.merged
Rscript Rscript_merge_blast_hits.R TE.Chanos_chanos.2.tblastn TE.Chanos_chanos.2.tblastn.merged


xargs samtools faidx N5.HOG0019032.Neoarius_graeffei.extended.1.fa < TE.Neoarius_graeffei.1.tblastn.merged > TE.Neoarius_graeffei.1.BEST.fa
blastx -query TE.Neoarius_graeffei.1.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Neoarius_graeffei.1.BEST.blastx -max_target_seqs 1
cut -f1 TE.Neoarius_graeffei.1.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Neoarius_graeffei.1.BEST.blastx >> temp  ; done ; mv temp TE.Neoarius_graeffei.1.BEST.blastx

xargs samtools faidx N5.HOG0019032.Neoarius_graeffei.extended.2.fa < TE.Neoarius_graeffei.2.tblastn.merged > TE.Neoarius_graeffei.2.BEST.fa
blastx -query TE.Neoarius_graeffei.2.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Neoarius_graeffei.2.BEST.blastx -max_target_seqs 1
cut -f1 TE.Neoarius_graeffei.2.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Neoarius_graeffei.2.BEST.blastx >> temp  ; done ; mv temp TE.Neoarius_graeffei.2.BEST.blastx


xargs samtools faidx N5.HOG0019032.Neoarius_graeffei.extended.3.fa < TE.Neoarius_graeffei.3.tblastn.merged > TE.Neoarius_graeffei.3.BEST.fa
blastx -query TE.Neoarius_graeffei.3.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Neoarius_graeffei.3.BEST.blastx -max_target_seqs 1
cut -f1 TE.Neoarius_graeffei.3.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Neoarius_graeffei.3.BEST.blastx >> temp  ; done ; mv temp TE.Neoarius_graeffei.3.BEST.blastx


xargs samtools faidx N5.HOG0019032.Chanos_chanos.extended.1.fa < TE.Chanos_chanos.1.tblastn.merged > TE.Chanos_chanos.1.BEST.fa
blastx -query TE.Chanos_chanos.1.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Chanos_chanos.1.BEST.blastx -max_target_seqs 1
cut -f1  TE.Chanos_chanos.1.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Chanos_chanos.1.BEST.blastx >> temp  ; done ; mv temp TE.Chanos_chanos.1.BEST.blastx


xargs samtools faidx  N5.HOG0019032.Chanos_chanos.extended.2.fa < TE.Chanos_chanos.2.tblastn.merged > TE.Chanos_chanos.2.BEST.fa
blastx -query TE.Chanos_chanos.2.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Chanos_chanos.2.BEST.blastx -max_target_seqs 1
cut -f1  TE.Chanos_chanos.2.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Chanos_chanos.2.BEST.blastx >> temp  ; done ; mv temp TE.Chanos_chanos.2.BEST.blastx


#Now find shared elements

cut -f2 TE.Neoarius_graeffei.1.BEST.blastx | sort | uniq > TE.Neoarius_graeffei.1.uniqTE
cut -f2 TE.Neoarius_graeffei.2.BEST.blastx | sort | uniq > TE.Neoarius_graeffei.2.uniqTE
cut -f2 TE.Neoarius_graeffei.3.BEST.blastx | sort | uniq > TE.Neoarius_graeffei.3.uniqTE

cut -f2 TE.Chanos_chanos.1.BEST.blastx | sort | uniq > TE.Chanos_chanos.1.uniqTE
cut -f2 TE.Chanos_chanos.2.BEST.blastx | sort | uniq > TE.Chanos_chanos.2.uniqTE

comm -12 TE.Neoarius_graeffei.1.uniqTE TE.Chanos_chanos.1.uniqTE
comm -12 TE.Neoarius_graeffei.1.uniqTE TE.Chanos_chanos.2.uniqTE

comm -12 TE.Neoarius_graeffei.2.uniqTE TE.Chanos_chanos.1.uniqTE
comm -12 TE.Neoarius_graeffei.2.uniqTE TE.Chanos_chanos.2.uniqTE

comm -12 TE.Neoarius_graeffei.3.uniqTE TE.Chanos_chanos.1.uniqTE
comm -12 TE.Neoarius_graeffei.3.uniqTE TE.Chanos_chanos.2.uniqTE

### 1 shared TEs

#Merlin-1_DR_1p:ClassII:TIR:Merlin:?


#### Now let's find every copy of these elements in the genome of both species + closely related species


samtools faidx Dfam_plus_Repbase.cdhit80.prot Merlin-1_DR_1p:ClassII:TIR:Merlin:? > Merlin-1.prot

#Neoarius_graeffei
tblastn -query Merlin-1.prot -db GCF_027579695.1_fNeoGra1.pri_genomic.fna -evalue 1e-1 -outfmt 6 -out Merlin-1.Neoarius_graeffei.tblastn -num_threads 8
sed -i 's/#//g' Merlin-1.Neoarius_graeffei.tblastn
Rscript Rscript_merge_blast_hits.R Merlin-1.Neoarius_graeffei.tblastn Merlin-1.Neoarius_graeffei.tblastn.merged
xargs samtools faidx GCF_027579695.1_fNeoGra1.pri_genomic.fna < Merlin-1.Neoarius_graeffei.tblastn.merged > Merlin-1.Neoarius_graeffei.BEST.fa
diamond blastx --query Merlin-1.Neoarius_graeffei.BEST.fa --db Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out Merlin-1.Neoarius_graeffei.BEST.blastx --max-target-seqs 1
grep "Merlin-1" Merlin-1.Neoarius_graeffei.BEST.blastx > Merlin-1.Neoarius_graeffei.list

#Chanos_chanos
tblastn -query Merlin-1.prot -db GCF_902362185.1_fChaCha1.1_genomic.fna -evalue 1e-1 -outfmt 6 -out Merlin-1.Chanos_chanos.tblastn -num_threads 8
sed -i 's/#//g' Merlin-1.Chanos_chanos.tblastn
Rscript Rscript_merge_blast_hits.R Merlin-1.Chanos_chanos.tblastn Merlin-1.Chanos_chanos.tblastn.merged
xargs samtools faidx GCF_902362185.1_fChaCha1.1_genomic.fna  < Merlin-1.Chanos_chanos.tblastn.merged > Merlin-1.Chanos_chanos.BEST.fa
diamond blastx --query Merlin-1.Chanos_chanos.BEST.fa --db Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out Merlin-1.Chanos_chanos.BEST.blastx --max-target-seqs 1
grep "Merlin-1" Merlin-1.Chanos_chanos.BEST.blastx > Merlin-1.Chanos_chanos.list


#Megalops_atlanticus
tblastn -query Merlin-1.prot -db GCA_019176425.1_MATL_1.0_genomic.fna -evalue 1e-1 -outfmt 6 -out Merlin-1.Megalops_atlanticus.tblastn -num_threads 8
sed -i 's/#//g' Merlin-1.Megalops_atlanticus.tblastn
Rscript Rscript_merge_blast_hits.R Merlin-1.Megalops_atlanticus.tblastn Merlin-1.Megalops_atlanticus.tblastn.merged
xargs samtools faidx GCA_019176425.1_MATL_1.0_genomic.fna  < Merlin-1.Megalops_atlanticus.tblastn.merged > Merlin-1.Megalops_atlanticus.BEST.fa
diamond blastx --query Merlin-1.Megalops_atlanticus.BEST.fa --db Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out Merlin-1.Megalops_atlanticus.BEST.blastx --max-target-seqs 1
grep "Merlin-1" Merlin-1.Megalops_atlanticus.BEST.blastx > Merlin-1.Megalops_atlanticus.list

#Anguilla_anguilla
tblastn -query Merlin-1.prot -db GCF_013347855.1_fAngAng1.pri_genomic.fna -evalue 1e-1 -outfmt 6 -out Merlin-1.Anguilla_anguilla.tblastn -num_threads 8
sed -i 's/#//g' Merlin-1.Anguilla_anguilla.tblastn
Rscript Rscript_merge_blast_hits.R Merlin-1.Anguilla_anguilla.tblastn Merlin-1.Anguilla_anguilla.tblastn.merged
xargs samtools faidx GCF_013347855.1_fAngAng1.pri_genomic.fna  < Merlin-1.Anguilla_anguilla.tblastn.merged > Merlin-1.Anguilla_anguilla.BEST.fa
diamond blastx --query Merlin-1.Anguilla_anguilla.BEST.fa --db Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out Merlin-1.Anguilla_anguilla.BEST.blastx --max-target-seqs 1
grep "Merlin-1" Merlin-1.Anguilla_anguilla.BEST.blastx > Merlin-1.Anguilla_anguilla.list


#Clupea_harengus
tblastn -query Merlin-1.prot -db GCF_900700415.2_Ch_v2.0.2_genomic.fna -evalue 1e-1 -outfmt 6 -out Merlin-1.Clupea_harengus.tblastn -num_threads 8
sed -i 's/#//g' Merlin-1.Clupea_harengus.tblastn
Rscript Rscript_merge_blast_hits.R Merlin-1.Clupea_harengus.tblastn Merlin-1.Clupea_harengus.tblastn.merged
xargs samtools faidx GCF_900700415.2_Ch_v2.0.2_genomic.fna  < Merlin-1.Clupea_harengus.tblastn.merged > Merlin-1.Clupea_harengus.BEST.fa
diamond blastx --query Merlin-1.Clupea_harengus.BEST.fa --db Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out Merlin-1.Clupea_harengus.BEST.blastx --max-target-seqs 1
grep "Merlin-1" Merlin-1.Clupea_harengus.BEST.blastx > Merlin-1.Clupea_harengus.list


#Ictalurus_punctatus
tblastn -query Merlin-1.prot -db GCF_001660625.3_Coco_2.0_genomic.fna -evalue 1e-1 -outfmt 6 -out Merlin-1.Ictalurus_punctatus.tblastn -num_threads 8
sed -i 's/#//g' Merlin-1.Ictalurus_punctatus.tblastn
Rscript Rscript_merge_blast_hits.R Merlin-1.Ictalurus_punctatus.tblastn Merlin-1.Ictalurus_punctatus.tblastn.merged
xargs samtools faidx GCF_001660625.3_Coco_2.0_genomic.fna  < Merlin-1.Ictalurus_punctatus.tblastn.merged > Merlin-1.Ictalurus_punctatus.BEST.fa
diamond blastx --query Merlin-1.Ictalurus_punctatus.BEST.fa --db Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out Merlin-1.Ictalurus_punctatus.BEST.blastx --max-target-seqs 1
grep "Merlin-1" Merlin-1.Ictalurus_punctatus.BEST.blastx > Merlin-1.Ictalurus_punctatus.list

#Silurus_meridionalis
tblastn -query Merlin-1.prot -db GCF_014805685.1_ASM1480568v1_genomic.fna -evalue 1e-1 -outfmt 6 -out Merlin-1.Silurus_meridionalis.tblastn -num_threads 8
sed -i 's/#//g' Merlin-1.Silurus_meridionalis.tblastn
Rscript Rscript_merge_blast_hits.R Merlin-1.Silurus_meridionalis.tblastn Merlin-1.Silurus_meridionalis.tblastn.merged
xargs samtools faidx GCF_014805685.1_ASM1480568v1_genomic.fna  < Merlin-1.Silurus_meridionalis.tblastn.merged > Merlin-1.Silurus_meridionalis.BEST.fa
diamond blastx --query Merlin-1.Silurus_meridionalis.BEST.fa --db Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out Merlin-1.Silurus_meridionalis.BEST.blastx --max-target-seqs 1
grep "Merlin-1" Merlin-1.Silurus_meridionalis.BEST.blastx > Merlin-1.Silurus_meridionalis.list

#Trichomycterus_rosablanca
tblastn -query Merlin-1.prot -db GCF_030014385.1_fTriRos1.hap1_genomic.fna -evalue 1e-1 -outfmt 6 -out Merlin-1.Trichomycterus_rosablanca.tblastn -num_threads 8
sed -i 's/#//g' Merlin-1.Trichomycterus_rosablanca.tblastn
Rscript Rscript_merge_blast_hits.R Merlin-1.Trichomycterus_rosablanca.tblastn Merlin-1.Trichomycterus_rosablanca.tblastn.merged
xargs samtools faidx GCF_030014385.1_fTriRos1.hap1_genomic.fna  < Merlin-1.Trichomycterus_rosablanca.tblastn.merged > Merlin-1.Trichomycterus_rosablanca.BEST.fa
diamond blastx --query Merlin-1.Trichomycterus_rosablanca.BEST.fa --db Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out Merlin-1.Trichomycterus_rosablanca.BEST.blastx --max-target-seqs 1
grep "Merlin-1" Merlin-1.Trichomycterus_rosablanca.BEST.blastx > Merlin-1.Trichomycterus_rosablanca.list


wc -l Merlin-1.Megalops_atlanticus.list
wc -l Merlin-1.Anguilla_anguilla.list
wc -l Merlin-1.Clupea_harengus.list
wc -l Merlin-1.Chanos_chanos.list
wc -l Merlin-1.Neoarius_graeffei.list
wc -l Merlin-1.Ictalurus_punctatus.list
wc -l Merlin-1.Silurus_meridionalis.list
wc -l Merlin-1.Trichomycterus_rosablanca.list


### In a FINAL step we will make another synteny plot more zoomed


## Neoarius_graeffei ## first scaffold

rm seq_clustered_infos_ogg.TE.txt ; rm clusters_ID_TE.txt


scaffold=`grep ">" N5.HOG0019032.Neoarius_graeffei.extended.1.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0019032.Neoarius_graeffei.extended.1.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0019032.Neoarius_graeffei.extended.1.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Neoarius_graeffei,$scaffold,$length" > clusters_ID_TE.txt

grep "N5_HOG0019032"  GFF3_N5_OGGs/Neoarius_graeffei.gff.simplified.sorted.OGG.tiret

grep "XM_060939486" GFF3_files_per_species/Neoarius_graeffei.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons


nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene1.$nbexon,N5_HOG0019032,+,Neoarius_graeffei" >> seq_clustered_infos_ogg.TE.txt
done

cut -f1 TE.Neoarius_graeffei.1.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Neoarius_graeffei.1.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0019032.Neoarius_graeffei.extended.1.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0019032.Neoarius_graeffei.extended.1.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold,/g" | sed "s/$/,$strand,Neoarius_graeffei/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done



## Neoarius_graeffei ## second region


scaffold=`grep ">" N5.HOG0019032.Neoarius_graeffei.extended.2.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0019032.Neoarius_graeffei.extended.2.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0019032.Neoarius_graeffei.extended.2.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Neoarius_graeffei,$scaffold.sec,$length" >> clusters_ID_TE.txt

grep "N5_HOG0019032"  GFF3_N5_OGGs/Neoarius_graeffei.gff.simplified.sorted.OGG.tiret

grep "XM_060938827" GFF3_files_per_species/Neoarius_graeffei.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons


nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold.sec,$real_exon_start,$real_exon_stop,gene2.$nbexon,N5_HOG0019032,-,Neoarius_graeffei" >> seq_clustered_infos_ogg.TE.txt
done


## Neoarius_graeffei ## third region


scaffold=`grep ">" N5.HOG0019032.Neoarius_graeffei.extended.3.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0019032.Neoarius_graeffei.extended.3.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0019032.Neoarius_graeffei.extended.3.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Neoarius_graeffei,$scaffold.third,$length" >> clusters_ID_TE.txt

grep "N5_HOG0019032"  GFF3_N5_OGGs/Neoarius_graeffei.gff.simplified.sorted.OGG.tiret

grep "XM_060938829" GFF3_files_per_species/Neoarius_graeffei.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons


nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold.third,$real_exon_start,$real_exon_stop,gene3.$nbexon,N5_HOG0019032,+,Neoarius_graeffei" >> seq_clustered_infos_ogg.TE.txt
done



## Chanos_chanos


scaffold=`grep ">" N5.HOG0019032.Chanos_chanos.extended.1.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0019032.Chanos_chanos.extended.1.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0019032.Chanos_chanos.extended.1.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Chanos_chanos,$scaffold,$length" >> clusters_ID_TE.txt

grep -A3 -B3 "N5_HOG0019032"  GFF3_N5_OGGs/Chanos_chanos.gff.simplified.sorted.OGG.tiret

grep "XM_030788611" GFF3_files_per_species/Chanos_chanos.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons
grep "XM_030788612" GFF3_files_per_species/Chanos_chanos.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_2.exons
grep "XM_030788614" GFF3_files_per_species/Chanos_chanos.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_3.exons
grep "XM_030788615" GFF3_files_per_species/Chanos_chanos.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_4.exons
grep "XM_030788616" GFF3_files_per_species/Chanos_chanos.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_5.exons


nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene1.$nbexon,N5_HOG0019032,+,Chanos_chanos" >> seq_clustered_infos_ogg.TE.txt
done

nbexon=0
for line in `cat gene_2.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene2.$nbexon,N5_HOG0019032,-,Chanos_chanos" >> seq_clustered_infos_ogg.TE.txt
done
nbexon=0
for line in `cat gene_3.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene3.$nbexon,N5_HOG0019032,+,Chanos_chanos" >> seq_clustered_infos_ogg.TE.txt
done
nbexon=0
for line in `cat gene_4.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene4.$nbexon,N5_HOG0019032,-,Chanos_chanos" >> seq_clustered_infos_ogg.TE.txt
done
nbexon=0
for line in `cat gene_5.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene5.$nbexon,N5_HOG0019032,+,Chanos_chanos" >> seq_clustered_infos_ogg.TE.txt
done

cut -f1 TE.Chanos_chanos.1.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Chanos_chanos.1.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0019032.Chanos_chanos.extended.1.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0019032.Chanos_chanos.extended.1.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold,/g" | sed "s/$/,$strand,Chanos_chanos/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done



#Second region of Chanos


scaffold=`grep ">" N5.HOG0019032.Chanos_chanos.extended.2.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0019032.Chanos_chanos.extended.2.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0019032.Chanos_chanos.extended.2.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Chanos_chanos,$scaffold,$length" >> clusters_ID_TE.txt

grep "N5_HOG0019032"  GFF3_N5_OGGs/Chanos_chanos.gff.simplified.sorted.OGG.tiret
grep "XM_030792951" GFF3_files_per_species/Chanos_chanos.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons



nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene6.$nbexon,N5_HOG0019032,-,Chanos_chanos" >> seq_clustered_infos_ogg.TE.txt
done



cut -f1 TE.Chanos_chanos.2.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Chanos_chanos.2.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0019032.Chanos_chanos.extended.2.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0019032.Chanos_chanos.extended.2.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold,/g" | sed "s/$/,$strand,Chanos_chanos/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done







#### Check the observed conserved non coding region =


samtools faidx N5.HOG0019032.Neoarius_graeffei.extended.1.fa.rev NC_083582.1-27011828-27049295:31990-32700 > region1.Neoarius_graeffei.fa
samtools faidx N5.HOG0019032.Chanos_chanos.extended.1.dotplot.fa NC_044506.1-22821135-22851623:10702-11442 > region1.Chanos_chanos.fa

sed -i 's/>.*/>N.graeffei/g' region1.Neoarius_graeffei.fa
sed -i 's/>.*/>C.chanos/g' region1.Chanos_chanos.fa

needle -asequence region1.Neoarius_graeffei.fa -bsequence region1.Chanos_chanos.fa -outfile region1.aln -gapopen 10.0 -gapextend 0.5

sbatch --qos=6hours -c 10 --mem=20G launch_blast_allgenomes.sh region1.Neoarius_graeffei.fa 1e-20 region1.blastn.tsv




##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##======================Search for the gene in other Neoarius non-annotated genomes + alternative assembly of C. chanos  =============
##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================

tblastn -query Neoarius_graeffei---rna-XM_060939486.1.prot -db GCA_025727735.1_ASM2572773v1_genomic.fna > Neoarius_leptaspis.tblastn
tblastn -query Neoarius_graeffei---rna-XM_060939486.1.prot -db GCA_027576835.1_ASM2757683v1_genomic.fna > Neoarius_berneyi.tblastn
tblastn -query Neoarius_graeffei---rna-XM_060939486.1.prot -db GCA_029814605.1_ASM2981460v1_genomic.fna > Neoarius_pectoralis.tblastn



exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs" Neoarius_graeffei---rna-XM_060939486.1.prot ../../Genomic_data/Neoarius_leptaspis/GCA_025727735.1_ASM2572773v1_genomic.fna > Neoarius_leptaspis.whole.fa
exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs" Neoarius_graeffei---rna-XM_060939486.1.prot ../../Genomic_data/Neoarius_berneyi/GCA_027576835.1_ASM2757683v1_genomic.fna > Neoarius_berneyi.whole.fa
exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs" Neoarius_graeffei---rna-XM_060939486.1.prot ../../Genomic_data/Neoarius_pectoralis/GCA_029814605.1_ASM2981460v1_genomic.fna > Neoarius_pectoralis.whole.fa


samtools faidx GCA_025727735.1_ASM2572773v1_genomic.fna JAODTU010000006.1:15288793-15503138 > Neoarius_leptaspis.PotentialRegions.fa
samtools faidx GCA_025727735.1_ASM2572773v1_genomic.fna JAODTU010000024.1:7650710-7664992 >> Neoarius_leptaspis.PotentialRegions.fa

samtools faidx GCA_027576835.1_ASM2757683v1_genomic.fna JAODIF010000006.1:15939796-15974075 > Neoarius_berneyi.PotentialRegions.fa

samtools faidx GCA_029814605.1_ASM2981460v1_genomic.fna JAPTHB010000006.1:15128827-15413399 > Neoarius_pectoralis.PotentialRegions.fa
samtools faidx GCA_029814605.1_ASM2981460v1_genomic.fna JAPTHB010000022.1:70174794-70186686 >> Neoarius_pectoralis.PotentialRegions.fa


exonerate  -E True --showtargetgff TRUE --model protein2genome --ryo "%tcs" Neoarius_graeffei---rna-XM_060939486.1.prot Neoarius_leptaspis.PotentialRegions.fa > Neoarius_leptaspis.region.exo
exonerate  -E True --showtargetgff TRUE --model protein2genome --ryo "%tcs" Neoarius_graeffei---rna-XM_060939486.1.prot Neoarius_berneyi.PotentialRegions.fa > Neoarius_berneyi.region.exo
exonerate  -E True --showtargetgff TRUE --model protein2genome --ryo "%tcs" Neoarius_graeffei---rna-XM_060939486.1.prot Neoarius_pectoralis.PotentialRegions.fa > Neoarius_pectoralis.region.exo


grep "Chanos_chanos\|Neoarius" Coding_sequences_alignments/N5.HOG0019032.prot | sed 's/>//g' 



samtools faidx Coding_sequences_alignments/N5.HOG0019032.prot Chanos_chanos_rna_XM_030788615_1 > Chanos_chanos_rna_XM_030788615_1.prot
tblastn -query Chanos_chanos_rna_XM_030788615_1.prot -db ../../Alternative_assemblies/Chanos_chanos/GCA_018691325.1_CCHAN_1.0_genomic.fna > Chanos_chanos.tblastn
exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs" Chanos_chanos_rna_XM_030788615_1.prot ../../Alternative_assemblies/Chanos_chanos/GCA_018691325.1_CCHAN_1.0_genomic.fna > Chanos_chanos.whole.exo

samtools faidx GCA_018691325.1_CCHAN_1.0_genomic.fna JAFFIB010000115.1:1827722-1847722 > Chanos_chanos.PotentialRegions.fa
samtools faidx GCA_018691325.1_CCHAN_1.0_genomic.fna JAFFIB010000115.1:1827722-1847722 > Chanos_chanos.PotentialRegions.fa




#Neoarius_spp.HOG0019032.fa
>Neoarius_leptaspis.GCA_025727735.1.JAODTU010000006.1-15288793-15503138
TTTTACACAGGACTCTGTGGCAGTGCCCTCGCTAGTAGACAAGAGATACTGAAATGGTTATGGGAACATA
TCAAAGaattactTCCATTTACGCAATACTCAGCCCAGGACCAAAAGGTGCTCTTAACTCTCGTGGAGCT
GAACAAATGTGTTCTGAACAAGGAAGACCCCACAAGTCTAATCACAGAAATGACTGTGCAGCTAGAAGga
aagctGGATAAACACGGGAGAGCAATGTTCATCTCTAAAGAAGAGTTACTATATGACCAGAAGCATCAGA
CACTGAAGGAATACACACTGATGATATGTAAAAATTTCCAAGATTGGTTTGAGGAtgctctggtttcccc
cacagttcaaagacatgcaggGAGGCGAAATGCAAGGATTATACTAAAGATCATTCAGCTTCTTAATGAC
TGGACATAGGGAACAACATACAATTTTCTCCACAAAATGACAGTGGAAGATGTGGACCCCTCTGTTCTGA
GGGAGGAGGAAAGAAATTATATAGATGCTGGGCTGCTGGACAACTCTCCTTATTTCCCAATGCTGAGAAA
AGAGAGACACATTGATATCTTCATCTCCCTTGACTTCAGTGCTGGAGACCCCATGGAGGTAATG
>Neoarius_leptaspis.GCA_025727735.1.JAODTU010000024.1-7650710-7664992
ATGTCTCAGAGCAAACCTATCCAACGTCCTGAAGTCAGGACTGGTCATTCTCTGAATGAAGGAGAGCAGG
ACCATGTTTCCAGGAGGGGAGAGACTGTTCTGCGGTGCTTAAAGAGACATGGCATCAGCTGCAATCAGGA
TACAATGCCCAACATTGCACTGCTGGGCTCTGGGGGAGGTGAGCGAGCCATGCTGGGGTTACTGGGATCA
CTGGTTCAACTAAAGAACTGTGATCTCCTGGACAGCATAATGTACCTGAGTGGAGTCTCAGGATCCACAT
GGTGCATGGCATGCTTGTATAAAGAGTCAGACTGATCCACCAAACTGGCCGCTGTAAAGGAGAATATCAT
TAAAAGGCTTGATGGGCCTGCAGTCAGCTGGTCGGATGCTTGGAAAAAACTGAACAAATACTATGACCAG
AAGGACAACTTCTCCCTGACAGATGTCTGGGCTGTAGTGTTTGTCACTATGATCGTGAAAGAGATTGATG
AAAACACTATAAGCGgtcagagagagaatcacagcaAAGACCCATATCCCATCTACACTGTGATTGACAA
GCAGTGCAAAAAGAAGAAACTGGACCGAGATGTCTGGTTTGAGGTTgcgcctcatgaagctggctacTCC
CTCACTGGTGCGTTTGTGGACTCTGTCTATTTTGGGAGTCAGTTTGAAGAAGGACGTTTAGTGAAAAAAC
AACCTGAGATTGACATGCTGTATCTGCAAGGTATCTGT
>Neoarius_berneyi.GCA_027576835.1.JAODIF010000006.1-15939796-15974075
ATGTCTCAGAGCAAACCTATCCAACGTACTGAAGTCAGGATTGGTCATTCTCTGAATGAAGGAGAGCAGG
ACCATGTTTCCAGGAGGAGAGAGACTGTTCTGCAGTGCTTAAAGAGACATGGCATCAGCTGCAGTCAGGA
TACAATGCCCAACATTGCACTGCTGGGCTCTGGGGGAGGTGAGCGAGCCATGCTGGGGTTACTGGGATCA
CTGGTTCAACTAAAGAACTGTGATCTCCTGGACAGCATAATGTACCTGAGTGGAGTCTCAGGATCCACAT
GGTGCATGGCATCCTTATATAAAGAGTCAGACTGGTCCACCAAACTGGCCGCTGTGAAGGAGAATATCAT
TAAAAGGCTTGATGGGCCTGCAGTCAGCTGGTCAGATGCTTGGAAAAAACTGAACAAATACTATGACCAG
AATGACAACTTCTCTCTGACAGACGTCTGGGCTGTAATGTTTGTCACTATGATCATGAAAGAGATTGATG
AAAACACTATAAGCGgtcagagagagaatcacagcaAAGACCCATATCCCATCTACACTGTGATTGACAA
GCAGTGCAAAAAGAAGAAACTGAACCGAGATGTCTGGTTTGAGGTTgcgcctcatgaagctggctacTCC
CTCACTGGTGCGTTTGTGGACTCTGCCTATTTTGGGAGTCGGTTTGAAGAAGGACGTTTAGTGAAAAAAC
AACCTGAGATTGACATGCTGTATCTGCAAGGTATCTGT
>Neoarius_pectoralis.GCA_029814605.1.JAPTHB010000022.1-70174794-70186686
ACTGAAGTCAGGATTGGTCATTCTCTGAATGAAGGAGAGCAGGACCATGTTTCCAGGAGGGGAGAGACTG
TTCTGCAGTGCTTAAAGAGACATGGCATCAACTGCAGTCAGGCTACAATGCCCAACATTGCACTGCTGGG
CTCTGGGGGAGGTGAGCGAGCCATGCTGGGGTTACTGGGATCACTGGTTCAACTAAAGAACTGTGATCTC
CTGGACAGCATAATGTACCTGAGTGGAGTCTCAGGATCCACATGGTGCATGACATCCTTATACAAAGAGT
CAGACTGGTCCACCAAACTGGCCGCTGTGAAGGAGAATATCATTAAAAGGCTTGATGGGTCTGCAGTCAG
CTGGTCGGATGCttggaaaaaaatgaacaaatactATGACCAGAATGACAACTTCTCCCTGACAGACGTC
TGGGCTGTAATGTTTGTCACTACGATCATGAAAGAGATTGATGAAAACACTATAAGCGgtcagagagaga
atcacagcAAAGACCCATATCCCATCTACACTGTGATTGACAAGCAGTGCAAAAAGAAGAAACTGAACCG
AGATGTCTGGTTTGAGGTTgcgcctcatgaagctggctacTCCCTCACTGGTGCGTTTGTGGACTCTGCC
TATTTTGGGAGTCGGTTTGAAGAAGGACGTTTAGTGAAAAAACAACCTGAGATTGACATGCTGTATCTGC
AAGGTATCTGT


## Chanos.alternative.fa
>C.chanos.GCA_018691325.1.JAFFIB010000115.1-1837168-1845440
ATGCCCAACATTGCACTGGTGGGCTCTGGAGGAGGCGAGCGAGCCATGCTGGGATTACTGGGATCACTGG
TTCAGTTGAAGAAATGTGATCTGCTGGACTGtatactgtatctgtgtggaaCTTCAGGATCCACATGGTG
CATGGCGTCCTTATACAAAGAACCAGACTGGTCCACTCATCTGAATGCTGTGCTGAAAGACATCATGGAA
AGGCTTGGGGCTAAAGTCAGCCTGCAGGATGCTTGGACTAAACTGCACAAATACCATGATGAGCCCAACT
TCTCCCTCACACATGTCTGGGCTGTAATGTTTGTCTCTATGATCGTGAAAGAGATTGATGAGAACACTAT
AAGTAATCAGAAAGGAAATCACAGCAAAGACCCATATCCCATCTACACTGTGATTGACAAGCAGTACAAA
CAGAAGAAACTGAACCAAGATGCTTGGTTTGAGATTACACCTCATGAAGCTGGATACTCTGTCACTGGGG
CTTTTGTGGACTCAGCCTATTTTGGGAGTCagtttaaaaatggaaatttgATTAAAAAGCAGCCTGAGAT
TGACATGCTGTATCTGCAAGGACTCTGTGGCAGTGCCCTTGCTGATAAAGAAAAGATACTGAAATGGTTT
TGGGAAAAGATCAGAGGAAAGTGGGACAAACCTAGGAAAGAAATGTTTGTCTCTGCACAGGGGTTACTCC
ATGAGCGCAGGCATGACATACTGAAGCAGCACACACTGACCATATTTGAAAACTTCGAAGATTGGTTTGA
GACCTCAGAGAAAGGTGCAAGCTTTATCAGTATTCTGAATGACATCATAAGGCTTCTTAAGACCTGGACA
TGGGGAACAACATACAATTTCCTCTACAAAATGGAAGaagaagatgtgCGCCTCTCTGTTCTGAGTGAGG
AGAAGACATATTATGAGGATGCTGGGTTGCTGAACAACTCTCCCTATTTCCCAATGctgagaaaggaaag
agacatCGATCTCATCATCTCCCTTGACTTCAGTGCTGGAGACCCAATGGagACTGTGATCAATACCTCT
GAAATGTGCAGGGACTTGAAGATTCCTTTCCCTGAGGTGAAATTGCCAAGAGGTGTTTCTGAGCCAGATG
ACTTCTATCTGTTTGAGGGGGAGAACGGAGCTCCCACTGTGATCCACATTCCTCTTTTCAACAAAGTCAA
CTGTGATGGACATATTAAGGAGTGGGTGAAAAGATATAGAACCTTCCAGAGCGCCTACAGCCATGACATG
ATGATCGATCTGATTGAGAAAGCAGGCGAgaatgtgaaaaacaacaagGCAAGACTAGTGACAGCGATTC
GGGAGGTCATTGAGAAAAAGAAGGCCAAAGCATGC
>C.chanos.GCA_018691325.1.JAFFIB010000115.1-1862497-1867291
TTACCCAAAATTGCACTGCTGGGCTCTGGGGGAGGTGAACGGGCCATGGTGGGATTACTGGGATCACTGG
TTCAGCTGGAGAAATGTGATCTACTGGACTGCATGCTGTATCTGAGTGGAGTCTCAGGATCCACATGGTG
CATGACTTCTTTATACAAAGAATCAAACTGGTCCCGTAAACTGGAGGCTATACAGAATGCCTTCATTAAA
AGGCTGGATAGAAGTGGCGTCAGCCTGTGGGGGCAGGTATTTAAGTTACTGGAATACTTCAAAAAGGACA
ACTGCTCGGCAACTGACACCTGGGCCGCAATCATTGTTTCGAACGTTATGAGAGAGATTAATGAAGACAC
CATAACAGATCACAGAGACAATCACGGCAAAGACCCATATCCTGTCTACACTGTGATCGACAAACAGTGC
AAATATGACAGACTATCCAAAGATGCCTGGTTTGAGATTACACCTCATGAAGCTGGATACTCCCTAACTG
GGGCTTTTGTGGACTCAGTCTGCTTTGGAAATCagtttgaaaatggaaaaataataaaagagcAGCCTGA
GATTGACATGCTGTATCTGCAAGGACTCTGTGGAAGTGCCCTAGCTGACATGGAGGAGATACTAAAG
TGGACCTGGGGAACAACATACAACTTCCTCCATGGTTTGACAGTGGAAGAGGTGCACCCCTCTGTTCtga
atgagaaggagagatatTATATAGACGCTGGACTGCTGATCAACTCTCCCTATTTCCCAAtgctgagaaa
ggagagagacatcGATCTCATCATCTCCCTTGACTTCAGTGCTGGAGACCCCATGGAGACAGTGCTAAGA
ACTGCTGATATGAGCAAGGCCTTGAAGATCCCTTTTCCAGAGGTGAAACTGCCAGATGATACTTCTGAGC
CAGACGACTTTTACGTTTTTGAAGGTCACACAAAAGCTCCAACTGTGATCCACATCCCCCTTTTCAACAA
AGTCAACTGTGATGGTCAAATAAAGGAGTGGGCAGAAAGATTTAAAACTTTTCAAGGGCCCTACACCCAT
GACATGACCATCGACCTGATTAAGAAAGCAGGCGAGAATGTGATAAACACCAAGGAAAAAATACTGAAAC
AGATTATGAAGCTTGTTGAGCAAAAGGAAGCCAAA
>C.chanos.GCA_018691325.1.JAFFIB010000115.1-1807781-1811600
ATGCCCAACATTGCACTGCTGggctctggaggaggagagcgaGCCATGCTGGGGTTACTGGGATCACTGG
TTCAGCTGAAGAAATGTGATCTGCTGGACTGCATGCTGTATCTGAGTGGAATATCAGGATCCACATGGTG
CATGGCTTCATTATACAGAGAGTCAGAATGGTCCACCAAACTGGATACGGTGCAGGATGCCATTGTTAAA
AGGCTTAATGGCCCTGAAGTCACCTGGACAGACTCTtggataaaactgaaaaaataccACAAGAGGGACA
acttctctctgacagactTCTGGGCTGTGATGTTTGTCTCCTCAATTGTTAAAGAGATTGACGAGGATGG
TGTAACGCATCAGAGGGATTATCACACCAAAGACCCGTATCCCGTCTATGCTGTGATTGATGACCAGTGC
AAACGGGAGAGTTTGGACTCTgATGCCTGGTTTGAGATTACACCTCATGAAGCTGGATACTCTCTCACTG
GGGCTTTTGTGGACTCAGCCTATTTAGGGAGTCAGTTTAAAAATGGAACCTCAGTTAAAAAGCAGCCTGA
GATTGACATGCTGTATCTGCAAGGGTTCTGTGGCAGTGCCCTTGCTgatgagacagagataaagaaatgG
CTATGGGAAGAGATGAAAAGAAAGCTGGACAAACATGGGAGAGCAATGTTTATCTCAGCACAGGAGTTAC
TCCATGAGCGCAAGAGTGCCATACTGGAGCAGTTCACACTAACCATCTGTGAAAACTTTGGGGATTGGTT
TAAGGTCCAGGACACTGGT
>C.chanos.GCA_018691325.1.JAFFIB010000077.1-462554-471452
ATTCCCAACATTGCTCTTCTGGGCTCTGGAGGAGGTGAGCGAGCCATGTTGGGGTTACTGGGATCACTGG
TTCAACTGCAGAAATGTGATCTACTGGACAGCATACTGTATCTGGGTGGAGTGTCAGGATCCACATGGTG
CATGGCGTCCTTATACAAAGAGCCAGACTGGTCCGCCCAACTGGATGCTGTGCAGGAGAAAATCATTAAA
AGGCTTGACGGTCCTGCAGTCAGCTGGTCGGATGCTTGGACAAAACTGGAAAAATACTATAACCAGAAGG
ACAATTTCTCCCTGACAGACGTCTGGGCTGCAATGTTTGTCACTACGATCAAGAAAGAGATTGATGAGAC
CACTGTAAGCAATCAGAGAGAGTATCACAGCAAAGACCCATACCCCATCTACACTGTGATTGACAAGCAG
TGCAAAAAGGACAAACTGAACCGAGATGCCTGGTTTGGGATTGTGCCTCATGAAGCTGGATACTCCCTCA
CTGGTGCTTTTGTGGACTCAACCAGTTTTGGGAGTCAGTTTAAAAATGGAAGTCTGATTAAAAAACAGCC
TGAGATGGACATGCTGTATCTGCAAGGACTCTGTGGCAGTGCCCTTGCTGATAGAGAAGAGATATCGAAA
TGGATTTGGGAACAAATCCCACGAAAGCTGGACAAACATGGGAGAGCGATGTTTATCTCTGCAGAGGAGG
TACTCTATGAGCCCAAGCATGACAGACTGAAGCAGTACACACTGAACATATGTGAAAACTTCGAAGATTG
GTTTGGGGACTCAGCAACAGgtgCAAGGGTCTCACTGATTATCCCAAAGATCATACAGCTTCTTTGCAAC
TGGACATGGGGAACAACATACAATTTCGTCCACAACATGAAAGTGCAAGATGTGAATCCCTCTGTtctga
gtgagaaggagagatattATGAAGATGCTGGGCTGCTGATCAACTCTCCCTATTTCCCAGTGCTGAGAAA
GGAAAGACACATTGATCTCATCATTTCCTTTGACTTCAGTGCTGGAGACCCCATGGAGacagtgaaaaaa
acctCTCAAATGTGCAAGGAGCTGAAGATCCCTTTCCCTGAAGTGCATGTCCCAGAGCCTGTTCAGGAAC
CAACAGACTTCTATTTGtttgaggggaaaaatgaaGCTCCAACTGTGATCCACATTCCCCTTTTCAACAA
AGTCAACTGTGGTGGTCATATCAAGGAGTGGGAGAGCAGATATAAGACCTTTCAGGGGCCCTACAGCCAT
GACATGATCATTGATCTGATTGAGAAAGCAGGCGAGAATGTGATAAACAACAAGGAAAAACTAGTGAAAA
TAATTCAAGCCATCATTGAGAAAAAGAAGGCCAAATCA




transeq Neoarius_spp.HOG0019032.fa Neoarius_spp.HOG0019032.prot ; sed -i 's/_1$//g' Neoarius_spp.HOG0019032.prot
transeq Chanos.alternative.fa Chanos.alternative.prot ; sed -i 's/_1$//g' Chanos.alternative.prot



grep "Chanos_chanos\|Neoarius" Coding_sequences_alignments/N5.HOG0019032.prot | sed 's/>//g' > curr_seq.id
xargs samtools faidx Coding_sequences_alignments/N5.HOG0019032.prot < curr_seq.id > curr_seq.prot 

cat Neoarius_spp.HOG0019032.prot Chanos.alternative.prot curr_seq.prot  > Neoarius_spp.ChanosALT.HOG0019032.combined.prot

muscle5.1 -align Neoarius_spp.ChanosALT.HOG0019032.combined.prot -output Neoarius_spp.ChanosALT.HOG0019032.combined.prot.aln
trimal -in Neoarius_spp.ChanosALT.HOG0019032.combined.prot.aln -gt 0.8 -out Neoarius_spp.ChanosALT.HOG0019032.combined.prot.aln.trimmed


grep ">" Neoarius_spp.ChanosALT.HOG0019032.combined.prot.aln | grep -v "GCA_025727735.1.JAODTU010000006.1-15288793-15503138" | sed 's/>//g' > seq_to_retain
xargs samtools faidx Neoarius_spp.ChanosALT.HOG0019032.combined.prot.aln.trimmed < seq_to_retain > temp.prot ; mv temp.prot Neoarius_spp.ChanosALT.HOG0019032.combined.prot.aln.trimmed
xargs samtools faidx Neoarius_spp.ChanosALT.HOG0019032.combined.prot.aln < seq_to_retain > temp.prot ; mv temp.prot Neoarius_spp.ChanosALT.HOG0019032.combined.prot.aln



iqtree -s Neoarius_spp.ChanosALT.HOG0019032.combined.prot.aln -st AA -nt 8 -bb 1000 --redo -m LG+F+G4


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
sbatch --qos=1day -c 4 --mem=10G -e error.fitMG4.out -o slurm.fitMG4.out launch_fitMG4.sh N5.HOG0019032


Recipient branches : 
Node17
Node19
Neoarius_graeffei_rna_XM_060938829_1
Neoarius_graeffei_rna_XM_060938827_1
Neoarius_graeffei_rna_XM_060939486_1


#Test positive selection and relaxed selection on receiver branch

sbatch --qos=1week -c 4 --mem=10G -e error.absrel.cand.out -o slurm.absrel.cand.out --job-name=HOG0019032 launch_absrel_cand.sh N5.HOG0019032


#No positive selection detected, no run of RELAX



