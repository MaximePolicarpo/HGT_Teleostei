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


samtools faidx Proteomes_BUSCO80/Hypomesus_transpacificus.fa Hypomesus_transpacificus---rna-XM_047023826.1 > Hypomesus_transpacificus---rna-XM_047023826.1.prot

#blastp against the fish proteome (50 best hits) + Uniprot database without any fish species (10 best hits)

blastp -query Hypomesus_transpacificus---rna-XM_047023826.1.prot -db concatenated_proteomes.fa -outfmt 6 -out Gene_vs_FishProteome.blastp -num_threads 8 -max_target_seqs 51
blastp -query Hypomesus_transpacificus---rna-XM_047023826.1.prot -db non_actino_uniprot.fa -outfmt 6 -out Gene_vs_Uniprot.blastp -num_threads 8  -max_target_seqs 10 -evalue 1e-5


cut -f2 Gene_vs_FishProteome.blastp | sort | uniq  > closest_fish_seq.id
cut -f2 Gene_vs_Uniprot.blastp | sort | uniq > closest_nonfish_seq.id 

xargs samtools faidx concatenated_proteomes.fa < closest_fish_seq.id > closest_fish_seq.fa
xargs samtools faidx non_actino_uniprot.fa  < closest_nonfish_seq.id > closest_nonfish_seq.fa

#Align with muscle and trim with trimal

cat closest_fish_seq.fa closest_nonfish_seq.fa > Uniprot_plus_closefish.fa

muscle5.1 -align closest_fish_seq.fa -output closest_fish_seq.aln
mafft --add closest_nonfish_seq.fa  --keeplength closest_fish_seq.aln > Uniprot_plus_closefish.aln
trimal -in Uniprot_plus_closefish.aln -automated1 -out Uniprot_plus_closefish.aln.trimmed

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


#N5.HOG0018901.clade_hgt.txt=> clade where the HGT is extracted from iTOL (from N5.HOG0018901.prot.aln.treefile)

cut -f1,2 -d "_" N5.HOG0018901.clade_hgt.txt | sort | uniq  > species_to_draw.clade1
#Choose close species without the gene .. 


echo "Engraulis_encrasicolus" > non_transfer_species_1
echo "Alosa_sapidissima" >> non_transfer_species_1
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

curr_OGG=N5.HOG0018901

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



#First add Sardina_pilchardus

curr_OGG=N5.HOG0018901
curr_sp=Sardina_pilchardus
ref_sp=Clupea_harengus


grep -A10 -B10 "N5_HOG0018901" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "NC_085015_1" > Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0018901" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "NC_084997_1" >> Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0040087" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0046053" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0016660" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0026571" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret  >> Syn_tables_dir/$curr_sp.synt.df



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




#Now add Alosa_sapidissima

curr_OGG=N5.HOG0018901
curr_sp=Alosa_sapidissima
ref_sp=Sardina_pilchardus

rm Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0018901" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "NC_055980_1" > Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0018901" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "NC_055973_1" >> Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0040087" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0046053" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0016660" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0026571" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0041470" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0024447" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "NC_055962_1" >> Syn_tables_dir/$curr_sp.synt.df


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

curr_OGG=N5.HOG0018901
curr_sp=Engraulis_encrasicolus
ref_sp=Clupea_harengus

rm Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0018901" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df
grep -A15 -B15 "N5_HOG0046053" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0016660" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0026571" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret  >> Syn_tables_dir/$curr_sp.synt.df



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





#Now add Hypomesus_transpacificus

curr_OGG=N5.HOG0018901
curr_sp=Osmerus_eperlanus
ref_sp=Hypomesus_transpacificus

grep -A10 -B10 "N5_HOG0047849" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0050182" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "NC_085039_1" >> Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0000034" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "NC_085020_1" >> Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0020354" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df


awk '{ if (NF < 5 || $5 == "") $5 = "no_OGG"; print $0;}' Syn_tables_dir/$curr_sp.synt.df > Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/	/,/g' Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/ no_OGG/,no_OGG/g' Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/ /,/g' Syn_tables_dir/$curr_sp.synt.final.df
grep -v ",,,," Syn_tables_dir/$curr_sp.synt.final.df > temp ; mv temp Syn_tables_dir/$curr_sp.synt.final.df


rm Syn_tables_dir/$curr_sp

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

curr_OGG=N5.HOG0018901
curr_sp=Borostomias_antarcticus
ref_sp=Hypomesus_transpacificus

grep -A10 -B10 "N5_HOG0047849" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0028263" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "OX465208_1"  >> Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0035654" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret  >> Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0020354" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df




awk '{ if (NF < 5 || $5 == "") $5 = "no_OGG"; print $0;}' Syn_tables_dir/$curr_sp.synt.df > Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/	/,/g' Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/ no_OGG/,no_OGG/g' Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/ /,/g' Syn_tables_dir/$curr_sp.synt.final.df
grep -v ",,,," Syn_tables_dir/$curr_sp.synt.final.df > temp ; mv temp Syn_tables_dir/$curr_sp.synt.final.df


rm Syn_tables_dir/$curr_sp

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




#Now add Gadus morhua

curr_OGG=N5.HOG0018901
curr_sp=Gadus_morhua
ref_sp=Hypomesus_transpacificus

grep -A10 -B10 "N5_HOG0047849" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0028263" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0035654" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret  >> Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0020354" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df



awk '{ if (NF < 5 || $5 == "") $5 = "no_OGG"; print $0;}' Syn_tables_dir/$curr_sp.synt.df > Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/	/,/g' Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/ no_OGG/,no_OGG/g' Syn_tables_dir/$curr_sp.synt.final.df
sed -i 's/ /,/g' Syn_tables_dir/$curr_sp.synt.final.df
grep -v ",,,," Syn_tables_dir/$curr_sp.synt.final.df > temp ; mv temp Syn_tables_dir/$curr_sp.synt.final.df


rm Syn_tables_dir/$curr_sp

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
##========================================== Now let's find TE in the region  ===================================================
##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================


#Hypomesus_transpacificus
grep -A1 -B1 "HOG0018901"  GFF3_N5_OGGs/Hypomesus_transpacificus.gff.simplified.sorted.OGG.tiret
samtools faidx GCF_021917145.1_fHypTra1_genomic.fna NC_061067.1:9519063-9530253 > N5.HOG0018901.Hypomesus_transpacificus.extended.1.fa
samtools faidx GCF_021917145.1_fHypTra1_genomic.fna NC_061068.1:15332589-15352913 > N5.HOG0018901.Hypomesus_transpacificus.extended.2.fa
samtools faidx GCF_021917145.1_fHypTra1_genomic.fna NC_061080.1:2837156-2852140 > N5.HOG0018901.Hypomesus_transpacificus.extended.3.fa
samtools faidx GCF_021917145.1_fHypTra1_genomic.fna NC_061080.1:3421351-3435168 > N5.HOG0018901.Hypomesus_transpacificus.extended.4.fa

sed -i 's/:/-/g' N5.HOG0018901.Hypomesus_transpacificus.extended.*
makeblastdb -in N5.HOG0018901.Hypomesus_transpacificus.extended.1.fa -dbtype nucl
makeblastdb -in N5.HOG0018901.Hypomesus_transpacificus.extended.2.fa -dbtype nucl
makeblastdb -in N5.HOG0018901.Hypomesus_transpacificus.extended.3.fa -dbtype nucl
makeblastdb -in N5.HOG0018901.Hypomesus_transpacificus.extended.4.fa -dbtype nucl


#Clupea_harengus
grep -A1 -B1 "HOG0018901"  GFF3_N5_OGGs/Clupea_harengus.gff.simplified.sorted.OGG.tiret
samtools faidx GCF_900700415.2_Ch_v2.0.2_genomic.fna NC_045160.1:5483282-5495013 > N5.HOG0018901.Clupea_harengus.extended.1.fa
samtools faidx GCF_900700415.2_Ch_v2.0.2_genomic.fna NC_045166.1:15821580-15885645 > N5.HOG0018901.Clupea_harengus.extended.2.fa
samtools faidx GCF_900700415.2_Ch_v2.0.2_genomic.fna NC_045177.1:8616226-8697916 > N5.HOG0018901.Clupea_harengus.extended.3.fa

sed -i 's/:/-/g' N5.HOG0018901.Clupea_harengus.extended.*
makeblastdb -in N5.HOG0018901.Clupea_harengus.extended.1.fa -dbtype nucl
makeblastdb -in N5.HOG0018901.Clupea_harengus.extended.2.fa -dbtype nucl
makeblastdb -in N5.HOG0018901.Clupea_harengus.extended.3.fa -dbtype nucl



#Now lets find TEs using the repbase and dfam databases

tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0018901.Hypomesus_transpacificus.extended.1.fa -evalue 1e-5 -outfmt 6 -out TE.Hypomesus_transpacificus.1.tblastn -num_threads 8
sed -i 's/#//g' TE.Hypomesus_transpacificus.1.tblastn

tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0018901.Hypomesus_transpacificus.extended.2.fa -evalue 1e-5 -outfmt 6 -out TE.Hypomesus_transpacificus.2.tblastn -num_threads 8
sed -i 's/#//g' TE.Hypomesus_transpacificus.2.tblastn

tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0018901.Hypomesus_transpacificus.extended.3.fa -evalue 1e-5 -outfmt 6 -out TE.Hypomesus_transpacificus.3.tblastn -num_threads 8
sed -i 's/#//g' TE.Hypomesus_transpacificus.3.tblastn

tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0018901.Hypomesus_transpacificus.extended.4.fa -evalue 1e-5 -outfmt 6 -out TE.Hypomesus_transpacificus.4.tblastn -num_threads 8
sed -i 's/#//g' TE.Hypomesus_transpacificus.4.tblastn


tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0018901.Clupea_harengus.extended.1.fa -evalue 1e-5 -outfmt 6 -out TE.Clupea_harengus.1.tblastn -num_threads 8
sed -i 's/#//g' TE.Clupea_harengus.1.tblastn

tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0018901.Clupea_harengus.extended.2.fa -evalue 1e-5 -outfmt 6 -out TE.Clupea_harengus.2.tblastn -num_threads 8
sed -i 's/#//g' TE.Clupea_harengus.2.tblastn

tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0018901.Clupea_harengus.extended.3.fa -evalue 1e-5 -outfmt 6 -out TE.Clupea_harengus.3.tblastn -num_threads 8
sed -i 's/#//g' TE.Clupea_harengus.3.tblastn


#merge tblastn hits and find the best TE match by doing a blastx

Rscript Rscript_merge_blast_hits.R TE.Hypomesus_transpacificus.1.tblastn TE.Hypomesus_transpacificus.1.tblastn.merged
Rscript Rscript_merge_blast_hits.R TE.Hypomesus_transpacificus.2.tblastn TE.Hypomesus_transpacificus.2.tblastn.merged
Rscript Rscript_merge_blast_hits.R TE.Hypomesus_transpacificus.3.tblastn TE.Hypomesus_transpacificus.3.tblastn.merged
Rscript Rscript_merge_blast_hits.R TE.Hypomesus_transpacificus.4.tblastn TE.Hypomesus_transpacificus.4.tblastn.merged

Rscript Rscript_merge_blast_hits.R TE.Clupea_harengus.1.tblastn TE.Clupea_harengus.1.tblastn.merged
Rscript Rscript_merge_blast_hits.R TE.Clupea_harengus.2.tblastn TE.Clupea_harengus.2.tblastn.merged
Rscript Rscript_merge_blast_hits.R TE.Clupea_harengus.3.tblastn TE.Clupea_harengus.3.tblastn.merged


xargs samtools faidx N5.HOG0018901.Hypomesus_transpacificus.extended.1.fa < TE.Hypomesus_transpacificus.1.tblastn.merged > TE.Hypomesus_transpacificus.1.BEST.fa
blastx -query TE.Hypomesus_transpacificus.1.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Hypomesus_transpacificus.1.BEST.blastx -max_target_seqs 1
cut -f1  TE.Hypomesus_transpacificus.1.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Hypomesus_transpacificus.1.BEST.blastx >> temp  ; done ; mv temp TE.Hypomesus_transpacificus.1.BEST.blastx

xargs samtools faidx N5.HOG0018901.Hypomesus_transpacificus.extended.2.fa < TE.Hypomesus_transpacificus.2.tblastn.merged > TE.Hypomesus_transpacificus.2.BEST.fa
blastx -query TE.Hypomesus_transpacificus.2.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Hypomesus_transpacificus.2.BEST.blastx -max_target_seqs 1
cut -f1  TE.Hypomesus_transpacificus.2.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Hypomesus_transpacificus.2.BEST.blastx >> temp  ; done ; mv temp TE.Hypomesus_transpacificus.2.BEST.blastx

xargs samtools faidx N5.HOG0018901.Hypomesus_transpacificus.extended.3.fa < TE.Hypomesus_transpacificus.3.tblastn.merged > TE.Hypomesus_transpacificus.3.BEST.fa
blastx -query TE.Hypomesus_transpacificus.3.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Hypomesus_transpacificus.3.BEST.blastx -max_target_seqs 1
cut -f1  TE.Hypomesus_transpacificus.3.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Hypomesus_transpacificus.3.BEST.blastx >> temp  ; done ; mv temp TE.Hypomesus_transpacificus.3.BEST.blastx


xargs samtools faidx N5.HOG0018901.Hypomesus_transpacificus.extended.4.fa < TE.Hypomesus_transpacificus.4.tblastn.merged > TE.Hypomesus_transpacificus.4.BEST.fa
blastx -query TE.Hypomesus_transpacificus.4.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Hypomesus_transpacificus.4.BEST.blastx -max_target_seqs 1
cut -f1  TE.Hypomesus_transpacificus.4.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Hypomesus_transpacificus.4.BEST.blastx >> temp  ; done ; mv temp TE.Hypomesus_transpacificus.4.BEST.blastx

xargs samtools faidx N5.HOG0018901.Clupea_harengus.extended.1.fa < TE.Clupea_harengus.1.tblastn.merged > TE.Clupea_harengus.1.BEST.fa
blastx -query TE.Clupea_harengus.1.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Clupea_harengus.1.BEST.blastx -max_target_seqs 1
cut -f1  TE.Clupea_harengus.1.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Clupea_harengus.1.BEST.blastx >> temp  ; done ; mv temp TE.Clupea_harengus.1.BEST.blastx

xargs samtools faidx N5.HOG0018901.Clupea_harengus.extended.2.fa < TE.Clupea_harengus.2.tblastn.merged > TE.Clupea_harengus.2.BEST.fa
blastx -query TE.Clupea_harengus.2.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Clupea_harengus.2.BEST.blastx -max_target_seqs 1
cut -f1  TE.Clupea_harengus.2.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Clupea_harengus.2.BEST.blastx >> temp  ; done ; mv temp TE.Clupea_harengus.2.BEST.blastx

xargs samtools faidx N5.HOG0018901.Clupea_harengus.extended.3.fa < TE.Clupea_harengus.3.tblastn.merged > TE.Clupea_harengus.3.BEST.fa
blastx -query TE.Clupea_harengus.3.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Clupea_harengus.3.BEST.blastx -max_target_seqs 1
cut -f1  TE.Clupea_harengus.3.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Clupea_harengus.3.BEST.blastx >> temp  ; done ; mv temp TE.Clupea_harengus.3.BEST.blastx




#Now find shared elements

cut -f2 TE.Hypomesus_transpacificus.1.BEST.blastx | sort | uniq > TE.Hypomesus_transpacificus.1.uniqTE
cut -f2 TE.Hypomesus_transpacificus.2.BEST.blastx | sort | uniq > TE.Hypomesus_transpacificus.2.uniqTE
cut -f2 TE.Hypomesus_transpacificus.3.BEST.blastx | sort | uniq > TE.Hypomesus_transpacificus.3.uniqTE
cut -f2 TE.Hypomesus_transpacificus.4.BEST.blastx | sort | uniq > TE.Hypomesus_transpacificus.4.uniqTE

cut -f2 TE.Clupea_harengus.1.BEST.blastx | sort | uniq > TE.Clupea_harengus.1.uniqTE
cut -f2 TE.Clupea_harengus.2.BEST.blastx | sort | uniq > TE.Clupea_harengus.2.uniqTE
cut -f2 TE.Clupea_harengus.3.BEST.blastx | sort | uniq > TE.Clupea_harengus.3.uniqTE

comm -12 TE.Hypomesus_transpacificus.1.uniqTE TE.Clupea_harengus.1.uniqTE
comm -12 TE.Hypomesus_transpacificus.2.uniqTE TE.Clupea_harengus.1.uniqTE
comm -12 TE.Hypomesus_transpacificus.3.uniqTE TE.Clupea_harengus.1.uniqTE
comm -12 TE.Hypomesus_transpacificus.4.uniqTE TE.Clupea_harengus.1.uniqTE

comm -12 TE.Hypomesus_transpacificus.1.uniqTE TE.Clupea_harengus.2.uniqTE
comm -12 TE.Hypomesus_transpacificus.2.uniqTE TE.Clupea_harengus.2.uniqTE
comm -12 TE.Hypomesus_transpacificus.3.uniqTE TE.Clupea_harengus.2.uniqTE
comm -12 TE.Hypomesus_transpacificus.4.uniqTE TE.Clupea_harengus.2.uniqTE

comm -12 TE.Hypomesus_transpacificus.1.uniqTE TE.Clupea_harengus.3.uniqTE
comm -12 TE.Hypomesus_transpacificus.2.uniqTE TE.Clupea_harengus.3.uniqTE
comm -12 TE.Hypomesus_transpacificus.3.uniqTE TE.Clupea_harengus.3.uniqTE
comm -12 TE.Hypomesus_transpacificus.4.uniqTE TE.Clupea_harengus.3.uniqTE



###ONE SHARED ELEMENT = L1-19_SSa_2p:ClassI:LINE:L1:L1


#### Now let's find every copy of this element in the genome of both species + closely related species


samtools faidx Dfam_plus_Repbase.cdhit80.prot L1-19_SSa_2p:ClassI:LINE:L1:L1 > L1-19.prot



#Hypomesus_transpacificus
tblastn -query L1-19.prot -db GCF_021917145.1_fHypTra1_genomic.fna -evalue 1e-1 -outfmt 6 -out L1-19.Hypomesus_transpacificus.tblastn -num_threads 8
sed -i 's/#//g' L1-19.Hypomesus_transpacificus.tblastn
Rscript Rscript_merge_blast_hits.R L1-19.Hypomesus_transpacificus.tblastn L1-19.Hypomesus_transpacificus.tblastn.merged
xargs samtools faidx GCF_021917145.1_fHypTra1_genomic.fna < L1-19.Hypomesus_transpacificus.tblastn.merged > L1-19.Hypomesus_transpacificus.BEST.fa
diamond blastx --query L1-19.Hypomesus_transpacificus.BEST.fa --db Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out L1-19.Hypomesus_transpacificus.BEST.blastx --max-target-seqs 1
echo "NC_061067.1:9528226-9528625	L1-19_SSa_2p:ClassI:LINE:L1:L1	84.5	148	23	0	544	101	244	391	1.6e-5	265" >> L1-19.Hypomesus_transpacificus.BEST.blastx
grep "L1-19" L1-19.Hypomesus_transpacificus.BEST.blastx > L1-19.Hypomesus_transpacificus.list

#Clupea_harengus
tblastn -query L1-19.prot -db GCF_900700415.2_Ch_v2.0.2_genomic.fna -evalue 1e-1 -outfmt 6 -out L1-19.Clupea_harengus.tblastn -num_threads 8
sed -i 's/#//g' L1-19.Clupea_harengus.tblastn
Rscript Rscript_merge_blast_hits.R L1-19.Clupea_harengus.tblastn L1-19.Clupea_harengus.tblastn.merged
xargs samtools faidx GCF_900700415.2_Ch_v2.0.2_genomic.fna  < L1-19.Clupea_harengus.tblastn.merged > L1-19.Clupea_harengus.BEST.fa
diamond blastx --query L1-19.Clupea_harengus.BEST.fa --db Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out L1-19.Clupea_harengus.BEST.blastx --max-target-seqs 1
grep "L1-19" L1-19.Clupea_harengus.BEST.blastx > L1-19.Clupea_harengus.list


#Sardina_pilchardus
tblastn -query L1-19.prot -db GCF_963854185.1_fSarPil1.1_genomic.fna -evalue 1e-1 -outfmt 6 -out L1-19.Sardina_pilchardus.tblastn -num_threads 8
sed -i 's/#//g' L1-19.Sardina_pilchardus.tblastn
Rscript Rscript_merge_blast_hits.R L1-19.Sardina_pilchardus.tblastn L1-19.Sardina_pilchardus.tblastn.merged
xargs samtools faidx GCF_963854185.1_fSarPil1.1_genomic.fna  < L1-19.Sardina_pilchardus.tblastn.merged > L1-19.Sardina_pilchardus.BEST.fa
diamond blastx --query L1-19.Sardina_pilchardus.BEST.fa --db Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out L1-19.Sardina_pilchardus.BEST.blastx --max-target-seqs 1
grep "L1-19" L1-19.Sardina_pilchardus.BEST.blastx > L1-19.Sardina_pilchardus.list

#Alosa_sapidissima
tblastn -query L1-19.prot -db GCF_018492685.1_fAloSap1.pri_genomic.fna -evalue 1e-1 -outfmt 6 -out L1-19.Alosa_sapidissima.tblastn -num_threads 8
sed -i 's/#//g' L1-19.Alosa_sapidissima.tblastn
Rscript Rscript_merge_blast_hits.R L1-19.Alosa_sapidissima.tblastn L1-19.Alosa_sapidissima.tblastn.merged
xargs samtools faidx GCF_018492685.1_fAloSap1.pri_genomic.fna  < L1-19.Alosa_sapidissima.tblastn.merged > L1-19.Alosa_sapidissima.BEST.fa
diamond blastx --query L1-19.Alosa_sapidissima.BEST.fa --db Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out L1-19.Alosa_sapidissima.BEST.blastx --max-target-seqs 1
grep "L1-19" L1-19.Alosa_sapidissima.BEST.blastx > L1-19.Alosa_sapidissima.list


#Engraulis_encrasicolus
tblastn -query L1-19.prot -db GCF_034702125.1_IST_EnEncr_1.0_genomic.fna -evalue 1e-1 -outfmt 6 -out L1-19.Engraulis_encrasicolus.tblastn -num_threads 8
sed -i 's/#//g' L1-19.Engraulis_encrasicolus.tblastn
Rscript Rscript_merge_blast_hits.R L1-19.Engraulis_encrasicolus.tblastn L1-19.Engraulis_encrasicolus.tblastn.merged
xargs samtools faidx GCF_034702125.1_IST_EnEncr_1.0_genomic.fna  < L1-19.Engraulis_encrasicolus.tblastn.merged > L1-19.Engraulis_encrasicolus.BEST.fa
diamond blastx --query L1-19.Engraulis_encrasicolus.BEST.fa --db Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out L1-19.Engraulis_encrasicolus.BEST.blastx --max-target-seqs 1
grep "L1-19" L1-19.Engraulis_encrasicolus.BEST.blastx > L1-19.Engraulis_encrasicolus.list


#Osmerus_eperlanus
tblastn -query L1-19.prot -db GCF_963692335.1_fOsmEpe2.1_genomic.fna -evalue 1e-1 -outfmt 6 -out L1-19.Osmerus_eperlanus.tblastn -num_threads 8
sed -i 's/#//g' L1-19.Osmerus_eperlanus.tblastn
Rscript Rscript_merge_blast_hits.R L1-19.Osmerus_eperlanus.tblastn L1-19.Osmerus_eperlanus.tblastn.merged
xargs samtools faidx GCF_963692335.1_fOsmEpe2.1_genomic.fna  < L1-19.Osmerus_eperlanus.tblastn.merged > L1-19.Osmerus_eperlanus.BEST.fa
diamond blastx --query L1-19.Osmerus_eperlanus.BEST.fa --db Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out L1-19.Osmerus_eperlanus.BEST.blastx --max-target-seqs 1
grep "L1-19" L1-19.Osmerus_eperlanus.BEST.blastx > L1-19.Osmerus_eperlanus.list

#Borostomias_antarcticus
tblastn -query L1-19.prot -db GCA_949987555.1_fBorAnt1.1_genomic.fna -evalue 1e-1 -outfmt 6 -out L1-19.Borostomias_antarcticus.tblastn -num_threads 8
sed -i 's/#//g' L1-19.Borostomias_antarcticus.tblastn
Rscript Rscript_merge_blast_hits.R L1-19.Borostomias_antarcticus.tblastn L1-19.Borostomias_antarcticus.tblastn.merged
xargs samtools faidx GCA_949987555.1_fBorAnt1.1_genomic.fna  < L1-19.Borostomias_antarcticus.tblastn.merged > L1-19.Borostomias_antarcticus.BEST.fa
diamond blastx --query L1-19.Borostomias_antarcticus.BEST.fa --db Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out L1-19.Borostomias_antarcticus.BEST.blastx --max-target-seqs 1
grep "L1-19" L1-19.Borostomias_antarcticus.BEST.blastx > L1-19.Borostomias_antarcticus.list

#Gadus_morhua
tblastn -query L1-19.prot -db GCF_902167405.1_gadMor3.0_genomic.fna -evalue 1e-1 -outfmt 6 -out L1-19.Gadus_morhua.tblastn -num_threads 8
sed -i 's/#//g' L1-19.Gadus_morhua.tblastn
Rscript Rscript_merge_blast_hits.R L1-19.Gadus_morhua.tblastn L1-19.Gadus_morhua.tblastn.merged
xargs samtools faidx GCF_902167405.1_gadMor3.0_genomic.fna  < L1-19.Gadus_morhua.tblastn.merged > L1-19.Gadus_morhua.BEST.fa
diamond blastx --query L1-19.Gadus_morhua.BEST.fa --db Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out L1-19.Gadus_morhua.BEST.blastx --max-target-seqs 1
grep "L1-19" L1-19.Gadus_morhua.BEST.blastx > L1-19.Gadus_morhua.list


wc -l L1-19.Engraulis_encrasicolus.list
wc -l L1-19.Alosa_sapidissima.list
wc -l L1-19.Sardina_pilchardus.list
wc -l L1-19.Clupea_harengus.list
wc -l L1-19.Hypomesus_transpacificus.list
wc -l L1-19.Osmerus_eperlanus.list 
wc -l L1-19.Borostomias_antarcticus.list
wc -l L1-19.Gadus_morhua.list


### Test the similarity between the TE found in both Clupea and Hypomesus


samtools faidx GCF_021917145.1_fHypTra1_genomic.fna NC_061067.1:9528226-9528625 > Hypomesus.L1-19.fa
revseq Hypomesus.L1-19.fa Hypomesus.L1-19.fa.rev 
header=`grep ">" Hypomesus.L1-19.fa | sed 's/>//g'` ; sed -i "s/>.*/>$header/g" Hypomesus.L1-19.fa.rev  

samtools faidx GCF_900700415.2_Ch_v2.0.2_genomic.fna NC_045166.1:15874321-15874792 > Clupea.L1-19.fa

sed -i 's/>.*/>Htranspacificus_L1-19_HGT/g' Hypomesus.L1-19.fa.rev  
sed -i 's/>.*/>Charengus_L1-19_HGT/g' Clupea.L1-19.fa

cat Hypomesus.L1-19.fa Clupea.L1-19.fa > L1-19_combined.fa

needle -asequence Hypomesus.L1-19.fa.rev -bsequence Clupea.L1-19.fa -outfile L1-19_combined.aln -gapopen 10.0 -gapextend 0.5



### In a FINAL step we will make another synteny plot, only between Clupea and Hypomesus


## Clupea harengus ## FIRST SCAFFOLD

rm seq_clustered_infos_ogg.TE.txt ; rm clusters_ID_TE.txt


scaffold=`grep ">" N5.HOG0018901.Clupea_harengus.extended.1.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0018901.Clupea_harengus.extended.1.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0018901.Clupea_harengus.extended.1.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Clupea_harengus,$scaffold,$length" > clusters_ID_TE.txt

grep "N5_HOG0018901"  GFF3_N5_OGGs/Clupea_harengus.gff.simplified.sorted.OGG.tiret

grep "XM_042708796" GFF3_files_per_species/Clupea_harengus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons


nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene1.$nbexon,N5_HOG0018901,+,Clupea_harengus" >> seq_clustered_infos_ogg.TE.txt
done


cut -f1 TE.Clupea_harengus.1.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Clupea_harengus.1.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0018901.Clupea_harengus.extended.1.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0018901.Clupea_harengus.extended.1.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold,/g" | sed "s/$/,$strand,Clupea_harengus/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done


## Clupea harengus ## SECOND SCAFFOLD


scaffold=`grep ">" N5.HOG0018901.Clupea_harengus.extended.2.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0018901.Clupea_harengus.extended.2.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0018901.Clupea_harengus.extended.2.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Clupea_harengus,$scaffold,$length" >> clusters_ID_TE.txt

grep "N5_HOG0018901"  GFF3_N5_OGGs/Clupea_harengus.gff.simplified.sorted.OGG.tiret

grep "XM_042709983" GFF3_files_per_species/Clupea_harengus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons
grep "XM_012826723" GFF3_files_per_species/Clupea_harengus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_2.exons
grep "XM_042709984" GFF3_files_per_species/Clupea_harengus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_3.exons


nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene2.$nbexon,N5_HOG0018901,-,Clupea_harengus" >> seq_clustered_infos_ogg.TE.txt
done
nbexon=0
for line in `cat gene_2.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene3.$nbexon,N5_HOG0018901,-,Clupea_harengus" >> seq_clustered_infos_ogg.TE.txt
done
nbexon=0
for line in `cat gene_3.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene4.$nbexon,N5_HOG0018901,-,Clupea_harengus" >> seq_clustered_infos_ogg.TE.txt
done

cut -f1 TE.Clupea_harengus.2.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Clupea_harengus.2.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0018901.Clupea_harengus.extended.2.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0018901.Clupea_harengus.extended.2.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold,/g" | sed "s/$/,$strand,Clupea_harengus/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done





## Clupea harengus ## THIRD SCAFFOLD


scaffold=`grep ">" N5.HOG0018901.Clupea_harengus.extended.3.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0018901.Clupea_harengus.extended.3.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0018901.Clupea_harengus.extended.3.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Clupea_harengus,$scaffold,$length" >> clusters_ID_TE.txt

grep "N5_HOG0018901"  GFF3_N5_OGGs/Clupea_harengus.gff.simplified.sorted.OGG.tiret

grep "XM_042703844" GFF3_files_per_species/Clupea_harengus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons
grep "XM_042703846" GFF3_files_per_species/Clupea_harengus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_2.exons


nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene5.$nbexon,N5_HOG0018901,-,Clupea_harengus" >> seq_clustered_infos_ogg.TE.txt
done
nbexon=0
for line in `cat gene_2.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene6.$nbexon,N5_HOG0018901,-,Clupea_harengus" >> seq_clustered_infos_ogg.TE.txt
done


cut -f1 TE.Clupea_harengus.3.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Clupea_harengus.3.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0018901.Clupea_harengus.extended.3.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0018901.Clupea_harengus.extended.3.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold,/g" | sed "s/$/,$strand,Clupea_harengus/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done




## Hypomesus_transpacificus  ### FIRST SCAFFOLD


grep ">" N5.HOG0018901.Hypomesus_transpacificus.extended.1.fa
scaffold=`grep ">" N5.HOG0018901.Hypomesus_transpacificus.extended.1.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0018901.Hypomesus_transpacificus.extended.1.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0018901.Hypomesus_transpacificus.extended.1.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Hypomesus_transpacificus,$scaffold,$length" >> clusters_ID_TE.txt


grep "N5_HOG0018901"  GFF3_N5_OGGs/Hypomesus_transpacificus.gff.simplified.sorted.OGG.tiret
grep "XM_047023826" GFF3_files_per_species/Hypomesus_transpacificus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons


nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene1.exon.$nbexon,N5_HOG0018901,+,Hypomesus_transpacificus" >> seq_clustered_infos_ogg.TE.txt
done


cut -f1 TE.Hypomesus_transpacificus.1.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Hypomesus_transpacificus.1.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0018901.Hypomesus_transpacificus.extended.1.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0018901.Hypomesus_transpacificus.extended.1.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold,/g" | sed "s/$/,$strand,Hypomesus_transpacificus/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done



## Hypomesus_transpacificus  ### SECOND SCAFFOLD


grep ">" N5.HOG0018901.Hypomesus_transpacificus.extended.2.fa
scaffold=`grep ">" N5.HOG0018901.Hypomesus_transpacificus.extended.2.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0018901.Hypomesus_transpacificus.extended.2.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0018901.Hypomesus_transpacificus.extended.2.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Hypomesus_transpacificus,$scaffold,$length" >> clusters_ID_TE.txt


grep "N5_HOG0018901"  GFF3_N5_OGGs/Hypomesus_transpacificus.gff.simplified.sorted.OGG.tiret
grep "XM_047025851" GFF3_files_per_species/Hypomesus_transpacificus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons


nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene2.exon.$nbexon,N5_HOG0018901,-,Hypomesus_transpacificus" >> seq_clustered_infos_ogg.TE.txt
done


cut -f1 TE.Hypomesus_transpacificus.2.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Hypomesus_transpacificus.2.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0018901.Hypomesus_transpacificus.extended.2.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0018901.Hypomesus_transpacificus.extended.2.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold,/g" | sed "s/$/,$strand,Hypomesus_transpacificus/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done


## Hypomesus_transpacificus  ### THIRD SCAFFOLD


grep ">" N5.HOG0018901.Hypomesus_transpacificus.extended.3.fa
scaffold=`grep ">" N5.HOG0018901.Hypomesus_transpacificus.extended.3.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0018901.Hypomesus_transpacificus.extended.3.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0018901.Hypomesus_transpacificus.extended.3.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Hypomesus_transpacificus,$scaffold,$length" >> clusters_ID_TE.txt


grep "N5_HOG0018901"  GFF3_N5_OGGs/Hypomesus_transpacificus.gff.simplified.sorted.OGG.tiret
grep "XM_047043933" GFF3_files_per_species/Hypomesus_transpacificus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons


nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene3.exon.$nbexon,N5_HOG0018901,-,Hypomesus_transpacificus" >> seq_clustered_infos_ogg.TE.txt
done


cut -f1 TE.Hypomesus_transpacificus.3.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Hypomesus_transpacificus.3.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0018901.Hypomesus_transpacificus.extended.3.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0018901.Hypomesus_transpacificus.extended.3.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold,/g" | sed "s/$/,$strand,Hypomesus_transpacificus/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done



## Hypomesus_transpacificus  ### FOURTH SCAFFOLD


grep ">" N5.HOG0018901.Hypomesus_transpacificus.extended.4.fa
scaffold=`grep ">" N5.HOG0018901.Hypomesus_transpacificus.extended.4.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0018901.Hypomesus_transpacificus.extended.4.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0018901.Hypomesus_transpacificus.extended.4.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Hypomesus_transpacificus,$scaffold.sec,$length" >> clusters_ID_TE.txt


grep "N5_HOG0018901"  GFF3_N5_OGGs/Hypomesus_transpacificus.gff.simplified.sorted.OGG.tiret
grep "XM_047044257" GFF3_files_per_species/Hypomesus_transpacificus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons


nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold.sec,$real_exon_start,$real_exon_stop,gene4.exon.$nbexon,N5_HOG0018901,+,Hypomesus_transpacificus" >> seq_clustered_infos_ogg.TE.txt
done


cut -f1 TE.Hypomesus_transpacificus.4.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Hypomesus_transpacificus.4.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0018901.Hypomesus_transpacificus.extended.4.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0018901.Hypomesus_transpacificus.extended.4.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold.sec,/g" | sed "s/$/,$strand,Hypomesus_transpacificus/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done





### Analyse conserved non-coding regions



samtools faidx N5.HOG0018901.Hypomesus_transpacificus.extended.1.rev NC_061067.1-9519063-9530253:500-1628 > region1.Hypomesus_transpacificus.fa
samtools faidx N5.HOG0018901.Hypomesus_transpacificus.extended.1.rev NC_061067.1-9519063-9530253:2500-5300 > region2.Hypomesus_transpacificus.fa
samtools faidx N5.HOG0018901.Hypomesus_transpacificus.extended.1.rev NC_061067.1-9519063-9530253:7000-7300 > region3.Hypomesus_transpacificus.fa

samtools faidx N5.HOG0018901.Clupea_harengus.extended.dotplot.fa NC_045166.1-15868305-15885905:5000-6000 > region1.Clupea_harengus.fa
samtools faidx N5.HOG0018901.Clupea_harengus.extended.dotplot.fa NC_045166.1-15868305-15885905:8500-11900 > region2.Clupea_harengus.fa
samtools faidx N5.HOG0018901.Clupea_harengus.extended.dotplot.fa NC_045166.1-15868305-15885905:3780-4600 > region3.Clupea_harengus.fa

sed -i 's/>.*/>H.transpacificus/g' region1.Hypomesus_transpacificus.fa
sed -i 's/>.*/>H.transpacificus/g' region2.Hypomesus_transpacificus.fa
sed -i 's/>.*/>H.transpacificus/g' region3.Hypomesus_transpacificus.fa

sed -i 's/>.*/>C.harengus/g' region1.Clupea_harengus.fa
sed -i 's/>.*/>C.harengus/g' region2.Clupea_harengus.fa
sed -i 's/>.*/>C.harengus/g' region3.Clupea_harengus.fa

needle -asequence region1.Hypomesus_transpacificus.fa -bsequence region1.Clupea_harengus.fa -outfile region1.aln -gapopen 10.0 -gapextend 0.5
needle -asequence region2.Hypomesus_transpacificus.fa -bsequence region2.Clupea_harengus.fa -outfile region2.aln -gapopen 10.0 -gapextend 0.5
needle -asequence region3.Hypomesus_transpacificus.fa -bsequence region3.Clupea_harengus.fa -outfile region3.aln -gapopen 10.0 -gapextend 0.5


sbatch --qos=6hours -c 10 --mem=20G launch_blast_allgenomes.sh region1.Hypomesus_transpacificus.fa 1e-20 region1.blastn.tsv
sbatch --qos=6hours -c 10 --mem=20G launch_blast_allgenomes.sh region2.Hypomesus_transpacificus.fa 1e-20 region2.blastn.tsv
sbatch --qos=6hours -c 10 --mem=20G launch_blast_allgenomes.sh region3.Hypomesus_transpacificus.fa 1e-20 region3.blastn.tsv 







##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##========================================== Finally, search if the gene is really absent in the other sp  ===========================
##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================


for scaffold in NC_085032.1 NC_085020.1 NC_085039.1 NC_085040.1 ; do

	scaffold_name=`echo "$scaffold" | sed 's/.1$/_1/g'`
	start=`grep "Osmerus" seq_clustered_infos_ogg.clade1.csv | grep "$scaffold_name" | sort -k2 -n -t "," | head -1 | cut -f2 -d ","`
	stop=`grep "Osmerus" seq_clustered_infos_ogg.clade1.csv | grep "$scaffold_name" | sort -k2 -n -t "," | tail -1 | cut -f3 -d ","`

	samtools faidx GCF_963692335.1_fOsmEpe2.1_genomic.fna $scaffold:$start-$stop >> Osmerus_eperlanus.whole.region.fa
done



exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs"  Hypomesus_transpacificus---rna-XM_047023826.1.prot Osmerus_eperlanus.whole.region.fa > region.Osmerus_eperlanus.exo
exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs"  Hypomesus_transpacificus---rna-XM_047023826.1.prot GCF_963692335.1_fOsmEpe2.1_genomic.fna > FULL_Osmerus_eperlanus.exo
exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs"  Hypomesus_transpacificus---rna-XM_047023826.1.prot GCF_902167405.1_gadMor3.0_genomic.fna > FULL_Gadus_morhua.exo
exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs"  Hypomesus_transpacificus---rna-XM_047023826.1.prot GCA_949987555.1_fBorAnt1.1_genomic.fna > FULL_Borostomias_antarcticus.exo



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

exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs"  Hypomesus_transpacificus---rna-XM_047023826.1.prot GCA_021870715.1_mHypTra1_genomic.fna > Hypomesus_transpacificus.alternative.exo


grep "Hypomesus" Coding_sequences_alignments/N5.HOG0018901.prot | sed 's/>//g' > Hypomesus.primary.id
xargs samtools faidx Coding_sequences_alignments/N5.HOG0018901.prot < Hypomesus.primary.id > Hypomesus.primary.prot 
transeq Hypomesus.alternative.fa Hypomesus.alternative.prot ; sed -i 's/_1$//g' Hypomesus.alternative.prot

cat Hypomesus.primary.prot Hypomesus.alternative.prot > Hypomesus.combined.prot

muscle5.1 -align Hypomesus.combined.prot -output Hypomesus.combined.prot.aln

trimal -in Hypomesus.combined.prot.aln -automated1 -out Hypomesus.combined.prot.aln.trimmed

iqtree -s Hypomesus.combined.prot.aln.trimmed -st AA -nt 8 -bb 1000 --redo 



#Hypomesus.alternative.fa
>Hypomesus_transpacificus.alternative.GCA_021870715.1.CM038953.1-10522339-10520683
ATGGCTCTCCGTATCGCAACACTCCTCTTTGCCCTGCACTGTCTGTTGCCGACTTCGTGTGCGGCGTATA
GATGTAAGCGTGTCTCTGGTGCCCTGGTGCAGATTGATGTAAGTGTTGGACAAGTGTTCGGGGTCAACAA
ACACGACTGCATCTTCACAAGGTATGGATCGTCATGGACACAGTTGCCAGGCAAACTGAAGCATGTCACC
GTGGGACCTGCAGGAGTCTGGGGAGCCAACAGAGCGAATCACATCTTTAAACTAGTTGGGGCTAACTGGG
TGAATGTACCTGGTCTCTTGAAGCAGATTGATGCTGGTGGAGACCAGTTTGTGGCAGGAGCCAATCATGG
AGATGCCATCTTCTGTCTTCCAAAGAAGAACACAGTTGGCTACAGTAAAAGAAACAGTGCTTTGAACTGG
CGCAATATCCCTGGCAGACTGAAGTACTACAGCTGTGGTCCTTACAGCTGCTGGGGTGTTAACAGTAATG
ATTACATCTTTGTGAGGAAGGGGGTGAGCTCCTTTAATTGTGAGGGGGAAGGTACCTGGCAGGGGGTTCC
AGGTCGCCTGTCTATGATCGAGGTGGGAAGTGACGGGTCTGTCTATGGAGTGAATTCTGTAGGAGACGTG
TACCGCAGGGATTCCACGTCTACCTGTAAACCTGAAGGTACCGGCTGGACCCGCATTCCTCTCTACAGCA
GACAGGTGAAACATGTGTCTTATGATCTGGGCCATCTCTGGCTCATCCTGAAGAATGATGCCATCTACGA
CTGCGCTGAGCAC
>Hypomesus_transpacificus.alternative.GCA_021870715.1.CM038954.1-17204857-17206425
ATGCCTCTCTGTATTGCAACACTCCTCTTTGCCCTACACTGTGTGTTGCTGACTTCATGTGTGGCATTTA
GATGTACACATCTCTCTGGTGCCCTGGTGCAGATTGATGTAGGTGTTGGACAAGTGATGGGGGTCAACAA
AGACGACAGCATCTTCACAAGGTTTGGATCGTCATGGACACAGTTGCCAGGCAAACTGAAGCATGTCACC
GTAGGACCTGCAGGAGTCTGGGGAGCCAACAGAGAGAATCTCATCTTTAAACTAGTTGGGACCGACTGGC
TGAACGTACCTGGTCTCTTGAAGCAGATTGATGCTGGTGGAGACCAATTTGTGGCAGGAGCCAATCATGT
CGATGCCCCCTTCTGTCTACCAAAGGAGAACACAGTTGGCTACAGTGGAAGTAATGGTGCTTTGAACTGG
CGCCAAATTCCTGGCAGCCTGAAGTACTACAGCTGTGGTCCTTATAGCTGCTGGGGTGTTAACAGTGGTG
ATCAAATCTTTGTGAGGAAGGGGGTGAACTCCTTTAATTGTGAGGGAGACGGTACCTGGCAGAAGGTTCC
AGGAAGCCTGTCTATGATCGAGGTGGGAAGTGACGGATCTGTCTATGGAGTGAATTCTGCAGGAGGCGTG
TTCCACAGGGATTCCACGTCTGCCTGTCAACCTGAAGGTACCGGCTGGACCCACCTTCCCCTCTACAGCG
GACAGATGAAACATGTGTCCTATGATCTGGGCCATCTCTGGCTCATCCTGAACAATGATGCCATCTACGA
CTGCACTGAGCAC
>Hypomesus_transpacificus.alternative.GCA_021870715.1.CM038966.1-4264671-4268452
ATGCCTCTCTGTATCGCAACACTCCTCTTTGCCCTGCACTGTCTGTTGCTGACTTCATGTtttgtcATGC
AGGATTTGGATGTACACATATCTCTGATGCAGATTGATGTAGGTGTTGGACAAGTTTTTGGGGTCAACCA
AGACGGCAGCATCTTCACAAGGTTTGGATCATCATGGACACAGTTGCCAGGCAAACTGAAGCATGTCACC
GTAGGACCTGCAGGAGTCTGGGGAGCCAACAGAGAGAATCTCATCTTTAAACTAGTTGGGACTGACTGGG
TGAATGTACCTGGTCTCTTGAAGCAGATTGATGCTGGTGGAGACCAGTTTGTGGCAGGAGCAAATCATGT
CGATGCAATCTACTGTCTTCCAAAGAAGAACACAGTTGGTTACAGTGGAAGTAATAGTGCTTTGAACTGG
CGCCAAATTCCTGGCGGCCTGAAGTACTACAGCTGTGGTCCTTATAGCTGCTGGGGTGTTAACAGTGCTG
ATCACATCTTTGTGAGGAAGGGGGTGAACTCCTTTAATTGTGAGGGGGACGGTACCTGGCAGAAGGTTCC
AGGAAGCCTGTCTATGATCGAGGTGGGAAGTGACGGATCTGTCTATGGAGTGAATTCTGTAGGAGACGTG
TTCCGCAGGGATTCCACCTCTGCTTGTCAACCTGAAGGTACCGGCTGGACCAACATTCCTCTCTACAGCG
GCCAGGTGAAACATGTGTCCTATGATCTGGGCCATCTTTGGCTCATCCTGAAGAATGATGCCATCTACGA
CTGCACTGAGCAC
>Hypomesus_transpacificus.alternative.GCA_021870715.1.CM038966.1-3704141-3700534
ATGCCTCTCTGTATCGCAACACTCCTCTTTGCCCTGCACTGTCTGTTGCTGACTTCATGTGTGGGATTTG
GATGTACACATATCTCTGGTGCCCTGATGCAGATTGATGTAGGTGTTGGACAAGTTTTTGGGGTCAACCA
AGACGGCAGCATCTTCACAAGGTTTGGATCATCATGGACACAGTTGCCAGGCAAACTGAAGCATGTCACC
GTAGGACCTGCAGGAGTCTGGGGAGCCAACAGAGAGAATCTCATCTTTAAACTAGTTGGGACTGACTGGG
TGAATGTACCTGGTCTCTTGAAGCAGATTGATGCTGGTGGAGACCAGTTTGTGGCAGGAGCAAATCATGT
CGATGCCATCTACTGTCTTCCAAAGAAGAACACAGTTGGCTACAGTGGAAGTAATAGTGCTTTGAACTGG
CGCCAAATTCCTGGCGGCCTGAAGTACTACAGCTGTGGTCCTTATAGCTGCTGGGGTGTTAACAGTGCTG
ATCACATCTTTGTGAGGAAGGGGGTGAACTCCTTTAATTGTGAGGGGGACGGTACCTGGCAGGAGGTTCC
AGGAAGCCTGTCTATGATCGAGGTGGGAAGTGACGGATCTGTCTATGGAGTGAATTCTGTAGGAGACGTG
TTCCGCAGGGATTCCACGTCTGCCTGTCAACCTGAAGGTACCGGCTGGACCAACATTCCTCTCTACAGCG
GACGGGTGAAACATGTGTCCTATGATCTGGGCCATCTCTGGCTCATCCTGAACAATGATGACATCTATGA
CTGCACTGAGCAC





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
sbatch --qos=1day -c 4 --mem=10G -e error.fitMG4.out -o slurm.fitMG4.out launch_fitMG4.sh N5.HOG0018901


Recipient branches : 
Hypomesus transpacificus rna XM 047023826 1
OROOT 
Node729
Hypomesus transpacificus rna XM 047025851 1
Hypomesus transpacificus rna XM 047044257 1
Hypomesus transpacificus rna XM 047043933 1



#Test positive selection and relaxed selection on receiver branch ==> put $HOG.prot.aln.treefile on https://phylotree.hyphy.org/ ==> exort to $HOG.prot.aln.treefile.SelecMarked

sbatch --qos=1week -c 4 --mem=10G -e error.absrel.cand.out -o slurm.absrel.cand.out --job-name=HOG0018901 launch_absrel_cand.sh N5.HOG0018901




##RELAX not launched because no branch detected under accelerated evolution by abSREL



#Extract dN/dS to table

grep "LB\":" N5.HOG0018901.cds.aln.FITTER.json | sed 's/.*://g' | sed 's/,.*//g' > curr_LB_values.txt
grep "MLE\":" N5.HOG0018901.cds.aln.FITTER.json | sed 's/.*://g' | sed 's/,.*//g' > curr_MLE_values.txt
grep "UB\":" N5.HOG0018901.cds.aln.FITTER.json | sed 's/.*://g' | sed 's/,.*//g' > curr_UB_values.txt
grep "\"dN\"" N5.HOG0018901.cds.aln.FITTER.json  | grep -v "{" | sed 's/.*://g' | sed 's/,//g' > curr_dN_values.txt
grep "\"dS\"" N5.HOG0018901.cds.aln.FITTER.json  | grep -v "{" | sed 's/.*://g' | sed 's/,//g' > curr_dS_values.txt
grep -B2 "LB\":" N5.HOG0018901.cds.aln.FITTER.json | grep -v "\-\-" | grep -v "Confidence" | grep -v "LB\":"  | sed 's/\"//g' | sed 's/:.*//g' | sed 's/^ *//g' > curr_labels

paste -d "," curr_labels curr_LB_values.txt curr_MLE_values.txt curr_UB_values.txt curr_dN_values.txt curr_dS_values.txt > N5.HOG0018901.dN_dS.csv










