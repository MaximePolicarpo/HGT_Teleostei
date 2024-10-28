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


samtools faidx Proteomes_BUSCO80/Osmerus_eperlanus.fa Osmerus_eperlanus---rna-XM_062453720.1 > Osmerus_eperlanus---rna-XM_062453720.1.prot

#blastp against the fish proteome (50 best hits) + Uniprot database without any fish species (10 best hits)

blastp -query Osmerus_eperlanus---rna-XM_062453720.1.prot -db concatenated_proteomes.fa -outfmt 6 -out Gene_vs_FishProteome.blastp -num_threads 8 -max_target_seqs 51
blastp -query Osmerus_eperlanus---rna-XM_062453720.1.prot -db non_actino_uniprot.fa -outfmt 6 -out Gene_vs_Uniprot.blastp -num_threads 8  -max_target_seqs 10 -evalue 1e-5


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


#N5.HOG0004480.clade_hgt.txt=> clade where the HGT is extracted from iTOL (from N5.HOG0004480.prot.aln.treefile)

cut -f1,2 -d "_" N5.HOG0004480.clade_hgt.txt | sort | uniq  > species_to_draw.clade1
#Choose close species without the gene .. 


echo "Engraulis_encrasicolus" > non_transfer_species_1
echo "Alosa_sapidissima" >> non_transfer_species_1
echo "Sardina_pilchardus" >> non_transfer_species_1
echo "Hypomesus_transpacificus" > non_transfer_species_2
echo "Borostomias_antarcticus" >> non_transfer_species_2
echo "Gadus_morhua" >> non_transfer_species_2



rm -r Syn_tables_dir ; mkdir Syn_tables_dir


#First order the species names 

rm species_to_draw.clade1.ordered

echo "Clupea_harengus" >> species_to_draw.clade1.ordered
echo "Osmerus_eperlanus" >> species_to_draw.clade1.ordered


#If there are more than 4 species from the same order, keep only 4 


#First extract genes around 

curr_OGG=N5.HOG0004480

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

curr_OGG=N5.HOG0004480
curr_sp=Sardina_pilchardus
ref_sp=Clupea_harengus

grep -A10 -B10 "N5_HOG0007348" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "	19" > Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0043664" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df



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

curr_OGG=N5.HOG0004480
curr_sp=Alosa_sapidissima
ref_sp=Clupea_harengus

rm Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0004480" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret  > Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0007348" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret  >> Syn_tables_dir/$curr_sp.synt.df



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

curr_OGG=N5.HOG0004480
curr_sp=Engraulis_encrasicolus
ref_sp=Clupea_harengus

rm Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0004480" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret  > Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0007348" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret  >> Syn_tables_dir/$curr_sp.synt.df



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

curr_OGG=N5.HOG0004480
curr_sp=Hypomesus_transpacificus
ref_sp=Osmerus_eperlanus

grep -A15 -B15 "N5_HOG0031985" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df


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

curr_OGG=N5.HOG0004480
curr_sp=Borostomias_antarcticus
ref_sp=Osmerus_eperlanus

grep -A15 -B15 "N5_HOG0031985" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df


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

curr_OGG=N5.HOG0004480
curr_sp=Gadus_morhua
ref_sp=Osmerus_eperlanus

grep -A10 -B10 "N5_HOG0022933" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0018394" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df


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



#Introduce our new annotation


sed -i 's/3801670,3807375,rna_XM_047042892_1,N5_HOG0004453,+,Hypomesus_transpacificus/3801702,3809259,custom_annot,N5_HOG0004480,+,Hypomesus_transpacificus/g' seq_clustered_infos_ogg.clade1.csv
grep -v "NC_061079_1,3807374,3809351" seq_clustered_infos_ogg.clade1.csv > temp ; mv temp seq_clustered_infos_ogg.clade1.csv

sed -i 's/8497910,8498954,rna_XM_062530891_1,N5_HOG0004476,-,Sardina_pilchardus/8497446,8547022,custom_annot,N5_HOG0004480,-,Sardina_pilchardus/g' seq_clustered_infos_ogg.clade1.csv
grep -v "NC_084996_1,8499944,8505144,rna_XM_062530892_1,N5_HOG0004453,-,Sardina_pilchardus" seq_clustered_infos_ogg.clade1.csv > temp ; mv temp seq_clustered_infos_ogg.clade1.csv
grep -v "NC_084996_1,8519853,8521771,rna_XM_062530893_1,N5_HOG0004453,-,Sardina_pilchardus" seq_clustered_infos_ogg.clade1.csv > temp ; mv temp seq_clustered_infos_ogg.clade1.csv



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




#Osmerus_eperlanus
grep -A2 -B2 "HOG0004480"  GFF3_N5_OGGs/Osmerus_eperlanus.gff.simplified.sorted.OGG.tiret
samtools faidx GCF_963692335.1_fOsmEpe2.1_genomic.fna NC_085044.1:6434631-6467776 > N5.HOG0004480.Osmerus_eperlanus.extended.1.fa
sed -i 's/:/-/g' N5.HOG0004480.Osmerus_eperlanus.extended.1.fa
makeblastdb -in N5.HOG0004480.Osmerus_eperlanus.extended.1.fa -dbtype nucl
revseq N5.HOG0004480.Osmerus_eperlanus.extended.1.fa  N5.HOG0004480.Osmerus_eperlanus.extended.1.fa.rev

grep "XM_062453720" GFF3_files_per_species/Osmerus_eperlanus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > curr_gene_1.exons
real_start=6434631
for line in `cat curr_gene_1.exons` ; do exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","` ; real_exon_start=$(( exon_start - real_start )) ; real_exon_stop=$(( exon_stop - real_start )) ; echo "$real_exon_start , $real_exon_stop" ; done




#Clupea_harengus
grep -A2 -B2 "HOG0004480"  GFF3_N5_OGGs/Clupea_harengus.gff.simplified.sorted.OGG.tiret
samtools faidx GCF_900700415.2_Ch_v2.0.2_genomic.fna NC_045152.1:4827654-4852619 > N5.HOG0004480.Clupea_harengus.extended.1.fa
sed -i 's/:/-/g' N5.HOG0004480.Clupea_harengus.extended.1.fa
makeblastdb -in N5.HOG0004480.Clupea_harengus.extended.1.fa -dbtype nucl



##One conserved non-coding regions between the two
header=`grep ">" N5.HOG0004480.Osmerus_eperlanus.extended.1.fa.rev | sed 's/>//g'`
samtools faidx N5.HOG0004480.Osmerus_eperlanus.extended.1.fa.rev NC_085044.1-6434631-6467776:26371-29900 > similarity_region_1.Osmerus_eperlanus.fa

header=`grep ">" N5.HOG0004480.Clupea_harengus.extended.1.fa | sed 's/>//g'`
samtools faidx N5.HOG0004480.Clupea_harengus.extended.1.fa $header:19000-23000 > similarity_region_1.Clupea_harengus.fa

needle -asequence similarity_region_1.Osmerus_eperlanus.fa -bsequence similarity_region_1.Clupea_harengus.fa -outfile similarity_region_1.aln -gapopen 10.0 -gapextend 0.5

cp ../launch_blast_allgenomes.sh ./

sbatch --qos=6hours -c 10 --mem=20G launch_blast_allgenomes.sh similarity_region_1.Osmerus_eperlanus.fa 1e-20 region1.blastn.tsv




#Now lets find TEs using the repbase and dfam databases

tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0004480.Osmerus_eperlanus.extended.1.fa -evalue 1e-5 -outfmt 6 -out TE.Osmerus_eperlanus.1.tblastn -num_threads 8
sed -i 's/#//g' TE.Osmerus_eperlanus.1.tblastn

tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0004480.Clupea_harengus.extended.1.fa -evalue 1e-5 -outfmt 6 -out TE.Clupea_harengus.1.tblastn -num_threads 8
sed -i 's/#//g' TE.Clupea_harengus.1.tblastn


#merge tblastn hits and find the best TE match by doing a blastx

cp ../Rscript_merge_blast_hits.R ./
Rscript Rscript_merge_blast_hits.R TE.Osmerus_eperlanus.1.tblastn TE.Osmerus_eperlanus.1.tblastn.merged
Rscript Rscript_merge_blast_hits.R TE.Clupea_harengus.1.tblastn TE.Clupea_harengus.1.tblastn.merged


xargs samtools faidx N5.HOG0004480.Osmerus_eperlanus.extended.1.fa < TE.Osmerus_eperlanus.1.tblastn.merged > TE.Osmerus_eperlanus.1.BEST.fa
blastx -query TE.Osmerus_eperlanus.1.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Osmerus_eperlanus.1.BEST.blastx -max_target_seqs 1
cut -f1  TE.Osmerus_eperlanus.1.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Osmerus_eperlanus.1.BEST.blastx >> temp  ; done ; mv temp TE.Osmerus_eperlanus.1.BEST.blastx

xargs samtools faidx N5.HOG0004480.Clupea_harengus.extended.1.fa < TE.Clupea_harengus.1.tblastn.merged > TE.Clupea_harengus.1.BEST.fa
blastx -query TE.Clupea_harengus.1.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Clupea_harengus.1.BEST.blastx -max_target_seqs 1
cut -f1  TE.Clupea_harengus.1.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Clupea_harengus.1.BEST.blastx >> temp  ; done ; mv temp TE.Clupea_harengus.1.BEST.blastx



#Now find shared elements

cut -f2 TE.Osmerus_eperlanus.1.BEST.blastx | sort | uniq > TE.Osmerus_eperlanus.1.uniqTE
cut -f2 TE.Clupea_harengus.1.BEST.blastx | sort | uniq > TE.Clupea_harengus.1.uniqTE

comm -12 TE.Osmerus_eperlanus.1.uniqTE TE.Clupea_harengus.1.uniqTE



###ONE SHARED ELEMENT = L2-3_EL_1p:ClassI:LINE:Jockey:L2


#### Now let's find every copy of this element in the genome of both species + closely related species


samtools faidx Dfam_plus_Repbase.cdhit80.prot L2-3_EL_1p:ClassI:LINE:Jockey:L2 > L2-3_EL.prot


#Hypomesus_transpacificus
tblastn -query L2-3_EL.prot -db GCF_021917145.1_fHypTra1_genomic.fna -evalue 1e-2 -outfmt 6 -out L2-3_EL.Hypomesus_transpacificus.tblastn -num_threads 8
sed -i 's/#//g' L2-3_EL.Hypomesus_transpacificus.tblastn
Rscript Rscript_merge_blast_hits.R L2-3_EL.Hypomesus_transpacificus.tblastn L2-3_EL.Hypomesus_transpacificus.tblastn.merged
xargs samtools faidx GCF_021917145.1_fHypTra1_genomic.fna < L2-3_EL.Hypomesus_transpacificus.tblastn.merged > L2-3_EL.Hypomesus_transpacificus.BEST.fa
diamond blastx --query L2-3_EL.Hypomesus_transpacificus.BEST.fa --db Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out L2-3_EL.Hypomesus_transpacificus.BEST.blastx --max-target-seqs 1
grep "L2-3_EL" L2-3_EL.Hypomesus_transpacificus.BEST.blastx > L2-3_EL.Hypomesus_transpacificus.list

#Clupea_harengus
tblastn -query L2-3_EL.prot -db GCF_900700415.2_Ch_v2.0.2_genomic.fna -evalue 1e-2 -outfmt 6 -out L2-3_EL.Clupea_harengus.tblastn -num_threads 8
sed -i 's/#//g' L2-3_EL.Clupea_harengus.tblastn
Rscript Rscript_merge_blast_hits.R L2-3_EL.Clupea_harengus.tblastn L2-3_EL.Clupea_harengus.tblastn.merged
xargs samtools faidx GCF_900700415.2_Ch_v2.0.2_genomic.fna  < L2-3_EL.Clupea_harengus.tblastn.merged > L2-3_EL.Clupea_harengus.BEST.fa
diamond blastx --query L2-3_EL.Clupea_harengus.BEST.fa --db Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out L2-3_EL.Clupea_harengus.BEST.blastx --max-target-seqs 1
grep "L2-3_EL" L2-3_EL.Clupea_harengus.BEST.blastx > L2-3_EL.Clupea_harengus.list


#Sardina_pilchardus
tblastn -query L2-3_EL.prot -db GCF_963854185.1_fSarPil1.1_genomic.fna -evalue 1e-2 -outfmt 6 -out L2-3_EL.Sardina_pilchardus.tblastn -num_threads 8
sed -i 's/#//g' L2-3_EL.Sardina_pilchardus.tblastn
Rscript Rscript_merge_blast_hits.R L2-3_EL.Sardina_pilchardus.tblastn L2-3_EL.Sardina_pilchardus.tblastn.merged
xargs samtools faidx GCF_963854185.1_fSarPil1.1_genomic.fna  < L2-3_EL.Sardina_pilchardus.tblastn.merged > L2-3_EL.Sardina_pilchardus.BEST.fa
diamond blastx --query L2-3_EL.Sardina_pilchardus.BEST.fa --db Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out L2-3_EL.Sardina_pilchardus.BEST.blastx --max-target-seqs 1
grep "L2-3_EL" L2-3_EL.Sardina_pilchardus.BEST.blastx > L2-3_EL.Sardina_pilchardus.list

#Alosa_sapidissima
tblastn -query L2-3_EL.prot -db GCF_018492685.1_fAloSap1.pri_genomic.fna -evalue 1e-2 -outfmt 6 -out L2-3_EL.Alosa_sapidissima.tblastn -num_threads 8
sed -i 's/#//g' L2-3_EL.Alosa_sapidissima.tblastn
Rscript Rscript_merge_blast_hits.R L2-3_EL.Alosa_sapidissima.tblastn L2-3_EL.Alosa_sapidissima.tblastn.merged
xargs samtools faidx GCF_018492685.1_fAloSap1.pri_genomic.fna  < L2-3_EL.Alosa_sapidissima.tblastn.merged > L2-3_EL.Alosa_sapidissima.BEST.fa
diamond blastx --query L2-3_EL.Alosa_sapidissima.BEST.fa --db Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out L2-3_EL.Alosa_sapidissima.BEST.blastx --max-target-seqs 1
grep "L2-3_EL" L2-3_EL.Alosa_sapidissima.BEST.blastx > L2-3_EL.Alosa_sapidissima.list


#Engraulis_encrasicolus
tblastn -query L2-3_EL.prot -db GCF_034702125.1_IST_EnEncr_1.0_genomic.fna -evalue 1e-2 -outfmt 6 -out L2-3_EL.Engraulis_encrasicolus.tblastn -num_threads 8
sed -i 's/#//g' L2-3_EL.Engraulis_encrasicolus.tblastn
Rscript Rscript_merge_blast_hits.R L2-3_EL.Engraulis_encrasicolus.tblastn L2-3_EL.Engraulis_encrasicolus.tblastn.merged
xargs samtools faidx GCF_034702125.1_IST_EnEncr_1.0_genomic.fna  < L2-3_EL.Engraulis_encrasicolus.tblastn.merged > L2-3_EL.Engraulis_encrasicolus.BEST.fa
diamond blastx --query L2-3_EL.Engraulis_encrasicolus.BEST.fa --db Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out L2-3_EL.Engraulis_encrasicolus.BEST.blastx --max-target-seqs 1
grep "L2-3_EL" L2-3_EL.Engraulis_encrasicolus.BEST.blastx > L2-3_EL.Engraulis_encrasicolus.list


#Osmerus_eperlanus
tblastn -query L2-3_EL.prot -db GCF_963692335.1_fOsmEpe2.1_genomic.fna -evalue 1e-2 -outfmt 6 -out L2-3_EL.Osmerus_eperlanus.tblastn -num_threads 8
sed -i 's/#//g' L2-3_EL.Osmerus_eperlanus.tblastn
Rscript Rscript_merge_blast_hits.R L2-3_EL.Osmerus_eperlanus.tblastn L2-3_EL.Osmerus_eperlanus.tblastn.merged
xargs samtools faidx GCF_963692335.1_fOsmEpe2.1_genomic.fna  < L2-3_EL.Osmerus_eperlanus.tblastn.merged > L2-3_EL.Osmerus_eperlanus.BEST.fa
diamond blastx --query L2-3_EL.Osmerus_eperlanus.BEST.fa --db Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out L2-3_EL.Osmerus_eperlanus.BEST.blastx --max-target-seqs 1
grep "L2-3_EL" L2-3_EL.Osmerus_eperlanus.BEST.blastx > L2-3_EL.Osmerus_eperlanus.list

#Borostomias_antarcticus
tblastn -query L2-3_EL.prot -db GCA_949987555.1_fBorAnt1.1_genomic.fna -evalue 1e-2 -outfmt 6 -out L2-3_EL.Borostomias_antarcticus.tblastn -num_threads 8
sed -i 's/#//g' L2-3_EL.Borostomias_antarcticus.tblastn
Rscript Rscript_merge_blast_hits.R L2-3_EL.Borostomias_antarcticus.tblastn L2-3_EL.Borostomias_antarcticus.tblastn.merged
xargs samtools faidx GCA_949987555.1_fBorAnt1.1_genomic.fna  < L2-3_EL.Borostomias_antarcticus.tblastn.merged > L2-3_EL.Borostomias_antarcticus.BEST.fa
diamond blastx --query L2-3_EL.Borostomias_antarcticus.BEST.fa --db Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out L2-3_EL.Borostomias_antarcticus.BEST.blastx --max-target-seqs 1
grep "L2-3_EL" L2-3_EL.Borostomias_antarcticus.BEST.blastx > L2-3_EL.Borostomias_antarcticus.list

#Gadus_morhua
tblastn -query L2-3_EL.prot -db GCF_902167405.1_gadMor3.0_genomic.fna -evalue 1e-2 -outfmt 6 -out L2-3_EL.Gadus_morhua.tblastn -num_threads 8
sed -i 's/#//g' L2-3_EL.Gadus_morhua.tblastn
Rscript Rscript_merge_blast_hits.R L2-3_EL.Gadus_morhua.tblastn L2-3_EL.Gadus_morhua.tblastn.merged
xargs samtools faidx GCF_902167405.1_gadMor3.0_genomic.fna  < L2-3_EL.Gadus_morhua.tblastn.merged > L2-3_EL.Gadus_morhua.BEST.fa
diamond blastx --query L2-3_EL.Gadus_morhua.BEST.fa --db Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out L2-3_EL.Gadus_morhua.BEST.blastx --max-target-seqs 1
grep "L2-3_EL" L2-3_EL.Gadus_morhua.BEST.blastx > L2-3_EL.Gadus_morhua.list


wc -l L2-3_EL.Engraulis_encrasicolus.list
wc -l L2-3_EL.Alosa_sapidissima.list
wc -l L2-3_EL.Sardina_pilchardus.list
wc -l L2-3_EL.Clupea_harengus.list
wc -l L2-3_EL.Hypomesus_transpacificus.list
wc -l L2-3_EL.Osmerus_eperlanus.list #+1
wc -l L2-3_EL.Borostomias_antarcticus.list
wc -l L2-3_EL.Gadus_morhua.list

### In a FINAL step we will make another synteny plot, only between Clupea and Hypomesus to show N5.HOG0001647 and neighboring regions 


## Clupea harengus 

rm seq_clustered_infos_ogg.TE.txt ; rm clusters_ID_TE.txt


scaffold=`grep ">" N5.HOG0004480.Clupea_harengus.extended.1.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0004480.Clupea_harengus.extended.1.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0004480.Clupea_harengus.extended.1.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Clupea_harengus,$scaffold,$length" > clusters_ID_TE.txt

grep "N5_HOG0004480"  GFF3_N5_OGGs/Clupea_harengus.gff.simplified.sorted.OGG.tiret

grep "XM_031565097" GFF3_files_per_species/Clupea_harengus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons


nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene1.$nbexon,N5_HOG0004480,+,Clupea_harengus" >> seq_clustered_infos_ogg.TE.txt
done


cut -f1 TE.Clupea_harengus.1.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Clupea_harengus.1.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0004480.Clupea_harengus.extended.1.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0004480.Clupea_harengus.extended.1.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold,/g" | sed "s/$/,$strand,Clupea_harengus/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done



## Osmerus_eperlanus 


grep ">" N5.HOG0004480.Osmerus_eperlanus.extended.1.fa
scaffold=`grep ">" N5.HOG0004480.Osmerus_eperlanus.extended.1.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0004480.Osmerus_eperlanus.extended.1.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0004480.Osmerus_eperlanus.extended.1.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Osmerus_eperlanus,$scaffold,$length" >> clusters_ID_TE.txt


grep "N5_HOG0004480"  GFF3_N5_OGGs/Osmerus_eperlanus.gff.simplified.sorted.OGG.tiret
grep "XM_062453720" GFF3_files_per_species/Osmerus_eperlanus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons


nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene1.exon.$nbexon,N5_HOG0004480,-,Osmerus_eperlanus" >> seq_clustered_infos_ogg.TE.txt
done



cut -f1 TE.Osmerus_eperlanus.1.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Osmerus_eperlanus.1.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0004480.Osmerus_eperlanus.extended.1.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0004480.Osmerus_eperlanus.extended.1.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold,/g" | sed "s/$/,$strand,Osmerus_eperlanus/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done



### Now, try to make an alignment with the RTEX-2 element which is just next to the gene and compute the identity percentage



#NC_045152.1,6129,7341,L2-3_EL_1p:ClassI:LINE:Jockey:L2,L2-3_EL_1p:ClassI:LINE:Jockey:L2,-,Clupea_harengus
#NC_085044.1,10839,11423,L2-3_EL_1p:ClassI:LINE:Jockey:L2,L2-3_EL_1p:ClassI:LINE:Jockey:L2,+,Osmerus_eperlanus


samtools faidx N5.HOG0004480.Osmerus_eperlanus.extended.1.fa NC_085044.1-6434631-6467776:6129-7341 > Osmerus.HGT_L2-3_EL.fa
samtools faidx N5.HOG0004480.Clupea_harengus.extended.1.fa NC_045152.1-4827654-4852619:10839-11423 > Clupea.HGT_L2-3_EL.fa


sed -i 's/>.*/>Oeperlanus_L2-3_EL_HGT/g' Osmerus.HGT_L2-3_EL.fa
sed -i 's/>.*/>Charengus_L2-3_EL_HGT/g' Clupea.HGT_L2-3_EL.fa

cat Osmerus.HGT_L2-3_EL.fa Clupea.HGT_L2-3_EL.fa > L2-3_EL_combined.fa

needle -asequence Osmerus.HGT_L2-3_EL.fa -bsequence Clupea.HGT_L2-3_EL.fa -outfile L2-3_EL_combined.aln -gapopen 10.0 -gapextend 0.5
needle -asequence Osmerus.HGT_L2-3_EL.rev -bsequence Clupea.HGT_L2-3_EL.fa -outfile L2-3_EL_combined.aln -gapopen 10.0 -gapextend 0.5

#Verify that the gene is absent from Hyomesus, Borostomia and Gadus and Sardina

scaffold=`grep "Hypomesus" seq_clustered_infos_ogg.clade1.csv | head -1 | cut -f1 -d "," | sed 's/_1$/.1/g'`
start=`grep "Hypomesus" seq_clustered_infos_ogg.clade1.csv | sort -k2 -n -t "," | head -1 | cut -f2 -d ","`
stop=`grep "Hypomesus" seq_clustered_infos_ogg.clade1.csv | sort -k2 -n -t "," | tail -1 | cut -f3 -d ","`
samtools faidx GCF_021917145.1_fHypTra1_genomic.fna $scaffold:$start-$stop > Hypomesus_transpacificus.whole.region.fa


scaffold=`grep "Borostomia" seq_clustered_infos_ogg.clade1.csv | head -1 | cut -f1 -d "," | sed 's/_1$/.1/g'`
start=`grep "Borostomia" seq_clustered_infos_ogg.clade1.csv | sort -k2 -n -t "," | head -1 | cut -f2 -d ","`
stop=`grep "Borostomia" seq_clustered_infos_ogg.clade1.csv | sort -k2 -n -t "," | tail -1 | cut -f3 -d ","`
samtools faidx GCA_949987555.1_fBorAnt1.1_genomic.fna $scaffold:$start-$stop > Borostomias_antarcticus.whole.region.fa


scaffold=`grep "Gadus" seq_clustered_infos_ogg.clade1.csv | head -1 | cut -f1 -d "," | sed 's/_1$/.1/g'`
start=`grep "Gadus" seq_clustered_infos_ogg.clade1.csv | sort -k2 -n -t "," | head -1 | cut -f2 -d ","`
stop=`grep "Gadus" seq_clustered_infos_ogg.clade1.csv | sort -k2 -n -t "," | tail -1 | cut -f3 -d ","`
samtools faidx GCF_902167405.1_gadMor3.0_genomic.fna $scaffold:$start-$stop > Gadus_morhua.whole.region.fa


scaffold=`grep "Sardina" seq_clustered_infos_ogg.clade1.csv | head -1 | cut -f1 -d "," | sed 's/_1$/.1/g'`
start=`grep "Sardina" seq_clustered_infos_ogg.clade1.csv | sort -k2 -n -t "," | head -1 | cut -f2 -d ","`
stop=`grep "Sardina" seq_clustered_infos_ogg.clade1.csv | sort -k2 -n -t "," | tail -1 | cut -f3 -d ","`
samtools faidx GCF_963854185.1_fSarPil1.1_genomic.fna $scaffold:$start-$stop > Sardina_pilchardus.whole.region.fa

exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs"  Osmerus_eperlanus---rna-XM_062453720.1.prot Hypomesus_transpacificus.whole.region.fa > Hypomesus_transpacificus.exo
exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs"  Osmerus_eperlanus---rna-XM_062453720.1.prot Borostomias_antarcticus.whole.region.fa > Borostomias_antarcticus.exo
exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs"  Osmerus_eperlanus---rna-XM_062453720.1.prot Gadus_morhua.whole.region.fa > Gadus_morhua.exo
exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs"  Osmerus_eperlanus---rna-XM_062453720.1.prot Sardina_pilchardus.whole.region.fa > Sardina_pilchardus.exo


exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs"  Osmerus_eperlanus---rna-XM_062453720.1.prot GCF_902167405.1_gadMor3.0_genomic.fna> FULL_Gadus_morhua.exo
exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs"  Osmerus_eperlanus---rna-XM_062453720.1.prot GCA_949987555.1_fBorAnt1.1_genomic.fna > FULL_Borostomias_antarcticus.exo


## GENE IS PRESENT IN HYPOMESUS AND SARDINA ===>



Hypomesus : NC_061079.1:3675897-4009674 125805 -> 133362 (true region = NC_061079.1:3801702-3809259)
Sardina:  NC_084996.1:8268216-20010399 278806 -> 229230 (true region = NC_084996.1:8497446-8547022)

transeq missing_annot.fa missing_annot.prot ; sed -i 's/_1$//g' missing_annot.prot



grep "Clupea_harengus\|Alosa_sapidissima\|Engraulis_encrasicolus\|Osmerus\|Chanos" Coding_sequences_alignments/N5.HOG0004480.prot | sed 's/>//g' > seq.id
xargs samtools faidx Coding_sequences_alignments/N5.HOG0004480.prot < seq.id > seq.prot

cat seq.prot missing_annot.prot > missing_annot_plus_annot.prot


muscle5.1 -align missing_annot_plus_annot.prot -output missing_annot_plus_annot.aln

iqtree -s missing_annot_plus_annot.aln -st AA -nt 8 -bb 1000 --redo -m LG+F+G4




##missing_annot.fa

>Sardina_pilchardus-N5.HOG0004480
ATGTCTGTGAGTGACACCATCGAGCTTGCAGGGTTGGGCCGCCCATTCCAGCTGGGAATGCTGTATGACT
GCAGGAGAGATGTTCTCATTCCAGGTATTACACTCTGGGATGCTGACATGCTGCAGAAGAATATCAATAT
ACACTCACAGCCAAACACTGAGTTCAAGATCATTGCCTCGGACTCTACAGAGGCGAAGGCAGAAGCTCTA
CATGTGTCAGCATCACTAGAGGCTAGCTTTATGGGTGGGTTAGTAAGTGTAAAAGGCTCTGCAGAGTTTC
TGAATGACAAGAAAAAGTCAAAAAATCAGTCCAGAGTTACTCTTCAGTATCGCACTACCACACGCTTTGA
GCAGCTTACCATGGAACACCTGGGAACTGGAAATTTCCAACATTGCAATGTCTTCCAGGAAGGATCTGCT
ACACATGTGGTGACTGCAATATTGTATGGTGCTCAAGCCTTCTTTGTATTTGATCGTGAAGTCTCTTCTA
ATGAAAATCATCAGTTAATTCAAGGAAGTCTGCAAGCAACAATCAAAGTGATTCCAATGGTATCTATTGA
GGGACAGGCATCCCTGAAGATGACTGAAAAAGAACAACAAGACGCAGACAAATTCAACTGTACATTTCAT
GGTGATTTTGCTTTGGAGAATAATCCTGTCAGTTTTCAAGATGCAATCAGAGTTTATTCCCATCTTCCAC
GACTGCTTGGGGAAAATGGTGAACATGCTGTGCCCATGACAGTGTGGTTGTATCCCCTCAAAAATCTGGA
CTCCACTGCAGCCCAACTTGTTAGGGAGATAAGTATAGGTCTGGTGCGCCAGGCTCAGCGCATCCTTGAT
GAGATAGACGAGGCTGACATTCAGTGCCAAGACTTGATGAAAGATGAGATAGCCATCCTCTTCTCTGAGA
TTACAACCaagcttcgaaaattcaaaggTCTGATTTCTGAGTATAAATTGGTGTTTCAAAAGAAGCTTTG
CAAATTACTTCCCACTATTCGCGGAGGAGGAACAGAAGAGCAAGAGTTGGTGACCACTCTGAACAGCAAA
GAGAGATCCCCCTTCCAGAGTGCACTAATCACAGAATATCTGAGTGACCGAGAACGTGAGATGAATGTTA
TCAGATCTTACCTTGAGATTATGAAAGAAATCAGGGTGACTAAATCCAGTAGTGAACTGGATAAAGTTGT
CTTACATCCAAGAAATGAATTTGTAGTGGCTTTTGTGTTCTCTTCCCTCAATGGTGATTGCAAATATCTA
CAAGATTTGGAAAACTACCTGAATGATGAGGCCACTGGCTCAGACACATTTTATGAAGCCAAAAGTGGCG
ACAGACAGGAGCATTGGTTTCATTCTGGGGAGGTGACTGCTCTTACCAGGCAAAGCATTCACCTGTTCCT
TGATTTCATGCAATCCAACAGGGACAGGGAAAACATTGAGTTTTGCATAGCATCAGTTCCAAACAAGTCC
ATCACAGCCTCATCAATTCATGTGTATGAGAAAGGGACACTTCTTAATGCTCAGTTTTTGCTGCCTAGCA
AGCCTCCAAGCCTGACTGTTTTGAAACTGGAACATGATCGTGTCCATTTGCTAATCAACCCTCCTACCAT
GGGAGAAAGTTTTGTCGTCTCCTATCGCATTCTGTACCAGACCTCTGAGGGCTCTGACTGGAAAGAAATA
GGCACTGATGAAAAAGCAACAGATGTTACTGTTAGTCGCTTACAGCCTCATACAGAGTACCGATTTAGCT
GTAAGGCAGTGTGTCGACCGGGAGTGAGTCTTGCTAGCGATGCAACAGAGTTCATCCGAACACGTCCATG
CAGTCCTCCAGGCCCGCCAAATGAGAAAAGTTCAGGGGCTGAAAGCATCACTATTACATGGGATATTCCC
ACAGCAGTTGGAGACGATGTGTCAGTCATCAGTTATGTTGTTGAATACCGACAGCACATGGAAGATGATA
CCAAGAGCAAACCATGGGAGTCTGAAAAAGCCTCTAGTAGGGAGTGTACATTGGAGGGATTGAAGATGAA
CACCACCTACAGCATCCGTGTTTTGGCCAATACAGGGAAATCAGGAAAGAGTCTGCCAAGTCCTGAAACA
AATATTACTACATCAGGACCTGATGCTAAATCTATGAAAAAGAAATCAGCACAGGGCAGGAGTGAGCAGT
TTGTAGAACAATCATACTGTCTGGAAAAGGGCAACCCATCTGTATATCAGATAGAGTTGGATCAAATATA
TGGAGAGAATGTTGATTTCTTTCAGTATGTGTTTGGTAGGAAGGTGGAAGATGTGAAAAACAAAGTCATA
CTACTACTGGGATCCACTGGTGCAGGAAAAACAACACTTGTAAATGCAATGATAAATTACATACTTGGAG
TCAAATGGGAGGACCAGTACCGCTACAAACTAATCCATGAGGTGACTAACCGAACACAAGCAGAAAGCCA
GACATCTGTGGTCACATCTTATGAGCTATATAACCAACCAGGCTTTCAGGTCCCATACTCACTCACAATT
ATCGACACACCTGGATTTGGGGACACCAGAGGAATGTCTCAAGATAAATTGATCACAGAGCAAGTGAAAA
GATTCCTTTGCAGCCCCTTAGGGGTTGAACATATTGATGCTGTGTGCTTTGTTGTTCAGGCATCTCTTGC
CCGTCTCAGTGCGAACCAGAAATACATATTTGATTCTATCCTGTCAATTTTTGGAAAAGACATTGGAGAA
AACATCATGGTTCTTGTAACATTTGCTGACAGTGACAACATCCCAGTTTTAGAAGCAATTAAGGCAGCAG
AGTTGCCATGTCAGAAGAACAAGAAAGGTCAACCCACCCACTTTAAGTTTAACAACTCTACTGTATATGT
CCAGAAAAAGGATGTACACAAAAAAACTGATGAAGATGAttcagatgaagatgaagatgatgatgacgat
gaggaGAAAATGAAGAGTATTGTCTGGTCAACAACTTTCAAACAGATGAAGCAGTTTTTCAAGGCTCTCG
AAAGCATAGAGAGCAAAGATCTCACTTTGACAAAGAGAGTTCTAGAGGAACGTGAACGTCTTGAGAATGC
GACAACCAGACTGACCCCGCAAATTACAGCTGGGCTTGCCAAGCTCAGTGAGATCAAAAGCATCAAACAG
TGCCTGGAGAATGAGGATGAGATGATGAAACAAAGCAACAACTTTGAGACAGAAGTAGAGGTACTGACAG
CTATTCGGAGCAATGTGAACTTCTTTGCAACAAACTGTAATAACTGTTTGTTCACATGCCATTCAGGTTG
TTTCCTTCCTGAGGGTGATGGTGTCCACACATGTGCTGTGATGGATGAGAGTGGGAAGTGTATAGTGTGT
CCAAAAAATTGCCTTTCTTCTATGCATATAAGAGAGAAAGCTCTATGGACATATGAGTCAAAGATGGAGA
AAAAAACCATTGAAGAGTTGAAGAAAAACTTCATGGATGCCCGGGGGAAATTTATGGATAAGAAACAGAT
GTTAGATACATTAGAGGAGGATTTTCATGAAATTGAGGACAAACTTATGTATTTGATCAAGCTATCCTCT
AACTGCTTGAAGAGGCTTAACGAGATTGCACTGAAGCCAAGTTCTATATCCACCTTCGAGTACATTGAAC
TACTTATCAGAAcagaagaggatgagagaaagcCGGGATTTGAGGATCGAATCATTGGATTGAGAAAAAT
GAAACATGAAGCTGAAATTCTGGATAAAATTGCAAGAGGAGTAAATGTCATGCCACAAGATCGTCGCATG
GTTAAACAAAAACTG
>Hypomesus_transpacificus-N5.HOG0004480
ATGCTGTATGACTGCAGGAAAGATGCCC
TCATCCCAGGTATTACACTCTGGGATGGCGACATGCTGCAGAAGAACATTAATATCCAaccacatccaaa
cactgagttCAAGATCATTGCCTCAGATTCCACAGAGAGCAAGGCAGAGGCTCTGAAAGTGTCTGCATCA
CTTGAGGCAAGCTTTATGACTGGGTTGGTGAGCGTAAAAGGATCTGCAGCCTTTCTGAATGACAAGAAAA
AATCGAAACAACAGTCCAGAGTCACTCTTCAGTATCACACCACCACACGTTATGAGCAGCTTACAATGGA
ACACTTAGGGGCTGGAAATTTGAAACACTGCGATGTCTTGAAGGAAGGATCTGCCACACATGTGGTGACT
GCAATATTGTATGGTGCTAATGCCTTCTTTGTGTTTGATCGTGAAGTCTCTTCGAATGAAAATCATCAGG
AAATACAAGGACATCTTGAGgtagaaataaaaaaaatcccatTTATATCTATTAAGGCAGAGGCAGATTT
GAAGATGAGtgaaacagaacaacaacaagcaGACAAATTCAACTGTACATTCCATGGTGATTTTGCATTG
GTGAATAATCCTGTCAGCTTTCAAGATGCAATCAAGGTGTACTCACAGCTTCCACAGCTGCTGGGAGAAA
ACGGGCAACATGCTGTACCCATGGAAGTGTGGCTGTATCCCCTTAAAAAACTGGACTCGACTGCGGCCCA
GCTTGTTAGGGAGATAAGTGTAGGTCTTGTTCGCCATTCTCAGCGCATCCTTGATGACATAGACGAGGCT
GACATTGATTGCCATGACCTGATGAAAGAGGAGATAGCCATCCACTTCTCTGAGATTACAACCAAGCTTA
GGAAATTCAAAGGCCTGATTTCTGAATACAAACTGGTATTTCAAAAGAAGCTTTGCAAATTACTTCCTAC
TATTCgcggaggggggacagaggagcaAGAGTTGGTGACTACCCTgaacagcagagagagatccCCATTC
CAGAGCGCATTAATCACAGAATATCTGAGTGACAGAGAACGTGAGATAAATGTTATCAGATCCTACCTTG
ACATTATGAAAGAAATCAGGGTGATTAAATCCAGTAGTGAACTGGACAAGGTTGTCCTCAATCCAAGTAA
TGACTTTGTAGTGGCTTGTGTGTTCTCATCCCTCAGTGGGGATACAGAATATCTACATGATTTGGAAAAC
TATCTGAAGGACGAGGCAAAGGCCATTAGGTCTGACACATCTTATGAGGTCCAAAGTGCTAACAAAAAGG
AACAATGGTTTCACTCTGGGGAGGTGACTGCTCTTACCAGGCAAAGCATTCGCCTGTTCCTTGATTTCAA
GCAGTCCAACCAAGACAGAGAAAACATTGAGTTTTGCATAGGATCAGTTCCAAACAAGTCCATCACAGCC
TCCTCAATTCATGTGTATGAGAAAGGGACACTTCTTGATTCTCAGTTTTTGCTgccctccacacctccta
gCCCGACCCTTTTGAGACAGGAACATGATCGTGTCCACTTGCAAATCCACCCCCCGGCCACAGGAGAACC
TTTTGTGGATTCTTATCGAATCCTGTACCAGACGTCTGAGTCTGACTGGAAGGAGATAGGCACTGATGCA
AAAGTAAGAGATTTTACTGTTAGTCGCCTACAGCCTCATACAGAGTACCGATTTAGCTGTAAGGCAGTGT
GTCGACCTGGAGTGAGCCTTGCTAGCGATGCAACAGACTTCATCCAAACACGTTCGTGCAGCCCTCCTGG
TCCACCAGAGGAGAAACGTTCAGGGTCTGAGAGCATTACTATTACATGGGATATCCCAACAGCAGTTGGA
GAAGATGTAACAGTCATCAGCTATGATGTTGAATACAGACAACACATGGAAGATGACACTAAGAGCAAAC
CATGGGAGTGTGTAAAAGCCACAACGAGGGAGTGCACTTTGGAAGGACTGAAGACGAGTACCACCTACAG
CATCCGTGTTTTGGCGAATACAGGGAAACCAGGAAAGAGTCTGCCAAGCCCTGAAAGTAATATTTCTACA
TTAGCACCTGACACCAAGCTCCCAAAAAAGAAATCGGCACAGGGCAAGAGTGAAAATTTTATGCAACCAT
CGAATCGTCTGGAAAAGGGCAACCCATCTGTATATCAGCTAGAGTTAGAAAGGAAATATGGAGAGAATGT
TGATTTCTTGCAATACATGTTTGGTAGGAAGGTGGAAGATGTGCAAAACAAGGTCATACTTCTTCTGGGA
TCTAGTGGCGCAGGAAAAACAACTCTTGTCAACGTGATGATGAATTACATTCTGGGGGTTAAATGGGAAG
ACCAGTACCGCTTCAAACTGATCCATGAAGTGACTAACCGGACACAAGCAGAAAGCCAGACATCGCTGGT
CACATCCTATGAGCTCTATAACCAACCAGGTTTTCAGATCCCATACTCACTTACAATTGTTGACACACCT
GGGTTTGGAGACACCAGAGGAATGGCGCAAGACAAACTGATTACAGAACAAGTGAAACATTTCCTGTGTA
ACCCTTTAGGGGTTGATCACATAGATGCTGTATGCTTTGTCGTTCAGGCATCTCTCGCCCGCCTCAGTGC
CAACCAGAAATACATCTTTGACTCGATTCTGTCTATTTTTGGAAAGGATATTGGCGAAAACATCATGGTT
CTTGTAACATTTGTTGACAGTGACAACATTCCAGTTTTAGAGGCAATTAAGGCTGCAGAGTTGCCATGTC
AGAAGAACAAGAAGGGCCAACCCACCCACTTTAAGTTTAACAACTCCACCTTATATGTGCAGAAAAAGGA
AATGGACGGAGGCACTGAGGAAGACAGTtcagatgaagatgatgatgaagaggaggaggagaaaatgaAG
AGCATTATCTGGTCAACAACCTTCAAACAGATGAAGTTGTTCTTCAAGGCCCTTGGTAGCACAGAGGGCA
AAGATCTCAAAATGACAAGGAAAGTTCTTGAAGAACGTGAACGTCTTGAGAATGCCACTGCAAGACTGAC
TCCTCAGATCACAGCTGGGCTTGCCAAGCTCAGTGAGATCAAAAACTTTAAACAATGCCTGGAGAATGAG
AATGAAATGATGGCTCAAAATAAGGACTTTGAGACGGAAGTAGAGGTACTGACAGCCAAACGGAACAACG
TGAACTTCTTTGCAATGAACTGTAACAGTTGTCTGTTTACATGTCATTCTGGTTGTTTCCTTCCTAAAGG
GGATAATCTCCTCACATGTGCTGTGATGGATGACGATGGGAACTGTGTGTTATGCCCAAATAACTGCTCT
TCCTTTGAGCATCTCAGAGAGCAAGCCTTGTGGACATATGAGACAAAGATGGAGAAAAAGACCATACAAG
AGCTGAAGGAGAACTTCATGAAAGCACAAGGAAAATTCATGGACACCAAGCAGATGCTTGAGGAATTGGA
GAATGAATTTCATGTGATTGAGGATAAACTGATGTATTTGATCAAGCTATCCTCTAACTGCATGAAGAGA
CTCAATGAGATTGCACTGAAGCCAAGTTCTATGTCAACCTTCGATTACATTGAAATACTTATccgaacag
aggaggaggagaaaaagccTGGGTTTGAGGATCGAATCATTGGGTTGAAGAAAATGAAACAGGAAGCGGA
AATTCTACATAAGATTGCAAGAGGAGAAGATGTACTGCCACAAGAGCGCCGCATGGTTAAACACAAACTG
CAC



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
sbatch --qos=1day -c 4 --mem=10G -e error.fitMG4.out -o slurm.fitMG4.out launch_fitMG4.sh N5.HOG0004480



Recipient branches : 
Osmerus_eperlanus_rna_XM_062453720_1



#Test positive selection and relaxed selection on receiver branch

sbatch --qos=1day -c 4 --mem=10G -e error.absrel.cand.out -o slurm.absrel.cand.out launch_absrel_cand.sh N5.HOG0004480

#Now test for relaxed selection
sbatch --qos=1week -c 4 --mem=10G -e error.relax.out -o slurm.relax.out --job-name=HOG0004480 launch_RELAX.sh N5.HOG0004480


#Extract dN/dS to table

grep "LB\":" N5.HOG0004480.cds.aln.FITTER.json | sed 's/.*://g' | sed 's/,.*//g' > curr_LB_values.txt
grep "MLE\":" N5.HOG0004480.cds.aln.FITTER.json | sed 's/.*://g' | sed 's/,.*//g' > curr_MLE_values.txt
grep "UB\":" N5.HOG0004480.cds.aln.FITTER.json | sed 's/.*://g' | sed 's/,.*//g' > curr_UB_values.txt
grep "\"dN\"" N5.HOG0004480.cds.aln.FITTER.json  | grep -v "{" | sed 's/.*://g' | sed 's/,//g' > curr_dN_values.txt
grep "\"dS\"" N5.HOG0004480.cds.aln.FITTER.json  | grep -v "{" | sed 's/.*://g' | sed 's/,//g' > curr_dS_values.txt
grep -B2 "LB\":" N5.HOG0004480.cds.aln.FITTER.json | grep -v "\-\-" | grep -v "Confidence" | grep -v "LB\":"  | sed 's/\"//g' | sed 's/:.*//g' | sed 's/^ *//g' > curr_labels

paste -d "," curr_labels curr_LB_values.txt curr_MLE_values.txt curr_UB_values.txt curr_dN_values.txt curr_dS_values.txt > N5.HOG0004480.dN_dS.csv



