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



samtools faidx Proteomes_BUSCO80/Hypomesus_transpacificus.fa Hypomesus_transpacificus---rna-XM_047045767.1 > Hypomesus_transpacificus---rna-XM_047045767.1.prot

#blastp against the fish proteome (50 best hits) + Uniprot database without any fish species (10 best hits)

blastp -query Hypomesus_transpacificus---rna-XM_047045767.1.prot -db concatenated_proteomes.fa -outfmt 6 -out Gene_vs_FishProteome.blastp -num_threads 8 -max_target_seqs 51
blastp -query Hypomesus_transpacificus---rna-XM_047045767.1.prot -db non_actino_uniprot.fa -outfmt 6 -out Gene_vs_Uniprot.blastp -num_threads 8  -max_target_seqs 10 -evalue 1e-5


cut -f2 Gene_vs_FishProteome.blastp | sort | uniq  > closest_fish_seq.id
cut -f2 Gene_vs_Uniprot.blastp | sort | uniq > closest_nonfish_seq.id 

xargs samtools faidx concatenated_proteomes.fa < closest_fish_seq.id > closest_fish_seq.fa
xargs samtools faidx non_actino_uniprot.fa  < closest_nonfish_seq.id > closest_nonfish_seq.fa

#Align with muscle and trim with trimal

cat closest_fish_seq.fa closest_nonfish_seq.fa > Uniprot_plus_closefish.fa

muscle5.1 -align closest_fish_seq.fa -output closest_fish_seq.aln
mafft --add closest_nonfish_seq.fa --keeplength closest_fish_seq.aln > Uniprot_plus_closefish.aln
trimal -in Uniprot_plus_closefish.aln -automated1 -out Uniprot_plus_closefish.aln.trimmed

#Make a ML tree with IQ-TREE
iqtree -s Uniprot_plus_closefish.aln -st AA -nt 8 -m TEST -mrate G4
iqtree -s Uniprot_plus_closefish.aln.trimmed -st AA -nt 8 -m TEST -mrate G4


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


#N5.HOG0013514.clade_hgt.txt=> clade where the HGT is extracted from iTOL (from N5.HOG0013514.prot.aln.treefile)

cut -f1,2 -d "_" N5.HOG0013514.clade_hgt.txt | sort | uniq  > species_to_draw.clade1
#Choose close species without the gene .. 


echo "Engraulis_encrasicolus" > non_transfer_species_1
echo "Alosa_alosa" >> non_transfer_species_1
echo "Sardina_pilchardus" >> non_transfer_species_1
echo "Borostomias_antarcticus" > non_transfer_species_2
echo "Gadus_morhua" >> non_transfer_species_2



rm -r Syn_tables_dir ; mkdir Syn_tables_dir


#First order the species names 

rm species_to_draw.clade1.ordered


echo "Clupea_harengus" > species_to_draw.clade1.ordered
echo "Hypomesus_transpacificus" >> species_to_draw.clade1.ordered
echo "Osmerus_eperlanus" >> species_to_draw.clade1.ordered


#If there are more than 4 species from the same order, keep only 4 


#First extract genes around 

curr_OGG=N5.HOG0013514

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




#Now add the non transfer speices



curr_OGG=N5.HOG0013514


for curr_sp in `cat non_transfer_species_1` ; do


	rm Syn_tables_dir/$curr_sp.synt.df

	ref_sp=Clupea_harengus #ref species with the HGT .. 
	
	i=0
	while [ $i -le 10 ]; do
	
		i=$((i + 1))
		close_up=`grep -A$i "$curr_OGG" GFF3_N5_OGGs/$ref_sp.gff.simplified.sorted.OGG.tiret | tail -1 | cut -f5`
		close_down=`grep -B$i "$curr_OGG" GFF3_N5_OGGs/$ref_sp.gff.simplified.sorted.OGG.tiret | head -1 | cut -f5`
	
	
		if grep -q "N5_HOG0025783" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret ; then 
	
			grep -A10 -B10 "N5_HOG0025783" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df
			i=11
	
		elif grep -q "N5_HOG0025783" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret ; then
	
			grep -A10 -B10 "N5_HOG0025783" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df
			i=11
	
		fi	
	
	done 


done




for curr_sp in `cat non_transfer_species_1` ; do

	ref_sp=Clupea_harengus #ref species with the HGT .. 
	
	i=0
	while [ $i -le 10 ]; do
	
		i=$((i + 1))
		close_up=`grep -A$i "$curr_OGG" GFF3_N5_OGGs/$ref_sp.gff.simplified.sorted.OGG.tiret | tail -1 | cut -f5`
		close_down=`grep -B$i "$curr_OGG" GFF3_N5_OGGs/$ref_sp.gff.simplified.sorted.OGG.tiret | head -1 | cut -f5`
	
	
		if grep -q "N5_HOG0027636" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret ; then 
	
			grep -A10 -B10 "N5_HOG0027636" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df
			i=11
	
		elif grep -q "N5_HOG0027636" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret ; then
	
			grep -A10 -B10 "N5_HOG0027636" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df
			i=11
	
		fi	
	
	done 

	awk '{ if (NF < 5 || $5 == "") $5 = "no_OGG"; print $0;}' Syn_tables_dir/$curr_sp.synt.df > Syn_tables_dir/$curr_sp.synt.final.df
	sed -i 's/	/,/g' Syn_tables_dir/$curr_sp.synt.final.df
	sed -i 's/ no_OGG/,no_OGG/g' Syn_tables_dir/$curr_sp.synt.final.df
	sed -i 's/ /,/g' Syn_tables_dir/$curr_sp.synt.final.df
	grep -v ",,,," Syn_tables_dir/$curr_sp.synt.final.df > temp ; mv temp Syn_tables_dir/$curr_sp.synt.final.df


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

done








curr_OGG=N5.HOG0013514

for curr_sp in `cat non_transfer_species_2` ; do


	rm Syn_tables_dir/$curr_sp.synt.df

	ref_sp=Hypomesus_transpacificus #ref species with the HGT .. 
	
	i=1
	while [ $i -le 10 ]; do
	
		i=$((i + 1))
		close_up=`grep -A$i "$curr_OGG" GFF3_N5_OGGs/$ref_sp.gff.simplified.sorted.OGG.tiret | tail -1 | cut -f5`
		close_down=`grep -B$i "$curr_OGG" GFF3_N5_OGGs/$ref_sp.gff.simplified.sorted.OGG.tiret | head -1 | cut -f5`
	
	
		if grep -q "$close_up" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret ; then 
	
			grep -A10 -B10 "$close_up" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df
			i=11
	
		elif grep -q "$close_down" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret ; then
	
			grep -A10 -B10 "$close_down" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df
			i=11
	
		fi	
	
	done 


	awk '{ if (NF < 5 || $5 == "") $5 = "no_OGG"; print $0;}' Syn_tables_dir/$curr_sp.synt.df > Syn_tables_dir/$curr_sp.synt.final.df
	sed -i 's/	/,/g' Syn_tables_dir/$curr_sp.synt.final.df
	sed -i 's/ no_OGG/,no_OGG/g' Syn_tables_dir/$curr_sp.synt.final.df
	sed -i 's/ /,/g' Syn_tables_dir/$curr_sp.synt.final.df
	grep -v ",,,," Syn_tables_dir/$curr_sp.synt.final.df > temp ; mv temp Syn_tables_dir/$curr_sp.synt.final.df

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


done

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


#We will take 5kb upstream and downstream of the gene OR gene cluster (or less if there is an other gene before)


#Hypomesus_transpacificus
grep -A2 -B2 "N5_HOG0013514"  GFF3_N5_OGGs/Hypomesus_transpacificus.gff.simplified.sorted.OGG.tiret
samtools faidx GCF_021917145.1_fHypTra1_genomic.fna NC_061062.1:6295541-6339487 > N5.HOG0013514.Hypomesus_transpacificus.extended.1.fa
sed -i 's/:/-/g' N5.HOG0013514.Hypomesus_transpacificus.extended.1.fa

makeblastdb -in N5.HOG0013514.Hypomesus_transpacificus.extended.1.fa -dbtype nucl


#Clupea_harengus
grep -A2 -B2 "N5_HOG0013514"  GFF3_N5_OGGs/Clupea_harengus.gff.simplified.sorted.OGG.tiret
samtools faidx GCF_900700415.2_Ch_v2.0.2_genomic.fna NC_045155.1:5717415-5788992 > N5.HOG0013514.Clupea_harengus.extended.1.fa
samtools faidx GCF_900700415.2_Ch_v2.0.2_genomic.fna NC_045155.1:1092148-1106954 > N5.HOG0013514.Clupea_harengus.extended.2.fa

sed -i 's/:/-/g' N5.HOG0013514.Clupea_harengus.extended.1.fa
makeblastdb -in N5.HOG0013514.Clupea_harengus.extended.1.fa -dbtype nucl

sed -i 's/:/-/g' N5.HOG0013514.Clupea_harengus.extended.2.fa
makeblastdb -in N5.HOG0013514.Clupea_harengus.extended.2.fa -dbtype nucl




#Now lets find TEs using the repbase and dfam databases

tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0013514.Hypomesus_transpacificus.extended.1.fa -evalue 1e-5 -outfmt 6 -out TE.Hypomesus_transpacificus.1.tblastn -num_threads 8
sed -i 's/#//g' TE.Hypomesus_transpacificus.1.tblastn


tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0013514.Clupea_harengus.extended.1.fa -evalue 1e-5 -outfmt 6 -out TE.Clupea_harengus.1.tblastn -num_threads 8
sed -i 's/#//g' TE.Clupea_harengus.1.tblastn

tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0013514.Clupea_harengus.extended.2.fa -evalue 1e-5 -outfmt 6 -out TE.Clupea_harengus.2.tblastn -num_threads 8
sed -i 's/#//g' TE.Clupea_harengus.2.tblastn


#merge tblastn hits and find the best TE match by doing a blastx

cp ../Rscript_merge_blast_hits.R ./
Rscript Rscript_merge_blast_hits.R TE.Hypomesus_transpacificus.1.tblastn TE.Hypomesus_transpacificus.1.tblastn.merged
Rscript Rscript_merge_blast_hits.R TE.Clupea_harengus.1.tblastn TE.Clupea_harengus.1.tblastn.merged
Rscript Rscript_merge_blast_hits.R TE.Clupea_harengus.2.tblastn TE.Clupea_harengus.2.tblastn.merged


xargs samtools faidx N5.HOG0013514.Hypomesus_transpacificus.extended.1.fa < TE.Hypomesus_transpacificus.1.tblastn.merged > TE.Hypomesus_transpacificus.1.BEST.fa
blastx -query TE.Hypomesus_transpacificus.1.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Hypomesus_transpacificus.1.BEST.blastx -max_target_seqs 1
cut -f1  TE.Hypomesus_transpacificus.1.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Hypomesus_transpacificus.1.BEST.blastx >> temp  ; done ; mv temp TE.Hypomesus_transpacificus.1.BEST.blastx

xargs samtools faidx N5.HOG0013514.Clupea_harengus.extended.1.fa < TE.Clupea_harengus.1.tblastn.merged > TE.Clupea_harengus.1.BEST.fa
blastx -query TE.Clupea_harengus.1.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Clupea_harengus.1.BEST.blastx -max_target_seqs 1
cut -f1  TE.Clupea_harengus.1.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Clupea_harengus.1.BEST.blastx >> temp  ; done ; mv temp TE.Clupea_harengus.1.BEST.blastx


xargs samtools faidx N5.HOG0013514.Clupea_harengus.extended.2.fa < TE.Clupea_harengus.2.tblastn.merged > TE.Clupea_harengus.2.BEST.fa
blastx -query TE.Clupea_harengus.2.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Clupea_harengus.2.BEST.blastx -max_target_seqs 1
cut -f1  TE.Clupea_harengus.2.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Clupea_harengus.2.BEST.blastx >> temp  ; done ; mv temp TE.Clupea_harengus.2.BEST.blastx





#Now find shared elements

cut -f2 TE.Hypomesus_transpacificus.1.BEST.blastx | sort | uniq > TE.Hypomesus_transpacificus.1.uniqTE
cut -f2 TE.Clupea_harengus.1.BEST.blastx | sort | uniq > TE.Clupea_harengus.1.uniqTE
cut -f2 TE.Clupea_harengus.2.BEST.blastx | sort | uniq > TE.Clupea_harengus.2.uniqTE

comm -12 TE.Hypomesus_transpacificus.1.uniqTE TE.Clupea_harengus.1.uniqTE
comm -12 TE.Hypomesus_transpacificus.1.uniqTE TE.Clupea_harengus.2.uniqTE



##ONE SHARED ELEMENT : HERO-1_AFC



samtools faidx Dfam_plus_Repbase.cdhit80.prot HERO-1_AFC_pol#LINE/R2-Hero > HERO-1.prot


#### Now let's find every copy of this element in the genome of both species + closely related species


#Hypomesus_transpacificus
tblastn -query HERO-1.prot -db GCF_021917145.1_fHypTra1_genomic.fna -outfmt 6 -out HERO-1.Hypomesus_transpacificus.tblastn -num_threads 8
sed -i 's/#//g' HERO-1.Hypomesus_transpacificus.tblastn
Rscript Rscript_merge_blast_hits.R HERO-1.Hypomesus_transpacificus.tblastn HERO-1.Hypomesus_transpacificus.tblastn.merged
xargs samtools faidx GCF_021917145.1_fHypTra1_genomic.fna < HERO-1.Hypomesus_transpacificus.tblastn.merged > HERO-1.Hypomesus_transpacificus.BEST.fa
diamond blastx --query HERO-1.Hypomesus_transpacificus.BEST.fa --db Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out HERO-1.Hypomesus_transpacificus.BEST.blastx --max-target-seqs 1
grep "HERO-1" HERO-1.Hypomesus_transpacificus.BEST.blastx > HERO-1.Hypomesus_transpacificus.list

#Clupea_harengus
tblastn -query HERO-1.prot -db GCF_900700415.2_Ch_v2.0.2_genomic.fna -outfmt 6 -out HERO-1.Clupea_harengus.tblastn -num_threads 8
sed -i 's/#//g' HERO-1.Clupea_harengus.tblastn
Rscript Rscript_merge_blast_hits.R HERO-1.Clupea_harengus.tblastn HERO-1.Clupea_harengus.tblastn.merged
xargs samtools faidx GCF_900700415.2_Ch_v2.0.2_genomic.fna  < HERO-1.Clupea_harengus.tblastn.merged > HERO-1.Clupea_harengus.BEST.fa
diamond blastx --query HERO-1.Clupea_harengus.BEST.fa --db Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out HERO-1.Clupea_harengus.BEST.blastx --max-target-seqs 1
grep "HERO-1" HERO-1.Clupea_harengus.BEST.blastx > HERO-1.Clupea_harengus.list


#Sardina_pilchardus
tblastn -query HERO-1.prot -db GCF_963854185.1_fSarPil1.1_genomic.fna -outfmt 6 -out HERO-1.Sardina_pilchardus.tblastn -num_threads 8
sed -i 's/#//g' HERO-1.Sardina_pilchardus.tblastn
Rscript Rscript_merge_blast_hits.R HERO-1.Sardina_pilchardus.tblastn HERO-1.Sardina_pilchardus.tblastn.merged
xargs samtools faidx GCF_963854185.1_fSarPil1.1_genomic.fna  < HERO-1.Sardina_pilchardus.tblastn.merged > HERO-1.Sardina_pilchardus.BEST.fa
diamond blastx --query HERO-1.Sardina_pilchardus.BEST.fa --db Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out HERO-1.Sardina_pilchardus.BEST.blastx --max-target-seqs 1
grep "HERO-1" HERO-1.Sardina_pilchardus.BEST.blastx > HERO-1.Sardina_pilchardus.list

#Alosa_alosa
tblastn -query HERO-1.prot -db GCF_017589495.1_AALO_Geno_1.1_genomic.fna -outfmt 6 -out HERO-1.Alosa_alosa.tblastn -num_threads 8
sed -i 's/#//g' HERO-1.Alosa_alosa.tblastn
Rscript Rscript_merge_blast_hits.R HERO-1.Alosa_alosa.tblastn HERO-1.Alosa_alosa.tblastn.merged
xargs samtools faidx GCF_017589495.1_AALO_Geno_1.1_genomic.fna  < HERO-1.Alosa_alosa.tblastn.merged > HERO-1.Alosa_alosa.BEST.fa
diamond blastx --query HERO-1.Alosa_alosa.BEST.fa --db Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out HERO-1.Alosa_alosa.BEST.blastx --max-target-seqs 1
grep "HERO-1" HERO-1.Alosa_alosa.BEST.blastx > HERO-1.Alosa_alosa.list


#Engraulis_encrasicolus
tblastn -query HERO-1.prot -db GCF_034702125.1_IST_EnEncr_1.0_genomic.fna -outfmt 6 -out HERO-1.Engraulis_encrasicolus.tblastn -num_threads 8
sed -i 's/#//g' HERO-1.Engraulis_encrasicolus.tblastn
Rscript Rscript_merge_blast_hits.R HERO-1.Engraulis_encrasicolus.tblastn HERO-1.Engraulis_encrasicolus.tblastn.merged
xargs samtools faidx GCF_034702125.1_IST_EnEncr_1.0_genomic.fna  < HERO-1.Engraulis_encrasicolus.tblastn.merged > HERO-1.Engraulis_encrasicolus.BEST.fa
diamond blastx --query HERO-1.Engraulis_encrasicolus.BEST.fa --db Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out HERO-1.Engraulis_encrasicolus.BEST.blastx --max-target-seqs 1
grep "HERO-1" HERO-1.Engraulis_encrasicolus.BEST.blastx > HERO-1.Engraulis_encrasicolus.list


#Osmerus_eperlanus
tblastn -query HERO-1.prot -db GCF_963692335.1_fOsmEpe2.1_genomic.fna -outfmt 6 -out HERO-1.Osmerus_eperlanus.tblastn -num_threads 8
sed -i 's/#//g' HERO-1.Osmerus_eperlanus.tblastn
Rscript Rscript_merge_blast_hits.R HERO-1.Osmerus_eperlanus.tblastn HERO-1.Osmerus_eperlanus.tblastn.merged
xargs samtools faidx GCF_963692335.1_fOsmEpe2.1_genomic.fna  < HERO-1.Osmerus_eperlanus.tblastn.merged > HERO-1.Osmerus_eperlanus.BEST.fa
diamond blastx --query HERO-1.Osmerus_eperlanus.BEST.fa --db Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out HERO-1.Osmerus_eperlanus.BEST.blastx --max-target-seqs 1
grep "HERO-1" HERO-1.Osmerus_eperlanus.BEST.blastx > HERO-1.Osmerus_eperlanus.list

#Borostomias_antarcticus
tblastn -query HERO-1.prot -db GCA_949987555.1_fBorAnt1.1_genomic.fna -outfmt 6 -out HERO-1.Borostomias_antarcticus.tblastn -num_threads 8
sed -i 's/#//g' HERO-1.Borostomias_antarcticus.tblastn
Rscript Rscript_merge_blast_hits.R HERO-1.Borostomias_antarcticus.tblastn HERO-1.Borostomias_antarcticus.tblastn.merged
xargs samtools faidx GCA_949987555.1_fBorAnt1.1_genomic.fna  < HERO-1.Borostomias_antarcticus.tblastn.merged > HERO-1.Borostomias_antarcticus.BEST.fa
diamond blastx --query HERO-1.Borostomias_antarcticus.BEST.fa --db Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out HERO-1.Borostomias_antarcticus.BEST.blastx --max-target-seqs 1
grep "HERO-1" HERO-1.Borostomias_antarcticus.BEST.blastx > HERO-1.Borostomias_antarcticus.list

#Gadus_morhua
tblastn -query HERO-1.prot -db GCF_902167405.1_gadMor3.0_genomic.fna -outfmt 6 -out HERO-1.Gadus_morhua.tblastn -num_threads 8
sed -i 's/#//g' HERO-1.Gadus_morhua.tblastn
Rscript Rscript_merge_blast_hits.R HERO-1.Gadus_morhua.tblastn HERO-1.Gadus_morhua.tblastn.merged
xargs samtools faidx GCF_902167405.1_gadMor3.0_genomic.fna  < HERO-1.Gadus_morhua.tblastn.merged > HERO-1.Gadus_morhua.BEST.fa
diamond blastx --query HERO-1.Gadus_morhua.BEST.fa --db Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out HERO-1.Gadus_morhua.BEST.blastx --max-target-seqs 1
grep "HERO-1" HERO-1.Gadus_morhua.BEST.blastx > HERO-1.Gadus_morhua.list


wc -l HERO-1.Engraulis_encrasicolus.list
wc -l HERO-1.Alosa_alosa.list
wc -l HERO-1.Sardina_pilchardus.list
wc -l HERO-1.Clupea_harengus.list
wc -l HERO-1.Hypomesus_transpacificus.list
wc -l HERO-1.Osmerus_eperlanus.list
wc -l HERO-1.Borostomias_antarcticus.list
wc -l HERO-1.Gadus_morhua.list


### In a FINAL step we will make another synteny plot, only between Clupea and Hypomesus to show N5.HOG0001647 and neighboring regions 


## Clupea harengus 

rm seq_clustered_infos_ogg.TE.txt ; rm clusters_ID_TE.txt


scaffold=`grep ">" N5.HOG0013514.Clupea_harengus.extended.1.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0013514.Clupea_harengus.extended.1.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0013514.Clupea_harengus.extended.1.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Clupea_harengus,$scaffold,$length" > clusters_ID_TE.txt

grep "N5_HOG0013514"  GFF3_N5_OGGs/Clupea_harengus.gff.simplified.sorted.OGG.tiret

grep "XM_042707511" GFF3_files_per_species/Clupea_harengus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons
grep "XM_031565905" GFF3_files_per_species/Clupea_harengus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_2.exons
grep "XM_031565911" GFF3_files_per_species/Clupea_harengus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_3.exons
grep "XM_031565910" GFF3_files_per_species/Clupea_harengus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_4.exons


nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene1.$nbexon,N5_HOG0013514,-,Clupea_harengus" >> seq_clustered_infos_ogg.TE.txt
done

nbexon=0
for line in `cat gene_2.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene2.$nbexon,N5_HOG0013514,-,Clupea_harengus" >> seq_clustered_infos_ogg.TE.txt
done
nbexon=0
for line in `cat gene_3.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene3.$nbexon,N5_HOG0013514,-,Clupea_harengus" >> seq_clustered_infos_ogg.TE.txt
done
nbexon=0
for line in `cat gene_4.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene4.$nbexon,N5_HOG0013514,-,Clupea_harengus" >> seq_clustered_infos_ogg.TE.txt
done


cut -f1 TE.Clupea_harengus.1.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Clupea_harengus.1.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0013514.Clupea_harengus.extended.1.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0013514.Clupea_harengus.extended.1.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold,/g" | sed "s/$/,$strand,Clupea_harengus/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done



#Second region of Clupea


scaffold=`grep ">" N5.HOG0013514.Clupea_harengus.extended.2.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0013514.Clupea_harengus.extended.2.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0013514.Clupea_harengus.extended.2.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Clupea_harengus,$scaffold.sec,$length" >> clusters_ID_TE.txt

grep "N5_HOG0013514"  GFF3_N5_OGGs/Clupea_harengus.gff.simplified.sorted.OGG.tiret
grep "XM_031565680" GFF3_files_per_species/Clupea_harengus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons

nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold.sec,$real_exon_start,$real_exon_stop,gene1.$nbexon,N5_HOG0013514,-,Clupea_harengus" >> seq_clustered_infos_ogg.TE.txt
done


cut -f1 TE.Clupea_harengus.2.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Clupea_harengus.2.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0013514.Clupea_harengus.extended.2.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0013514.Clupea_harengus.extended.2.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold.sec,/g" | sed "s/$/,$strand,Clupea_harengus/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done




## Hypomesus_transpacificus 


grep ">" N5.HOG0013514.Hypomesus_transpacificus.extended.1.fa
scaffold=`grep ">" N5.HOG0013514.Hypomesus_transpacificus.extended.1.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0013514.Hypomesus_transpacificus.extended.1.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0013514.Hypomesus_transpacificus.extended.1.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Hypomesus_transpacificus,$scaffold,$length" >> clusters_ID_TE.txt


grep "N5_HOG0013514"  GFF3_N5_OGGs/Hypomesus_transpacificus.gff.simplified.sorted.OGG.tiret
grep "XM_047045786" GFF3_files_per_species/Hypomesus_transpacificus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons
grep "XM_047045767" GFF3_files_per_species/Hypomesus_transpacificus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_2.exons


nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene1.exon.$nbexon,N5_HOG0013514,+,Hypomesus_transpacificus" >> seq_clustered_infos_ogg.TE.txt
done

nbexon=0
for line in `cat gene_2.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene2.exon.$nbexon,N5_HOG0013514,+,Hypomesus_transpacificus" >> seq_clustered_infos_ogg.TE.txt
done


cut -f1 TE.Hypomesus_transpacificus.1.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Hypomesus_transpacificus.1.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0013514.Hypomesus_transpacificus.extended.1.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0013514.Hypomesus_transpacificus.extended.1.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold,/g" | sed "s/$/,$strand,Hypomesus_transpacificus/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done




##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##========================================== Extract full gene length on Clupea  ===================================================
##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================



grep ">Hypomesus\|Clupea" Coding_sequences_alignments/N5.HOG0013514.prot | sed 's/>//g' > toetxract.id

exonerate -E True --showtargetgff TRUE --model protein2genome --ryo "%tcs" Hypomesus_transpacificus---rna-XM_047045767.1.prot N5.HOG0013514.Clupea_harengus.extended.1.fa > gene.clupea.1.exo
exonerate -E True --showtargetgff TRUE --model protein2genome --ryo "%tcs" Hypomesus_transpacificus---rna-XM_047045767.1.prot N5.HOG0013514.Clupea_harengus.extended.2.fa > gene.clupea.2.exo

transeq Full_length_Clupea.fa Full_length_Clupea.prot ; sed -i 's/_1$//g' Full_length_Clupea.prot


xargs samtools faidx Coding_sequences_alignments/N5.HOG0013514.prot < toetxract.id > toetxract.prot 
cat Full_length_Clupea.prot >> toetxract.prot 

xargs samtools faidx Coding_sequences_alignments/N5.HOG0013514.cds < toetxract.id > toetxract.fa

/scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Opsins_analysis_2024/muscle5.1.linux_intel64 -align toetxract.prot  -output Full_length_clupea.aln

iqtree -s Full_length_clupea.aln -st AA -nt 8 -m TEST -mrate G4 


#Recompute the stats with the full length alignment

cat Full_length_Clupea.fa toetxract.fa > ALL.fa

trimal -in Full_length_clupea.aln -cons 100 -backtrans ALL.fa -out Full_length_clupea.cds.aln

Rscript Compute_aligned_CDS_stats_FAST.R Full_length_clupea.cds.aln Full_length.stats




### Check the conserved regions between Hypomesus and Clupea


samtools faidx N5.HOG0013514.Hypomesus_transpacificus.extended.1.fa NC_061062.1-6295541-6339487:12352-13352 > region1.Hypomesus_transpacificus.fa
samtools faidx N5.HOG0013514.Hypomesus_transpacificus.extended.1.fa NC_061062.1-6295541-6339487:39823-40823 > region2.Hypomesus_transpacificus.fa

samtools faidx N5.HOG0013514.Clupea_harengus.extended.1.fa.rev NC_045155.1-5717415-5788992:13636-14636 > region1.Clupea_harengus.fa
samtools faidx N5.HOG0013514.Clupea_harengus.extended.1.fa.rev NC_045155.1-5717415-5788992:39960-40960 > region2.Clupea_harengus.fa
samtools faidx N5.HOG0013514.Clupea_harengus.extended.1.fa.rev NC_045155.1-5717415-5788992:51540-52540 > region3.Clupea_harengus.fa
samtools faidx N5.HOG0013514.Clupea_harengus.extended.1.fa.rev NC_045155.1-5717415-5788992:65989-66989 > region4.Clupea_harengus.fa

sed -i 's/>.*/>H.transpacificus/g' region1.Hypomesus_transpacificus.fa
sed -i 's/>.*/>H.transpacificus/g' region2.Hypomesus_transpacificus.fa

sed -i 's/>.*/>C.harengus/g' region1.Clupea_harengus.fa
sed -i 's/>.*/>C.harengus/g' region2.Clupea_harengus.fa
sed -i 's/>.*/>C.harengus/g' region3.Clupea_harengus.fa
sed -i 's/>.*/>C.harengus/g' region4.Clupea_harengus.fa

needle -asequence region1.Hypomesus_transpacificus.fa -bsequence region1.Clupea_harengus.fa -outfile region1.aln -gapopen 10.0 -gapextend 0.5
needle -asequence region1.Hypomesus_transpacificus.fa -bsequence region2.Clupea_harengus.fa -outfile region2.aln -gapopen 10.0 -gapextend 0.5
needle -asequence region1.Hypomesus_transpacificus.fa -bsequence region3.Clupea_harengus.fa -outfile region3.aln -gapopen 10.0 -gapextend 0.5
needle -asequence region1.Hypomesus_transpacificus.fa -bsequence region4.Clupea_harengus.fa -outfile region4.aln -gapopen 10.0 -gapextend 0.5

needle -asequence region2.Hypomesus_transpacificus.fa -bsequence region1.Clupea_harengus.fa -outfile region5.aln -gapopen 10.0 -gapextend 0.5
needle -asequence region2.Hypomesus_transpacificus.fa -bsequence region2.Clupea_harengus.fa -outfile region6.aln -gapopen 10.0 -gapextend 0.5
needle -asequence region2.Hypomesus_transpacificus.fa -bsequence region3.Clupea_harengus.fa -outfile region7.aln -gapopen 10.0 -gapextend 0.5
needle -asequence region2.Hypomesus_transpacificus.fa -bsequence region4.Clupea_harengus.fa -outfile region8.aln -gapopen 10.0 -gapextend 0.5


samtools faidx N5.HOG0013514.Clupea_harengus.extended.2.fa.rev NC_045155.1-1092148-1106954:13564-14364 > region5.Clupea_harengus.fa


needle -asequence region1.Hypomesus_transpacificus.fa -bsequence region5.Clupea_harengus.fa -outfile region9.aln -gapopen 10.0 -gapextend 0.5
needle -asequence region2.Hypomesus_transpacificus.fa -bsequence region5.Clupea_harengus.fa -outfile region10.aln -gapopen 10.0 -gapextend 0.5



## blastn the conserved region against all genomes


blastn -query region1.Hypomesus_transpacificus.fa -db ../../Concatenated_Genomes_assemblies/Concatenated_assemblies.fa -outfmt 6 -out regionA.blastn.tsv -num_threads 10 -evalue 1e-20




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
sbatch --qos=1day -c 4 --mem=10G -e error.fitMG4.out -o slurm.fitMG4.out launch_fitMG4.sh N5.HOG0013514


Recipient branches : 
Node360
Hypomesus transpacificus rna XM 047045767 1
Hypomesus transpacificus rna XM 047045786 1
Osmerus_eperlanus_rna_XM_062468699_1
Node371
Osmerus eperlanus rna XM 062468715 1
Osmerus eperlanus rna XM 062468701 1




#Test positive selection and relaxed selection on receiver branch ==> put $HOG.prot.aln.treefile on https://phylotree.hyphy.org/ ==> exort to $HOG.prot.aln.treefile.SelecMarked

sbatch --qos=1day -c 4 --mem=10G -e error.absrel.cand.out -o slurm.absrel.cand.out --job-name=HOG0013514 launch_absrel_cand.sh N5.HOG0013514



#Launch RELAX

### Adaptive branch site random effects likelihood test 
#Likelihood ratio test for episodic diversifying positive selection at Holm-Bonferroni corrected _p =   0.0500_ found **2** branches under selection among **7** tested.
#* Node371, p-value =  0.00086
#* Hypomesus_transpacificus_rna_XM_047045786_1, p-value =  0.01603

===> $HOG.prot.aln.treefile.SelecMarked.relax => select branches detected with accelerated evolution
sbatch --qos=1week -c 6 --mem=10G -e error.relax.out -o slurm.relax.out --job-name=HOG0013514 launch_RELAX.sh N5.HOG0013514




#Extract dN/dS to table

grep "LB\":" N5.HOG0013514.cds.aln.FITTER.json | sed 's/.*://g' | sed 's/,.*//g' > curr_LB_values.txt
grep "MLE\":" N5.HOG0013514.cds.aln.FITTER.json | sed 's/.*://g' | sed 's/,.*//g' > curr_MLE_values.txt
grep "UB\":" N5.HOG0013514.cds.aln.FITTER.json | sed 's/.*://g' | sed 's/,.*//g' > curr_UB_values.txt
grep "\"dN\"" N5.HOG0013514.cds.aln.FITTER.json  | grep -v "{" | sed 's/.*://g' | sed 's/,//g' > curr_dN_values.txt
grep "\"dS\"" N5.HOG0013514.cds.aln.FITTER.json  | grep -v "{" | sed 's/.*://g' | sed 's/,//g' > curr_dS_values.txt
grep -B2 "LB\":" N5.HOG0013514.cds.aln.FITTER.json | grep -v "\-\-" | grep -v "Confidence" | grep -v "LB\":"  | sed 's/\"//g' | sed 's/:.*//g' | sed 's/^ *//g' > curr_labels

paste -d "," curr_labels curr_LB_values.txt curr_MLE_values.txt curr_UB_values.txt curr_dN_values.txt curr_dS_values.txt > N5.HOG0013514.dN_dS.csv




