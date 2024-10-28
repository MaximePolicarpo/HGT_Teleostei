conda activate miniprot


LMOD_DISABLE_SAME_NAME_AUTOSWAP=no

module purge
module load R/4.2.0-foss-2021a
module load BLAST/2.12.0-Linux_x86_64
module load SAMtools/1.15-GCC-10.3.0




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


samtools faidx Proteomes_BUSCO80/Hypomesus_transpacificus.fa Hypomesus_transpacificus---rna-XM_047033889.1 > Hypomesus_transpacificus---rna-XM_047033889.1.prot

#blastp against the fish proteome (50 best hits) + Uniprot database without any fish species (10 best hits)

blastp -query Hypomesus_transpacificus---rna-XM_047033889.1.prot -db concatenated_proteomes.fa -outfmt 6 -out Gene_vs_FishProteome.blastp -num_threads 8 -max_target_seqs 51
blastp -query Hypomesus_transpacificus---rna-XM_047033889.1.prot -db non_actino_uniprot.fa -outfmt 6 -out Gene_vs_Uniprot.blastp -num_threads 8  -max_target_seqs 10 -evalue 1e-5


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
##========================================== Now make a synteny plot  ================================================================
##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================


#N5.HOG0047987.clade_hgt.txt=> clade where the HGT is extracted from iTOL (from N5.HOG0047987.prot.aln.treefile)

cut -f1,2 -d "_" N5.HOG0047987.clade_hgt.txt | sort | uniq  > species_to_draw.clade1
#Choose close species without the gene .. 



echo "Engraulis_encrasicolus" > non_transfer_species_1
echo "Borostomias_antarcticus" > non_transfer_species_2
echo "Gadus_morhua" >> non_transfer_species_2



rm -r Syn_tables_dir ; mkdir Syn_tables_dir


#First order the species names 

rm species_to_draw.clade1.ordered

echo "Alosa_alosa" > species_to_draw.clade1.ordered
echo "Sardina_pilchardus" >> species_to_draw.clade1.ordered
echo "Clupea_harengus" >> species_to_draw.clade1.ordered
echo "Hypomesus_transpacificus" >> species_to_draw.clade1.ordered
echo "Osmerus_eperlanus" >> species_to_draw.clade1.ordered


#If there are more than 4 species from the same order, keep only 4 


#First extract genes around 

curr_OGG=N5.HOG0047987

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




#Now add Engraulis_encrasicolus

curr_OGG=N5.HOG0047987
curr_sp=Engraulis_encrasicolus
ref_sp=Alosa_alosa

rm Syn_tables_dir/$curr_sp.synt.df
#grep -A10 -B10 "N5_HOG0025113" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "NC_085875_1" > Syn_tables_dir/$curr_sp.synt.df

grep -A10 -B10 "N5_HOG0023665" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df
grep -A20 -B10 "N5_HOG0004644" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "NC_085875_1" >> Syn_tables_dir/$curr_sp.synt.df


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


#Now add the second scaffold of Osmerus

curr_OGG=N5.HOG0047987
curr_sp=Osmerus_eperlanus
ref_sp=Hypomesus_transpacificus

grep -A10 -B10 "N5_HOG0035715" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret  >> Syn_tables_dir/$curr_sp.synt.df


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

curr_OGG=N5.HOG0047987
curr_sp=Borostomias_antarcticus
ref_sp=Hypomesus_transpacificus

grep -A10 -B10 "N5_HOG0035715" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret  > Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0043866" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df


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



#Now add Gadus_morhua



curr_OGG=N5.HOG0047987
curr_sp=Gadus_morhua
ref_sp=Hypomesus_transpacificus

grep -A10 -B10 "N5_HOG0035715" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret  > Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0043866" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df


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

#We will take 5kb upstream and downstream of the gene OR gene cluster (or less if there is an other gene before)


#Hypomesus_transpacificus
grep -A2 -B2 "N5_HOG0047987"  GFF3_N5_OGGs/Hypomesus_transpacificus.gff.simplified.sorted.OGG.tiret
samtools faidx GCF_021917145.1_fHypTra1_genomic.fna NC_061064.1:5562537-5571517 > N5.HOG0047987.Hypomesus_transpacificus.extended.1.fa
samtools faidx GCF_021917145.1_fHypTra1_genomic.fna NC_061073.1:4640619-4650388 > N5.HOG0047987.Hypomesus_transpacificus.extended.2.fa

sed -i 's/:/-/g' N5.HOG0047987.Hypomesus_transpacificus.extended.1.fa
sed -i 's/:/-/g' N5.HOG0047987.Hypomesus_transpacificus.extended.2.fa

makeblastdb -in N5.HOG0047987.Hypomesus_transpacificus.extended.1.fa -dbtype nucl
makeblastdb -in N5.HOG0047987.Hypomesus_transpacificus.extended.2.fa -dbtype nucl


#Clupea_harengus
grep -A2 -B2 "N5_HOG0047987"  GFF3_N5_OGGs/Clupea_harengus.gff.simplified.sorted.OGG.tiret
samtools faidx GCF_900700415.2_Ch_v2.0.2_genomic.fna NC_045165.1:18002350-18021184 > N5.HOG0047987.Clupea_harengus.extended.1.fa
sed -i 's/:/-/g' N5.HOG0047987.Clupea_harengus.extended.1.fa
makeblastdb -in N5.HOG0047987.Clupea_harengus.extended.1.fa -dbtype nucl


#Clupea_harengus
grep -A2 -B2 "N5_HOG0047987"  GFF3_N5_OGGs/Osmerus_eperlanus.gff.simplified.sorted.OGG.tiret
samtools faidx GCF_963692335.1_fOsmEpe2.1_genomic.fna NC_085036.1:6379548-6386677 > N5.HOG0047987.Osmerus_eperlanus.extended.1.fa
sed -i 's/:/-/g' N5.HOG0047987.Osmerus_eperlanus.extended.1.fa
makeblastdb -in N5.HOG0047987.Osmerus_eperlanus.extended.1.fa -dbtype nucl





#Now lets find TEs using the repbase and dfam databases

tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0047987.Hypomesus_transpacificus.extended.1.fa -evalue 1e-5 -outfmt 6 -out TE.Hypomesus_transpacificus.1.tblastn -num_threads 8
sed -i 's/#//g' TE.Hypomesus_transpacificus.1.tblastn

tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0047987.Hypomesus_transpacificus.extended.2.fa -evalue 1e-5 -outfmt 6 -out TE.Hypomesus_transpacificus.2.tblastn -num_threads 8
sed -i 's/#//g' TE.Hypomesus_transpacificus.2.tblastn

tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0047987.Clupea_harengus.extended.1.fa -evalue 1e-5 -outfmt 6 -out TE.Clupea_harengus.1.tblastn -num_threads 8
sed -i 's/#//g' TE.Clupea_harengus.1.tblastn

tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0047987.Osmerus_eperlanus.extended.1.fa -evalue 1e-5 -outfmt 6 -out TE.Osmerus_eperlanus.1.tblastn -num_threads 8
sed -i 's/#//g' TE.Osmerus_eperlanus.1.tblastn


#merge tblastn hits and find the best TE match by doing a blastx

Rscript Rscript_merge_blast_hits.R TE.Hypomesus_transpacificus.1.tblastn TE.Hypomesus_transpacificus.1.tblastn.merged
Rscript Rscript_merge_blast_hits.R TE.Hypomesus_transpacificus.2.tblastn TE.Hypomesus_transpacificus.2.tblastn.merged
Rscript Rscript_merge_blast_hits.R TE.Clupea_harengus.1.tblastn TE.Clupea_harengus.1.tblastn.merged
Rscript Rscript_merge_blast_hits.R TE.Osmerus_eperlanus.1.tblastn TE.Osmerus_eperlanus.1.tblastn.merged


xargs samtools faidx N5.HOG0047987.Hypomesus_transpacificus.extended.1.fa < TE.Hypomesus_transpacificus.1.tblastn.merged > TE.Hypomesus_transpacificus.1.BEST.fa
blastx -query TE.Hypomesus_transpacificus.1.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Hypomesus_transpacificus.1.BEST.blastx -max_target_seqs 1
cut -f1  TE.Hypomesus_transpacificus.1.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Hypomesus_transpacificus.1.BEST.blastx >> temp  ; done ; mv temp TE.Hypomesus_transpacificus.1.BEST.blastx

xargs samtools faidx N5.HOG0047987.Hypomesus_transpacificus.extended.2.fa < TE.Hypomesus_transpacificus.2.tblastn.merged > TE.Hypomesus_transpacificus.2.BEST.fa
blastx -query TE.Hypomesus_transpacificus.2.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Hypomesus_transpacificus.2.BEST.blastx -max_target_seqs 1
cut -f1  TE.Hypomesus_transpacificus.2.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Hypomesus_transpacificus.2.BEST.blastx >> temp  ; done ; mv temp TE.Hypomesus_transpacificus.2.BEST.blastx


xargs samtools faidx N5.HOG0047987.Clupea_harengus.extended.1.fa < TE.Clupea_harengus.1.tblastn.merged > TE.Clupea_harengus.1.BEST.fa
blastx -query TE.Clupea_harengus.1.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Clupea_harengus.1.BEST.blastx -max_target_seqs 1
cut -f1  TE.Clupea_harengus.1.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Clupea_harengus.1.BEST.blastx >> temp  ; done ; mv temp TE.Clupea_harengus.1.BEST.blastx

xargs samtools faidx N5.HOG0047987.Osmerus_eperlanus.extended.1.fa < TE.Osmerus_eperlanus.1.tblastn.merged > TE.Osmerus_eperlanus.1.BEST.fa
blastx -query TE.Osmerus_eperlanus.1.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Osmerus_eperlanus.1.BEST.blastx -max_target_seqs 1
cut -f1  TE.Osmerus_eperlanus.1.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Osmerus_eperlanus.1.BEST.blastx >> temp  ; done ; mv temp TE.Osmerus_eperlanus.1.BEST.blastx





#Now find shared elements

cut -f2 TE.Hypomesus_transpacificus.1.BEST.blastx | sort | uniq > TE.Hypomesus_transpacificus.1.uniqTE
cut -f2 TE.Hypomesus_transpacificus.2.BEST.blastx | sort | uniq > TE.Hypomesus_transpacificus.2.uniqTE
cut -f2 TE.Clupea_harengus.1.BEST.blastx | sort | uniq > TE.Clupea_harengus.1.uniqTE

comm -12 TE.Hypomesus_transpacificus.1.uniqTE TE.Clupea_harengus.1.uniqTE
comm -12 TE.Hypomesus_transpacificus.2.uniqTE TE.Clupea_harengus.1.uniqTE


#There is one shared element : Togen-1_DR_1p:ClassI:?:?:?

samtools faidx Dfam_plus_Repbase.cdhit80.prot Togen-1_DR_1p:ClassI:?:?:? > Togen-1.prot


#### Now let's find every copy of this element in the genome of both species + closely related species


#Hypomesus_transpacificus
tblastn -query Togen-1.prot -db GCF_021917145.1_fHypTra1_genomic.fna -evalue 1e-5 -outfmt 6 -out Togen-1.Hypomesus_transpacificus.tblastn -num_threads 8
sed -i 's/#//g' Togen-1.Hypomesus_transpacificus.tblastn
Rscript Rscript_merge_blast_hits.R Togen-1.Hypomesus_transpacificus.tblastn Togen-1.Hypomesus_transpacificus.tblastn.merged
xargs samtools faidx GCF_021917145.1_fHypTra1_genomic.fna < Togen-1.Hypomesus_transpacificus.tblastn.merged > Togen-1.Hypomesus_transpacificus.BEST.fa
diamond blastx --query Togen-1.Hypomesus_transpacificus.BEST.fa --db Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out Togen-1.Hypomesus_transpacificus.BEST.blastx --max-target-seqs 1
grep "Togen-1" Togen-1.Hypomesus_transpacificus.BEST.blastx > Togen-1.Hypomesus_transpacificus.list

#Clupea_harengus
tblastn -query Togen-1.prot -db GCF_900700415.2_Ch_v2.0.2_genomic.fna -evalue 1e-5 -outfmt 6 -out Togen-1.Clupea_harengus.tblastn -num_threads 8
sed -i 's/#//g' Togen-1.Clupea_harengus.tblastn
Rscript Rscript_merge_blast_hits.R Togen-1.Clupea_harengus.tblastn Togen-1.Clupea_harengus.tblastn.merged
xargs samtools faidx GCF_900700415.2_Ch_v2.0.2_genomic.fna  < Togen-1.Clupea_harengus.tblastn.merged > Togen-1.Clupea_harengus.BEST.fa
diamond blastx --query Togen-1.Clupea_harengus.BEST.fa --db Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out Togen-1.Clupea_harengus.BEST.blastx --max-target-seqs 1
grep "Togen-1" Togen-1.Clupea_harengus.BEST.blastx > Togen-1.Clupea_harengus.list


#Sardina_pilchardus
tblastn -query Togen-1.prot -db GCF_963854185.1_fSarPil1.1_genomic.fna -evalue 1e-5 -outfmt 6 -out Togen-1.Sardina_pilchardus.tblastn -num_threads 8
sed -i 's/#//g' Togen-1.Sardina_pilchardus.tblastn
Rscript Rscript_merge_blast_hits.R Togen-1.Sardina_pilchardus.tblastn Togen-1.Sardina_pilchardus.tblastn.merged
xargs samtools faidx GCF_963854185.1_fSarPil1.1_genomic.fna  < Togen-1.Sardina_pilchardus.tblastn.merged > Togen-1.Sardina_pilchardus.BEST.fa
diamond blastx --query Togen-1.Sardina_pilchardus.BEST.fa --db Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out Togen-1.Sardina_pilchardus.BEST.blastx --max-target-seqs 1
grep "Togen-1" Togen-1.Sardina_pilchardus.BEST.blastx > Togen-1.Sardina_pilchardus.list

#Alosa_alosa
tblastn -query Togen-1.prot -db GCF_017589495.1_AALO_Geno_1.1_genomic.fna -evalue 1e-5 -outfmt 6 -out Togen-1.Alosa_alosa.tblastn -num_threads 8
sed -i 's/#//g' Togen-1.Alosa_alosa.tblastn
Rscript Rscript_merge_blast_hits.R Togen-1.Alosa_alosa.tblastn Togen-1.Alosa_alosa.tblastn.merged
xargs samtools faidx GCF_017589495.1_AALO_Geno_1.1_genomic.fna  < Togen-1.Alosa_alosa.tblastn.merged > Togen-1.Alosa_alosa.BEST.fa
diamond blastx --query Togen-1.Alosa_alosa.BEST.fa --db Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out Togen-1.Alosa_alosa.BEST.blastx --max-target-seqs 1
grep "Togen-1" Togen-1.Alosa_alosa.BEST.blastx > Togen-1.Alosa_alosa.list


#Engraulis_encrasicolus
tblastn -query Togen-1.prot -db GCF_034702125.1_IST_EnEncr_1.0_genomic.fna -evalue 1e-5 -outfmt 6 -out Togen-1.Engraulis_encrasicolus.tblastn -num_threads 8
sed -i 's/#//g' Togen-1.Engraulis_encrasicolus.tblastn
Rscript Rscript_merge_blast_hits.R Togen-1.Engraulis_encrasicolus.tblastn Togen-1.Engraulis_encrasicolus.tblastn.merged
xargs samtools faidx GCF_034702125.1_IST_EnEncr_1.0_genomic.fna  < Togen-1.Engraulis_encrasicolus.tblastn.merged > Togen-1.Engraulis_encrasicolus.BEST.fa
diamond blastx --query Togen-1.Engraulis_encrasicolus.BEST.fa --db Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out Togen-1.Engraulis_encrasicolus.BEST.blastx --max-target-seqs 1
grep "Togen-1" Togen-1.Engraulis_encrasicolus.BEST.blastx > Togen-1.Engraulis_encrasicolus.list


#Osmerus_eperlanus
tblastn -query Togen-1.prot -db GCF_963692335.1_fOsmEpe2.1_genomic.fna -evalue 1e-5 -outfmt 6 -out Togen-1.Osmerus_eperlanus.tblastn -num_threads 8
sed -i 's/#//g' Togen-1.Osmerus_eperlanus.tblastn
Rscript Rscript_merge_blast_hits.R Togen-1.Osmerus_eperlanus.tblastn Togen-1.Osmerus_eperlanus.tblastn.merged
xargs samtools faidx GCF_963692335.1_fOsmEpe2.1_genomic.fna  < Togen-1.Osmerus_eperlanus.tblastn.merged > Togen-1.Osmerus_eperlanus.BEST.fa
diamond blastx --query Togen-1.Osmerus_eperlanus.BEST.fa --db Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out Togen-1.Osmerus_eperlanus.BEST.blastx --max-target-seqs 1
grep "Togen-1" Togen-1.Osmerus_eperlanus.BEST.blastx > Togen-1.Osmerus_eperlanus.list

#Borostomias_antarcticus
tblastn -query Togen-1.prot -db ../../Genomic_data/Borostomias_antarcticus/GCA_949987555.1_fBorAnt1.1_genomic.fna -evalue 1e-5 -outfmt 6 -out Togen-1.Borostomias_antarcticus.tblastn -num_threads 8
sed -i 's/#//g' Togen-1.Borostomias_antarcticus.tblastn
Rscript Rscript_merge_blast_hits.R Togen-1.Borostomias_antarcticus.tblastn Togen-1.Borostomias_antarcticus.tblastn.merged
xargs samtools faidx ../../Genomic_data/Borostomias_antarcticus/GCA_949987555.1_fBorAnt1.1_genomic.fna  < Togen-1.Borostomias_antarcticus.tblastn.merged > Togen-1.Borostomias_antarcticus.BEST.fa
diamond blastx --query Togen-1.Borostomias_antarcticus.BEST.fa --db Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out Togen-1.Borostomias_antarcticus.BEST.blastx --max-target-seqs 1
grep "Togen-1" Togen-1.Borostomias_antarcticus.BEST.blastx > Togen-1.Borostomias_antarcticus.list

#Gadus_morhua
tblastn -query Togen-1.prot -db ../../Genomic_data/Gadus_morhua/GCF_902167405.1_gadMor3.0_genomic.fna -evalue 1e-5 -outfmt 6 -out Togen-1.Gadus_morhua.tblastn -num_threads 8
sed -i 's/#//g' Togen-1.Gadus_morhua.tblastn
Rscript Rscript_merge_blast_hits.R Togen-1.Gadus_morhua.tblastn Togen-1.Gadus_morhua.tblastn.merged
xargs samtools faidx ../../Genomic_data/Gadus_morhua/GCF_902167405.1_gadMor3.0_genomic.fna  < Togen-1.Gadus_morhua.tblastn.merged > Togen-1.Gadus_morhua.BEST.fa
diamond blastx --query Togen-1.Gadus_morhua.BEST.fa --db Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out Togen-1.Gadus_morhua.BEST.blastx --max-target-seqs 1
grep "Togen-1" Togen-1.Gadus_morhua.BEST.blastx > Togen-1.Gadus_morhua.list


wc -l Togen-1.Engraulis_encrasicolus.list
wc -l Togen-1.Alosa_alosa.list
wc -l Togen-1.Sardina_pilchardus.list
wc -l Togen-1.Clupea_harengus.list
wc -l Togen-1.Hypomesus_transpacificus.list
wc -l Togen-1.Osmerus_eperlanus.list
wc -l Togen-1.Borostomias_antarcticus.list
wc -l Togen-1.Gadus_morhua.list




### In a FINAL step we will make another synteny plot, only between Clupea and Hypomesus to show N5.HOG0001647 and neighboring regions 


## Clupea harengus 

rm seq_clustered_infos_ogg.TE.txt ; rm clusters_ID_TE.txt


scaffold=`grep ">" N5.HOG0047987.Clupea_harengus.extended.1.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0047987.Clupea_harengus.extended.1.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0047987.Clupea_harengus.extended.1.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Clupea_harengus,$scaffold,$length" > clusters_ID_TE.txt

grep "N5_HOG0047987"  GFF3_N5_OGGs/Clupea_harengus.gff.simplified.sorted.OGG.tiret
grep "XM_031580274" GFF3_files_per_species/Clupea_harengus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons

nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene1.$nbexon,N5_HOG0047987,+,Clupea_harengus" >> seq_clustered_infos_ogg.TE.txt
done


cut -f1 TE.Clupea_harengus.1.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Clupea_harengus.1.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0047987.Clupea_harengus.extended.1.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0047987.Clupea_harengus.extended.1.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold,/g" | sed "s/$/,$strand,Clupea_harengus/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done






## Hypomesus_transpacificus 


grep ">" N5.HOG0047987.Hypomesus_transpacificus.extended.1.fa
scaffold=`grep ">" N5.HOG0047987.Hypomesus_transpacificus.extended.1.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0047987.Hypomesus_transpacificus.extended.1.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0047987.Hypomesus_transpacificus.extended.1.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Hypomesus_transpacificus,$scaffold,$length" >> clusters_ID_TE.txt


grep "N5_HOG0047987"  GFF3_N5_OGGs/Hypomesus_transpacificus.gff.simplified.sorted.OGG.tiret
grep "XM_047020147" GFF3_files_per_species/Hypomesus_transpacificus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons

nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,exon.$nbexon,N5_HOG0047987,+,Hypomesus_transpacificus" >> seq_clustered_infos_ogg.TE.txt
done


cut -f1 TE.Hypomesus_transpacificus.1.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Hypomesus_transpacificus.1.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0047987.Hypomesus_transpacificus.extended.1.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0047987.Hypomesus_transpacificus.extended.1.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold,/g" | sed "s/$/,$strand,Hypomesus_transpacificus/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done





## Second scaffold

grep ">" N5.HOG0047987.Hypomesus_transpacificus.extended.2.fa
scaffold=`grep ">" N5.HOG0047987.Hypomesus_transpacificus.extended.2.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0047987.Hypomesus_transpacificus.extended.2.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0047987.Hypomesus_transpacificus.extended.2.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Hypomesus_transpacificus,$scaffold,$length" >> clusters_ID_TE.txt


grep "N5_HOG0047987"  GFF3_N5_OGGs/Hypomesus_transpacificus.gff.simplified.sorted.OGG.tiret
grep "XM_047033889" GFF3_files_per_species/Hypomesus_transpacificus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons

nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,exon.$nbexon,N5_HOG0047987,-,Hypomesus_transpacificus" >> seq_clustered_infos_ogg.TE.txt
done


cut -f1 TE.Hypomesus_transpacificus.2.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Hypomesus_transpacificus.2.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0047987.Hypomesus_transpacificus.extended.2.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0047987.Hypomesus_transpacificus.extended.2.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold,/g" | sed "s/$/,$strand,Hypomesus_transpacificus/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done



### Check the conserved regions



samtools faidx N5.HOG0047987.Hypomesus_transpacificus.extended.2.fa.rev NC_061073.1-4640619-4650388:6703-7600 > region1.Hypomesus_transpacificus.fa
samtools faidx N5.HOG0047987.Clupea_harengus.extended.1.fa NC_045165.1-18002350-18021184:3790-5300 > region1.Clupea_harengus.fa
samtools faidx N5.HOG0047987.Clupea_harengus.extended.1.fa NC_045165.1-18002350-18021184:12390-13900 > region2.Clupea_harengus.fa

sed -i 's/>.*/>H.transpacificus/g' region1.Hypomesus_transpacificus.fa

sed -i 's/>.*/>C.harengus/g' region1.Clupea_harengus.fa
sed -i 's/>.*/>C.harengus/g' region2.Clupea_harengus.fa

needle -asequence region1.Hypomesus_transpacificus.fa -bsequence region1.Clupea_harengus.fa -outfile region1.aln -gapopen 10.0 -gapextend 0.5
needle -asequence region1.Hypomesus_transpacificus.fa -bsequence region2.Clupea_harengus.fa -outfile region2.aln -gapopen 10.0 -gapextend 0.5

## blastn the conserved region against all genomes

blastn -query region1.Hypomesus_transpacificus.fa -db Concatenated_assemblies.fa -outfmt 6 -out region1.blastn.tsv -num_threads 10 -evalue 1e-20






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



#Receiver_sequence.txt: 
Hypomesus_transpacificus_rna_XM_047020147_1
Osmerus_eperlanus_rna_XM_062485952_1
Hypomesus_transpacificus_rna_XM_047033889_1

#Compute the dN/dS on every branches
sbatch --qos=1day -c 4 --mem=10G -e error.fitMG4.out -o slurm.fitMG4.out launch_fitMG4.sh N5.HOG0047987


Recipient branches : 
Node105
Node107 
Hypomesus_transpacificus_rna_XM_047033889_1
Hypomesus_transpacificus_rna_XM_047020147_1
Osmerus_eperlanus_rna_XM_062485952_1


#Test positive selection and relaxed selection on receiver branch ==> put $HOG.prot.aln.treefile on https://phylotree.hyphy.org/ ==> exort to $HOG.prot.aln.treefile.SelecMarked

sbatch --qos=1day -c 4 --mem=10G -e error.absrel.cand.out -o slurm.absrel.cand.out --job-name=HOG0047987 launch_absrel_cand.sh N5.HOG0047987



#Launch RELAX

### Adaptive branch site random effects likelihood test 
#Likelihood ratio test for episodic diversifying positive selection at Holm-Bonferroni corrected _p =   0.0500_ found **2** branches under selection among **5** tested.
#* Osmerus_eperlanus_rna_XM_062485952_1, p-value =  0.00002
#* Hypomesus_transpacificus_rna_XM_047033889_1, p-value =  0.00005

===> $HOG.prot.aln.treefile.SelecMarked.relax => select branches detected with accelerated evolution
sbatch --qos=1week -c 6 --mem=10G -e error.relax.out -o slurm.relax.out --job-name=HOG0047987 launch_RELAX.sh N5.HOG0047987


#Extract dN/dS to table

grep "LB\":" N5.HOG0047987.cds.aln.FITTER.json | sed 's/.*://g' | sed 's/,.*//g' > curr_LB_values.txt
grep "MLE\":" N5.HOG0047987.cds.aln.FITTER.json | sed 's/.*://g' | sed 's/,.*//g' > curr_MLE_values.txt
grep "UB\":" N5.HOG0047987.cds.aln.FITTER.json | sed 's/.*://g' | sed 's/,.*//g' > curr_UB_values.txt
grep "\"dN\"" N5.HOG0047987.cds.aln.FITTER.json  | grep -v "{" | sed 's/.*://g' | sed 's/,//g' > curr_dN_values.txt
grep "\"dS\"" N5.HOG0047987.cds.aln.FITTER.json  | grep -v "{" | sed 's/.*://g' | sed 's/,//g' > curr_dS_values.txt
grep -B2 "LB\":" N5.HOG0047987.cds.aln.FITTER.json | grep -v "\-\-" | grep -v "Confidence" | grep -v "LB\":"  | sed 's/\"//g' | sed 's/:.*//g' | sed 's/^ *//g' > curr_labels

paste -d "," curr_labels curr_LB_values.txt curr_MLE_values.txt curr_UB_values.txt curr_dN_values.txt curr_dS_values.txt > N5.HOG0047987.dN_dS.csv









