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



samtools faidx Proteomes_BUSCO80/Hypomesus_transpacificus.fa Hypomesus_transpacificus---rna-XM_047031800.1 > Hypomesus_transpacificus---rna-XM_047031800.1.prot

grep -v "Hypomesus_transpacificus---rna-XM_047021690.1" Gene_vs_FishProteome.blastp  > temp ; mv temp Gene_vs_FishProteome.blastp
#blastp against the fish proteome (50 best hits) + Uniprot database without any fish species (10 best hits)

blastp -query Hypomesus_transpacificus---rna-XM_047031800.1.prot -db concatenated_proteomes.fa -outfmt 6 -out Gene_vs_FishProteome.blastp -num_threads 8 -max_target_seqs 51
blastp -query Hypomesus_transpacificus---rna-XM_047031800.1.prot -db non_actino_uniprot.fa -outfmt 6 -out Gene_vs_Uniprot.blastp -num_threads 8  -max_target_seqs 10 -evalue 1e-5


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


#N5.HOG0028899.clade_hgt.txt=> clade where the HGT is extracted from iTOL (from N5.HOG0028899.prot.aln.treefile)

cut -f1,2 -d "_" N5.HOG0028899.clade_hgt.txt | sort | uniq  > species_to_draw.clade1
#Choose close species without the gene .. 



echo "Engraulis_encrasicolus" > non_transfer_species_1
echo "Alosa_alosa" >> non_transfer_species_1
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

curr_OGG=N5.HOG0028899

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


#Add Sardina because the location of the OGG is not in synteny with Clupea


curr_OGG=N5.HOG0028899
curr_sp=Sardina_pilchardus
ref_sp=Clupea_harengus

rm Syn_tables_dir/$curr_sp.synt.df

grep -A10 -B10 "N5_HOG0030030" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0028899" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0044165" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df



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



#Add Alosa_saloa because the location of the OGG is not in synteny with Clupea


curr_OGG=N5.HOG0028899
curr_sp=Alosa_alosa
ref_sp=Clupea_harengus

rm Syn_tables_dir/$curr_sp.synt.df

grep -A10 -B10 "N5_HOG0028899" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df



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

curr_OGG=N5.HOG0028899
curr_sp=Engraulis_encrasicolus
ref_sp=Alosa_alosa

rm Syn_tables_dir/$curr_sp.synt.df

grep -A10 -B10 "N5_HOG0018734" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0028899" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0018734" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret >> Syn_tables_dir/$curr_sp.synt.df


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


#Now add Osmerus eperlanus

curr_OGG=N5.HOG0028899
curr_sp=Osmerus_eperlanus
ref_sp=Hypomesus_transpacificus

grep -A10 -B10 "N5_HOG0014375" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df




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


echo "NC_085029_1,16932918,16937143,custom_annot,N5_HOG0028899,+,Osmerus_eperlanus" >> Syn_tables_dir/$curr_sp.synt.final.df




#Now add Borostomias_antarcticus

curr_OGG=N5.HOG0028899
curr_sp=Borostomias_antarcticus
ref_sp=Hypomesus_transpacificus

grep -A10 -B10 "N5_HOG0039668" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret  > Syn_tables_dir/$curr_sp.synt.df


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



curr_OGG=N5.HOG0028899
curr_sp=Gadus_morhua
ref_sp=Hypomesus_transpacificus

grep -A10 -B10 "N5_HOG0014375" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret  > Syn_tables_dir/$curr_sp.synt.df


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
##========================================== Find the gene in Osmerus  =============================================================
##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================


grep "N5.HOG0014375"  GFF3_N5_OGGs/Osmerus_eperlanus.gff.simplified.sorted.OGG.tiret
grep "N5.HOG0039668"  GFF3_N5_OGGs/Osmerus_eperlanus.gff.simplified.sorted.OGG.tiret
samtools faidx GCF_963692335.1_fOsmEpe2.1_genomic.fna NC_085029.1:16907226-16971149 > Osmerus_eperlanus.region1.fa
revseq Osmerus_eperlanus.region1.fa Osmerus_eperlanus.region1.fa.rev

/scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Database/Scripts_opsins/exonerate-2.2.0-x86_64/bin/exonerate --bestn 1 --showtargetgff TRUE --model protein2genome --ryo "%tcs" --maxintron 5000 Hypomesus_transpacificus---rna-XM_047031800.1.prot Osmerus_eperlanus.region1.fa > gene.Osmerus_eperlanus.exo




#Osmerus_found_gene.fa
>Osmerus_eperlanus.N5.HOG0028899
ATGAATAGGCCTTTCACCATGATCCTACTGGTGGTGTGTGCTTTGGCTGCTGTTTCTTTGGCTCTGCCTC
CGGgaagtgacatcactgaAGCTCGGCCTATTTCAGTCAGTGCTTTCCCTACCATATTTGCAGTGTGCTG
GACAAATTGGTACGATCGTGACGATCCTTCCGGTACCGGTGACTGGGAGGTGCTGTCGGACCTGCGGAGT
GAGAACCCTGGTGAGATCTGTGCCAACCCGATTGCCATTGAGAGCAAGACTGTGGACACAGGGACTGCTG
CCGCCCTGACTGGCCAAAACTTCCTCCACTTCTCTCCGACAATTGGGTTCGTATGCAAGAACGACCCCAA
ACAATACTGCCGGGACTACCAAGTCCGGTTTCAGTGTCCTTGTGGGCCTGGTGGGTCTGGTGGGCCTCGT
GTGCCT

grep "Clupea\|Alosa\|Sardina\|Engraulis" Coding_sequences_alignments/N5.HOG0028899.prot | sed 's/>//g' | grep -v "XM_031582176\|XM_031564016" > toextract.id
xargs samtools faidx  Coding_sequences_alignments/N5.HOG0028899.prot < toextract.id > toextract.prot

transeq Osmerus_found_gene.fa Osmerus_found_gene.prot ; sed -i 's/_1$//g' Osmerus_found_gene.prot

cat Hypomesus_transpacificus---rna-XM_047031800.1.prot Osmerus_found_gene.prot toextract.prot > Osmerus_phylogeny.prot


muscle5.1 -align Osmerus_phylogeny.prot -output Osmerus_phylogeny.aln

trimal -in Osmerus_phylogeny.aln -gt 0.7 -out Osmerus_phylogeny.aln.trimmed

iqtree -s Osmerus_phylogeny.aln -st AA -nt 8 -m TEST -mrate G4 --redo



#Likely a pseudogene. Much shorter than the Clupea gene version



samtools faidx Proteomes_BUSCO80/Clupea_harengus.fa Clupea_harengus---rna-XM_042709989.1 > Clupea_harengus---rna-XM_042709989.1.prot
samtools faidx Proteomes_BUSCO80/Clupea_harengus.fa Clupea_harengus---rna-XM_031581960.1 > Clupea_harengus---rna-XM_031581960.1.prot

exonerate --bestn 1 --showtargetgff TRUE --model protein2genome --ryo "%tcs" --maxintron 5000 Clupea_harengus---rna-XM_031581960.1.prot Hypomesus_transpacificus.region1.fa
exonerate --bestn 1 --showtargetgff TRUE --model protein2genome --ryo "%tcs" --maxintron 5000 Clupea_harengus---rna-XM_042709989.1.prot Osmerus_eperlanus.region1.fa



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
grep -A2 -B2 "N5_HOG0028899"  GFF3_N5_OGGs/Hypomesus_transpacificus.gff.simplified.sorted.OGG.tiret
samtools faidx GCF_021917145.1_fHypTra1_genomic.fna NC_061072.1:7354826-7366858 > N5.HOG0028899.Hypomesus_transpacificus.extended.1.fa
sed -i 's/:/-/g' N5.HOG0028899.Hypomesus_transpacificus.extended.1.fa

makeblastdb -in N5.HOG0028899.Hypomesus_transpacificus.extended.1.fa -dbtype nucl


#Clupea_harengus
grep -A2 -B2 "N5_HOG0028899"  GFF3_N5_OGGs/Clupea_harengus.gff.simplified.sorted.OGG.tiret
samtools faidx GCF_900700415.2_Ch_v2.0.2_genomic.fna NC_045166.1:16191554-16269185 > N5.HOG0028899.Clupea_harengus.extended.1.fa
samtools faidx GCF_900700415.2_Ch_v2.0.2_genomic.fna NC_045177.1:8606891-8619407 > N5.HOG0028899.Clupea_harengus.extended.2.fa

sed -i 's/:/-/g' N5.HOG0028899.Clupea_harengus.extended.1.fa
makeblastdb -in N5.HOG0028899.Clupea_harengus.extended.1.fa -dbtype nucl

sed -i 's/:/-/g' N5.HOG0028899.Clupea_harengus.extended.2.fa
makeblastdb -in N5.HOG0028899.Clupea_harengus.extended.2.fa -dbtype nucl




#Now lets find TEs using the repbase and dfam databases

tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0028899.Hypomesus_transpacificus.extended.1.fa -evalue 1e-5 -outfmt 6 -out TE.Hypomesus_transpacificus.1.tblastn -num_threads 8
sed -i 's/#//g' TE.Hypomesus_transpacificus.1.tblastn


tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0028899.Clupea_harengus.extended.1.fa -evalue 1e-5 -outfmt 6 -out TE.Clupea_harengus.1.tblastn -num_threads 8
sed -i 's/#//g' TE.Clupea_harengus.1.tblastn

tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0028899.Clupea_harengus.extended.2.fa -evalue 1e-5 -outfmt 6 -out TE.Clupea_harengus.2.tblastn -num_threads 8
sed -i 's/#//g' TE.Clupea_harengus.2.tblastn


#merge tblastn hits and find the best TE match by doing a blastx

Rscript Rscript_merge_blast_hits.R TE.Hypomesus_transpacificus.1.tblastn TE.Hypomesus_transpacificus.1.tblastn.merged
Rscript Rscript_merge_blast_hits.R TE.Clupea_harengus.1.tblastn TE.Clupea_harengus.1.tblastn.merged
Rscript Rscript_merge_blast_hits.R TE.Clupea_harengus.2.tblastn TE.Clupea_harengus.2.tblastn.merged


xargs samtools faidx N5.HOG0028899.Hypomesus_transpacificus.extended.1.fa < TE.Hypomesus_transpacificus.1.tblastn.merged > TE.Hypomesus_transpacificus.1.BEST.fa
blastx -query TE.Hypomesus_transpacificus.1.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Hypomesus_transpacificus.1.BEST.blastx -max_target_seqs 1
cut -f1  TE.Hypomesus_transpacificus.1.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Hypomesus_transpacificus.1.BEST.blastx >> temp  ; done ; mv temp TE.Hypomesus_transpacificus.1.BEST.blastx

xargs samtools faidx N5.HOG0028899.Clupea_harengus.extended.1.fa < TE.Clupea_harengus.1.tblastn.merged > TE.Clupea_harengus.1.BEST.fa
blastx -query TE.Clupea_harengus.1.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Clupea_harengus.1.BEST.blastx -max_target_seqs 1
cut -f1  TE.Clupea_harengus.1.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Clupea_harengus.1.BEST.blastx >> temp  ; done ; mv temp TE.Clupea_harengus.1.BEST.blastx


xargs samtools faidx N5.HOG0028899.Clupea_harengus.extended.2.fa < TE.Clupea_harengus.2.tblastn.merged > TE.Clupea_harengus.2.BEST.fa
blastx -query TE.Clupea_harengus.2.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Clupea_harengus.2.BEST.blastx -max_target_seqs 1
cut -f1  TE.Clupea_harengus.2.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Clupea_harengus.2.BEST.blastx >> temp  ; done ; mv temp TE.Clupea_harengus.2.BEST.blastx


#Now find shared elements

cut -f2 TE.Hypomesus_transpacificus.1.BEST.blastx | sort | uniq > TE.Hypomesus_transpacificus.1.uniqTE
cut -f2 TE.Clupea_harengus.1.BEST.blastx | sort | uniq > TE.Clupea_harengus.1.uniqTE
cut -f2 TE.Clupea_harengus.2.BEST.blastx | sort | uniq > TE.Clupea_harengus.2.uniqTE

comm -12 TE.Hypomesus_transpacificus.1.uniqTE TE.Clupea_harengus.1.uniqTE
comm -12 TE.Hypomesus_transpacificus.1.uniqTE TE.Clupea_harengus.2.uniqTE



## No common TE found


#### Lets draw a ZOOM in the regions of Clupea and Hypomesus



## Clupea harengus 

rm seq_clustered_infos_ogg.TE.txt ; rm clusters_ID_TE.txt


scaffold=`grep ">" N5.HOG0028899.Clupea_harengus.extended.1.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0028899.Clupea_harengus.extended.1.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0028899.Clupea_harengus.extended.1.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Clupea_harengus,$scaffold,$length" > clusters_ID_TE.txt

grep "N5_HOG0028899"  GFF3_N5_OGGs/Clupea_harengus.gff.simplified.sorted.OGG.tiret

grep "XM_031582176" GFF3_files_per_species/Clupea_harengus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons
grep "XM_042709989" GFF3_files_per_species/Clupea_harengus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_2.exons
grep "XM_031581960" GFF3_files_per_species/Clupea_harengus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_3.exons


nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene1.$nbexon,N5_HOG0028899,-,Clupea_harengus" >> seq_clustered_infos_ogg.TE.txt
done

nbexon=0
for line in `cat gene_2.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene2.$nbexon,N5_HOG0028899,-,Clupea_harengus" >> seq_clustered_infos_ogg.TE.txt
done
nbexon=0
for line in `cat gene_3.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene3.$nbexon,N5_HOG0028899,+,Clupea_harengus" >> seq_clustered_infos_ogg.TE.txt
done


cut -f1 TE.Clupea_harengus.1.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Clupea_harengus.1.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0028899.Clupea_harengus.extended.1.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0028899.Clupea_harengus.extended.1.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold,/g" | sed "s/$/,$strand,Clupea_harengus/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done



#Second region of Clupea


scaffold=`grep ">" N5.HOG0028899.Clupea_harengus.extended.2.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0028899.Clupea_harengus.extended.2.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0028899.Clupea_harengus.extended.2.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Clupea_harengus,$scaffold.sec,$length" >> clusters_ID_TE.txt

grep "N5_HOG0028899"  GFF3_N5_OGGs/Clupea_harengus.gff.simplified.sorted.OGG.tiret
grep "XM_031564016" GFF3_files_per_species/Clupea_harengus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons

nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold.sec,$real_exon_start,$real_exon_stop,gene4.$nbexon,N5_HOG0028899,-,Clupea_harengus" >> seq_clustered_infos_ogg.TE.txt
done


cut -f1 TE.Clupea_harengus.2.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Clupea_harengus.2.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0028899.Clupea_harengus.extended.2.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0028899.Clupea_harengus.extended.2.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold.sec,/g" | sed "s/$/,$strand,Clupea_harengus/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done




## Hypomesus_transpacificus 


grep ">" N5.HOG0028899.Hypomesus_transpacificus.extended.1.fa
scaffold=`grep ">" N5.HOG0028899.Hypomesus_transpacificus.extended.1.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0028899.Hypomesus_transpacificus.extended.1.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0028899.Hypomesus_transpacificus.extended.1.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Hypomesus_transpacificus,$scaffold,$length" >> clusters_ID_TE.txt


grep "N5_HOG0028899"  GFF3_N5_OGGs/Hypomesus_transpacificus.gff.simplified.sorted.OGG.tiret
grep "XM_047031800" GFF3_files_per_species/Hypomesus_transpacificus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons


nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene1.exon.$nbexon,N5_HOG0028899,-,Hypomesus_transpacificus" >> seq_clustered_infos_ogg.TE.txt
done


cut -f1 TE.Hypomesus_transpacificus.1.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Hypomesus_transpacificus.1.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0028899.Hypomesus_transpacificus.extended.1.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0028899.Hypomesus_transpacificus.extended.1.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold,/g" | sed "s/$/,$strand,Hypomesus_transpacificus/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done



## Check in more detail the conserved non coding region found


header=`grep ">" N5.HOG0028899.Hypomesus_transpacificus.extended.1.fa | sed 's/>//'`
samtools faidx N5.HOG0028899.Hypomesus_transpacificus.extended.1.fa $header:5050-5800 > similarity_region_1.Hypomesus_transpacificus.fa

header=`grep ">" N5.HOG0028899.Clupea_harengus.extended.1.fa | sed 's/>//'`
samtools faidx N5.HOG0028899.Clupea_harengus.extended.1.fa $header:10000-12000 > similarity_region_1.Clupea_harengus.1.fa


header=`grep ">" N5.HOG0028899.Clupea_harengus.extended.1.fa | sed 's/>//'`
samtools faidx N5.HOG0028899.Clupea_harengus.extended.1.fa $header:40000-42000 > similarity_region_1.Clupea_harengus.2.fa


header=`grep ">" N5.HOG0028899.Clupea_harengus.extended.2.fa | sed 's/>//'`
samtools faidx N5.HOG0028899.Clupea_harengus.extended.2.fa $header:8000-9200 > similarity_region_1.Clupea_harengus.3.fa

sed 's/>/>Hypomesus_transpacificus-/g' similarity_region_1.Hypomesus_transpacificus.fa > similarity_region_1.combined.fa

sed 's/>/>Clupea_harengus-/g' similarity_region_1.Clupea_harengus.1.fa  >> similarity_region_1.combined.fa
sed 's/>/>Clupea_harengus-/g' similarity_region_1.Clupea_harengus.2.fa  >> similarity_region_1.combined.fa
sed 's/>/>Clupea_harengus-/g' similarity_region_1.Clupea_harengus.3.fa >> similarity_region_1.combined.fa


trimal -in similarity_region_1.combined.aln -out similarity_region_1.combined.aln.trimmed -gt 0.9

muscle5.1 -align similarity_region_1.combined.fa -output similarity_region_1.combined.aln


iqtree -s similarity_region_1.combined.aln -st DNA -nt 8 -m TEST --redo

iqtree -s similarity_region_1.combined.aln.trimmed -st DNA -nt 8 -m TEST --redo


blastx -query similarity_region_1.Hypomesus_transpacificus.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out similarity_region_1.Hypomesus_transpacificus.TE.blastx 


sed -i 's/:/-/g' similarity_region_1.combined.fa 
samtools faidx similarity_region_1.combined.fa  Hypomesus_transpacificus-NC_061072.1-7354826-7366858-5050-5800 > test1.fa
samtools faidx similarity_region_1.combined.fa  Clupea_harengus-NC_045166.1-16191554-16269185-10000-12000 > test2.fa
samtools faidx similarity_region_1.combined.fa  Clupea_harengus-NC_045166.1-16191554-16269185-40000-42000 > test3.fa
samtools faidx similarity_region_1.combined.fa  Clupea_harengus-NC_045177.1-8606891-8619407-8000-9200 > test4.fa

needle -asequence test1.fa -bsequence test2.fa -outfile test2.aln -gapopen 10.0 -gapextend 0.5 #=> 91.3
needle -asequence test1.fa -bsequence test3.fa -outfile test3.aln -gapopen 10.0 -gapextend 0.5 #=> 90.5
needle -asequence test1.fa -bsequence test4.fa -outfile test4.aln -gapopen 10.0 -gapextend 0.5 #=> 87.3

sbatch --qos=6hours -c 10 --mem=20G launch_blast_allgenomes.sh test1.fa 1e-20 region1.blastn.tsv




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
sbatch --qos=1day -c 4 --mem=10G -e error.fitMG4.out -o slurm.fitMG4.out launch_fitMG4.sh N5.HOG0028899


Recipient branches : 
Hypomesus_transpacificus_rna_XM_047031800_1


#Test positive selection and relaxed selection on receiver branch

sbatch --qos=1day -c 4 --mem=10G -e error.absrel.cand.out -o slurm.absrel.cand.out launch_absrel_cand.sh N5.HOG0028899

#Now, test for relaxed selection
sbatch --qos=1week -c 4 --mem=10G -e error.relax.out -o slurm.relax.out --job-name=HOG0028899 launch_RELAX.sh N5.HOG0028899


#Extract dN/dS to table

grep "LB\":" N5.HOG0028899.cds.aln.FITTER.json | sed 's/.*://g' | sed 's/,.*//g' > curr_LB_values.txt
grep "MLE\":" N5.HOG0028899.cds.aln.FITTER.json | sed 's/.*://g' | sed 's/,.*//g' > curr_MLE_values.txt
grep "UB\":" N5.HOG0028899.cds.aln.FITTER.json | sed 's/.*://g' | sed 's/,.*//g' > curr_UB_values.txt
grep "\"dN\"" N5.HOG0028899.cds.aln.FITTER.json  | grep -v "{" | sed 's/.*://g' | sed 's/,//g' > curr_dN_values.txt
grep "\"dS\"" N5.HOG0028899.cds.aln.FITTER.json  | grep -v "{" | sed 's/.*://g' | sed 's/,//g' > curr_dS_values.txt
grep -B2 "LB\":" N5.HOG0028899.cds.aln.FITTER.json | grep -v "\-\-" | grep -v "Confidence" | grep -v "LB\":"  | sed 's/\"//g' | sed 's/:.*//g' | sed 's/^ *//g' > curr_labels

paste -d "," curr_labels curr_LB_values.txt curr_MLE_values.txt curr_UB_values.txt curr_dN_values.txt curr_dS_values.txt > N5.HOG0028899.dN_dS.csv









