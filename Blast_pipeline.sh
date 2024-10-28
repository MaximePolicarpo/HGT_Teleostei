#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
################## Blast All vs All  ##################################################################################################
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================


mkdir Coding_sequences_files_per_species


#Make a BLAST database for each species

for species in `cat all_annotated_species.txt` ; do
	makeblastdb -in Coding_sequences_files_per_species/$species.cds -dbtype nucl
done



#Launch a blast of every sequence from one species against the sequences of all other species


mkdir Concatenated_CDS_files
mkdir Blastn_results


for species in `cat all_annotated_species.txt` ; do

	./run_blast_proc.sh $species

done




##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
################## Extract BLAST outliers with the tree strategies  ##################################################################
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================


grep -v "Polypterus_senegalus\|Erpetoichthys_calabaricus\|Amia_calva\|Atractosteus_spatula\|Lepisosteus_oculatus\|Acipenser_oxyrinchus_oxyrinchus\|Acipenser_ruthenus\|Huso_huso\|Polyodon_spathula" all_annotated_species.txt > teleost_annotated_species.txt


mkdir Outliers_BLAST_stratA
mkdir Outliers_BLAST_stratB
mkdir Outliers_BLAST_stratC


for curr_species in `cat teleost_annotated_species.txt` ; do 
	Rscript Parse_blast_results.R $curr_species Outliers_BLAST_stratA/$curr_species.csv Outliers_BLAST_stratB/$curr_species.csv Outliers_BLAST_stratC/$curr_species.csv
done



head -1 Outliers_BLAST_stratA/Danio_rerio.csv > All_blast_outliers_table.A.csv
head -1 Outliers_BLAST_stratB/Danio_rerio.csv > All_blast_outliers_table.B.csv
head -1 Outliers_BLAST_stratC/Danio_rerio.csv > All_blast_outliers_table.C.csv

for curr_species in `cat teleost_annotated_species.txt` ; do tail -n+2 Outliers_BLAST_stratA/$curr_species.csv >> All_blast_outliers_table.A.csv ; done
for curr_species in `cat teleost_annotated_species.txt` ; do tail -n+2 Outliers_BLAST_stratB/$curr_species.csv >> All_blast_outliers_table.B.csv ; done
for curr_species in `cat teleost_annotated_species.txt` ; do tail -n+2 Outliers_BLAST_stratC/$curr_species.csv >> All_blast_outliers_table.C.csv ; done


#ADD the OGG information ..

ortho_file=N5.tsv 
cp $ortho_file ./N5_copy.tsv
sed -i 's/, /	/g' N5_copy.tsv


cut -f1 -d "," All_blast_outliers_table.A.csv | tail -n+2 > c1_genes.A
cut -f2 -d "," All_blast_outliers_table.A.csv | tail -n+2 > c2_genes.A
awk 'FNR==NR{genes[$1]; next} {for (i=2; i<=NF; i++) if ($i in genes) {print $i "\t" $1}}' c1_genes.A N5_copy.tsv > c1_genes.OGG.A.txt
awk 'FNR==NR{genes[$1]; next} {for (i=2; i<=NF; i++) if ($i in genes) {print $i "\t" $1}}' c2_genes.A N5_copy.tsv > c2_genes.OGG.A.txt


cut -f1 -d "," All_blast_outliers_table.B.csv | tail -n+2 > c1_genes.B
cut -f2 -d "," All_blast_outliers_table.B.csv | tail -n+2 > c2_genes.B
awk 'FNR==NR{genes[$1]; next} {for (i=2; i<=NF; i++) if ($i in genes) {print $i "\t" $1}}' c1_genes.B N5_copy.tsv > c1_genes.OGG.B.txt
awk 'FNR==NR{genes[$1]; next} {for (i=2; i<=NF; i++) if ($i in genes) {print $i "\t" $1}}' c2_genes.B N5_copy.tsv > c2_genes.OGG.B.txt



cut -f1 -d "," All_blast_outliers_table.C.csv | tail -n+2 > c1_genes.C
cut -f2 -d "," All_blast_outliers_table.C.csv | tail -n+2 > c2_genes.C
awk 'FNR==NR{genes[$1]; next} {for (i=2; i<=NF; i++) if ($i in genes) {print $i "\t" $1}}' c1_genes.C N5_copy.tsv > c1_genes.OGG.C.txt
awk 'FNR==NR{genes[$1]; next} {for (i=2; i<=NF; i++) if ($i in genes) {print $i "\t" $1}}' c2_genes.C N5_copy.tsv > c2_genes.OGG.C.txt


##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
################## Extract a micro synteny score only for pairwise gene which pass the blast filters #################################
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================


tail -n+2 blast_outliers_df_filtered.A.csv > blast_outliers_df_filtered.noheader.A.csv #file generated using the script HGT_analysis.R
tail -n+2 blast_outliers_df_filtered.B.csv > blast_outliers_df_filtered.noheader.B.csv #file generated using the script HGT_analysis.R
tail -n+2 blast_outliers_df_filtered.C.csv > blast_outliers_df_filtered.noheader.C.csv #file generated using the script HGT_analysis.R

#cp MicroSynteny_scores_table.csv coffre.MicroSynteny_scores_table.csv

rm MicroSynteny_scores_table.csv
rm blast_outliers_df_filtered.noheader.A.csv.MicroSynteny_scores_table.csv
rm blast_outliers_df_filtered.noheader.B.csv.MicroSynteny_scores_table.csv
rm blast_outliers_df_filtered.noheader.C.csv.MicroSynteny_scores_table.csv

./extract_synt.blast.AB.sh blast_outliers_df_filtered.noheader.A.csv
./extract_synt.blast.AB.sh blast_outliers_df_filtered.noheader.B.csv
./extract_synt.blast.C.sh blast_outliers_df_filtered.noheader.C.csv



cat blast_outliers_df_filtered.noheader.A.csv.MicroSynteny_scores_table.csv  blast_outliers_df_filtered.noheader.B.csv.MicroSynteny_scores_table.csv blast_outliers_df_filtered.noheader.C.csv.MicroSynteny_scores_table.csv > MicroSynteny_scores_table.blast.csv



================== extract_synt.blast.AB.sh =====



#!/bin/bash


#SBATCH --job-name=extract_synt  # Job name

curr_file=$1 
nb_line=`wc -l < $curr_file`

for line_nb in `seq 1 $nb_line` ; do 

	line=`head -$line_nb $curr_file | tail -1`
	seq1_full=`echo "$line" | cut -f1 -d ","`
	seq2_full=`echo "$line" | cut -f2 -d ","`

	sp1=`echo "$seq1_full" | sed 's/---.*//g'`
	sp2=`echo "$seq2_full" | sed 's/---.*//g'`

	seq1_name=`echo "$seq1_full" | sed "s/$sp1---//g" |  sed 's/gene-/gene_/g' | sed 's/rna-/rna_/g' | sed 's/\.1/_1/g' | sed 's/\.2/_2/g' | sed 's/\.3/_3/g' | sed 's/\.4/_4/g'`
	seq2_name=`echo "$seq2_full" | sed "s/$sp2---//g" |  sed 's/gene-/gene_/g' | sed 's/rna-/rna_/g' | sed 's/\.1/_1/g' | sed 's/\.2/_2/g' | sed 's/\.3/_3/g' | sed 's/\.4/_4/g'`


	#echo "$seq1_name,$seq2_name"

	if [ -f All_pairwise_files_synteny/$sp1.$sp2.SyntenyMeasures.synt ] ; then

		grep "$seq1_name.*$seq2_name" All_pairwise_files_synteny/$sp1.$sp2.SyntenyMeasures.synt | sed "s/$seq1_name/$seq1_full/g" | sed "s/$seq2_name/$seq2_full/g" >> $curr_file.MicroSynteny_scores_table.csv
	
	elif grep -q "$seq2_name.*$seq1_name" All_pairwise_files_synteny/$sp2.$sp1.SyntenyMeasures.synt ; then

		grep "$seq2_name.*$seq1_name" All_pairwise_files_synteny/$sp2.$sp1.SyntenyMeasures.synt | sed "s/$seq1_name/$seq1_full/g" | sed "s/$seq2_name/$seq2_full/g" >> $curr_file.MicroSynteny_scores_table.csv

	else 

		echo "not found"
	fi

done 






================== extract_synt.blast.C.sh =====



#!/bin/bash


#SBATCH --job-name=extract_synt  # Job name

curr_file=$1 
nb_line=`wc -l < $curr_file`

for line_nb in `seq 1 $nb_line` ; do 

	line=`head -$line_nb $curr_file | tail -1`
	seq1_full=`echo "$line" | cut -f1 -d ","`
	seq2_full=`echo "$line" | cut -f6 -d ","`

	sp1=`echo "$seq1_full" | sed 's/---.*//g'`
	sp2=`echo "$seq2_full" | sed 's/---.*//g'`

	seq1_name=`echo "$seq1_full" | sed "s/$sp1---//g" |  sed 's/gene-/gene_/g' | sed 's/rna-/rna_/g' | sed 's/\.1/_1/g' | sed 's/\.2/_2/g' | sed 's/\.3/_3/g' | sed 's/\.4/_4/g'`
	seq2_name=`echo "$seq2_full" | sed "s/$sp2---//g" |  sed 's/gene-/gene_/g' | sed 's/rna-/rna_/g' | sed 's/\.1/_1/g' | sed 's/\.2/_2/g' | sed 's/\.3/_3/g' | sed 's/\.4/_4/g'`


	#echo "$seq1_name,$seq2_name"

	if [ -f ../BetweenActino_HGT/All_pairwise_files_synteny/$sp1.$sp2.SyntenyMeasures.synt ] ; then

		grep "$seq1_name.*$seq2_name" All_pairwise_files_synteny/$sp1.$sp2.SyntenyMeasures.synt | sed "s/$seq1_name/$seq1_full/g" | sed "s/$seq2_name/$seq2_full/g" >> $curr_file.MicroSynteny_scores_table.csv
	
	elif grep -q "$seq2_name.*$seq1_name" All_pairwise_files_synteny/$sp2.$sp1.SyntenyMeasures.synt ; then

		grep "$seq2_name.*$seq1_name" All_pairwise_files_synteny/$sp2.$sp1.SyntenyMeasures.synt | sed "s/$seq1_name/$seq1_full/g" | sed "s/$seq2_name/$seq2_full/g" >> $curr_file.MicroSynteny_scores_table.csv

	else 

		echo "not found"
	fi

done 



##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
################## Add KS informations to retained genes  ############################################################################
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================


IFS=$'\n'
tail -n+2 blast_outliers_df_A_micro.tsv > blast_outliers_df_A_micro.tsv.noheader
header_line=`head -1 blast_outliers_df_A_micro.tsv` ; echo "$header_line,ks,ka,cds_identity" > blast_outliers_df_A_micro.kaks.tsv
for line in `cat blast_outliers_df_A_micro.tsv.noheader` ; do
	curr_sp=`echo "$line" | cut -f16`
	curr_qseq=`echo "$line" | cut -f1`
	curr_sseq=`echo "$line" | cut -f2`

	seq1_name=`echo "$curr_qseq" | sed 's/---/_/g' | sed 's/gene-/gene_/g' | sed 's/rna-/rna_/g' | sed 's/\.1/_1/g' | sed 's/\.2/_2/g' | sed 's/\.3/_3/g' | sed 's/\.4/_4/g' | sed 's/.model/_model/g'`
	seq2_name=`echo "$curr_sseq" | sed 's/---/_/g' | sed 's/gene-/gene_/g' | sed 's/rna-/rna_/g' | sed 's/\.1/_1/g' | sed 's/\.2/_2/g' | sed 's/\.3/_3/g' | sed 's/\.4/_4/g' | sed 's/.model/_model/g'`

	if grep -q "$seq1_name.*$seq2_name" HGT_stats_Results_per_sp/$curr_sp.cds.stats ; then
		ks_line=`grep "$seq1_name.*$seq2_name" HGT_stats_Results_per_sp/$curr_sp.cds.stats`
	elif grep -q "$seq2_name.*$seq1_name" HGT_stats_Results_per_sp/$curr_sp.cds.stats ; then
		ks_line=`grep "$seq2_name.*$seq1_name" HGT_stats_Results_per_sp/$curr_sp.cds.stats`
	else
		echo "not found"
		echo "$line"

	fi

	curr_ks=`echo "$ks_line" | cut -f4 -d ","`
	curr_ka=`echo "$ks_line" | cut -f5 -d ","`
	curr_ident=`echo "$ks_line" | cut -f10 -d ","`


	echo "$line,$curr_ks,$curr_ka,$curr_ident" >> blast_outliers_df_A_micro.kaks.tsv

done

cat blast_outliers_df_A_micro.kaks.tsv | tr ',' "\t"  > temp 
mv temp blast_outliers_df_A_micro.kaks.tsv



IFS=$'\n'
tail -n+2 blast_outliers_df_B_micro.tsv > blast_outliers_df_B_micro.tsv.noheader
header_line=`head -1 blast_outliers_df_B_micro.tsv` ; echo "$header_line,ks,ka,cds_identity" > blast_outliers_df_B_micro.kaks.tsv
for line in `cat blast_outliers_df_B_micro.tsv.noheader` ; do
	curr_sp=`echo "$line" | cut -f16`
	curr_qseq=`echo "$line" | cut -f1`
	curr_sseq=`echo "$line" | cut -f2`

	seq1_name=`echo "$curr_qseq" | sed 's/---/_/g' | sed 's/gene-/gene_/g' | sed 's/rna-/rna_/g' | sed 's/\.1/_1/g' | sed 's/\.2/_2/g' | sed 's/\.3/_3/g' | sed 's/\.4/_4/g' | sed 's/.model/_model/g'`
	seq2_name=`echo "$curr_sseq" | sed 's/---/_/g' | sed 's/gene-/gene_/g' | sed 's/rna-/rna_/g' | sed 's/\.1/_1/g' | sed 's/\.2/_2/g' | sed 's/\.3/_3/g' | sed 's/\.4/_4/g' | sed 's/.model/_model/g'`

	if grep -q "$seq1_name.*$seq2_name" HGT_stats_Results_per_sp/$curr_sp.cds.stats ; then
		ks_line=`grep "$seq1_name.*$seq2_name" HGT_stats_Results_per_sp/$curr_sp.cds.stats`
	elif grep -q "$seq2_name.*$seq1_name" HGT_stats_Results_per_sp/$curr_sp.cds.stats ; then
		ks_line=`grep "$seq2_name.*$seq1_name" HGT_stats_Results_per_sp/$curr_sp.cds.stats`
	else
		echo "not found"
		echo "$line"

	fi

	curr_ks=`echo "$ks_line" | cut -f4 -d ","`
	curr_ka=`echo "$ks_line" | cut -f5 -d ","`
	curr_ident=`echo "$ks_line" | cut -f10 -d ","`


	echo "$line,$curr_ks,$curr_ka,$curr_ident" >> blast_outliers_df_B_micro.kaks.tsv

done

cat blast_outliers_df_B_micro.kaks.tsv | tr ',' "\t"  > temp 
mv temp blast_outliers_df_B_micro.kaks.tsv


IFS=$'\n'
tail -n+2 blast_outliers_df_C_micro.tsv > blast_outliers_df_C_micro.tsv.noheader
header_line=`head -1 blast_outliers_df_C_micro.tsv` ; echo "$header_line,ks,ka,cds_identity" > blast_outliers_df_C_micro.kaks.tsv
for line in `cat blast_outliers_df_C_micro.tsv.noheader` ; do
	curr_sp=`echo "$line" | cut -f16`
	curr_qseq=`echo "$line" | cut -f1`
	curr_sseq=`echo "$line" | cut -f2`

	seq1_name=`echo "$curr_qseq" | sed 's/---/_/g' | sed 's/gene-/gene_/g' | sed 's/rna-/rna_/g' | sed 's/\.1/_1/g' | sed 's/\.2/_2/g' | sed 's/\.3/_3/g' | sed 's/\.4/_4/g' | sed 's/.model/_model/g'`
	seq2_name=`echo "$curr_sseq" | sed 's/---/_/g' | sed 's/gene-/gene_/g' | sed 's/rna-/rna_/g' | sed 's/\.1/_1/g' | sed 's/\.2/_2/g' | sed 's/\.3/_3/g' | sed 's/\.4/_4/g' | sed 's/.model/_model/g'`

	if grep -q "$seq1_name.*$seq2_name" HGT_stats_Results_per_sp/$curr_sp.cds.stats ; then
		ks_line=`grep "$seq1_name.*$seq2_name" HGT_stats_Results_per_sp/$curr_sp.cds.stats`
	elif grep -q "$seq2_name.*$seq1_name" HGT_stats_Results_per_sp/$curr_sp.cds.stats ; then
		ks_line=`grep "$seq2_name.*$seq1_name" HGT_stats_Results_per_sp/$curr_sp.cds.stats`
	else
		echo "not found"
		echo "$line"

	fi

	curr_ks=`echo "$ks_line" | cut -f4 -d ","`
	curr_ka=`echo "$ks_line" | cut -f5 -d ","`
	curr_ident=`echo "$ks_line" | cut -f10 -d ","`


	echo "$line,$curr_ks,$curr_ka,$curr_ident" >> blast_outliers_df_C_micro.kaks.tsv

done

cat blast_outliers_df_C_micro.kaks.tsv | tr ',' "\t"  > temp 
mv temp blast_outliers_df_C_micro.kaks.tsv





