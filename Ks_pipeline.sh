##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
################## Launch OrthoFinder ################################################################################################
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================

orthofinder -t 100 -a 50 -f Proteomes_BUSCO80/ -s AMAS_concatenated_alignment_BUSCO.fa.timetree.nwk


##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
################## Align each orthogroup coding sequences ############################################################################
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================



mkdir Coding_sequences_alignments
mkdir Log_files

ortho_file=N5.tsv 

cut -f1 $ortho_file | tail -n+2 > list_orthogroups.txt


#Launch coding sequences alignment 
for curr_OGG in `cat list_orthogroups.txt` ; do sbatch --qos=6hours --time=2:00:00 -c 10 --mem=20G -e Log_files/error.$curr_OGG.out -o Log_files/slurm.$curr_OGG.out --job-name=HGT_stats measure_Ks_stats_HGT.sh $curr_OGG $ortho_file ; done


### USE THE SAME PROCEDURE FOR HOG GENES AND BUSCO GENES


##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
################## Extract Ks Statistics #############################################################################################
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================


mkdir HGT_stats_Results_per_sp

for curr_species in `cat all_annotated_species.txt` ; do

	for curr_OGG in `cat list_orthogroups.txt` ; do 
		grep "^$curr_species" Coding_sequences_alignments/$curr_OGG.cds.stats | sed "s/^/$curr_OGG,/g" >> HGT_stats_Results_per_sp/$curr_species.cds.stats
		grep "$curr_species" Intron_sequences_alignments/$curr_OGG.fa.stats | sed "s/^/$curr_OGG	/g"  >> HGT_stats_Results_per_sp/$curr_species.introns.stats
		grep "$curr_species" Intergenic_sequences_alignments/$curr_OGG.downstream.fa.stats | sed "s/^/$curr_OGG	/g" >> HGT_stats_Results_per_sp/$curr_species.downstream.stats
		grep "$curr_species" Intergenic_sequences_alignments/$curr_OGG.upstream.fa.stats | sed "s/^/$curr_OGG	/g" >> HGT_stats_Results_per_sp/$curr_species.upstream.stats
	done


done





##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
################## Compute a micro synteny score only for pairwise species ###########################################################
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================



#We will use the GFF3 files ending with "simplified.sorted.OGG.prez.f145.clustered.tiret"
#Those are clustered GFF3 files. If two flanking genes belong to the same HOG, these are concatenated so that recent duplications don't artifically decrease
#or increase the number of shared HOG computed between species.


mkdir All_pairwise_files_synteny


line_nb=0
for master_line in `cat species_pairs.csv` ; do line_nb=$(( line_nb + 1)) ; sbatch -qos=1day -c 4 --mem=8G Compute_synteny_score.sh $master_line $line_nb ; done



#Make summary tables

rm Summary_SyntenyBlocks_Pairwise.csv
for line in `cat species_pairs.csv` ; do 
	species1=`echo "$line" | cut -f1 -d ","`
	species2=`echo "$line" | cut -f2 -d ","`
	#only keep lines with at-least 20 genes. This will be used to compute stats and to have the same number of genes compared over all comparisons
	awk -F',' '$4 >= 20 && $5 >= 20' All_pairwise_files_synteny/$species1.$species2.SyntenyMeasures.synt |  cut -f3 -d "," | sort | uniq -c | sed 's/^ *//g' | sed 's/ /,/g' | sed "s/^/$species1,$species2,/g" >> Summary_SyntenyBlocks_Pairwise.csv
done




##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
################## For each gene, compute the number of other gene in the same scaffold ##############################################
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================



### Count the number of genes on the same scaffold as each gene ... SO that I can remove genes if they are on a isolated scaffold 

for curr_species in `cat all_annotated_species.txt` ; do

	#Only keep genes present in the proteome in the GFF3 file

	grep ">" Proteomes_BUSCO80/$curr_species.fa | sed 's/>//g' | sed 's/.*---//g' | sed 's/_1$//g' > GFF3_files_per_species/$curr_species.pid.id
	awk 'NR==FNR {good_genes[$1]; next} $4 in good_genes' GFF3_files_per_species/$curr_species.pid.id GFF3_files_per_species/$curr_species.gff.simplified.sorted > GFF3_files_per_species/$curr_species.gff.simplified.sorted.inproteome

	#Count the number of genes per scaffold
	awk -F'\t' '{count[$1]++} END {for (scaffold in count) print scaffold, count[scaffold]}' GFF3_files_per_species/$curr_species.gff.simplified.sorted.inproteome > GFF3_files_per_species/$curr_species.scaffold_counts.txt
	sed -i 's/ /	/g' GFF3_files_per_species/$curr_species.scaffold_counts.txt

	#Add the count to the GFF3 file

	awk -F'\t' 'NR==FNR {counts[$1]=$2; next} {print $1, $2, $3, $4, counts[$1]}' GFF3_files_per_species/$curr_species.scaffold_counts.txt GFF3_files_per_species/$curr_species.gff.simplified.sorted.inproteome > GFF3_files_per_species/$curr_species.gff.simplified.sorted.inproteome.counts
	rm GFF3_files_per_species/$curr_species.scaffold_counts.txt ; rm GFF3_files_per_species/$curr_species.pid.id

done


rm ALL_SCAFFS.counts
for curr_species in `cat all_annotated_species.txt` ; do
	cut -f1 -d " " GFF3_files_per_species/$curr_species.gff.simplified.sorted.inproteome.counts >> ALL_SCAFFS.counts
done

rm ALL_GENES.counts
for curr_species in `cat all_annotated_species.txt` ; do
	cut -f4 -d " " GFF3_files_per_species/$curr_species.gff.simplified.sorted.inproteome.counts | sed "s/^/$curr_species\_/g" >> ALL_GENES.counts
done

rm ALL_COUNTS.counts
for curr_species in `cat all_annotated_species.txt` ; do
	cut -f5 -d " " GFF3_files_per_species/$curr_species.gff.simplified.sorted.inproteome.counts  >> ALL_COUNTS.counts
done

rm ALL_SPECIES.counts
for curr_species in `cat all_annotated_species.txt` ; do
	cut -f6 -d " " GFF3_files_per_species/$curr_species.gff.simplified.sorted.inproteome.counts  | sed "s/^/$curr_species/g"  >> ALL_SPECIES.counts
done


paste -d "," ALL_SPECIES.counts ALL_SCAFFS.counts ALL_GENES.counts ALL_COUNTS.counts > ALL_species.inproteome.counts



##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
################## Identify genes corresponding to transposable elements #############################################################
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================


#RepeatPeps.lib => RepeatMasker Database
#repbase22.05_aaSeq_cleaned_TE.fa => RepBase database

module load DIAMOND
module load BLAST

#Build diamond database

diamond makedb --in RepeatPeps.lib -d RepeatPeps --threads 20
diamond makedb --in repbase22.05_aaSeq_cleaned_TE.fa -d RepBase --threads 20


#Build BLAST database

makeblastdb -in RepeatPeps.lib -dbtype prot 
makeblastdb -in repbase22.05_aaSeq_cleaned_TE.fa --dbtype prot 




#Launch diamond and blast for both database


mkdir Transposon_search_blast


for curr_species in `cat all_annotated_species.txt` ; do
	diamond blastp -p 20 --query Proteomes_BUSCO80/$curr_species.fa --max-target-seqs 5 --out Transposon_search_blast/$curr_species.tsv.diamond --db RepeatPeps --outfmt 6
	blastp -query Proteomes_BUSCO80/$curr_species.fa -db RepeatPeps.lib -outfmt 6 -out Transposon_search_blast/$curr_species.BLAST.tsv -max_target_seqs 5 -num_threads 20
done


mkdir Transposon_search_blast_Repbase

for curr_species in `cat all_annotated_species.txt` ; do
	diamond blastp -p 20 --query Proteomes_BUSCO80/$curr_species.fa --max-target-seqs 5 --out Transposon_search_blast/$curr_species.tsv.diamond --db RepBase --outfmt 6
	blastp -query Proteomes_BUSCO80/$curr_species.fa -db repbase22.05_aaSeq_cleaned_TE.fa -outfmt 6 -out Transposon_search_blast/$curr_species.BLAST.tsv -max_target_seqs 5 -num_threads 20 
done




#Launch interproscan to find genes with TE associated domains

for curr_species in `cat all_annotated_species.txt` ; do
	interproscan.sh -i Proteomes_BUSCO80/$curr_species.fa -d out_interpro/ -t p -appl Pfam -f TSV -cpu 20 -T temp.$curr_species/
	grep "$element" out_interpro/$curr_species.fa.tsv | cut -f1 >> Transposable_elements_candidates_PFAM.txt
done
sort Transposable_elements_candidates_PFAM.txt | uniq > temp ; mv temp Transposable_elements_candidates_PFAM.txt



#Extract results

cut -f1 Transposon_search_blast/*.diamond | sort | uniq > Transposable_elements_candidates.repeatpeps.diamond.txt
cut -f1 Transposon_search_blast/*.BLAST.tsv | sort | uniq > Transposable_elements_candidates.repeatpeps.blast.txt

cut -f1 Transposon_search_blast_Repbase/*.diamond | sort | uniq > Transposable_elements_candidates.repbase.diamond.txt
cut -f1 Transposon_search_blast_Repbase/*.BLAST.tsv | sort | uniq > Transposable_elements_candidates.repbase.blast.txt



#Combine results 


cat Transposable_elements_candidates.repeatpeps.diamond.txt Transposable_elements_candidates.repeatpeps.blast.txt Transposable_elements_candidates.repbase.diamond.txt Transposable_elements_candidates.repbase.blast.txt Transposable_elements_candidates_PFAM.txt | sort | uniq > Transposable_elements_candidates_PFAM_BLAST_DIAMOND.txt


#Extract OGG corresponding to genes matchig to TEs. Consider a OGG as a 'TE OGG' if it contains >20% of TE


ortho_file=N5.tsv 
cp N5.tsv  N5_copy.tsv
sed -i 's/, /	/g' N5_copy.tsv
sed -i $'s/[^[:print:]\t]//g' N5_copy.tsv

awk 'FNR==NR{genes[$1]; next} {for (i=2; i<=NF; i++) if ($i in genes) {print $i "\t" $1}}' Transposable_elements_candidates_PFAM_BLAST_DIAMOND.txt N5_copy.tsv > N5_OGG_Transposable_elements_candidates_PFAM_BLAST_DIAMOND_atleastOne.txt
cut -f2 N5_OGG_Transposable_elements_candidates_PFAM_BLAST_DIAMOND_atleastOne.txt | sort | uniq > temp ; mv temp N5_OGG_Transposable_elements_candidates_PFAM_BLAST_DIAMOND_atleastOne.txt


for curr_OGG in `cat N5_OGG_Transposable_elements_candidates_PFAM_BLAST_DIAMOND_atleastOne.txt` ; do 
	grep "$curr_OGG	" $ortho_file | sed 's/, /,/g' | cut -f4- | tr ',' '\n' | tr '\t' '\n' | grep -v "$curr_OGG" | sed '/^$/d' | sed 's/_1$//g' > $curr_OGG.list
	sed $'s/[^[:print:]\t]//g' $curr_OGG.list | sort | uniq > $curr_OGG.list.reformat ; rm $curr_OGG.list

	nb_TE_genes=`comm -12  Transposable_elements_candidates_PFAM_BLAST.txt $curr_OGG.list.reformat | wc -l`
	nb_total_genes=`wc -l < $curr_OGG.list.reformat`
	rm $curr_OGG.list.reformat

	perc20_gene=`echo $(( nb_total_genes*20/100 ))`

	if [ $nb_TE_genes -ge $perc20_gene ] ; then echo "$curr_OGG" >> N5_OGG_Transposable_elements_candidates_PFAM_BLAST_20perc.txt ; fi

done


cp N5_copy.tsv N5_copy.test.tsv
sed -i 's/\t\{2,\}/\t/g' N5_copy.test.tsv
awk -F'\t' '{print $1, NF-3}' N5_copy.test.tsv > N5.OGG_copynumber.tsv
tail -n+2 N5.OGG_copynumber.tsv  > temp ; sed -i 's/ /,/g' temp ; mv temp N5.OGG_copynumber.tsv

awk 'FNR==NR{genes[$1]; next} {for (i=2; i<=NF; i++) if ($i in genes) {print $i "\t" $1}}' Transposable_elements_candidates_PFAM_BLAST_DIAMOND.txt N5_copy.tsv | cut -f2 > temp
sort temp | uniq -c > TE_ogg.copynumber.csv

sed 's/^ *//g' TE_ogg.copynumber.csv | sed 's/ /,/g' > temp
cut -f2 -d "," temp > c1.txt
cut -f1 -d "," temp > c2.txt
paste -d "," c1.txt c2.txt > temp
mv temp TE_ogg.copynumber.csv


join -t, -1 1 -2 1 TE_ogg.copynumber.csv N5.OGG_copynumber.tsv > merged.csv
awk -F, '($2 >= $3 * 0.2) {print $1}' merged.csv >  N5_OGG_Transposable_elements_candidates_PFAM_BLAST_20perc.diamond_blast_pfam.txt 








##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
################## Reduce the list of candidate HGTs ################################################################################
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================



#Create file with a list of  teleost species 
grep -v "Polypterus_senegalus\|Erpetoichthys_calabaricus\|Amia_calva\|Atractosteus_spatula\|Lepisosteus_oculatus\|Acipenser_oxyrinchus_oxyrinchus\|Acipenser_ruthenus\|Huso_huso\|Polyodon_spathula" all_annotated_species.txt > teleost_annotated_species.txt


#First, Define the ks distributions and their quantile + mean + median between each pair of species

mkdir KS_Stats_per_species
for curr_species in `cat teleost_annotated_species.txt` ; do
	Rscript Extract_ks_stats_OGG.R $curr_species KS_Stats_per_species/$curr_species.csv
done

head -1 KS_Stats_per_species/Danio_rerio.csv > All_OGG_Ks_stats.csv
for curr_species in `cat teleost_annotated_species.txt` ; do tail -n+2 KS_Stats_per_species/$curr_species.csv >> All_OGG_Ks_stats.csv ; done


#Now extract only pairwise comparisons wich fall below the quantile 0.5% (keep maximum between the quantile values computed with BUSCO or with HOG genes)

mkdir Candidate_HGT_tables

for curr_species in `cat teleost_annotated_species.txt` ; do
	Rscript Extract_candidate_HGT.R $curr_species Candidate_HGT_tables/$curr_species.csv
done

head -1 Candidate_HGT_tables/Danio_rerio.csv > HGT_candidates_dataframe.csv
for curr_species in `cat teleost_annotated_species.txt` ; do tail -n+2 Candidate_HGT_tables/$curr_species.csv >> HGT_candidates_dataframe.csv ; done











##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
################## Add Synt score to the dataframe of HGT candidates #################################################################
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================


mkdir Gene_LogFiles
mkdir Temp_dir_synt
tail -n+2 HGT_candidates_dataframe.filtered.csv > HGT_candidates_dataframe.filtered.noheader.csv ##HGT_candidates_dataframe.filtered.csv > generated using the Rscript "HGT_analysis.R"
split -l 10000 --numeric-suffixes HGT_candidates_dataframe.filtered.noheader.csv Temp_dir_synt/splitted.HGT_candidates_dataframe.filtered. 
ls -l Temp_dir_synt/ | sed 's/.* //g' | grep "splitted" > list_splitted_df.txt


rm MicroSynteny_scores_table.csv #remove file if exist


sbatch submit_array_synt.sh


rm -rf Gene_LogFiles 

cat Temp_dir_synt/*.csv  > MicroSynteny_scores_table.csv




================== extract_synt.sh =====



#!/bin/bash


#SBATCH --job-name=extract_synt  # Job name

curr_file=$1 
nb_line=`wc -l < Temp_dir_synt/$curr_file`

for line_nb in `seq 1 $nb_line` ; do 

	line=`head -$line_nb Temp_dir_synt/$curr_file | tail -1`
	seq1_full=`echo "$line" | cut -f4 -d ","`
	seq2_full=`echo "$line" | cut -f5 -d ","`

	sp1=`echo "$line" | cut -f1 -d ","`
	sp2=`echo "$line" | cut -f2 -d ","`

	seq1_name=`echo "$seq1_full" | sed "s/$sp1\_//g"`
	seq2_name=`echo "$seq2_full" | sed "s/$sp2\_//g"`

	#echo "$seq1_name,$seq2_name"

	if [ -f All_pairwise_files_synteny/$sp1.$sp2.SyntenyMeasures.synt ] ; then

		grep "$seq1_name.*$seq2_name" All_pairwise_files_synteny/$sp1.$sp2.SyntenyMeasures.synt | sed "s/$seq1_name/$seq1_full/g" | sed "s/$seq2_name/$seq2_full/g" >> Temp_dir_synt/$curr_file.csv
	
	elif grep -q "$seq2_name.*$seq1_name" All_pairwise_files_synteny/$sp2.$sp1.SyntenyMeasures.synt ; then

		grep "$seq2_name.*$seq1_name" All_pairwise_files_synteny/$sp2.$sp1.SyntenyMeasures.synt | sed "s/$seq1_name/$seq1_full/g" | sed "s/$seq2_name/$seq2_full/g" >> Temp_dir_synt/$curr_file.csv

	else 

		echo "not found"
	fi

done 


================== submit_array_synt.sh =====

#!/bin/bash

#SBATCH --job-name=synt_score
#SBATCH --qos=6hours
#SBATCH -c 2
#SBATCH --mem=5G
#SBATCH --error=Gene_LogFiles/error.%A_%a.out
#SBATCH --output=Gene_LogFiles/slurm.%A_%a.out
#SBATCH --array=1-XXX%500 => replace X by the total number of lines in the file list_splitted_df.txt

#Remove existing log files
rm -f Gene_LogFiles/error.*.out
rm -f Gene_LogFiles/slurm.*.out

# Extract the gene corresponding to the array task ID
file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" list_splitted_df.txt)

# Call the script to process each gene
./extract_synt.sh $file







