#!/bin/bash


#SBATCH --job-name=array_synt   # Job name


master_line=$1
number_line=$2

species1=`echo "$master_line" | cut -f1 -d ','`
species2=`echo "$master_line" | cut -f2 -d ','`

for line in `cat All_pairwise_files_synteny/$species1.$species2.csv` ; do 

	gene_sp1=`echo "$line" |cut -f2 -d "," | sed "s/$species1\_//g"`
	gene_sp2=`echo "$line" | cut -f3 -d "," | sed "s/$species2\_//g"`
	
	scaff_sp1=`grep -m1 ",$gene_sp1$\|,$gene_sp1," GFF3_N5_OGGs/$species1\.gff.simplified.sorted.OGG.prez.f145.clustered.tiret  | cut -f1 -d ","`
	scaff_sp2=`grep -m1 ",$gene_sp2$\|,$gene_sp2," GFF3_N5_OGGs/$species2\.gff.simplified.sorted.OGG.prez.f145.clustered.tiret  | cut -f1 -d ","`

	grep -B10 -A10 -m1 ",$gene_sp1$\|,$gene_sp1," GFF3_N5_OGGs/$species1\.gff.simplified.sorted.OGG.prez.f145.clustered.tiret | grep "$scaff_sp1" | grep -v ",$gene_sp1$\|,$gene_sp1," | sed 's/.*N5_HOG/N5.HOG/g' | sed 's/,.*//g' | grep "[A-Z]"  > All_pairwise_files_synteny/curr_sp1_HOGS.$number_line.txt
	grep -B10 -A10 -m1 ",$gene_sp2$\|,$gene_sp2," GFF3_N5_OGGs/$species2\.gff.simplified.sorted.OGG.prez.f145.clustered.tiret | grep "$scaff_sp2" | grep -v ",$gene_sp2$\|,$gene_sp2," | sed 's/.*N5_HOG/N5.HOG/g' | sed 's/,.*//g' | grep "[A-Z]"  > All_pairwise_files_synteny/curr_sp2_HOGS.$number_line.txt

	common_HOG_nb=$(comm -12 <(sort All_pairwise_files_synteny/curr_sp1_HOGS.$number_line.txt) <(sort All_pairwise_files_synteny/curr_sp2_HOGS.$number_line.txt) | wc -l)

	nb_OGG_first_sp=`wc -l < All_pairwise_files_synteny/curr_sp1_HOGS.$number_line.txt`
	nb_OGG_second_sp=`wc -l < All_pairwise_files_synteny/curr_sp2_HOGS.$number_line.txt`

	echo "$gene_sp1,$gene_sp2,$common_HOG_nb,$nb_OGG_first_sp,$nb_OGG_second_sp" >> All_pairwise_files_synteny/$species1.$species2.SyntenyMeasures.synt

	rm -f All_pairwise_files_synteny/curr_sp1_HOGS.$number_line.txt ; rm -f All_pairwise_files_synteny/curr_sp2_HOGS.$number_line.txt

done

