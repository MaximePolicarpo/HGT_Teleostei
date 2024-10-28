#!/bin/bash


#SBATCH --job-name=HGT_stats   # Job name

LMOD_DISABLE_SAME_NAME_AUTOSWAP=no

module purge
module load R/4.3.0-foss-2021a
module load SAMtools/1.15-GCC-10.3.0
module load MAFFT/7.467-GCCcore-7.3.0-with-extensions
module load FASTX-Toolkit/0.0.14-goolf-1.7.20
module load HyPhy
module load trimAl
module load Biopython
module load EMBOSS


curr_OGG=$1 
ortho_file=$2

IFS=$'\n'



dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt"

echo "Computing stats on $curr_OGG"


#Extract the list of sequences present in the OGG
grep "$curr_OGG	" $ortho_file | sed 's/, /,/g' | cut -f4- | tr ',' '\n' | tr '\t' '\n' | grep -v "$curr_OGG" | sed '/^$/d' | sed 's/_1$//g' > $curr_OGG.list
sed $'s/[^[:print:]\t]//g' $curr_OGG.list > $curr_OGG.list.reformat ; rm $curr_OGG.list

#Extract the corresponding coding sequences 
xargs samtools faidx concatenated_nucleotides.fa < $curr_OGG.list.reformat > Coding_sequences_alignments/$curr_OGG.cds

#Remove stop codons
hyphy cln Universal Coding_sequences_alignments/$curr_OGG.cds "No/No" Coding_sequences_alignments/$curr_OGG.clean.cds ; sed -i 's/?//g' Coding_sequences_alignments/$curr_OGG.clean.cds
mv Coding_sequences_alignments/$curr_OGG.clean.cds Coding_sequences_alignments/$curr_OGG.cds
sed -i 's/---$//g' Coding_sequences_alignments/$curr_OGG.cds
sed -i 's/-$//g' Coding_sequences_alignments/$curr_OGG.cds
sed -i 's/--$//g' Coding_sequences_alignments/$curr_OGG.cds

#Translate coding sequences into proteins
transeq Coding_sequences_alignments/$curr_OGG.cds Coding_sequences_alignments/$curr_OGG.prot ; sed -i 's/_1$//g' Coding_sequences_alignments/$curr_OGG.prot

#Remove sequences with X characters -- create problems when performing alignment
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Coding_sequences_alignments/$curr_OGG.prot | sed 's/\*$//g' | awk -F '\t'  '!($2 ~ /\X/)' | tr "\t" "\n" > Coding_sequences_alignments/$curr_OGG.prot.clean
mv Coding_sequences_alignments/$curr_OGG.prot.clean Coding_sequences_alignments/$curr_OGG.prot
grep ">" Coding_sequences_alignments/$curr_OGG.prot | sed 's/>//g' > Coding_sequences_alignments/$curr_OGG.id
xargs samtools faidx Coding_sequences_alignments/$curr_OGG.cds < Coding_sequences_alignments/$curr_OGG.id  > Coding_sequences_alignments/$curr_OGG.cds.cleaned ; mv Coding_sequences_alignments/$curr_OGG.cds.cleaned Coding_sequences_alignments/$curr_OGG.cds
rm Coding_sequences_alignments/$curr_OGG.cds.fai

awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Coding_sequences_alignments/$curr_OGG.cds | sed 's/\*$//g' | awk -F '\t'  '!($2 ~ /\-/)' | tr "\t" "\n" > Coding_sequences_alignments/$curr_OGG.cds.cleaned
mv Coding_sequences_alignments/$curr_OGG.cds.cleaned Coding_sequences_alignments/$curr_OGG.cds
grep ">" Coding_sequences_alignments/$curr_OGG.cds | sed 's/>//g' > Coding_sequences_alignments/$curr_OGG.id
xargs samtools faidx Coding_sequences_alignments/$curr_OGG.prot < Coding_sequences_alignments/$curr_OGG.id  > Coding_sequences_alignments/$curr_OGG.prot.cleaned ; mv Coding_sequences_alignments/$curr_OGG.prot.cleaned Coding_sequences_alignments/$curr_OGG.prot
rm Coding_sequences_alignments/$curr_OGG.prot.fai

#Count the number of sequences remaining and the number of different species included in the sequence file

grep ">" Coding_sequences_alignments/$curr_OGG.prot | sed 's/>//g' > $curr_OGG.list.reformat
nb_seq=`grep -c ">" Coding_sequences_alignments/$curr_OGG.cds`

grep ">" Coding_sequences_alignments/$curr_OGG.cds | sed 's/---.*//g' | uniq | sed 's/>//g' > Coding_sequences_alignments/$curr_OGG.temp

nb_uniq_sp=0
for curr_sp in `cat all_annotated_species.txt` ; do 
	curr_species_name=`echo "$curr_sp-tiret" | sed 's/-tiret/_/g'`
	if grep -q "^$curr_species_name" Coding_sequences_alignments/$curr_OGG.temp ; then 

		nb_uniq_sp=$((nb_uniq_sp + 1))
		echo "$curr_sp" >> Coding_sequences_alignments/$curr_OGG.uniqSp
	fi
done
rm Coding_sequences_alignments/$curr_OGG.temp


#Lets make an anignment if there are at-least two species in the orthogroup
#If more than 300 sequences, then use muscle super5 algorithm

if [ $nb_uniq_sp -ge 2 ] ; then 

	if [ $nb_seq -ge 300 ] ; then 
		muscle5 -super5 Coding_sequences_alignments/$curr_OGG.prot -output Coding_sequences_alignments/$curr_OGG.prot.aln
	else 
		muscle5 -align Coding_sequences_alignments/$curr_OGG.prot -output Coding_sequences_alignments/$curr_OGG.prot.aln
	fi


	#trim the alignment
	trimal -in Coding_sequences_alignments/$curr_OGG.prot.aln -cons 100 -backtrans Coding_sequences_alignments/$curr_OGG.cds -out Coding_sequences_alignments/$curr_OGG.cds.aln
	trimal -in Coding_sequences_alignments/$curr_OGG.prot.aln -backtrans Coding_sequences_alignments/$curr_OGG.cds -automated1 -out Coding_sequences_alignments/$curr_OGG.cds.trimmed
	trimal -in Coding_sequences_alignments/$curr_OGG.prot.aln -automated1 -out Coding_sequences_alignments/$curr_OGG.prot.trimmed.aln

	#Compute pairwise statistics. Because we are interested in pairwise comparison statistics, we will compute those on non trimmed alignments.
	#trimmed alignments will be used to compute ML phylogenies (which will also be performed on non-trimmed aligmnents as a control)
	Rscript Compute_aligned_CDS_stats_FAST.R Coding_sequences_alignments/$curr_OGG.cds.aln Coding_sequences_alignments/$curr_OGG.cds.stats

fi



echo "HGT stats computed correctly !"

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt"

