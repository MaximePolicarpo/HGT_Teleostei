#####====================================================================================================================================
#####====================================================================================================================================
####====================================================================================================================================
#####====================================================================================================================================
#####====================================================================================================================================
#####========================================== Extract mitochondrial genome from assembly =============================================
#####====================================================================================================================================
#####====================================================================================================================================
####====================================================================================================================================
#####====================================================================================================================================
#####====================================================================================================================================

#Extract/Download mitochondrial genomes

#Clupea_harengus
mt_scaff=`grep "^MT	" GCF_900700415.2_Ch_v2.0.2_assembly_report.txt | cut -f7`
samtools faidx GCF_900700415.2_Ch_v2.0.2_genomic.fna $mt_scaff > Mitochondrial_Genomes/Clupea_harengus.mt.fa

#Hypomesus_transpacificus
https://www.ncbi.nlm.nih.gov/nuccore/NC_072210.1

#Anoplopoma_fimbria
https://www.ncbi.nlm.nih.gov/nuccore/NC_018119

#Channa_striata
https://www.ncbi.nlm.nih.gov/nuccore/NC_032037


#Electrona_antarctica
mt_scaff=`grep "^scaffold_MT_1" GCA_951216825.1_fEleAnt2.1_assembly_report.txt | cut -f5`
samtools faidx .GCA_951216825.1_fEleAnt2.1_genomic.fna $mt_scaff > Mitochondrial_Genomes/Electrona_antarctica.mt.fa

#Oncorhynchus_mykiss
mt_scaff=`grep "^MT	" GCF_013265735.2_USDA_OmykA_1.1_assembly_report.txt | cut -f7`
samtools faidx GCF_013265735.2_USDA_OmykA_1.1_genomic.fna $mt_scaff > Mitochondrial_Genomes/Oncorhynchus_mykiss.mt.fa


#Oryzias_melastigma

mt_scaff=`grep "^MT	" GCF_002922805.2_ASM292280v2_assembly_report.txt | cut -f7`
samtools faidx GCF_002922805.2_ASM292280v2_genomic.fna  $mt_scaff > Mitochondrial_Genomes/Oryzias_melastigma.mt.fa


#Salarias_fasciatus

mt_scaff=`grep "^MT	" GCF_902148845.1_fSalaFa1.1_assembly_report.txt | cut -f7`
samtools faidx GCF_902148845.1_fSalaFa1.1_genomic.fna  $mt_scaff > Mitochondrial_Genomes/Salarias_fasciatus.mt.fa


#Borostomias_antarcticus
mt_scaff=`grep "^scaffold_MT_1" GCA_949987555.1_fBorAnt1.1_assembly_report.txt | cut -f5`
samtools faidx GCA_949987555.1_fBorAnt1.1_genomic.fna $mt_scaff > Mitochondrial_Genomes/Borostomias_antarcticus.mt.fa



#Siniperca_chuatsi

mt_scaff=`grep "^MT	" GCF_020085105.1_ASM2008510v1_assembly_report.txt | cut -f7`
samtools faidx GCF_020085105.1_ASM2008510v1_genomic.fna  $mt_scaff > Mitochondrial_Genomes/Siniperca_chuatsi.mt.fa



#Takifugu_flavidus

mt_scaff=`grep "^MT	" GCF_003711565.1_ASM371156v2_assembly_report.txt  | cut -f7`
samtools faidx GCF_003711565.1_ASM371156v2_genomic.fna $mt_scaff > Mitochondrial_Genomes/Takifugu_flavidus.mt.fa



#Xiphophorus_couchianus

mt_scaff=`grep "^MT	" GCF_001444195.1_X_couchianus-1.0_assembly_report.txt | cut -f7`
samtools faidx GCF_001444195.1_X_couchianus-1.0_genomic.fna $mt_scaff > Mitochondrial_Genomes/Xiphophorus_couchianus.mt.fa




#Now create bwamem index for every mitochondrial genomes


module purge ; module load bwa-mem2
bwa-mem2 index Mitochondrial_Genomes/Clupea_harengus.mt.fa
bwa-mem2 index Mitochondrial_Genomes/Hypomesus_transpacificus.mt.fa
bwa-mem2 index Mitochondrial_Genomes/Anoplopoma_fimbria.mt.fa
bwa-mem2 index Mitochondrial_Genomes/Borostomias_antarcticus.mt.fa
bwa-mem2 index Mitochondrial_Genomes/Channa_striata.mt.fa
bwa-mem2 index Mitochondrial_Genomes/Electrona_antarctica.mt.fa
bwa-mem2 index Mitochondrial_Genomes/Oncorhynchus_mykiss.mt.fa
bwa-mem2 index Mitochondrial_Genomes/Oryzias_melastigma.mt.fa
bwa-mem2 index Mitochondrial_Genomes/Salarias_fasciatus.mt.fa
bwa-mem2 index Mitochondrial_Genomes/Siniperca_chuatsi.mt.fa
bwa-mem2 index Mitochondrial_Genomes/Takifugu_flavidus.mt.fa
bwa-mem2 index Mitochondrial_Genomes/Xiphophorus_couchianus.mt.fa

#####====================================================================================================================================
#####====================================================================================================================================
####====================================================================================================================================
#####====================================================================================================================================
#####====================================================================================================================================
#####========================================== Download and map reads to mitochondrial genomes =========================================
#####====================================================================================================================================
#####====================================================================================================================================
####====================================================================================================================================
#####====================================================================================================================================
#####====================================================================================================================================



./download_and_launch_map.clupea.sh Hypomesus_transpacificus
./download_and_launch_map.clupea.sh Anoplopoma_fimbria
./download_and_launch_map.clupea.sh Borostomias_antarcticus
./download_and_launch_map.clupea.sh Channa_striata
./download_and_launch_map.clupea.sh Electrona_antarctica
./download_and_launch_map.clupea.sh Oncorhynchus_mykiss
./download_and_launch_map.clupea.sh Oryzias_melastigma
./download_and_launch_map.clupea.sh Salarias_fasciatus
./download_and_launch_map.clupea.sh Siniperca_chuatsi
./download_and_launch_map.clupea.sh Takifugu_flavidus
./download_and_launch_map.clupea.sh Xiphophorus_couchianus


sed 's/Average =  //g' Mapping.*.csv > Summary_mapping.csv




#####====================================================================================================================================
#####====================================================================================================================================
####====================================================================================================================================
#####====================================================================================================================================
#####====================================================================================================================================
#####========================================== Make a phylogeny of mitochondrial genomes ===============================================
#####====================================================================================================================================
#####====================================================================================================================================
####====================================================================================================================================
#####====================================================================================================================================
#####====================================================================================================================================

cd Mitochondrial_genes/


#Extract mito genes from already annotated mitochondrial genomes

grep ">" Clupea_harengus.fa  | sed 's/.*gene=//g' | sed 's/].*//g'  > list_genes_mito.txt

mkdir Alignments
for gene in `cat list_genes_mito.txt` ; do
	for file in *.fa ; do
		species_name=`echo "$file" | cut -f1 -d "."`
		curr_name=`grep "gene=$gene]" $file | sed 's/>//g' | cut -f1 -d " "`
		samtools faidx $file $curr_name  | sed "s/>.*/>$species_name/g" >> Alignments/$gene.fa
		samtools faidx $file $curr_name  | sed "s/>.*/>$gene.$species_name/g" >> Mito_DB.fa
	done
done

transeq Mito_DB.fa Mito_DB.prot -table 2 ; sed -i 's/_1$//g' Mito_DB.prot
makeblastdb -in Mito_DB.prot -dbtype prot


#Annotate mito genes for non-annotated mitochondria 

getorf -sequence ../Mitochondrial_Genomes/Electrona_antarctica.mt.fa -outseq Electrona_antarctica.orf -minsize 168 -find 3 -table 2
getorf -sequence ../Mitochondrial_Genomes/Borostomias_antarcticus.mt.fa -outseq Borostomias_antarcticus.orf -minsize 168 -find 3 -table 2

transeq Electrona_antarctica.orf Electrona_antarctica.orf.prot -table 2 ; sed -i 's/_1 / /g' Electrona_antarctica.orf.prot
transeq Borostomias_antarcticus.orf Borostomias_antarcticus.orf.prot -table 2 ; sed -i 's/_1 / /g' Borostomias_antarcticus.orf.prot


blastp -query Electrona_antarctica.orf.prot -db Mito_DB.prot -outfmt 6 -max_target_seqs 1 -num_threads 8 -out Electrona_antarctica.vs.MitoDb.blastp
blastp -query Borostomias_antarcticus.orf.prot -db Mito_DB.prot -outfmt 6 -max_target_seqs 1 -num_threads 8 -out Borostomias_antarcticus.vs.MitoDb.blastp

for gene in `cat list_genes_mito.txt` ; do
	good_ID=`grep "$gene\." Electrona_antarctica.vs.MitoDb.blastp | sort -k12 | tail -1 | cut -f1` 
	samtools faidx Electrona_antarctica.orf $good_ID | sed "s/>.*/>Electrona_antarctica/g" >> Alignments/$gene.fa

	good_ID_sec=`grep "$gene\." Borostomias_antarcticus.vs.MitoDb.blastp | sort -k12 | tail -1 | cut -f1` 
	samtools faidx Borostomias_antarcticus.orf $good_ID_sec | sed "s/>.*/>Borostomias_antarcticus/g" >> Alignments/$gene.fa
done


#Finally, Add Lepisosteus mito genes to be used outgroup


cd Alignments/

for file in *.fa ; do transeq $file $file.prot -table 2 ; sed -i 's/_1$//g' $file.prot  ; done
for file in *.prot ; do muscle5.1.linux_intel64 -align $file -output $file.aln ; done
for file in *.aln ; do trimal -in $file -gt 0.8 -cons 70 -out $file.trimal ; done

cd ../

python3 AMAS.py concat -f fasta -d aa -i Alignments/*.trimal --part-format nexus
sed -i 's/?/-/g' concatenated.out
mv concatenated.out Concatenated_mitochondrial_genes.aln

iqtree -s Concatenated_mitochondrial_genes.aln --seqtype AA -m mtVer+F+G4 -nt 8 -bb 1000

#####====================================================================================================================================
#####====================================================================================================================================
####====================================================================================================================================
#####====================================================================================================================================
#####====================================================================================================================================
#####========================================== Accessory scripts ======================================================================
#####====================================================================================================================================
#####====================================================================================================================================
####====================================================================================================================================
#####====================================================================================================================================
#####====================================================================================================================================


================================================================================
================================================================================
========== download_and_launch_map.clupea.sh  ==================================
================================================================================
================================================================================

#!/bin/bash

module load SRA-Toolkit

curr_species=$1

IFS=$'\n'

grep "$curr_species" Table_reads.tsv > $curr_species.reads.tsv

for line in `cat $curr_species.reads.tsv` ; do
	SRA_accession=`echo "$line" | cut -f2`
	read_type=`echo "$line" | cut -f3`

	fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip $SRA_accession

	if [[ $curr_species == "Clupea_harengus" ]] ; then 
		sbatch --qos=1week -c 16 --mem=50G map_reads.clupea.clupea.sh $curr_species $SRA_accession $read_type
	else
		sbatch --qos=1week -c 16 --mem=50G map_reads.clupea.sh $curr_species $SRA_accession $read_type

	fi

done

================================================================================
================================================================================
================================================================================
================================================================================




================================================================================
================================================================================
========== map_reads.clupea.sh  ==================================
================================================================================
================================================================================


#!/bin/bash


#SBATCH --job-name=map_mit   # Job name

eval "$(conda shell.bash hook)"

curr_species=$1
SRA_accession=$2
read_type=$3

if [[ $read_type == "short" ]] ; then 

	conda activate miniprot
	module load fastp

	fastp -i "${SRA_accession}_pass_1.fastq.gz" -I "${SRA_accession}_pass_2.fastq.gz" -o "${SRA_accession}_pass_1.fastp.fastq.gz" -O "${SRA_accession}_pass_2.fastp.fastq.gz" --thread 16
	rm "${SRA_accession}_pass_1.fastq.gz" "${SRA_accession}_pass_2.fastq.gz"

	module purge ; module load bwa-mem2

	bwa-mem2 mem -t 16 Mitochondrial_Genomes/Clupea_harengus.mt.fa "${SRA_accession}_pass_1.fastp.fastq.gz" "${SRA_accession}_pass_2.fastp.fastq.gz" > "${SRA_accession}.Clupea_harengus.sam"
	samtools sort "${SRA_accession}.Clupea_harengus.sam" -o "${SRA_accession}.Clupea_harengus.sorted.bam" -@ 16  ; samtools index "${SRA_accession}.Clupea_harengus.sorted.bam" -@ 16 ; rm "${SRA_accession}.Clupea_harengus.sam"

	total_reads=`samtools view -c $SRA_accession.Clupea_harengus.sorted.bam`
	mapped_reads=`samtools view -c -F 260 $SRA_accession.Clupea_harengus.sorted.bam`
	mean_depth=`samtools depth -a $SRA_accession.Clupea_harengus.sorted.bam -q 0 | awk '{sum+=$3} END { print "Average = ",sum/NR}'`
	echo "Clupea_harengus,$SRA_accession,$total_reads,$mapped_reads,$mean_depth" > Mapping.$SRA_accession.Clupea_harengus.csv
	rm $SRA_accession.Clupea_harengus*

	bwa-mem2 mem -t 16 Mitochondrial_Genomes/$curr_species.mt.fa "${SRA_accession}_pass_1.fastp.fastq.gz" "${SRA_accession}_pass_2.fastp.fastq.gz" > $SRA_accession.$curr_species.sam
	samtools sort $SRA_accession.$curr_species.sam -o $SRA_accession.$curr_species.sorted.bam -@ 16  ; samtools index $SRA_accession.$curr_species.sorted.bam -@ 16 ; rm $SRA_accession.$curr_species.sam

	total_reads=`samtools view -c $SRA_accession.$curr_species.sorted.bam`
	mapped_reads=`samtools view -c -F 260 $SRA_accession.$curr_species.sorted.bam`
	mean_depth=`samtools depth -a $SRA_accession.$curr_species.sorted.bam -q 0 | awk '{sum+=$3} END { print "Average = ",sum/NR}'`
	echo "$curr_species,$SRA_accession,$total_reads,$mapped_reads,$mean_depth" > Mapping.$SRA_accession.$curr_species.csv
	rm $SRA_accession.$curr_species*


	rm "${SRA_accession}_pass_1.fastp.fastq.gz" "${SRA_accession}_pass_2.fastp.fastq.gz"



else 

	conda activate align_long_env
	module purge ; module load minimap2

	fastplong -i "${SRA_accession}_pass.fastq.gz" -o "${SRA_accession}_pass.fastp.fastq.gz" --thread 16
	rm "${SRA_accession}_pass.fastq.gz"

	minimap2 -t 16 -a Mitochondrial_Genomes/Clupea_harengus.mt.fa "${SRA_accession}_pass.fastp.fastq.gz" > $SRA_accession.Clupea_harengus.sam 
	minimap2 -t 16 -a Mitochondrial_Genomes/$curr_species.mt.fa "${SRA_accession}_pass.fastp.fastq.gz" > $SRA_accession.$curr_species.sam 

	conda deactivate ; conda activate miniprot

	samtools sort $SRA_accession.Clupea_harengus.sam -o $SRA_accession.Clupea_harengus.sorted.bam -@ 16  ; samtools index $SRA_accession.Clupea_harengus.sorted.bam -@ 16 ; rm $SRA_accession.Clupea_harengus.sam
	samtools sort $SRA_accession.$curr_species.sam -o $SRA_accession.$curr_species.sorted.bam -@ 16  ; samtools index $SRA_accession.$curr_species.sorted.bam -@ 16 ; rm $SRA_accession.$curr_species.sam


	total_reads=`samtools view -c $SRA_accession.Clupea_harengus.sorted.bam`
	mapped_reads=`samtools view -c -F 260 $SRA_accession.Clupea_harengus.sorted.bam`
	mean_depth=`samtools depth -a $SRA_accession.Clupea_harengus.sorted.bam -q 0 | awk '{sum+=$3} END { print "Average = ",sum/NR}'`
	echo "Clupea_harengus,$SRA_accession,$total_reads,$mapped_reads,$mean_depth" > Mapping.$SRA_accession.Clupea_harengus.csv
	rm $SRA_accession.Clupea_harengus*


	total_reads=`samtools view -c $SRA_accession.$curr_species.sorted.bam`
	mapped_reads=`samtools view -c -F 260 $SRA_accession.$curr_species.sorted.bam`
	mean_depth=`samtools depth -a $SRA_accession.$curr_species.sorted.bam -q 0 | awk '{sum+=$3} END { print "Average = ",sum/NR}'`
	echo "$curr_species,$SRA_accession,$total_reads,$mapped_reads,$mean_depth" > Mapping.$SRA_accession.$curr_species.csv
	rm $SRA_accession.$curr_species*


	rm "${SRA_accession}_pass.fastp.fastq.gz"



fi


================================================================================
================================================================================
================================================================================
================================================================================

