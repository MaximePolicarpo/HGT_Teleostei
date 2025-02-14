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



#Extract or download mitochondrial genomes

#Clupea_harengus
mt_scaff=`grep "^MT	" Clupea_harengus/GCF_900700415.2_Ch_v2.0.2_assembly_report.txt | cut -f7`
samtools faidx Clupea_harengus/GCF_900700415.2_Ch_v2.0.2_genomic.fna $mt_scaff > Mitochondrial_Genomes/Clupea_harengus.mt.fa

#Hypomesus_transpacificus
https://www.ncbi.nlm.nih.gov/nuccore/NC_072210.1


#Anoplopoma_fimbria

https://www.ncbi.nlm.nih.gov/nuccore/NC_018119



#Channa_striata
https://www.ncbi.nlm.nih.gov/nuccore/NC_032037


#Oncorhynchus_mykiss
mt_scaff=`grep "^MT	" Oncorhynchus_mykiss/GCF_013265735.2_USDA_OmykA_1.1_assembly_report.txt | cut -f7`
samtools faidx Oncorhynchus_mykiss/GCF_013265735.2_USDA_OmykA_1.1_genomic.fna $mt_scaff > Mitochondrial_Genomes/Oncorhynchus_mykiss.mt.fa


#Oryzias_melastigma

mt_scaff=`grep "^MT	" Oryzias_melastigma/GCF_002922805.2_ASM292280v2_assembly_report.txt | cut -f7`
samtools faidx Oryzias_melastigma/GCF_002922805.2_ASM292280v2_genomic.fna  $mt_scaff > Mitochondrial_Genomes/Oryzias_melastigma.mt.fa


#Salarias_fasciatus

mt_scaff=`grep "^MT	" Salarias_fasciatus/GCF_902148845.1_fSalaFa1.1_assembly_report.txt | cut -f7`
samtools faidx Salarias_fasciatus/GCF_902148845.1_fSalaFa1.1_genomic.fna  $mt_scaff > Mitochondrial_Genomes/Salarias_fasciatus.mt.fa


#Borostomias_antarcticus
mt_scaff=`grep "^scaffold_MT_1" Borostomias_antarcticus/GCA_949987555.1_fBorAnt1.1_assembly_report.txt | cut -f5`
samtools faidx Borostomias_antarcticus/GCA_949987555.1_fBorAnt1.1_genomic.fna $mt_scaff > Mitochondrial_Genomes/Borostomias_antarcticus.mt.fa



#Siniperca_chuatsi

mt_scaff=`grep "^MT	" Siniperca_chuatsi/GCF_020085105.1_ASM2008510v1_assembly_report.txt | cut -f7`
samtools faidx Siniperca_chuatsi/GCF_020085105.1_ASM2008510v1_genomic.fna  $mt_scaff > Mitochondrial_Genomes/Siniperca_chuatsi.mt.fa



#Takifugu_flavidus

mt_scaff=`grep "^MT	" Takifugu_flavidus/GCF_003711565.1_ASM371156v2_assembly_report.txt  | cut -f7`
samtools faidx Takifugu_flavidus/GCF_003711565.1_ASM371156v2_genomic.fna $mt_scaff > Mitochondrial_Genomes/Takifugu_flavidus.mt.fa



#Xiphophorus_couchianus

mt_scaff=`grep "^MT	" Xiphophorus_couchianus/GCF_001444195.1_X_couchianus-1.0_assembly_report.txt | cut -f7`
samtools faidx Xiphophorus_couchianus/GCF_001444195.1_X_couchianus-1.0_genomic.fna $mt_scaff > Mitochondrial_Genomes/Xiphophorus_couchianus.mt.fa




## Now create bwamem index for every mitochondrial genomes


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


cd Mitochondrial_Genomes/

for file in *.fa ; do 
	species=`echo "$file" | sed 's/.mt.fa//g'`
	sed "s/>.*/>$species/g" $file >> Combined_all_species.fa
done
bwa-mem2 index Combined_all_species.fa

cd ../



cut -f1 Mitochondrial_Genomes/Combined_all_species.fa.fai > names
cut -f2 Mitochondrial_Genomes/Combined_all_species.fa.fai | sed 's/^/1,/g' | tr ',' '\t' > lengths
paste -d '\t' names lengths > mitochondrial_genomes.bed

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




module purge ; module load SRA-Toolkit ; fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip SRR26322996 ; sbatch --qos=6hours -c 16 --mem=50G map_reads.mitochondria.sh SRR26322996 long
module purge ; module load SRA-Toolkit ; fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip SRR26322995 ; sbatch --qos=6hours -c 16 --mem=50G map_reads.mitochondria.sh SRR26322995 long
module purge ; module load SRA-Toolkit ; fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip SRR30009981 ; sbatch --qos=6hours -c 16 --mem=50G map_reads.mitochondria.sh SRR30009981 short
module purge ; module load SRA-Toolkit ; fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip SRR30009984 ; sbatch --qos=6hours -c 16 --mem=50G map_reads.mitochondria.sh SRR30009984 short
module purge ; module load SRA-Toolkit ; fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip SRR17921806 ; sbatch --qos=6hours -c 16 --mem=50G map_reads.mitochondria.sh SRR17921806 long
module purge ; module load SRA-Toolkit ; fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip SRR17921807 ; sbatch --qos=6hours -c 16 --mem=50G map_reads.mitochondria.sh SRR17921807 short
module purge ; module load SRA-Toolkit ; fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip SRR17799534 ; sbatch --qos=6hours -c 16 --mem=50G map_reads.mitochondria.sh SRR17799534 long
module purge ; module load SRA-Toolkit ; fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip SRR17799535 ; sbatch --qos=6hours -c 16 --mem=50G map_reads.mitochondria.sh SRR17799535 short
module purge ; module load SRA-Toolkit ; fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip SRR27482199 ; sbatch --qos=6hours -c 16 --mem=50G map_reads.mitochondria.sh SRR27482199 short
module purge ; module load SRA-Toolkit ; fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip SRR27482200 ; sbatch --qos=6hours -c 16 --mem=50G map_reads.mitochondria.sh SRR27482200 short
module purge ; module load SRA-Toolkit ; fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip ERR3566168 ; sbatch --qos=6hours -c 16 --mem=50G map_reads.mitochondria.sh ERR3566168 short
module purge ; module load SRA-Toolkit ; fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip ERR3566167 ; sbatch --qos=6hours -c 16 --mem=50G map_reads.mitochondria.sh ERR3566167 short
module purge ; module load SRA-Toolkit ; fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip ERR3629070 ; sbatch --qos=6hours -c 16 --mem=50G map_reads.mitochondria.sh ERR3629070 long
module purge ; module load SRA-Toolkit ; fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip SRR2131027 ; sbatch --qos=6hours -c 16 --mem=50G map_reads.mitochondria.sh SRR2131027 short
module purge ; module load SRA-Toolkit ; fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip SRR2131118 ; sbatch --qos=6hours -c 16 --mem=50G map_reads.mitochondria.sh SRR2131118 short
module purge ; module load SRA-Toolkit ; fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip SRR7464412 ; sbatch --qos=6hours -c 16 --mem=50G map_reads.mitochondria.sh SRR7464412 long
module purge ; module load SRA-Toolkit ; fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip SRR7464416 ; sbatch --qos=6hours -c 16 --mem=50G map_reads.mitochondria.sh SRR7464416 long
module purge ; module load SRA-Toolkit ; fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip SRR29199351 ; sbatch --qos=6hours -c 16 --mem=50G map_reads.mitochondria.sh SRR29199351 short
module purge ; module load SRA-Toolkit ; fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip SRR29199350 ; sbatch --qos=6hours -c 16 --mem=50G map_reads.mitochondria.sh SRR29199350 short
module purge ; module load SRA-Toolkit ; fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip SRR29199348 ; sbatch --qos=6hours -c 16 --mem=50G map_reads.mitochondria.sh SRR29199348 long
module purge ; module load SRA-Toolkit ; fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip SRR15111159 ; sbatch --qos=6hours -c 16 --mem=50G map_reads.mitochondria.sh SRR15111159 short
module purge ; module load SRA-Toolkit ; fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip SRR7881547 ; sbatch --qos=6hours -c 16 --mem=50G map_reads.mitochondria.sh SRR7881547 long
module purge ; module load SRA-Toolkit ; fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip ERR10936409 ; sbatch --qos=6hours -c 16 --mem=50G map_reads.mitochondria.sh ERR10936409 short
module purge ; module load SRA-Toolkit ; fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip ERR10934070 ; sbatch --qos=6hours -c 16 --mem=50G map_reads.mitochondria.sh ERR10934070 long
module purge ; module load SRA-Toolkit ; fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip ERR3629069 ; sbatch --qos=6hours -c 16 --mem=50G map_reads.mitochondria.sh ERR3629069 long
module purge ; module load SRA-Toolkit ; fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip ERR3629068 ; sbatch --qos=6hours -c 16 --mem=50G map_reads.mitochondria.sh ERR3629068 long
module purge ; module load SRA-Toolkit ; fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip SRR6002358 ; sbatch --qos=6hours -c 16 --mem=50G map_reads.mitochondria.sh SRR6002358 short
module purge ; module load SRA-Toolkit ; fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip SRR6002359 ; sbatch --qos=6hours -c 16 --mem=50G map_reads.mitochondria.sh SRR6002359 short
module purge ; module load SRA-Toolkit ; fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip SRR6002363 ; sbatch --qos=6hours -c 16 --mem=50G map_reads.mitochondria.sh SRR6002363 short
module purge ; module load SRA-Toolkit ; fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip SRR2127225 ; sbatch --qos=6hours -c 16 --mem=50G map_reads.mitochondria.sh SRR2127225 short
module purge ; module load SRA-Toolkit ; fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip SRR7881551 ; sbatch --qos=6hours -c 16 --mem=50G map_reads.mitochondria.sh SRR7881551 short
module purge ; module load SRA-Toolkit ; fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip SRR7881549 ; sbatch --qos=6hours -c 16 --mem=50G map_reads.mitochondria.sh SRR7881549 long
module purge ; module load SRA-Toolkit ; fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip SRR7881546 ; sbatch --qos=6hours -c 16 --mem=50G map_reads.mitochondria.sh SRR7881546 long
module purge ; module load SRA-Toolkit ; fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip SRR21199007 ; sbatch --qos=6hours -c 16 --mem=50G map_reads.mitochondria.sh SRR21199007 short
module purge ; module load SRA-Toolkit ; fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip SRR21199008 ; sbatch --qos=6hours -c 16 --mem=50G map_reads.mitochondria.sh SRR21199008 long
module purge ; module load SRA-Toolkit ; fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip SRR21199004 ; sbatch --qos=6hours -c 16 --mem=50G map_reads.mitochondria.sh SRR21199004 long
module purge ; module load SRA-Toolkit ; fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip ERR10936410 ; sbatch --qos=6hours -c 16 --mem=50G map_reads.mitochondria.sh ERR10936410 short
module purge ; module load SRA-Toolkit ; fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip ERR10934071 ; sbatch --qos=6hours -c 16 --mem=50G map_reads.mitochondria.sh ERR10934071 long
module purge ; module load SRA-Toolkit ; fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip ERR10934069 ; sbatch --qos=6hours -c 16 --mem=50G map_reads.mitochondria.sh ERR10934069 long




rm Compet_mapping_results.tsv ; for file in *.idxstats ; do SRA_ID=`echo "$file" | cut -f1 -d "."` ; sed "s/^/$SRA_ID,/g" $file | tr ',' '\t' | grep -v "^\*"  >> Compet_mapping_results.tsv ; done

#sed 's/Average =  //g' Mapping.*.csv > Summary_mapping.csv







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


#to annotate:

getorf -sequence ../Mitochondrial_Genomes/Borostomias_antarcticus.mt.fa -outseq Borostomias_antarcticus.orf -minsize 168 -find 3 -table 2

transeq Borostomias_antarcticus.orf Borostomias_antarcticus.orf.prot -table 2 ; sed -i 's/_1 / /g' Borostomias_antarcticus.orf.prot


blastp -query Borostomias_antarcticus.orf.prot -db Mito_DB.prot -outfmt 6 -max_target_seqs 1 -num_threads 8 -out Borostomias_antarcticus.vs.MitoDb.blastp

for gene in `cat list_genes_mito.txt` ; do
	good_ID_sec=`grep "$gene\." Borostomias_antarcticus.vs.MitoDb.blastp | sort -k12 | tail -1 | cut -f1` 
	samtools faidx Borostomias_antarcticus.orf $good_ID_sec | sed "s/>.*/>Borostomias_antarcticus/g" >> Alignments/$gene.fa
done

#Finally, add lepisosteus as outgroup





cd Alignments/

for file in *.fa ; do transeq $file $file.prot -table 2 ; sed -i 's/_1$//g' $file.prot  ; done
for file in *.prot ; do muscle5.1.linux_intel64 -align $file -output $file.aln ; done
for file in *.aln ; do trimal -in $file -gt 0.8 -cons 70 -out $file.trimal ; done

cd ../

python3 AMAS/amas/AMAS.py concat -f fasta -d aa -i Alignments/*.trimal --part-format nexus
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
========== map_reads.mitochondria.sh  ==================================
================================================================================
================================================================================



#!/bin/bash


#SBATCH --job-name=map_mit   # Job name

eval "$(conda shell.bash hook)"

SRA_accession=$1
read_type=$2

if [[ $read_type == "short" ]] ; then 

	conda activate miniprot
	module load fastp

	fastp -i "${SRA_accession}_pass_1.fastq.gz" -I "${SRA_accession}_pass_2.fastq.gz" -o "${SRA_accession}_pass_1.fastp.fastq.gz" -O "${SRA_accession}_pass_2.fastp.fastq.gz" --thread 16
	rm "${SRA_accession}_pass_1.fastq.gz" "${SRA_accession}_pass_2.fastq.gz"

	module purge ; module load bwa-mem2

	bwa-mem2 mem -t 16 Mitochondrial_Genomes/Combined_all_species.fa "${SRA_accession}_pass_1.fastp.fastq.gz" "${SRA_accession}_pass_2.fastp.fastq.gz" > $SRA_accession.sam
	rm "${SRA_accession}_pass_1.fastp.fastq.gz" "${SRA_accession}_pass_2.fastp.fastq.gz"

	samtools sort $SRA_accession.sam  -o $SRA_accession.sorted.bam  -@ 16  ; samtools index $SRA_accession.sorted.bam -@ 16
	rm $SRA_accession.sam

	samtools idxstats $SRA_accession.sorted.bam > $SRA_accession.idxstats
	for species in `cat species_list.txt` ; do samtools depth -a $SRA_accession.sorted.bam -q 0 -r $species | awk '{sum+=$3} END { print "Average = ",sum/NR}' | sed "s/^/$species,/g" ; done > $SRA_accession.depth

	rm $SRA_accession.sorted.bam
	

else 

	conda activate align_long_env
	module purge ; module load minimap2

	fastplong -i "${SRA_accession}_pass.fastq.gz" -o "${SRA_accession}_pass.fastp.fastq.gz" --thread 16
	rm "${SRA_accession}_pass.fastq.gz"

	minimap2 -t 16 -a Mitochondrial_Genomes/Combined_all_species.fa "${SRA_accession}_pass.fastp.fastq.gz" > $SRA_accession.sam
	rm "${SRA_accession}_pass.fastp.fastq.gz"

	conda deactivate ; conda activate miniprot

	samtools sort $SRA_accession.sam  -o $SRA_accession.sorted.bam  -@ 16  ; samtools index $SRA_accession.sorted.bam -@ 16
	rm $SRA_accession.sam
	
	samtools idxstats $SRA_accession.sorted.bam > $SRA_accession.idxstats
	for species in `cat species_list.txt` ; do samtools depth -a $SRA_accession.sorted.bam -q 0 -r $species | awk '{sum+=$3} END { print "Average = ",sum/NR}' | sed "s/^/$species,/g" ; done > $SRA_accession.depth

	rm $SRA_accession.sorted.bam
	



fi
