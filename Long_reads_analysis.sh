#####====================================================================================================================================
#####====================================================================================================================================
####====================================================================================================================================
#####====================================================================================================================================
#####====================================================================================================================================
#####========================================== Download Clupea and Hypomesus and Borstoma (control) reads ==============================
#####====================================================================================================================================
#####====================================================================================================================================
####====================================================================================================================================
#####====================================================================================================================================
#####====================================================================================================================================


module purge ; module load SRA-Toolkit

fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip SRR26322996
fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip SRR26322995
fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip SRR17921806
fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip SRR17799534
fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip SRR32203436
fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip ERR10934070
fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip ERR10934071
fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip ERR10934069

#####====================================================================================================================================
#####====================================================================================================================================
#####========================================== Map reads to the hypomesus + clupea genomes =============================================
#####====================================================================================================================================
#####====================================================================================================================================


module purge ; module load bwa-mem2

bwa-mem2 index GCF_021917145.1_fHypTra1_genomic.fna
bwa-mem2 index GCF_900700415.2_Ch_v2.0.2_genomic.fna


sbatch --qos=1week -c 16 --mem=50G filter_map_reads.sh SRR26322996 long
sbatch --qos=1week -c 16 --mem=50G filter_map_reads.sh SRR26322995 long
sbatch --qos=1week -c 16 --mem=50G filter_map_reads.sh SRR17921806 long
sbatch --qos=1week -c 16 --mem=50G filter_map_reads.sh SRR17799534 long
sbatch --qos=1week -c 16 --mem=50G filter_map_reads.sh SRR32203436 long
sbatch --qos=1week -c 16 --mem=50G filter_map_reads.sh ERR10934070 long
sbatch --qos=1week -c 16 --mem=50G filter_map_reads.sh ERR10934071 long
sbatch --qos=1week -c 16 --mem=50G filter_map_reads.sh ERR10934069 long


#####====================================================================================================================================
#####====================================================================================================================================
#####========================================== Compute read depth on HGT regions =======================================================
#####====================================================================================================================================
#####====================================================================================================================================

#Extract scaffolds with at-least one horizontal transfer gene on it in both Clupea and Hypomesus

grep "Clupea_harengus" table_genes_loc.tsv  | cut -f4 | sed 's/:.*//g' | sort | uniq  > list_Clupea_scaffolds.txt
grep "Hypomesus_transpacificus" table_genes_loc.tsv  | cut -f4 | sed 's/:.*//g' | sort | uniq  > list_Hypomesus_scaffolds.txt

for Scaffold_name in `cat list_Clupea_scaffolds.txt` ; do
	length_scaffold=`grep "$Scaffold_name	" GCF_900700415.2_Ch_v2.0.2_genomic.fna.fai | cut -f2`
	echo "$Scaffold_name,1,$length_scaffold" | tr ',' '\t' > $Scaffold_name.bed
done


for Scaffold_name in `cat list_Hypomesus_scaffolds.txt` ; do
	length_scaffold=`grep "$Scaffold_name	" GCF_021917145.1_fHypTra1_genomic.fna.fai | cut -f2`
	echo "$Scaffold_name,1,$length_scaffold" | tr ',' '\t' > $Scaffold_name.bed
done


#Launch read depth calculation on Clupea scaffolds

rm -r Depth_files.Clupea_harengus ; mkdir Depth_files.Clupea_harengus
for SRA_accession in `cat list_SRA_reads.txt` ; do
	for Scaffold_name in `cat list_Clupea_scaffolds.txt` ; do
		sbatch --qos=30min -c 3 --mem=8G compute_depth.sh $SRA_accession $Scaffold_name Clupea_harengus
	done
done


cd Depth_files.Clupea_harengus
ls -l | grep ".txt" | sed 's/.* //g' | cut -f3,4 -d '.' | sort | uniq > list_scaffold

for curr_scaff in `cat list_scaffold` ; do


	for file in Depth.*$curr_scaff*.txt ; do
		SRR_accession=`echo "$file" | cut -f2 -d "."`
		sed -i "s/^/$SRR_accession,/g" $file
	done


	cat Depth.*$curr_scaff*.txt | tr ',' '\t' > $curr_scaff.all.depth

	rm Depth.*$curr_scaff*
done


#Launch read depth calculation on Hypomesus scaffolds

rm -rf Depth_files.Hypomesus_transpacificus ; mkdir Depth_files.Hypomesus_transpacificus
for SRA_accession in `cat list_SRA_reads.txt` ; do
	for Scaffold_name in `cat list_Hypomesus_scaffolds.txt` ; do
		sbatch --qos=30min -c 3 --mem=8G compute_depth.sh $SRA_accession $Scaffold_name Hypomesus_transpacificus
	done
done




#####====================================================================================================================================
#####====================================================================================================================================
#####========================================== Count the number of long reads spanning entire HGT genes ==============================
#####====================================================================================================================================
#####====================================================================================================================================

IFS=$'\n'
rm -rf Reads_spanning_regions_tables/ ; mkdir Reads_spanning_regions_tables/
for line in `cat Table_longreads_count.tsv` ; do

	#Define variables
	species=`echo "$line" | cut -f2`
	genome_fai=`ls -l Genomic_data/$species/ | grep "genomic.fna.fai" | sed 's/.* //g'`
	
	SRA_accession=`echo "$line" | cut -f5`
	Region=`echo "$line" | cut -f4`

	curr_scaffold=`echo "$Region" | cut -f1 -d ":"`
	scaffold_length=`grep "$curr_scaffold	" Genomic_data/$species/$genome_fai | cut -f2`
	start=`echo "$Region" | cut -f2 -d ":" | cut -f1 -d "-"`
	stop=`echo "$Region" | cut -f2 -d ":" | cut -f2 -d "-"`
	


	start_ext=$((start - 1000000))
	if [ "$start_ext" -le 0 ]; then start_ext=1 ; fi
	stop_ext=$((stop + 1000000))
	if [ "$stop_ext" -ge "$scaffold_length" ] ; then stop_ext=$((scaffold_length - 10)) ; fi
	Region_extended=`echo "$curr_scaffold:$start_ext-$stop_ext"`




	#Extract the exact position of each read mapped 
	samtools view $SRA_accession.$species.sorted.bam $Region_extended -@ 8 | awk '{
    cigar = $6;
    scaffold = $3;
    start = $4;
    len = 0;
    
    while (match(cigar, /([0-9]+)([MDN])/)) {
        len += substr(cigar, RSTART, RLENGTH-1);
        cigar = substr(cigar, RSTART+RLENGTH);
    }
    
    end = start + len - 1;
    print $1, scaffold, start, end;
	}' | tr ' ' '\t' > curr_reads_positions.tsv



	awk -v start="$start" -v stop="$stop" '$3 <= start && $4 >= stop' curr_reads_positions.tsv | sed "s/^/$Region,/g" | tr ',' '\t' > curr_reads_positions.named.tsv

	nb_reads_full_span=`wc -l < curr_reads_positions.named.tsv`

	cat curr_reads_positions.named.tsv >> Reads_spanning_regions_tables/$species.$Region
	echo "$nb_reads_full_span"

	rm curr_reads_positions.named.tsv ; rm curr_reads_positions.tsv

done


#####====================================================================================================================================
#####====================================================================================================================================
#####========================================== Align long reads for each HOG  =========================================================
#####====================================================================================================================================
#####====================================================================================================================================

awk '$2 == "Clupea_harengus"' Table_longreads_count.final.tsv > Table_longreads_count.clupea.txt

cut -f1 Table_longreads_count.clupea.txt | sort | uniq > list_HOG.txt

for curr_HOG in `cat list_HOG.txt` ; do

	echo "$curr_HOG"

	hypo_line_nb=`grep "$curr_HOG" Table_longreads_count.clupea.txt | awk '$6 == "Hypomesus_transpacificus"' | awk '$7 > 0' | awk '$7 < 350' | wc -l`

	if [ $hypo_line_nb -ge 1 ] ; then

		rm -rf $curr_HOG ; mkdir $curr_HOG

		curr_XM=`grep "$curr_HOG" Table_longreads_count.clupea.txt | awk '$6 == "Hypomesus_transpacificus"' | awk '$7 > 0' | awk '$7 < 350' | head -1 | cut -f3`
		full_XM_name=`grep "$curr_XM" Clupea_harengus.prot.nostop | sed 's/>//g'`
		samtools faidx Clupea_harengus.prot.nostop $full_XM_name > $curr_HOG/$curr_XM.prot

		Hypomesus_SRR=`grep "$curr_HOG" Table_longreads_count.clupea.txt | grep "$curr_XM" | awk '$6 == "Hypomesus_transpacificus"' | awk '$7 > 0' | head -1 | cut -f5`
		Clupea_SRR=`grep "$curr_HOG" Table_longreads_count.clupea.txt | grep "$curr_XM" | awk '$6 == "Clupea_harengus"' | awk '$7 > 0' | head -1 | cut -f5`


		curr_scaffold=`grep "$curr_XM" Table_longreads_count.clupea.txt | cut -f4 | sort | uniq | cut -f1 -d ":"`
		start=`grep "$curr_XM" Table_longreads_count.clupea.txt | cut -f4 | sort | uniq | cut -f2 -d ":" | cut -f1 -d "-"`
		stop=`grep "$curr_XM" Table_longreads_count.clupea.txt | cut -f4 | sort | uniq | cut -f2 -d ":" | cut -f2 -d "-"`
		scaffold_length=`grep "$curr_scaffold	" GCF_900700415.2_Ch_v2.0.2_genomic.fna.fai | cut -f2`
		start_ext=$((start - 1000000))
		if [ "$start_ext" -le 0 ]; then start_ext=1 ; fi
		stop_ext=$((stop + 1000000))
		if [ "$stop_ext" -ge "$scaffold_length" ] ; then stop_ext=$((scaffold_length - 10)) ; fi
		Region_extended=`echo "$curr_scaffold:$start_ext-$stop_ext"`


		Hypomesus_read_name=`grep "$Hypomesus_SRR" ../Reads_spanning_regions_tables/Clupea_harengus.$curr_scaffold\:$start-$stop | grep -v "split-by-adapter" | head -1 | cut -f2`
		Clupea_read_name=`grep "$Clupea_SRR" ../Reads_spanning_regions_tables/Clupea_harengus.$curr_scaffold\:$start-$stop | grep -v "split-by-adapter" | head -1 | cut -f2`

		echo ">$Hypomesus_read_name" > $curr_HOG/Hypomesus_transpacificus.fa
		samtools view ../$Hypomesus_SRR.Clupea_harengus.sorted.bam $Region_extended -@ 8 | grep "$Hypomesus_read_name" | cut -f10 | awk '{print $0,length}' | tr ' ' ','  | sort -k2 -n -t ',' | tail -1 | cut -f1 -d ","  >> $curr_HOG/Hypomesus_transpacificus.fa
		
		echo ">$Clupea_read_name" > $curr_HOG/Clupea_harengus.fa
		samtools view ../$Clupea_SRR.Clupea_harengus.sorted.bam $Region_extended -@ 8 | grep "$Clupea_read_name" | cut -f10 | awk '{print $0,length}' | tr ' ' ','  | sort -k2 -n -t ',' | tail -1 | cut -f1 -d "," >> $curr_HOG/Clupea_harengus.fa

		strand_hypo=`exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/$curr_XM.prot $curr_HOG/Hypomesus_transpacificus.fa | grep "vulgar" | cut -f9 -d " "`
		if [ $strand_hypo == "-" ] ; then revseq $curr_HOG/Hypomesus_transpacificus.fa $curr_HOG/Hypomesus_transpacificus.rev ; sed -i 's/ Reversed://g' $curr_HOG/Hypomesus_transpacificus.rev ; mv $curr_HOG/Hypomesus_transpacificus.rev $curr_HOG/Hypomesus_transpacificus.fa ; fi

		strand_clupea=`exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/$curr_XM.prot $curr_HOG/Clupea_harengus.fa | grep "vulgar" | cut -f9 -d " "`
		if [ $strand_clupea == "-" ] ; then revseq $curr_HOG/Clupea_harengus.fa $curr_HOG/Clupea_harengus.rev ; sed -i 's/ Reversed://g' $curr_HOG/Clupea_harengus.rev ; mv $curr_HOG/Clupea_harengus.rev $curr_HOG/Clupea_harengus.fa ; fi

		exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/$curr_XM.prot $curr_HOG/Hypomesus_transpacificus.fa > $curr_HOG/Hypomesus_transpacificus.exo
		exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/$curr_XM.prot $curr_HOG/Clupea_harengus.fa > $curr_HOG/Clupea_harengus.exo


	fi

done



### Now do the same but on Hypomesus genome for N5.HOG0047987 


awk '$2 == "Hypomesus_transpacificus"' Table_longreads_count.final.tsv > Table_longreads_count.hypomesus.txt

curr_HOG=N5.HOG0047987

rm -rf $curr_HOG ; mkdir $curr_HOG

curr_XM=XM_047020147
full_XM_name=`grep "$curr_XM" Hypomesus_transpacificus.prot.nostop | sed 's/>//g'`
samtools faidx Hypomesus_transpacificus.prot.nostop $full_XM_name > $curr_HOG/$curr_XM.prot


grep "$curr_HOG" Table_longreads_count.hypomesus.txt | grep "$curr_XM"

Hypomesus_SRR=SRR17921806
Clupea_SRR=SRR26322996

curr_scaffold=`grep "$curr_XM" Table_longreads_count.hypomesus.txt | cut -f4 | sort | uniq | cut -f1 -d ":"`
start=`grep "$curr_XM" Table_longreads_count.hypomesus.txt | cut -f4 | sort | uniq | cut -f2 -d ":" | cut -f1 -d "-"`
stop=`grep "$curr_XM" Table_longreads_count.hypomesus.txt | cut -f4 | sort | uniq | cut -f2 -d ":" | cut -f2 -d "-"`
scaffold_length=`grep "$curr_scaffold	" GCF_021917145.1_fHypTra1_genomic.fna.fai | cut -f2`
start_ext=$((start - 1000000))
if [ "$start_ext" -le 0 ]; then start_ext=1 ; fi
stop_ext=$((stop + 1000000))
if [ "$stop_ext" -ge "$scaffold_length" ] ; then stop_ext=$((scaffold_length - 10)) ; fi
Region_extended=`echo "$curr_scaffold:$start_ext-$stop_ext"`


grep "$Hypomesus_SRR" ../Reads_spanning_regions_tables/Hypomesus_transpacificus.$curr_scaffold\:$start-$stop
grep "$Clupea_SRR" ../Reads_spanning_regions_tables/Hypomesus_transpacificus.$curr_scaffold\:$start-$stop

Hypomesus_read_name=SRR17921806.669654
Clupea_read_name=SRR26322996.533176

echo ">$Hypomesus_read_name" > $curr_HOG/Hypomesus_transpacificus.fa
samtools view ../$Hypomesus_SRR.Hypomesus_transpacificus.sorted.bam $Region_extended -@ 8 | grep "$Hypomesus_read_name" | cut -f10 | awk '{print $0,length}' | tr ' ' ','  | sort -k2 -n -t ',' | tail -1 | cut -f1 -d ","  >> $curr_HOG/Hypomesus_transpacificus.fa

echo ">$Clupea_read_name" > $curr_HOG/Clupea_harengus.fa
samtools view ../$Clupea_SRR.Hypomesus_transpacificus.sorted.bam $Region_extended -@ 8 | grep "$Clupea_read_name" | cut -f10 | awk '{print $0,length}' | tr ' ' ','  | sort -k2 -n -t ',' | tail -1 | cut -f1 -d "," >> $curr_HOG/Clupea_harengus.fa


strand_hypo=`exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/$curr_XM.prot $curr_HOG/Hypomesus_transpacificus.fa | grep "vulgar" | cut -f9 -d " "`
if [ $strand_hypo == "-" ] ; then revseq $curr_HOG/Hypomesus_transpacificus.fa $curr_HOG/Hypomesus_transpacificus.rev ; sed -i 's/ Reversed://g' $curr_HOG/Hypomesus_transpacificus.rev ; mv $curr_HOG/Hypomesus_transpacificus.rev $curr_HOG/Hypomesus_transpacificus.fa ; fi

strand_clupea=`exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/$curr_XM.prot $curr_HOG/Clupea_harengus.fa | grep "vulgar" | cut -f9 -d " "`
if [ $strand_clupea == "-" ] ; then revseq $curr_HOG/Clupea_harengus.fa $curr_HOG/Clupea_harengus.rev ; sed -i 's/ Reversed://g' $curr_HOG/Clupea_harengus.rev ; mv $curr_HOG/Clupea_harengus.rev $curr_HOG/Clupea_harengus.fa ; fi

exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/$curr_XM.prot $curr_HOG/Hypomesus_transpacificus.fa > $curr_HOG/Hypomesus_transpacificus.exo
exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/$curr_XM.prot $curr_HOG/Clupea_harengus.fa > $curr_HOG/Clupea_harengus.exo



### Now do the same but on Hypomesus genome for N5.HOG0004480 


curr_HOG=N5.HOG0004480

rm -rf $curr_HOG ; mkdir $curr_HOG

Hypomesus_SRR=SRR17799534
Region_extended=NC_061079.1:3701703-3909259
grep "$Hypomesus_SRR" ../Reads_spanning_regions_tables/Hypomesus_transpacificus.NC_061079.1\:3801703-3809259
Hypomesus_read_name=SRR17799534.753713
echo ">$Hypomesus_read_name" > $curr_HOG/Hypomesus_transpacificus.fa
samtools view ../$Hypomesus_SRR.Hypomesus_transpacificus.sorted.bam $Region_extended -@ 8 | grep "$Hypomesus_read_name" | cut -f10 | awk '{print $0,length}' | tr ' ' ','  | sort -k2 -n -t ',' | tail -1 | cut -f1 -d ","  >> $curr_HOG/Hypomesus_transpacificus.fa
strand_hypo=`exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/XM_031565097.prot $curr_HOG/Hypomesus_transpacificus.fa | grep "vulgar" | cut -f9 -d " "`
if [ $strand_hypo == "-" ] ; then revseq $curr_HOG/Hypomesus_transpacificus.fa $curr_HOG/Hypomesus_transpacificus.rev ; sed -i 's/ Reversed://g' $curr_HOG/Hypomesus_transpacificus.rev ; mv $curr_HOG/Hypomesus_transpacificus.rev $curr_HOG/Hypomesus_transpacificus.fa ; fi
exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/XM_031565097.prot $curr_HOG/Hypomesus_transpacificus.fa > $curr_HOG/Hypomesus_transpacificus.exo



Clupea_SRR=SRR26322996
Region_extended=NC_045152.1:4732716-4946865
grep "$Clupea_SRR" ../Reads_spanning_regions_tables/Clupea_harengus.NC_045152.1\:4832716-4846865
Clupea_read_name=SRR26322996.158256
echo ">$Clupea_read_name" > $curr_HOG/Clupea_harengus.fa
samtools view ../$Clupea_SRR.Clupea_harengus.sorted.bam $Region_extended -@ 8 | grep "$Clupea_read_name" | cut -f10 | awk '{print $0,length}' | tr ' ' ','  | sort -k2 -n -t ',' | tail -1 | cut -f1 -d ","  >> $curr_HOG/Clupea_harengus.fa
strand_hypo=`exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/XM_031565097.prot $curr_HOG/Clupea_harengus.fa | grep "vulgar" | cut -f9 -d " "`
if [ $strand_hypo == "-" ] ; then revseq $curr_HOG/Clupea_harengus.fa $curr_HOG/Clupea_harengus.rev ; sed -i 's/ Reversed://g' $curr_HOG/Clupea_harengus.rev ; mv $curr_HOG/Clupea_harengus.rev $curr_HOG/Clupea_harengus.fa ; fi
exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/XM_031565097.prot $curr_HOG/Clupea_harengus.fa > $curr_HOG/Clupea_harengus.exo



awk '
    BEGIN {start=1; end=10500; count=0}
    /^>/ {print; next} 
    {
        line = ""
        for (i = 1; i <= length($0); i++) {
            count++
            if (count >= start && count <= end) {
                line = line "N"
            } else {
                line = line substr($0, i, 1)
            }
        }
        print line
    }
' Hypomesus_transpacificus.fa > Hypomesus_transpacificus.modif.fa
mv Hypomesus_transpacificus.modif.fa  Hypomesus_transpacificus.fa




strand_hypo=`exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/$curr_XM.prot $curr_HOG/Hypomesus_transpacificus.fa | grep "vulgar" | cut -f9 -d " "`
if [ $strand_hypo == "-" ] ; then revseq $curr_HOG/Hypomesus_transpacificus.fa $curr_HOG/Hypomesus_transpacificus.rev ; sed -i 's/ Reversed://g' $curr_HOG/Hypomesus_transpacificus.rev ; mv $curr_HOG/Hypomesus_transpacificus.rev $curr_HOG/Hypomesus_transpacificus.fa ; fi

strand_clupea=`exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/$curr_XM.prot $curr_HOG/Clupea_harengus.fa | grep "vulgar" | cut -f9 -d " "`
if [ $strand_clupea == "-" ] ; then revseq $curr_HOG/Clupea_harengus.fa $curr_HOG/Clupea_harengus.rev ; sed -i 's/ Reversed://g' $curr_HOG/Clupea_harengus.rev ; mv $curr_HOG/Clupea_harengus.rev $curr_HOG/Clupea_harengus.fa ; fi

exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/$curr_XM.prot $curr_HOG/Hypomesus_transpacificus.fa > $curr_HOG/Hypomesus_transpacificus.exo
exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/$curr_XM.prot $curr_HOG/Clupea_harengus.fa > $curr_HOG/Clupea_harengus.exo






#####====================================================================================================================================
#####====================================================================================================================================
#####========================================== Accessory scripts =======================================================================
#####====================================================================================================================================
#####====================================================================================================================================



================================================================================
================================================================================
========== compute_depth.sh  ==================================
================================================================================
================================================================================


#!/bin/bash


#SBATCH --job-name=depth   # Job name

module load SAMtools

SRA_accession=$1
Scaffold_name=$2
Species_name=$3


samtools depth -a $SRA_accession.$Species_name.sorted.bam -q 0 -b $Scaffold_name.bed --min-MQ 30 --min-BQ 30 -g SECONDARY,SUPPLEMENTARY | awk '$3 < 600' > Depth_files.$Species_name/Depth.$SRA_accession.$Scaffold_name.$Species_name.txt


================================================================================
================================================================================
================================================================================
================================================================================




================================================================================
================================================================================
========== filter_map_reads.sh  ==================================
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

	bwa-mem2 mem -t 16 GCF_021917145.1_fHypTra1_genomic.fna "${SRA_accession}_pass_1.fastp.fastq.gz" "${SRA_accession}_pass_2.fastp.fastq.gz" > "${SRA_accession}.Hypomesus_transpacificus.sam"
	bwa-mem2 mem -t 16 GCF_900700415.2_Ch_v2.0.2_genomic.fna "${SRA_accession}_pass_1.fastp.fastq.gz" "${SRA_accession}_pass_2.fastp.fastq.gz" > "${SRA_accession}.Clupea_harengus.sam"

	samtools sort "${SRA_accession}.Hypomesus_transpacificus.sam" -o "${SRA_accession}.Hypomesus_transpacificus.sorted.bam" -@ 16  ; samtools index "${SRA_accession}.Hypomesus_transpacificus.sorted.bam" -@ 16 ; rm "${SRA_accession}.Hypomesus_transpacificus.sam"
	samtools sort "${SRA_accession}.Clupea_harengus.sam" -o "${SRA_accession}.Clupea_harengus.sorted.bam" -@ 16  ; samtools index "${SRA_accession}.Clupea_harengus.sorted.bam" -@ 16 ; rm "${SRA_accession}.Clupea_harengus.sam"


	rm "${SRA_accession}_pass_1.fastp.fastq.gz" "${SRA_accession}_pass_2.fastp.fastq.gz"



else 

	conda activate align_long_env
	module purge ; module load minimap2

	fastplong -i "${SRA_accession}_pass.fastq.gz" -o "${SRA_accession}_pass.fastp.fastq.gz" --thread 16
	rm "${SRA_accession}_pass.fastq.gz"

	minimap2 -t 16 -a GCF_021917145.1_fHypTra1_genomic.fna "${SRA_accession}_pass.fastp.fastq.gz" > $SRA_accession.Hypomesus_transpacificus.sam 
	minimap2 -t 16 -a GCF_900700415.2_Ch_v2.0.2_genomic.fna "${SRA_accession}_pass.fastp.fastq.gz" > $SRA_accession.Clupea_harengus.sam 

	conda deactivate ; conda activate miniprot

	samtools sort $SRA_accession.Clupea_harengus.sam -o $SRA_accession.Clupea_harengus.sorted.bam -@ 16  ; samtools index $SRA_accession.Clupea_harengus.sorted.bam -@ 16 ; rm $SRA_accession.Clupea_harengus.sam
	samtools sort $SRA_accession.Hypomesus_transpacificus.sam -o $SRA_accession.Hypomesus_transpacificus.sorted.bam -@ 16  ; samtools index $SRA_accession.Hypomesus_transpacificus.sorted.bam -@ 16 ; rm $SRA_accession.Hypomesus_transpacificus.sam

	rm "${SRA_accession}_pass.fastp.fastq.gz"


fi

================================================================================
================================================================================
================================================================================
================================================================================


