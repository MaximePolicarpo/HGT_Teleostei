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


cd ~/Horizontal_transfer_project/Hypomesus_and_Clupea_mapping

module purge ; module load SRA-Toolkit


fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip SRR26322996
fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip SRR26322995
fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip SRR17921806
fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip SRR17921807
fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip SRR17799534
fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip SRR17799535
fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip SRR32203436
fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip SRR30009981
fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip SRR30009984
fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip ERR10936409
fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip ERR10936410
fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip ERR10934070
fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip ERR10934071
fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip ERR10934069



grep "$gene" GCF_900700415.2_Ch_v2.0.2_genomic.gff | cut -f1,4,5,9 | grep "rna-XM" | sed 's/;.*//g' | sed 's/ID=rna-//g' | awk '{print $0, $3 - $2}' | tr ' ' '\t' > All_genes_pos_length.clupea.tsv


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
sbatch --qos=1week -c 16 --mem=50G filter_map_reads.sh SRR17921807 short
sbatch --qos=1week -c 16 --mem=50G filter_map_reads.sh SRR17799534 long
sbatch --qos=1week -c 16 --mem=50G filter_map_reads.sh SRR17799535 short
sbatch --qos=1week -c 16 --mem=50G filter_map_reads.sh SRR32203436 long
sbatch --qos=1week -c 16 --mem=50G filter_map_reads.sh SRR30009981 short
sbatch --qos=1week -c 16 --mem=50G filter_map_reads.sh SRR30009984 short
sbatch --qos=1week -c 16 --mem=50G filter_map_reads.sh ERR10936409 short
sbatch --qos=1week -c 16 --mem=50G filter_map_reads.sh ERR10936410 short
sbatch --qos=1week -c 16 --mem=50G filter_map_reads.sh ERR10934070 long
sbatch --qos=1week -c 16 --mem=50G filter_map_reads.sh ERR10934071 long
sbatch --qos=1week -c 16 --mem=50G filter_map_reads.sh ERR10934069 long



#####====================================================================================================================================
#####====================================================================================================================================
#####========================================== Count the number of long reads spanning entire HGT genes ==============================
#####====================================================================================================================================
#####====================================================================================================================================

#samtools view SRR26322996.Clupea_harengus.sorted.bam NC_045165.1:18004543-18016184 -@ 8 | awk '{if($4 <= 18004543 && $4 + length($10) - 1 >= 18016184) print}' | wc -l

IFS=$'\n'
rm -rf Reads_spanning_regions_tables/ ; mkdir Reads_spanning_regions_tables/
for line in `cat Table_regions.tsv` ; do

	#Define variables
	species=`echo "$line" | cut -f2`
	genome_fai=`ls -l /scicore/home/salzburg/polica0000/Horizontal_transfer_project/Genomic_data/$species/ | grep "genomic.fna.fai" | sed 's/.* //g'`
	
	Region=`echo "$line" | cut -f4`

	curr_scaffold=`echo "$Region" | cut -f1 -d ":"`
	scaffold_length=`grep "$curr_scaffold	" /scicore/home/salzburg/polica0000/Horizontal_transfer_project/Genomic_data/$species/$genome_fai | cut -f2`
	start=`echo "$Region" | cut -f2 -d ":" | cut -f1 -d "-"`
	stop=`echo "$Region" | cut -f2 -d ":" | cut -f2 -d "-"`
	


	start_ext=$((start - 1000000))
	if [ "$start_ext" -le 0 ]; then start_ext=1 ; fi
	stop_ext=$((stop + 1000000))
	if [ "$stop_ext" -ge "$scaffold_length" ] ; then stop_ext=$((scaffold_length - 10)) ; fi
	Region_extended=`echo "$curr_scaffold:$start_ext-$stop_ext"`


	for SRA_read_line in `cat list_SRA_reads.txt` ; do


		SRA_accession=`echo "$SRA_read_line" | cut -f2`
		species_SRA=`echo "$SRA_read_line" | cut -f1`

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
	
		rm curr_reads_positions.named.tsv ; rm curr_reads_positions.tsv

		echo "$line,$SRA_accession,$species_SRA,$nb_reads_full_span" | tr "," "\t"


	done

done





#####====================================================================================================================================
#####====================================================================================================================================
#####================================================== N5.HOG0001647  ==================================================================
#####====================================================================================================================================
#####====================================================================================================================================

cd /scicore/home/salzburg/polica0000/Horizontal_transfer_project/Hypomesus_and_Clupea_mapping/Second_Analysis_reads/
awk '$2 == "Hypomesus_transpacificus"' Table_longreads_count.final.tsv > Table_longreads_count.hypomesus.txt
awk '$2 == "Clupea_harengus"' Table_longreads_count.final.tsv > Table_longreads_count.clupea.txt


curr_HOG=N5.HOG0001647

rm -rf $curr_HOG ; mkdir $curr_HOG


curr_XM_clupea=XM_031563501
full_XM_name=`grep "$curr_XM_clupea" Clupea_harengus.prot.nostop | sed 's/>//g'`
samtools faidx Clupea_harengus.prot.nostop $full_XM_name > $curr_HOG/$curr_XM_clupea.prot
Clupea_SRR=SRR26322996
curr_scaffold=`grep "$curr_XM_clupea" Table_longreads_count.clupea.txt | cut -f4 | sort | uniq | cut -f1 -d ":"`
start=`grep "$curr_XM_clupea" Table_longreads_count.clupea.txt | cut -f4 | sort | uniq | cut -f2 -d ":" | cut -f1 -d "-"`
stop=`grep "$curr_XM_clupea" Table_longreads_count.clupea.txt | cut -f4 | sort | uniq | cut -f2 -d ":" | cut -f2 -d "-"`
scaffold_length=`grep "$curr_scaffold	" GCF_900700415.2_Ch_v2.0.2_genomic.fna.fai | cut -f2`
start_ext=$((start - 1000000))
if [ "$start_ext" -le 0 ]; then start_ext=1 ; fi
stop_ext=$((stop + 1000000))
if [ "$stop_ext" -ge "$scaffold_length" ] ; then stop_ext=$((scaffold_length - 10)) ; fi
Region_extended=`echo "$curr_scaffold:$start_ext-$stop_ext"`
Clupea_read_name=`grep "$Clupea_SRR" ../Reads_spanning_regions_tables/Clupea_harengus.$curr_scaffold\:$start-$stop | awk '{print $0, $5 - $4}' | tr ',' '\t' | sort -k6 -n | tail -1 | cut -f2`
echo ">$Clupea_read_name" > $curr_HOG/Clupea_harengus.fa
samtools view ../$Clupea_SRR.Clupea_harengus.sorted.bam $Region_extended -@ 8 | grep "$Clupea_read_name" | cut -f10 | awk '{print $0,length}' | tr ' ' ','  | sort -k2 -n -t ',' | tail -1 | cut -f1 -d ","  >> $curr_HOG/Clupea_harengus.fa
strand_clupea=`exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/$curr_XM_clupea.prot $curr_HOG/Clupea_harengus.fa | grep "vulgar" | cut -f9 -d " "`
if [ $strand_clupea == "-" ] ; then revseq $curr_HOG/Clupea_harengus.fa $curr_HOG/Clupea_harengus.rev ; sed -i 's/ Reversed://g' $curr_HOG/Clupea_harengus.rev ; mv $curr_HOG/Clupea_harengus.rev $curr_HOG/Clupea_harengus.fa ; fi
exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/$curr_XM_clupea.prot $curr_HOG/Clupea_harengus.fa > $curr_HOG/Clupea_harengus.exo

curr_XM_hypo=XM_047032478
full_XM_name=`grep "$curr_XM_hypo" Hypomesus_transpacificus.prot.nostop | sed 's/>//g'`
samtools faidx Hypomesus_transpacificus.prot.nostop $full_XM_name > $curr_HOG/$curr_XM_hypo.prot
Hypomesus_SRR=SRR17921806
curr_scaffold=`grep "$curr_XM_hypo" Table_longreads_count.hypomesus.txt | cut -f4 | sort | uniq | cut -f1 -d ":"`
start=`grep "$curr_XM_hypo" Table_longreads_count.hypomesus.txt | cut -f4 | sort | uniq | cut -f2 -d ":" | cut -f1 -d "-"`
stop=`grep "$curr_XM_hypo" Table_longreads_count.hypomesus.txt | cut -f4 | sort | uniq | cut -f2 -d ":" | cut -f2 -d "-"`
scaffold_length=`grep "$curr_scaffold	" GCF_021917145.1_fHypTra1_genomic.fna.fai | cut -f2`
start_ext=$((start - 1000000))
if [ "$start_ext" -le 0 ]; then start_ext=1 ; fi
stop_ext=$((stop + 1000000))
if [ "$stop_ext" -ge "$scaffold_length" ] ; then stop_ext=$((scaffold_length - 10)) ; fi
Region_extended=`echo "$curr_scaffold:$start_ext-$stop_ext"`
Hypomesus_read_name=`grep "$Hypomesus_SRR" ../Reads_spanning_regions_tables/Hypomesus_transpacificus.$curr_scaffold\:$start-$stop | awk '{print $0, $5 - $4}' | tr ',' '\t' | sort -k6 -n | tail -1 | cut -f2`
echo ">$Hypomesus_read_name" > $curr_HOG/Hypomesus_transpacificus.fa
samtools view ../$Hypomesus_SRR.Hypomesus_transpacificus.sorted.bam $Region_extended -@ 8 | grep "$Hypomesus_read_name" | cut -f10 | awk '{print $0,length}' | tr ' ' ','  | sort -k2 -n -t ',' | tail -1 | cut -f1 -d ","  >> $curr_HOG/Hypomesus_transpacificus.fa
strand_hypo=`exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/$curr_XM_hypo.prot $curr_HOG/Hypomesus_transpacificus.fa | grep "vulgar" | cut -f9 -d " "`
if [ $strand_hypo == "-" ] ; then revseq $curr_HOG/Hypomesus_transpacificus.fa $curr_HOG/Hypomesus_transpacificus.rev ; sed -i 's/ Reversed://g' $curr_HOG/Hypomesus_transpacificus.rev ; mv $curr_HOG/Hypomesus_transpacificus.rev $curr_HOG/Hypomesus_transpacificus.fa ; fi
exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/$curr_XM_hypo.prot $curr_HOG/Hypomesus_transpacificus.fa > $curr_HOG/Hypomesus_transpacificus.exo


cd $curr_HOG ; for i in *.fa ; do samtools faidx $i ; done
for i in *.fai ; do echo "$i" ; cat $i ; done

grep "	exon	" Clupea_harengus.exo | cut -f4,5 | sed 's/^/Clupea     /g'
grep "	exon	" Hypomesus_transpacificus.exo | cut -f4,5 | sed 's/^/Hypomesus     /g'




#####====================================================================================================================================
#####====================================================================================================================================
#####================================================== N5.HOG0008885  ==================================================================
#####====================================================================================================================================
#####====================================================================================================================================


curr_HOG=N5.HOG0008885

rm -rf $curr_HOG ; mkdir $curr_HOG


curr_XM_clupea=XM_031562957
full_XM_name=`grep "$curr_XM_clupea" Clupea_harengus.prot.nostop | sed 's/>//g'`
samtools faidx Clupea_harengus.prot.nostop $full_XM_name > $curr_HOG/$curr_XM_clupea.prot
Clupea_SRR=SRR26322995
curr_scaffold=`grep "$curr_XM_clupea" Table_longreads_count.clupea.txt | cut -f4 | sort | uniq | cut -f1 -d ":"`
start=`grep "$curr_XM_clupea" Table_longreads_count.clupea.txt | cut -f4 | sort | uniq | cut -f2 -d ":" | cut -f1 -d "-"`
stop=`grep "$curr_XM_clupea" Table_longreads_count.clupea.txt | cut -f4 | sort | uniq | cut -f2 -d ":" | cut -f2 -d "-"`
scaffold_length=`grep "$curr_scaffold	" GCF_900700415.2_Ch_v2.0.2_genomic.fna.fai | cut -f2`
start_ext=$((start - 1000000))
if [ "$start_ext" -le 0 ]; then start_ext=1 ; fi
stop_ext=$((stop + 1000000))
if [ "$stop_ext" -ge "$scaffold_length" ] ; then stop_ext=$((scaffold_length - 10)) ; fi
Region_extended=`echo "$curr_scaffold:$start_ext-$stop_ext"`
Clupea_read_name=`grep "$Clupea_SRR" ../Reads_spanning_regions_tables/Clupea_harengus.$curr_scaffold\:$start-$stop | awk '{print $0, $5 - $4}' | tr ',' '\t' | sort -k6 -n | tail -1 | cut -f2`
echo ">$Clupea_read_name" > $curr_HOG/Clupea_harengus.fa
samtools view ../$Clupea_SRR.Clupea_harengus.sorted.bam $Region_extended -@ 8 | grep "$Clupea_read_name" | cut -f10 | awk '{print $0,length}' | tr ' ' ','  | sort -k2 -n -t ',' | tail -1 | cut -f1 -d ","  >> $curr_HOG/Clupea_harengus.fa
strand_clupea=`exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/$curr_XM_clupea.prot $curr_HOG/Clupea_harengus.fa | grep "vulgar" | cut -f9 -d " "`
if [ $strand_clupea == "-" ] ; then revseq $curr_HOG/Clupea_harengus.fa $curr_HOG/Clupea_harengus.rev ; sed -i 's/ Reversed://g' $curr_HOG/Clupea_harengus.rev ; mv $curr_HOG/Clupea_harengus.rev $curr_HOG/Clupea_harengus.fa ; fi
exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/$curr_XM_clupea.prot $curr_HOG/Clupea_harengus.fa > $curr_HOG/Clupea_harengus.exo

curr_XM_hypo=XM_047048755
full_XM_name=`grep "$curr_XM_hypo" Hypomesus_transpacificus.prot.nostop | sed 's/>//g'`
samtools faidx Hypomesus_transpacificus.prot.nostop $full_XM_name > $curr_HOG/$curr_XM_hypo.prot
Hypomesus_SRR=SRR17921806
curr_scaffold=`grep "$curr_XM_hypo" Table_longreads_count.hypomesus.txt | cut -f4 | sort | uniq | cut -f1 -d ":"`
start=`grep "$curr_XM_hypo" Table_longreads_count.hypomesus.txt | cut -f4 | sort | uniq | cut -f2 -d ":" | cut -f1 -d "-"`
stop=`grep "$curr_XM_hypo" Table_longreads_count.hypomesus.txt | cut -f4 | sort | uniq | cut -f2 -d ":" | cut -f2 -d "-"`
scaffold_length=`grep "$curr_scaffold	" GCF_021917145.1_fHypTra1_genomic.fna.fai | cut -f2`
start_ext=$((start - 1000000))
if [ "$start_ext" -le 0 ]; then start_ext=1 ; fi
stop_ext=$((stop + 1000000))
if [ "$stop_ext" -ge "$scaffold_length" ] ; then stop_ext=$((scaffold_length - 10)) ; fi
Region_extended=`echo "$curr_scaffold:$start_ext-$stop_ext"`
Hypomesus_read_name=`grep "$Hypomesus_SRR" ../Reads_spanning_regions_tables/Hypomesus_transpacificus.$curr_scaffold\:$start-$stop | awk '{print $0, $5 - $4}' | tr ',' '\t' | sort -k6 -n | tail -1 | cut -f2`
echo ">$Hypomesus_read_name" > $curr_HOG/Hypomesus_transpacificus.fa
samtools view ../$Hypomesus_SRR.Hypomesus_transpacificus.sorted.bam $Region_extended -@ 8 | grep "$Hypomesus_read_name" | cut -f10 | awk '{print $0,length}' | tr ' ' ','  | sort -k2 -n -t ',' | tail -1 | cut -f1 -d ","  >> $curr_HOG/Hypomesus_transpacificus.fa
strand_hypo=`exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/$curr_XM_hypo.prot $curr_HOG/Hypomesus_transpacificus.fa | grep "vulgar" | cut -f9 -d " "`
if [ $strand_hypo == "-" ] ; then revseq $curr_HOG/Hypomesus_transpacificus.fa $curr_HOG/Hypomesus_transpacificus.rev ; sed -i 's/ Reversed://g' $curr_HOG/Hypomesus_transpacificus.rev ; mv $curr_HOG/Hypomesus_transpacificus.rev $curr_HOG/Hypomesus_transpacificus.fa ; fi
exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/$curr_XM_hypo.prot $curr_HOG/Hypomesus_transpacificus.fa > $curr_HOG/Hypomesus_transpacificus.exo


cd $curr_HOG ; for i in *.fa ; do samtools faidx $i ; done
for i in *.fai ; do echo "$i" ; cat $i ; done

grep "	exon	" Clupea_harengus.exo | cut -f4,5 | sed 's/^/Clupea	/g'
grep "	exon	" Hypomesus_transpacificus.exo | cut -f4,5 | sed 's/^/Hypomesus	/g'






#####====================================================================================================================================
#####====================================================================================================================================
#####================================================== N5.HOG0013514  ==================================================================
#####====================================================================================================================================
#####====================================================================================================================================


curr_HOG=N5.HOG0013514

rm -rf $curr_HOG ; mkdir $curr_HOG


curr_XM_clupea=XM_031565680
full_XM_name=`grep "$curr_XM_clupea" Clupea_harengus.prot.nostop | sed 's/>//g'`
samtools faidx Clupea_harengus.prot.nostop $full_XM_name > $curr_HOG/$curr_XM_clupea.prot
Clupea_SRR=SRR26322996
curr_scaffold=`grep "$curr_XM_clupea" Table_longreads_count.clupea.txt | cut -f4 | sort | uniq | cut -f1 -d ":"`
start=`grep "$curr_XM_clupea" Table_longreads_count.clupea.txt | cut -f4 | sort | uniq | cut -f2 -d ":" | cut -f1 -d "-"`
stop=`grep "$curr_XM_clupea" Table_longreads_count.clupea.txt | cut -f4 | sort | uniq | cut -f2 -d ":" | cut -f2 -d "-"`
scaffold_length=`grep "$curr_scaffold	" GCF_900700415.2_Ch_v2.0.2_genomic.fna.fai | cut -f2`
start_ext=$((start - 1000000))
if [ "$start_ext" -le 0 ]; then start_ext=1 ; fi
stop_ext=$((stop + 1000000))
if [ "$stop_ext" -ge "$scaffold_length" ] ; then stop_ext=$((scaffold_length - 10)) ; fi
Region_extended=`echo "$curr_scaffold:$start_ext-$stop_ext"`
Clupea_read_name=`grep "$Clupea_SRR" ../Reads_spanning_regions_tables/Clupea_harengus.$curr_scaffold\:$start-$stop | awk '{print $0, $5 - $4}' | tr ',' '\t' | sort -k6 -n | tail -1 | cut -f2`
echo ">$Clupea_read_name" > $curr_HOG/Clupea_harengus.fa
samtools view ../$Clupea_SRR.Clupea_harengus.sorted.bam $Region_extended -@ 8 | grep "$Clupea_read_name" | cut -f10 | awk '{print $0,length}' | tr ' ' ','  | sort -k2 -n -t ',' | tail -1 | cut -f1 -d ","  >> $curr_HOG/Clupea_harengus.fa
strand_clupea=`exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/$curr_XM_clupea.prot $curr_HOG/Clupea_harengus.fa | grep "vulgar" | cut -f9 -d " "`
if [ $strand_clupea == "-" ] ; then revseq $curr_HOG/Clupea_harengus.fa $curr_HOG/Clupea_harengus.rev ; sed -i 's/ Reversed://g' $curr_HOG/Clupea_harengus.rev ; mv $curr_HOG/Clupea_harengus.rev $curr_HOG/Clupea_harengus.fa ; fi
exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/$curr_XM_hypo.prot $curr_HOG/Clupea_harengus.fa > $curr_HOG/Clupea_harengus.exo

curr_XM_hypo=XM_047045767
full_XM_name=`grep "$curr_XM_hypo" Hypomesus_transpacificus.prot.nostop | sed 's/>//g'`
samtools faidx Hypomesus_transpacificus.prot.nostop $full_XM_name > $curr_HOG/$curr_XM_hypo.prot
Hypomesus_SRR=SRR17921806
curr_scaffold=`grep "$curr_XM_hypo" Table_longreads_count.hypomesus.txt | cut -f4 | sort | uniq | cut -f1 -d ":"`
start=`grep "$curr_XM_hypo" Table_longreads_count.hypomesus.txt | cut -f4 | sort | uniq | cut -f2 -d ":" | cut -f1 -d "-"`
stop=`grep "$curr_XM_hypo" Table_longreads_count.hypomesus.txt | cut -f4 | sort | uniq | cut -f2 -d ":" | cut -f2 -d "-"`
scaffold_length=`grep "$curr_scaffold	" GCF_021917145.1_fHypTra1_genomic.fna.fai | cut -f2`
start_ext=$((start - 1000000))
if [ "$start_ext" -le 0 ]; then start_ext=1 ; fi
stop_ext=$((stop + 1000000))
if [ "$stop_ext" -ge "$scaffold_length" ] ; then stop_ext=$((scaffold_length - 10)) ; fi
Region_extended=`echo "$curr_scaffold:$start_ext-$stop_ext"`
Hypomesus_read_name=`grep "$Hypomesus_SRR" ../Reads_spanning_regions_tables/Hypomesus_transpacificus.$curr_scaffold\:$start-$stop | awk '{print $0, $5 - $4}' | tr ',' '\t' | sort -k6 -n | tail -1 | cut -f2`
echo ">$Hypomesus_read_name" > $curr_HOG/Hypomesus_transpacificus.fa
samtools view ../$Hypomesus_SRR.Hypomesus_transpacificus.sorted.bam $Region_extended -@ 8 | grep "$Hypomesus_read_name" | cut -f10 | awk '{print $0,length}' | tr ' ' ','  | sort -k2 -n -t ',' | tail -1 | cut -f1 -d ","  >> $curr_HOG/Hypomesus_transpacificus.fa
strand_hypo=`exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/$curr_XM_hypo.prot $curr_HOG/Hypomesus_transpacificus.fa | grep "vulgar" | cut -f9 -d " "`
if [ $strand_hypo == "-" ] ; then revseq $curr_HOG/Hypomesus_transpacificus.fa $curr_HOG/Hypomesus_transpacificus.rev ; sed -i 's/ Reversed://g' $curr_HOG/Hypomesus_transpacificus.rev ; mv $curr_HOG/Hypomesus_transpacificus.rev $curr_HOG/Hypomesus_transpacificus.fa ; fi
exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/$curr_XM_hypo.prot $curr_HOG/Hypomesus_transpacificus.fa > $curr_HOG/Hypomesus_transpacificus.exo


cd $curr_HOG ; for i in *.fa ; do samtools faidx $i ; done
for i in *.fai ; do echo "$i" ; cat $i ; done

grep "	exon	" Clupea_harengus.exo | cut -f4,5 | sed 's/^/Clupea	/g'
grep "	exon	" Hypomesus_transpacificus.exo | cut -f4,5 | sed 's/^/Hypomesus	/g'



#####====================================================================================================================================
#####====================================================================================================================================
#####================================================== N5.HOG0018901  ==================================================================
#####====================================================================================================================================
#####====================================================================================================================================


curr_HOG=N5.HOG0018901

rm -rf $curr_HOG ; mkdir $curr_HOG


curr_XM_clupea=XM_012826723
full_XM_name=`grep "$curr_XM_clupea" Clupea_harengus.prot.nostop | sed 's/>//g'`
samtools faidx Clupea_harengus.prot.nostop $full_XM_name > $curr_HOG/$curr_XM_clupea.prot
Clupea_SRR=SRR26322996
curr_scaffold=`grep "$curr_XM_clupea" Table_longreads_count.clupea.txt | cut -f4 | sort | uniq | cut -f1 -d ":"`
start=`grep "$curr_XM_clupea" Table_longreads_count.clupea.txt | cut -f4 | sort | uniq | cut -f2 -d ":" | cut -f1 -d "-"`
stop=`grep "$curr_XM_clupea" Table_longreads_count.clupea.txt | cut -f4 | sort | uniq | cut -f2 -d ":" | cut -f2 -d "-"`
scaffold_length=`grep "$curr_scaffold	" GCF_900700415.2_Ch_v2.0.2_genomic.fna.fai | cut -f2`
start_ext=$((start - 1000000))
if [ "$start_ext" -le 0 ]; then start_ext=1 ; fi
stop_ext=$((stop + 1000000))
if [ "$stop_ext" -ge "$scaffold_length" ] ; then stop_ext=$((scaffold_length - 10)) ; fi
Region_extended=`echo "$curr_scaffold:$start_ext-$stop_ext"`
Clupea_read_name=`grep "$Clupea_SRR" ../Reads_spanning_regions_tables/Clupea_harengus.$curr_scaffold\:$start-$stop | awk '{print $0, $5 - $4}' | tr ',' '\t' | sort -k6 -n | tail -1 | cut -f2`
echo ">$Clupea_read_name" > $curr_HOG/Clupea_harengus.fa
samtools view ../$Clupea_SRR.Clupea_harengus.sorted.bam $Region_extended -@ 8 | grep "$Clupea_read_name" | cut -f10 | awk '{print $0,length}' | tr ' ' ','  | sort -k2 -n -t ',' | tail -1 | cut -f1 -d ","  >> $curr_HOG/Clupea_harengus.fa
strand_clupea=`exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/$curr_XM_clupea.prot $curr_HOG/Clupea_harengus.fa | grep "vulgar" | cut -f9 -d " "`
if [ $strand_clupea == "-" ] ; then revseq $curr_HOG/Clupea_harengus.fa $curr_HOG/Clupea_harengus.rev ; sed -i 's/ Reversed://g' $curr_HOG/Clupea_harengus.rev ; mv $curr_HOG/Clupea_harengus.rev $curr_HOG/Clupea_harengus.fa ; fi
exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/$curr_XM_clupea.prot $curr_HOG/Clupea_harengus.fa > $curr_HOG/Clupea_harengus.exo

curr_XM_hypo=XM_047044257
full_XM_name=`grep "$curr_XM_hypo" Hypomesus_transpacificus.prot.nostop | sed 's/>//g'`
samtools faidx Hypomesus_transpacificus.prot.nostop $full_XM_name > $curr_HOG/$curr_XM_hypo.prot
Hypomesus_SRR=SRR17921806
curr_scaffold=`grep "$curr_XM_hypo" Table_longreads_count.hypomesus.txt | cut -f4 | sort | uniq | cut -f1 -d ":"`
start=`grep "$curr_XM_hypo" Table_longreads_count.hypomesus.txt | cut -f4 | sort | uniq | cut -f2 -d ":" | cut -f1 -d "-"`
stop=`grep "$curr_XM_hypo" Table_longreads_count.hypomesus.txt | cut -f4 | sort | uniq | cut -f2 -d ":" | cut -f2 -d "-"`
scaffold_length=`grep "$curr_scaffold	" GCF_021917145.1_fHypTra1_genomic.fna.fai | cut -f2`
start_ext=$((start - 1000000))
if [ "$start_ext" -le 0 ]; then start_ext=1 ; fi
stop_ext=$((stop + 1000000))
if [ "$stop_ext" -ge "$scaffold_length" ] ; then stop_ext=$((scaffold_length - 10)) ; fi
Region_extended=`echo "$curr_scaffold:$start_ext-$stop_ext"`
Hypomesus_read_name=`grep "$Hypomesus_SRR" ../Reads_spanning_regions_tables/Hypomesus_transpacificus.$curr_scaffold\:$start-$stop | awk '{print $0, $5 - $4}' | tr ',' '\t' | sort -k6 -n | tail -1 | cut -f2`
echo ">$Hypomesus_read_name" > $curr_HOG/Hypomesus_transpacificus.fa
samtools view ../$Hypomesus_SRR.Hypomesus_transpacificus.sorted.bam $Region_extended -@ 8 | grep "$Hypomesus_read_name" | cut -f10 | awk '{print $0,length}' | tr ' ' ','  | sort -k2 -n -t ',' | tail -1 | cut -f1 -d ","  >> $curr_HOG/Hypomesus_transpacificus.fa
strand_hypo=`exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/$curr_XM_hypo.prot $curr_HOG/Hypomesus_transpacificus.fa | grep "vulgar" | cut -f9 -d " "`
if [ $strand_hypo == "-" ] ; then revseq $curr_HOG/Hypomesus_transpacificus.fa $curr_HOG/Hypomesus_transpacificus.rev ; sed -i 's/ Reversed://g' $curr_HOG/Hypomesus_transpacificus.rev ; mv $curr_HOG/Hypomesus_transpacificus.rev $curr_HOG/Hypomesus_transpacificus.fa ; fi
exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/$curr_XM_hypo.prot $curr_HOG/Hypomesus_transpacificus.fa > $curr_HOG/Hypomesus_transpacificus.exo


cd $curr_HOG ; for i in *.fa ; do samtools faidx $i ; done
for i in *.fai ; do echo "$i" ; cat $i ; done

grep "	exon	" Clupea_harengus.exo | cut -f4,5 | sed 's/^/Clupea	/g'
grep "	exon	" Hypomesus_transpacificus.exo | cut -f4,5 | sed 's/^/Hypomesus	/g'



#####====================================================================================================================================
#####====================================================================================================================================
#####================================================== N5.HOG0019111  ==================================================================
#####====================================================================================================================================
#####====================================================================================================================================


curr_HOG=N5.HOG0019111

rm -rf $curr_HOG ; mkdir $curr_HOG


curr_XM_clupea=XM_031577666
full_XM_name=`grep "$curr_XM_clupea" Clupea_harengus.prot.nostop | sed 's/>//g'`
samtools faidx Clupea_harengus.prot.nostop $full_XM_name > $curr_HOG/$curr_XM_clupea.prot
Clupea_SRR=SRR26322995
curr_scaffold=`grep "$curr_XM_clupea" Table_longreads_count.clupea.txt | cut -f4 | sort | uniq | cut -f1 -d ":"`
start=`grep "$curr_XM_clupea" Table_longreads_count.clupea.txt | cut -f4 | sort | uniq | cut -f2 -d ":" | cut -f1 -d "-"`
stop=`grep "$curr_XM_clupea" Table_longreads_count.clupea.txt | cut -f4 | sort | uniq | cut -f2 -d ":" | cut -f2 -d "-"`
scaffold_length=`grep "$curr_scaffold	" GCF_900700415.2_Ch_v2.0.2_genomic.fna.fai | cut -f2`
start_ext=$((start - 1000000))
if [ "$start_ext" -le 0 ]; then start_ext=1 ; fi
stop_ext=$((stop + 1000000))
if [ "$stop_ext" -ge "$scaffold_length" ] ; then stop_ext=$((scaffold_length - 10)) ; fi
Region_extended=`echo "$curr_scaffold:$start_ext-$stop_ext"`
Clupea_read_name=`grep "$Clupea_SRR" ../Reads_spanning_regions_tables/Clupea_harengus.$curr_scaffold\:$start-$stop | awk '{print $0, $5 - $4}' | tr ',' '\t' | sort -k6 -n | tail -1 | cut -f2`
echo ">$Clupea_read_name" > $curr_HOG/Clupea_harengus.fa
samtools view ../$Clupea_SRR.Clupea_harengus.sorted.bam $Region_extended -@ 8 | grep "$Clupea_read_name" | cut -f10 | awk '{print $0,length}' | tr ' ' ','  | sort -k2 -n -t ',' | tail -1 | cut -f1 -d ","  >> $curr_HOG/Clupea_harengus.fa
strand_clupea=`exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/$curr_XM_clupea.prot $curr_HOG/Clupea_harengus.fa | grep "vulgar" | cut -f9 -d " "`
if [ $strand_clupea == "-" ] ; then revseq $curr_HOG/Clupea_harengus.fa $curr_HOG/Clupea_harengus.rev ; sed -i 's/ Reversed://g' $curr_HOG/Clupea_harengus.rev ; mv $curr_HOG/Clupea_harengus.rev $curr_HOG/Clupea_harengus.fa ; fi
exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/$curr_XM_clupea.prot $curr_HOG/Clupea_harengus.fa > $curr_HOG/Clupea_harengus.exo

curr_XM_hypo=XM_047049897
full_XM_name=`grep "$curr_XM_hypo" Hypomesus_transpacificus.prot.nostop | sed 's/>//g'`
samtools faidx Hypomesus_transpacificus.prot.nostop $full_XM_name > $curr_HOG/$curr_XM_hypo.prot
Hypomesus_SRR=SRR17799534
curr_scaffold=`grep "$curr_XM_hypo" Table_longreads_count.hypomesus.txt | cut -f4 | sort | uniq | cut -f1 -d ":"`
start=`grep "$curr_XM_hypo" Table_longreads_count.hypomesus.txt | cut -f4 | sort | uniq | cut -f2 -d ":" | cut -f1 -d "-"`
stop=`grep "$curr_XM_hypo" Table_longreads_count.hypomesus.txt | cut -f4 | sort | uniq | cut -f2 -d ":" | cut -f2 -d "-"`
scaffold_length=`grep "$curr_scaffold	" GCF_021917145.1_fHypTra1_genomic.fna.fai | cut -f2`
start_ext=$((start - 1000000))
if [ "$start_ext" -le 0 ]; then start_ext=1 ; fi
stop_ext=$((stop + 1000000))
if [ "$stop_ext" -ge "$scaffold_length" ] ; then stop_ext=$((scaffold_length - 10)) ; fi
Region_extended=`echo "$curr_scaffold:$start_ext-$stop_ext"`
Hypomesus_read_name=`grep "$Hypomesus_SRR" ../Reads_spanning_regions_tables/Hypomesus_transpacificus.$curr_scaffold\:$start-$stop | awk '{print $0, $5 - $4}' | tr ',' '\t' | sort -k6 -n | tail -1 | cut -f2`
echo ">$Hypomesus_read_name" > $curr_HOG/Hypomesus_transpacificus.fa
samtools view ../$Hypomesus_SRR.Hypomesus_transpacificus.sorted.bam $Region_extended -@ 8 | grep "$Hypomesus_read_name" | cut -f10 | awk '{print $0,length}' | tr ' ' ','  | sort -k2 -n -t ',' | tail -1 | cut -f1 -d ","  >> $curr_HOG/Hypomesus_transpacificus.fa
strand_hypo=`exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/$curr_XM_hypo.prot $curr_HOG/Hypomesus_transpacificus.fa | grep "vulgar" | cut -f9 -d " "`
if [ $strand_hypo == "-" ] ; then revseq $curr_HOG/Hypomesus_transpacificus.fa $curr_HOG/Hypomesus_transpacificus.rev ; sed -i 's/ Reversed://g' $curr_HOG/Hypomesus_transpacificus.rev ; mv $curr_HOG/Hypomesus_transpacificus.rev $curr_HOG/Hypomesus_transpacificus.fa ; fi
exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/$curr_XM_hypo.prot $curr_HOG/Hypomesus_transpacificus.fa > $curr_HOG/Hypomesus_transpacificus.exo


cd $curr_HOG ; for i in *.fa ; do samtools faidx $i ; done
for i in *.fai ; do echo "$i" ; cat $i ; done

grep "	exon	" Clupea_harengus.exo | cut -f4,5 | sed 's/^/Clupea	/g'
grep "	exon	" Hypomesus_transpacificus.exo | cut -f4,5 | sed 's/^/Hypomesus	/g'






#####====================================================================================================================================
#####====================================================================================================================================
#####================================================== N5.HOG0028899  ==================================================================
#####====================================================================================================================================
#####====================================================================================================================================

cd /scicore/home/salzburg/polica0000/Horizontal_transfer_project/Hypomesus_and_Clupea_mapping/Second_Analysis_reads/


awk '$2 == "Hypomesus_transpacificus"' Table_longreads_count.final.tsv > Table_longreads_count.hypomesus.txt
awk '$2 == "Clupea_harengus"' Table_longreads_count.final.tsv > Table_longreads_count.clupea.txt

curr_HOG=N5.HOG0028899

rm -rf $curr_HOG ; mkdir $curr_HOG


curr_XM_clupea=XM_031581960
full_XM_name=`grep "$curr_XM_clupea" Clupea_harengus.prot.nostop | sed 's/>//g'`
samtools faidx Clupea_harengus.prot.nostop $full_XM_name > $curr_HOG/$curr_XM_clupea.prot
Clupea_SRR=SRR32203436
curr_scaffold=`grep "$curr_XM_clupea" Table_longreads_count.clupea.txt | cut -f4 | sort | uniq | cut -f1 -d ":"`
start=`grep "$curr_XM_clupea" Table_longreads_count.clupea.txt | cut -f4 | sort | uniq | cut -f2 -d ":" | cut -f1 -d "-"`
stop=`grep "$curr_XM_clupea" Table_longreads_count.clupea.txt | cut -f4 | sort | uniq | cut -f2 -d ":" | cut -f2 -d "-"`
scaffold_length=`grep "$curr_scaffold	" GCF_900700415.2_Ch_v2.0.2_genomic.fna.fai | cut -f2`
start_ext=$((start - 1000000))
if [ "$start_ext" -le 0 ]; then start_ext=1 ; fi
stop_ext=$((stop + 1000000))
if [ "$stop_ext" -ge "$scaffold_length" ] ; then stop_ext=$((scaffold_length - 10)) ; fi
Region_extended=`echo "$curr_scaffold:$start_ext-$stop_ext"`
Clupea_read_name=`grep "$Clupea_SRR" ../Reads_spanning_regions_tables/Clupea_harengus.$curr_scaffold\:$start-$stop | awk '{print $0, $5 - $4}' | tr ',' '\t' | sort -k6 -n | tail -1 | cut -f2`
echo ">$Clupea_read_name" > $curr_HOG/Clupea_harengus.fa
samtools view ../$Clupea_SRR.Clupea_harengus.sorted.bam $Region_extended -@ 8 | grep "$Clupea_read_name" | cut -f10 | awk '{print $0,length}' | tr ' ' ','  | sort -k2 -n -t ',' | tail -1 | cut -f1 -d ","  >> $curr_HOG/Clupea_harengus.fa
strand_clupea=`exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/$curr_XM_clupea.prot $curr_HOG/Clupea_harengus.fa | grep "vulgar" | cut -f9 -d " "`
if [ $strand_clupea == "-" ] ; then revseq $curr_HOG/Clupea_harengus.fa $curr_HOG/Clupea_harengus.rev ; sed -i 's/ Reversed://g' $curr_HOG/Clupea_harengus.rev ; mv $curr_HOG/Clupea_harengus.rev $curr_HOG/Clupea_harengus.fa ; fi
exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/$curr_XM_clupea.prot $curr_HOG/Clupea_harengus.fa > $curr_HOG/Clupea_harengus.exo

curr_XM_hypo=XM_047031800
full_XM_name=`grep "$curr_XM_hypo" Hypomesus_transpacificus.prot.nostop | sed 's/>//g'`
samtools faidx Hypomesus_transpacificus.prot.nostop $full_XM_name > $curr_HOG/$curr_XM_hypo.prot
Hypomesus_SRR=SRR17799534
curr_scaffold=`grep "$curr_XM_hypo" Table_longreads_count.hypomesus.txt | cut -f4 | sort | uniq | cut -f1 -d ":"`
start=`grep "$curr_XM_hypo" Table_longreads_count.hypomesus.txt | cut -f4 | sort | uniq | cut -f2 -d ":" | cut -f1 -d "-"`
stop=`grep "$curr_XM_hypo" Table_longreads_count.hypomesus.txt | cut -f4 | sort | uniq | cut -f2 -d ":" | cut -f2 -d "-"`
scaffold_length=`grep "$curr_scaffold	" GCF_021917145.1_fHypTra1_genomic.fna.fai | cut -f2`
start_ext=$((start - 1000000))
if [ "$start_ext" -le 0 ]; then start_ext=1 ; fi
stop_ext=$((stop + 1000000))
if [ "$stop_ext" -ge "$scaffold_length" ] ; then stop_ext=$((scaffold_length - 10)) ; fi
Region_extended=`echo "$curr_scaffold:$start_ext-$stop_ext"`
Hypomesus_read_name=`grep "$Hypomesus_SRR" ../Reads_spanning_regions_tables/Hypomesus_transpacificus.$curr_scaffold\:$start-$stop | awk '{print $0, $5 - $4}' | tr ',' '\t' | sort -k6 -n | tail -1 | cut -f2`
echo ">$Hypomesus_read_name" > $curr_HOG/Hypomesus_transpacificus.fa
samtools view ../$Hypomesus_SRR.Hypomesus_transpacificus.sorted.bam $Region_extended -@ 8 | grep "$Hypomesus_read_name" | cut -f10 | awk '{print $0,length}' | tr ' ' ','  | sort -k2 -n -t ',' | tail -1 | cut -f1 -d ","  >> $curr_HOG/Hypomesus_transpacificus.fa
strand_hypo=`exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/$curr_XM_hypo.prot $curr_HOG/Hypomesus_transpacificus.fa | grep "vulgar" | cut -f9 -d " "`
if [ $strand_hypo == "-" ] ; then revseq $curr_HOG/Hypomesus_transpacificus.fa $curr_HOG/Hypomesus_transpacificus.rev ; sed -i 's/ Reversed://g' $curr_HOG/Hypomesus_transpacificus.rev ; mv $curr_HOG/Hypomesus_transpacificus.rev $curr_HOG/Hypomesus_transpacificus.fa ; fi
exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/$curr_XM_hypo.prot $curr_HOG/Hypomesus_transpacificus.fa > $curr_HOG/Hypomesus_transpacificus.exo


cd $curr_HOG ; for i in *.fa ; do samtools faidx $i ; done
for i in *.fai ; do echo "$i" ; cat $i ; done

grep "	exon	" Clupea_harengus.exo | cut -f4,5 | sed 's/^/Clupea	/g'
grep "	exon	" Hypomesus_transpacificus.exo | cut -f4,5 | sed 's/^/Hypomesus	/g'





#####====================================================================================================================================
#####====================================================================================================================================
#####================================================== N5.HOG0047987  ==================================================================
#####====================================================================================================================================
#####====================================================================================================================================


curr_HOG=N5.HOG0047987

rm -rf $curr_HOG ; mkdir $curr_HOG


curr_XM_clupea=XM_031580274
full_XM_name=`grep "$curr_XM_clupea" Clupea_harengus.prot.nostop | sed 's/>//g'`
samtools faidx Clupea_harengus.prot.nostop $full_XM_name > $curr_HOG/$curr_XM_clupea.prot
Clupea_SRR=SRR26322995
curr_scaffold=`grep "$curr_XM_clupea" Table_longreads_count.clupea.txt | cut -f4 | sort | uniq | cut -f1 -d ":"`
start=`grep "$curr_XM_clupea" Table_longreads_count.clupea.txt | cut -f4 | sort | uniq | cut -f2 -d ":" | cut -f1 -d "-"`
stop=`grep "$curr_XM_clupea" Table_longreads_count.clupea.txt | cut -f4 | sort | uniq | cut -f2 -d ":" | cut -f2 -d "-"`
scaffold_length=`grep "$curr_scaffold	" GCF_900700415.2_Ch_v2.0.2_genomic.fna.fai | cut -f2`
start_ext=$((start - 1000000))
if [ "$start_ext" -le 0 ]; then start_ext=1 ; fi
stop_ext=$((stop + 1000000))
if [ "$stop_ext" -ge "$scaffold_length" ] ; then stop_ext=$((scaffold_length - 10)) ; fi
Region_extended=`echo "$curr_scaffold:$start_ext-$stop_ext"`
Clupea_read_name=`grep "$Clupea_SRR" ../Reads_spanning_regions_tables/Clupea_harengus.$curr_scaffold\:$start-$stop | awk '{print $0, $5 - $4}' | tr ',' '\t' | sort -k6 -n | tail -1 | cut -f2`
echo ">$Clupea_read_name" > $curr_HOG/Clupea_harengus.fa
samtools view ../$Clupea_SRR.Clupea_harengus.sorted.bam $Region_extended -@ 8 | grep "$Clupea_read_name" | cut -f10 | awk '{print $0,length}' | tr ' ' ','  | sort -k2 -n -t ',' | tail -1 | cut -f1 -d ","  >> $curr_HOG/Clupea_harengus.fa
strand_clupea=`exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/$curr_XM_clupea.prot $curr_HOG/Clupea_harengus.fa | grep "vulgar" | cut -f9 -d " "`
if [ $strand_clupea == "-" ] ; then revseq $curr_HOG/Clupea_harengus.fa $curr_HOG/Clupea_harengus.rev ; sed -i 's/ Reversed://g' $curr_HOG/Clupea_harengus.rev ; mv $curr_HOG/Clupea_harengus.rev $curr_HOG/Clupea_harengus.fa ; fi
exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/$curr_XM_clupea.prot $curr_HOG/Clupea_harengus.fa > $curr_HOG/Clupea_harengus.exo

curr_XM_hypo=XM_047020147
full_XM_name=`grep "$curr_XM_hypo" Hypomesus_transpacificus.prot.nostop | sed 's/>//g'`
samtools faidx Hypomesus_transpacificus.prot.nostop $full_XM_name > $curr_HOG/$curr_XM_hypo.prot
Hypomesus_SRR=SRR17921806
curr_scaffold=`grep "$curr_XM_hypo" Table_longreads_count.hypomesus.txt | cut -f4 | sort | uniq | cut -f1 -d ":"`
start=`grep "$curr_XM_hypo" Table_longreads_count.hypomesus.txt | cut -f4 | sort | uniq | cut -f2 -d ":" | cut -f1 -d "-"`
stop=`grep "$curr_XM_hypo" Table_longreads_count.hypomesus.txt | cut -f4 | sort | uniq | cut -f2 -d ":" | cut -f2 -d "-"`
scaffold_length=`grep "$curr_scaffold	" GCF_021917145.1_fHypTra1_genomic.fna.fai | cut -f2`
start_ext=$((start - 1000000))
if [ "$start_ext" -le 0 ]; then start_ext=1 ; fi
stop_ext=$((stop + 1000000))
if [ "$stop_ext" -ge "$scaffold_length" ] ; then stop_ext=$((scaffold_length - 10)) ; fi
Region_extended=`echo "$curr_scaffold:$start_ext-$stop_ext"`
Hypomesus_read_name=`grep "$Hypomesus_SRR" ../Reads_spanning_regions_tables/Hypomesus_transpacificus.$curr_scaffold\:$start-$stop | awk '{print $0, $5 - $4}' | tr ',' '\t' | sort -k6 -n | tail -1 | cut -f2`
echo ">$Hypomesus_read_name" > $curr_HOG/Hypomesus_transpacificus.fa
samtools view ../$Hypomesus_SRR.Hypomesus_transpacificus.sorted.bam $Region_extended -@ 8 | grep "$Hypomesus_read_name" | cut -f10 | awk '{print $0,length}' | tr ' ' ','  | sort -k2 -n -t ',' | tail -1 | cut -f1 -d ","  >> $curr_HOG/Hypomesus_transpacificus.fa
strand_hypo=`exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/$curr_XM_hypo.prot $curr_HOG/Hypomesus_transpacificus.fa | grep "vulgar" | cut -f9 -d " "`
if [ $strand_hypo == "-" ] ; then revseq $curr_HOG/Hypomesus_transpacificus.fa $curr_HOG/Hypomesus_transpacificus.rev ; sed -i 's/ Reversed://g' $curr_HOG/Hypomesus_transpacificus.rev ; mv $curr_HOG/Hypomesus_transpacificus.rev $curr_HOG/Hypomesus_transpacificus.fa ; fi
exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/$curr_XM_hypo.prot $curr_HOG/Hypomesus_transpacificus.fa > $curr_HOG/Hypomesus_transpacificus.exo


cd $curr_HOG ; for i in *.fa ; do samtools faidx $i ; done
for i in *.fai ; do echo "$i" ; cat $i ; done

grep "	exon	" Clupea_harengus.exo | cut -f4,5 | sed 's/^/Clupea	/g'
grep "	exon	" Hypomesus_transpacificus.exo | cut -f4,5 | sed 's/^/Hypomesus	/g'




#####====================================================================================================================================
#####====================================================================================================================================
#####================================================== N5.HOG0004480  ==================================================================
#####====================================================================================================================================
#####====================================================================================================================================


curr_HOG=N5.HOG0004480

rm -rf $curr_HOG ; mkdir $curr_HOG


curr_XM_clupea=XM_031565097
full_XM_name=`grep "$curr_XM_clupea" Clupea_harengus.prot.nostop | sed 's/>//g'`
samtools faidx Clupea_harengus.prot.nostop $full_XM_name > $curr_HOG/$curr_XM_clupea.prot
Clupea_SRR=SRR26322996
curr_scaffold=`grep "$curr_XM_clupea" Table_longreads_count.clupea.txt | cut -f4 | sort | uniq | cut -f1 -d ":"`
start=`grep "$curr_XM_clupea" Table_longreads_count.clupea.txt | cut -f4 | sort | uniq | cut -f2 -d ":" | cut -f1 -d "-"`
stop=`grep "$curr_XM_clupea" Table_longreads_count.clupea.txt | cut -f4 | sort | uniq | cut -f2 -d ":" | cut -f2 -d "-"`
scaffold_length=`grep "$curr_scaffold	" GCF_900700415.2_Ch_v2.0.2_genomic.fna.fai | cut -f2`
start_ext=$((start - 1000000))
if [ "$start_ext" -le 0 ]; then start_ext=1 ; fi
stop_ext=$((stop + 1000000))
if [ "$stop_ext" -ge "$scaffold_length" ] ; then stop_ext=$((scaffold_length - 10)) ; fi
Region_extended=`echo "$curr_scaffold:$start_ext-$stop_ext"`
Clupea_read_name=`grep "$Clupea_SRR" ../Reads_spanning_regions_tables/Clupea_harengus.$curr_scaffold\:$start-$stop | awk '{print $0, $5 - $4}' | tr ',' '\t' | sort -k6 -n | tail -1 | cut -f2`
echo ">$Clupea_read_name" > $curr_HOG/Clupea_harengus.fa
samtools view ../$Clupea_SRR.Clupea_harengus.sorted.bam $Region_extended -@ 8 | grep "$Clupea_read_name" | cut -f10 | awk '{print $0,length}' | tr ' ' ','  | sort -k2 -n -t ',' | tail -1 | cut -f1 -d ","  >> $curr_HOG/Clupea_harengus.fa
strand_clupea=`exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/$curr_XM_clupea.prot $curr_HOG/Clupea_harengus.fa | grep "vulgar" | cut -f9 -d " "`
if [ $strand_clupea == "-" ] ; then revseq $curr_HOG/Clupea_harengus.fa $curr_HOG/Clupea_harengus.rev ; sed -i 's/ Reversed://g' $curr_HOG/Clupea_harengus.rev ; mv $curr_HOG/Clupea_harengus.rev $curr_HOG/Clupea_harengus.fa ; fi
exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/$curr_XM_clupea.prot $curr_HOG/Clupea_harengus.fa > $curr_HOG/Clupea_harengus.exo

#curr_XM_hypo=XM_047020147
#full_XM_name=`grep "$curr_XM_hypo" Hypomesus_transpacificus.prot.nostop | sed 's/>//g'`
#samtools faidx Hypomesus_transpacificus.prot.nostop $full_XM_name > $curr_HOG/$curr_XM_hypo.prot
Hypomesus_SRR=SRR17799534
curr_scaffold=NC_061079.1
start=3801703
stop=3809259
scaffold_length=`grep "$curr_scaffold	" GCF_021917145.1_fHypTra1_genomic.fna.fai | cut -f2`
start_ext=$((start - 1000000))
if [ "$start_ext" -le 0 ]; then start_ext=1 ; fi
stop_ext=$((stop + 1000000))
if [ "$stop_ext" -ge "$scaffold_length" ] ; then stop_ext=$((scaffold_length - 10)) ; fi
Region_extended=`echo "$curr_scaffold:$start_ext-$stop_ext"`
Hypomesus_read_name=`grep "$Hypomesus_SRR" ../Reads_spanning_regions_tables/Hypomesus_transpacificus.$curr_scaffold\:$start-$stop | awk '{print $0, $5 - $4}' | tr ',' '\t' | sort -k6 -n | tail -1 | cut -f2`
echo ">$Hypomesus_read_name" > $curr_HOG/Hypomesus_transpacificus.fa
samtools view ../$Hypomesus_SRR.Hypomesus_transpacificus.sorted.bam $Region_extended -@ 8 | grep "$Hypomesus_read_name" | cut -f10 | awk '{print $0,length}' | tr ' ' ','  | sort -k2 -n -t ',' | tail -1 | cut -f1 -d ","  >> $curr_HOG/Hypomesus_transpacificus.fa
strand_hypo=`exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/$curr_XM_clupea.prot $curr_HOG/Hypomesus_transpacificus.fa | grep "vulgar" | cut -f9 -d " "`
if [ $strand_hypo == "-" ] ; then revseq $curr_HOG/Hypomesus_transpacificus.fa $curr_HOG/Hypomesus_transpacificus.rev ; sed -i 's/ Reversed://g' $curr_HOG/Hypomesus_transpacificus.rev ; mv $curr_HOG/Hypomesus_transpacificus.rev $curr_HOG/Hypomesus_transpacificus.fa ; fi
exonerate --showtargetgff TRUE --model protein2genome --ryo "%tcs" --bestn 1 $curr_HOG/$curr_XM_clupea.prot $curr_HOG/Hypomesus_transpacificus.fa > $curr_HOG/Hypomesus_transpacificus.exo


cd $curr_HOG ; for i in *.fa ; do samtools faidx $i ; done
for i in *.fai ; do echo "$i" ; cat $i ; done

grep "	exon	" Clupea_harengus.exo | cut -f4,5 | sed 's/^/Clupea	/g'
grep "	exon	" Hypomesus_transpacificus.exo | cut -f4,5 | sed 's/^/Hypomesus	/g'






#####====================================================================================================================================
#####====================================================================================================================================
#####========================================== Accessory scripts =======================================================================
#####====================================================================================================================================
#####====================================================================================================================================


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

