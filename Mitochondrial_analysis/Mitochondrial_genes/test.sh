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

#Finally, add lepisosteus as outgroup





cd Alignments/

for file in *.fa ; do transeq $file $file.prot -table 2 ; sed -i 's/_1$//g' $file.prot  ; done
for file in *.prot ; do /scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Opsins_analysis_2024/muscle5.1.linux_intel64 -align $file -output $file.aln ; done
for file in *.aln ; do trimal -in $file -gt 0.8 -cons 70 -out $file.trimal ; done

cd ../

python3 AMAS/amas/AMAS.py concat -f fasta -d aa -i Alignments/*.trimal --part-format nexus
sed -i 's/?/-/g' concatenated.out
mv concatenated.out Concatenated_mitochondrial_genes.aln


iqtree -s Concatenated_mitochondrial_genes.aln --seqtype AA -m mtVer+F+G4 -nt 8 -bb 1000
