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


samtools faidx Proteomes_BUSCO80/Paramormyrops_kingsleyae.fa Paramormyrops_kingsleyae---rna-XM_023822687.1 > Paramormyrops_kingsleyae---rna-XM_023822687.1.prot

#blastp against the fish proteome (50 best hits) + Uniprot database without any fish species (10 best hits). Here I will only take 20 catffishes otherwise there are way too much of them

blastp -query Paramormyrops_kingsleyae---rna-XM_023822687.1.prot -db concatenated_proteomes.fa -outfmt 6 -out Gene_vs_FishProteome.blastp -num_threads 8 -max_target_seqs 51
blastp -query Paramormyrops_kingsleyae---rna-XM_023822687.1.prot -db non_actino_uniprot.fa -outfmt 6 -out Gene_vs_Uniprot.blastp -num_threads 8  -max_target_seqs 10 -evalue 1e-5


cut -f2 Gene_vs_FishProteome.blastp | sort | uniq  > closest_fish_seq.id
cut -f2 Gene_vs_Uniprot.blastp | sort | uniq > closest_nonfish_seq.id 

xargs samtools faidx concatenated_proteomes.fa < closest_fish_seq.id > closest_fish_seq.fa
xargs samtools faidx non_actino_uniprot.fa  < closest_nonfish_seq.id > closest_nonfish_seq.fa

#Align with muscle and trim with trimal

cat closest_fish_seq.fa closest_nonfish_seq.fa > Uniprot_plus_closefish.fa

muscle5.1 -align closest_fish_seq.fa -output closest_fish_seq.aln
mafft --add closest_nonfish_seq.fa  --keeplength closest_fish_seq.aln > Uniprot_plus_closefish.aln
trimal -in Uniprot_plus_closefish.aln -automated1 -out Uniprot_plus_closefish.aln.trimmed

#Make a ML tree with IQ-TREE
iqtree -s Uniprot_plus_closefish.aln -st AA -nt 8 -m TEST -mrate G4 --redo
iqtree -s Uniprot_plus_closefish.aln.trimmed -st AA -nt 8 -m TEST -mrate G4 --redo


#====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##========================================== Correct the annotation and search in other genomes ===================================
##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================

#First, predict the good complete sequence on the current genome of Paramoryrops kinglsaye

tblastn -query Paramormyrops_kingsleyae---rna-XM_023822687.1.prot -db GCF_002872115.1_PKINGS_0.1_genomic.fna
samtools faidx GCF_002872115.1_PKINGS_0.1_genomic.fna NW_019712620.1:1-30000 > test.fa
exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs" Ictalurus_punctatus_rna_XM_053685012_1.prot test.fa

TRUE location = NW_019712620.1:11591-20856

##Full_length_Paramormyrops_kingsleyae.XM_023822687.1.fa
##>Paramormyrops_kingsleyae.Full_length_sequence.XM_023822687.1
##ATGGCCAGCAAGGGCAAACAGTGGGTTCTCCTGGCAGCAGGTTCTAATGGTTGGACAGACTACGGTGTCC
##AGGCCAATGTGTGCCATGCTTACCAGGTCGCTCACCATAATGGTGTTCCAGATGAGCAGATTGTGGTGAT
##GGTCTATGATGATATTGCTTATAATCAACAAAATCCCATTCCCGGAAACATAATCAATGTACCAAATGGT
##CCAAATGTTTACCCTGGGGTTCTAAAGGACTACACTGGAGAGGAGGTATCAGCTGGAAACTTTCTGGCTG
##TACTGTGCGGAGATTTGGCTGCTGTCAAGAAAACAGGACCAAAAAAAGTAATACAGAGtgGTGAAAATGA
##TTCTATCTTCATCTACCTATCAGACCATGGAAATAAAGGAATCTTTCACTTTCCTAATTCCACACTTTAC
##GCACATGATCTCATTAACACAGTGAAGGAGATGTCAAAGAGTCACAAATTTTCAAAGATGGTCATTTATT
##TGGAAGCGTGTCATGCTGGATCAATGCTGAATCAGCTGTCTGACAGCAACGtgtaTGCAGTTTCTTCATG
##CAACCCTGATGAATATACCTATGCTTGTTTCTTTGACAAGAAAAGAAATACCTTTCTATCTGATATCTTC
##AGCTTTAACTGGCTACATCACATGGACACAATAAAACTTACTGCGACTTCATTTGGAGAGCAGTTTTCCT
##ACCTGGAGAAAAACGTGAGTAAGGATGCGAGGAAGGCAGGAGTGACTATGACACCATGTAACTATGGGGA
##CAAGATCATGCTGAAGTTAATGTTGAGTGAATTTCTTGGAGAATCCCCGGCCTCTGTCAGAGAGACCTAC
##ATGTTCCAGCTCCAAGTGTCTGATGTGGTCGACACTACAGAAGTTCCACTGATAAtccaaaaaaacaaga
##taaaaaatgaacaggacCCAGAAAAAAGACAGACCCTGGAAAAACAATATGATGATCTTAAAAAAAAGAG
##GAAGACAGTGGATGAAGCACTGCAGAAGATTGCTGAGCGCACAAACACTTCAGGAGCTCTGACTGAAAAA
##CGTGAGGTGACTCGCACATACGAGCTCAAAGTCGTGGCTGAGCATTTTAGGAAAAATCTCTTTAACTGGG
##AGGAGGAGCcgtttgTAGTCACTCGTTCACACCTGCAGGTTTTAGTGAATCTGTGTGAATGTGGGTTGGA
##GCTTGAGAGCATCACTGACGCCATCACTCACGTGAGCCAGGAAATTACTTTC



#Sarch the gene on an alternative assembly of Paramormyrops kingslaye + Paramormyrops_hopkinsi + Paramormyrops sp. MAG

makeblastdb -in GCA_025177845.1_MSU160_NANOPORE_genomic.fna -dbtype nucl
makeblastdb -in GCA_026780785.1_ASM2678078v1_genomic.fna -dbtype nucl
makeblastdb -in GCA_027576875.1_ASM2757687v1_genomic.fna -dbtype nucl

samtools faidx Paramormyrops_kingsleyae.cds rna-XM_023822687.1 > Paramormyrops_kingsleyae---rna-XM_023822687.1.fa 
samtools faidx Proteomes_BUSCO80/Ictalurus_punctatus.fa Ictalurus_punctatus_rna_XM_053685012_1 > Ictalurus_punctatus_rna_XM_053685012_1.prot

tblastn -query Paramormyrops_kingsleyae---rna-XM_023822687.1.prot -db GCA_025177845.1_MSU160_NANOPORE_genomic.fna
samtools faidx GCA_025177845.1_MSU160_NANOPORE_genomic.fna JAMQYA010000385.1:395703-435703 > Paramormyrops_kingsleyae.alternative.fa
exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs"  Ictalurus_punctatus_rna_XM_053685012_1.prot Paramormyrops_kingsleyae.alternative.fa > Paramormyrops_kingsleyae.alternative.exo

tblastn -query Ictalurus_punctatus_rna_XM_053685012_1.prot -db GCA_026780785.1_ASM2678078v1_genomic.fna
samtools faidx GCA_026780785.1_ASM2678078v1_genomic.fna JAODKI010006440.1:1-24432 > Paramormyrops_hopkinsi.fa
exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs"  Ictalurus_punctatus_rna_XM_053685012_1.prot Paramormyrops_hopkinsi.fa > Paramormyrops_hopkinsi.exo

tblastn -query Ictalurus_punctatus_rna_XM_053685012_1.prot -db GCA_027576875.1_ASM2757687v1_genomic.fna
samtools faidx GCA_027576875.1_ASM2757687v1_genomic.fna JAODKJ010001283.1:1-24432 > Paramormyrops_sp_MAG.fa
exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs"  Ictalurus_punctatus_rna_XM_053685012_1.prot Paramormyrops_sp_MAG.fa > Paramormyrops_sp_MAG.exo






###other_paramoryrops_gene.fa
##>Paramormyrops_kingsleyae.Alternative_assembly.GCA_025177845.1
##ATGGCCAGCAAGGGCAAACAGTGGGTTCTCCTGGCAGCAGGTTCTAATGGTTGGACAGACTACGGTGTCC
##AGGCCAATGTGTGCCATGCTTACCAGGTCGCTCACCATAATGGTGTTCCAGATGAGCAGATTGTGGTGAT
##GGTCTATGATGATATTGCTTATAATCAACAAAATCCCATTCCTGGAAACATAATCAATGTACCAAATGGT
##CCAAATGTTTACCCGGGGGTTCTAAAGGACTACACTGGAGAGGAGGTATCAGCTGGAAACTTTCTGGCTG
##TACTGTGCGGAGATTTGGCTGCTGTCAAGAAAACAGGACCAAAAAGTAATACAGAGagtgGTGAAAATGA
##TTCTATCTTCATCTACTTATCAGACCATGGAAATAAAGGAATCTTTCACTTTCCTAATTCCACACTTTAC
##GCACATGATCTCATTAACACAGTGAAGGAGATGTCAAAGAGTCACAAATTTTCAAAGATGGTCATTTATT
##TGGAAGCGTGTCATGCTGGATCAATGCTGAATCAGCTGTCTGACAGCAACGtgtaTGCAGTTTCTTCATG
##CAACCCTGATGAATATACCTATGCTTGTTTCTTTGACAAGAAAAGAAATACCTTTCTATCTGATATCTTC
##AGCTTTAACTGGCTACATCACATGGACACAATAAAACTTACTGCGACTTCATTTGGAGAGCAGTTTTCCT
##ACCTGGAGAAAGTGAGTAAGGATGCGAGGAAGGCAGGAGTGACTATGACACCATGTAACTATGGGGACAA
##GATCATGCTGAAGTTAATGTTGAGTGAATTTCTTGAGAATCCCGGCCTCGTCAGAGAGACCTACATGTTC
##CAGCTCCAAGTGTCTGATGTGGTCGACACTACAGAAGTTCCACTGATAATCcaaaaacaagataaaaatG
##AACAGGACCCAGAAAAGACAGACCCTGGAAAACAATATGATGATCTTAAAAACAAGAGGAAGACAGTGGA
##TGAAGCACTGCAGAAGATTGCTGAGCGCACAAACACTTCAGGAGCTCTGACTGAAAAACGTGAGGTGACT
##ACATACGAGCTCAAAGTCGTGGCTGAGCATTTTAGGAATCTCTTTAACTGGGAGGAGGAGCcgtttgTAG
##TCACTCGTTCACACCTGCAGGTTTTAGTGAATCTGTGTGAATGTGGGTTGGAGCTTGAGAGCATCACTGA
##CGCCATCACTCACGTGAGCCAGGAAATTACTTTC
##>Paramormyrops_hopkinsi
##ATGGCCAGCAAGGGCAAACAGTGGGTTCTCCTGGCAGCAGGTTCTAATGGTTGGACAGACTACGGTGTCC
##AGGCCAATGTGTGCCATGCTTACCAGGTCGCTCACCATAACGGTGTTCCAGATGAGCAGATTGTGGTGAT
##GGTCTATGATGATATTGCTTATAATCAACAAAATCCCATTCCTGGAAACATAATCAATGTACCAAATGGT
##CCAAATGTTTACCCGGGGGTTCTAAAGGACTACACTGGAGAGGAGGTATCAGCTGGAAACTTTCTGGCTG
##TACTGTGTGGAGATTTGGCTGCTGTCAAGAAAACAGGACCAAAAAAAGTAATACAGAGtggTGAAAATGA
##TTCTATCTTCATCTACTTATCAGACCATGGAAATAAAGGAATCTTTCACTTTCCTAATTCCACACTTTAC
##GCACATGATCTCATTAACACAGTGAAGGAGATGTCAAAGAGTCACAAATTTTCAAAGATGGTCATTTATT
##TGGAAGCGTGTCATGCTGGATCAATGCTGAATCAGCTGTCTGACAGCAACGtgtaTGCAGTTTCTTCATG
##CAACCCTGATGAATATACCTATGCTTGTTTCTTTGACAAGAAAAGAAATACCTTTCTATCTGATATCTTC
##AGCTTTAACTGGCTACATCACATGGACACAATAAAACTTACTGCGACTTCATTTGGAGAGCAGTTTTCCT
##ACCTGGAGAAAAACGTGAGTAAGGATGCGAGGAAGGCAGGAGTGACTATGACACCATGTAACTATGGGGA
##CAAGATCATGCTGAAGTTAATGTTGAGTGAATTTCTTGGAGAATCCCCGGCCTCTGTCAGAGAGACATAC
##ATGTTCCAGCTCCAAGTGTCTGATGTGGTCGACACTACAGAAGTTCCACTGATAAtccaaaaaaacaaga
##taaaaaatgaacaggacCCAGAAAAAAGACAGACCCTGGAAAAACAATATGATGATCTTAAAAAAAAGAG
##GAAGACAGTGGATGAAGCACTGCAGAAGATTGCTGAGCGCACAAACACTTCAGGAGCTCTGACTGAAAAA
##CGTGAGGTGACTCGCACATATGAGCTCAAAGTCGTGGCTGAGCATTTTAGGAAAAATCTCTTTAACTGGG
##AGGAGGAGCCgtttgTAGTCACTCGTTCACACCTGCAGGTTTTAGTGAATCTGTGTGAATGTGGGTTGGA
##GCTTGAGAGCATCACTGACGCCATCACTCACGTGAGCCAGGAAATTACTTTC

transeq other_paramoryrops_gene.fa other_paramoryrops_gene.prot ; sed -i 's/_1$//g' other_paramoryrops_gene.prot
transeq Full_length_Paramormyrops_kingsleyae.XM_023822687.1.fa Full_length_Paramormyrops_kingsleyae.XM_023822687.1.prot ; sed -i 's/_1$//g' Full_length_Paramormyrops_kingsleyae.XM_023822687.1.prot

grep "Ictalurus_punctatus\|Pangasianodon_hypophthalmus\|Silurus_meridionalis\|Xyrauchen_texanus\|Colossoma_macropomum\|Electrophorus_electricus" all.id > curr_id
xargs samtools faidx Coding_sequences_alignments/N5.HOG0029633.prot < curr_id > curr_id.prot 
cat Paramormyrops_kingsleyae---rna-XM_023822687.1.prot curr_id.prot other_paramoryrops_gene.prot Full_length_Paramormyrops_kingsleyae.XM_023822687.1.prot > OGG_plus_ManuallyFound.prot


muscle5.1 -align OGG_plus_ManuallyFound.prot -output OGG_plus_ManuallyFound.aln
trimal -gt 0.7 -cons 60 -in OGG_plus_ManuallyFound.aln -out OGG_plus_ManuallyFound.aln.trimmed


iqtree -s OGG_plus_ManuallyFound.aln.trimmed -st AA -nt 8 -bb 1000 --redo -m LG+F+G4


#Recompute the Ks with the full length gene

xargs samtools faidx Coding_sequences_alignments/N5.HOG0029633.cds < curr_id > curr_id.fa
cat other_paramoryrops_gene.fa Full_length_Paramormyrops_kingsleyae.XM_023822687.1.fa curr_id.fa >  OGG_plus_ManuallyFound.fa
samtools faidx Coding_sequences_alignments/N5.HOG0029633.cds Paramormyrops_kingsleyae_rna_XM_023822687_1 | sed "s/>.*/>Paramormyrops_kingsleyae---rna-XM_023822687.1/g" >>  OGG_plus_ManuallyFound.fa


trimal -in OGG_plus_ManuallyFound.aln -cons 100 -backtrans OGG_plus_ManuallyFound.fa -out OGG_plus_ManuallyFound.fa.aln

Rscript Compute_aligned_CDS_stats_FAST.R OGG_plus_ManuallyFound.fa.aln Full_length.stats

grep "Paramormyrops_kingsleyae.Full_length_sequence.*Ictalurus_punct" Full_length.stats > Full_length.stats.1
grep "Paramormyrops_kingsleyae.Full_length_sequence.*Silurus_mer" Full_length.stats > Full_length.stats.2


#====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##========================================== Search for HOG0046344 in the alternative assembly and other paramormyrops ===============
##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================



tblastn -query Paramormyrops_kingsleyae_rna_XM_023822671_1.prot -db GCA_025177845.1_MSU160_NANOPORE_genomic.fna
samtools faidx GCA_025177845.1_MSU160_NANOPORE_genomic.fna JAMQYA010000385.1:395703-435703 > Paramormyrops_kingsleyae.HOG0046344.region.alternative.fa
exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs"  Paramormyrops_kingsleyae_rna_XM_023822671_1.prot Paramormyrops_kingsleyae.HOG0046344.region.alternative.fa > Paramormyrops_kingsleyae.HOG0046344.region.alternative.exo


tblastn -query Paramormyrops_kingsleyae_rna_XM_023822671_1.prot -db GCA_026780785.1_ASM2678078v1_genomic.fna
samtools faidx GCA_026780785.1_ASM2678078v1_genomic.fna JAODKI010006440.1:1-24432 > Paramormyrops_hopkinsi.HOG0046344.region.fa
exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs"  Paramormyrops_kingsleyae_rna_XM_023822671_1.prot Paramormyrops_hopkinsi.HOG0046344.region.fa > Paramormyrops_hopkinsi.HOG0046344.region.exp

#HOG0046344.otherPara.fa
>Paramormyrops_kingsleyae.Alternative_assembly.GCA_025177845.1
ATGATGGTGCATTTGACCTTAAACCAGCAGAAGAGGCGGacacacaccttTTCACTGTTATTACATTTCC
CAGTTACTAGTATTTTACTCAAAGTGATCAAGATGAAGAGAAAAGCAGACCTTGATGGATGTGAACCAAG
AAAGCGTGTATGTAAAGATGGAGCCAGAACCAAGGCCAAGATAACACGTGTCTCAGAAGAAATGTGCGCT
CAGAGTCCATATGCAGGGCTAGAGCGACACATGCATTTAGAAGCAAAACAAGGAACGTTTACATTTGGAG
CTGCAGCTACAGCTAAGATTAGTGCAGGACGAGTCAGAGCAAAGCGGAGCATCGTTTCAGCTGAAGCTAA
CAGTTCAGATGCTTCTGAGATAGTTGTGACAGAGATTGATGCTGAGGGAAAGTATGCAGAGGGAGGAGCA
TATGtgcatttagaaacaaaacaaggaaCATTTAGATTTGGAGCTGGAGCCACAGCTAAGGTTAGTGCAG
GACGAGTCAGAGCAAAGCGGAGCATCGTTTCAGCTGAAGCTAACAGTTCCTCCATAGTTGATATCGAAGG
ACCATATTCTGAGGCAGAGGCAAACATGAATGTTGAAAGAACACAAAGAAATCTAGAATTTTTAGCTAAG
GCCAGAGCTAAAGTTAGTGCAGGACGGATAAGAGTAACAGAGGATGACGTTTCAGCTGAAGCTAGTGGTC
CAAATGCTTCTGCCATGTTTGTGGCTAATAATAAAGAAGCACAGATCATGGTCGGAGCTGAGATTGCAAG
TGCCTCAGTTTCTGCAGGTCCAGTTGGAATGAAGGTTGGTCTTGGGGTTGAAACCGGTGCTCGCGTTGGT
GATGATGGCGCAGAAGTGAAAGTCCTGGGCTGTGGCATGTCAATCGGCAGAACGATGGGAGCttcttttt
gg
>Paramormyrops_hopkinsi
ATGACGGTGCATTTGACCTGAAACCAGCAGAAGAGGCGGacacacaccttTTCACTGTTATTACATTTCC
CAGTTACTAGTATTTTACTCAAAGTGATCAAGATGAAGAGAAAAGCAGACCTTGATGGATGTGAACCAAG
AAAGCGTGTATGTAAAGATGGAGCCAGAACCAAGGCCAAGATAACACGTGTCTCAGAAGAAATGTGCGCT
CAGAGTCCATATGCAGGGCTAGAGCGACACATGCATTTAGAAGCAAAACAAGGAACGTTTACATTTGGAG
CTGCAGCTACAGCTAAGATTAGTGCAGGACGAGTCAGAGCAAAGCGGAGCATCGTTTCAGCTGAAGCTAA
CAGTTCAGATGCTTCTGAGATGGTTGTGACAGAGATTGATGCTGAGGGAAAGTATGCAGAGGGAGGAGCA
TATGTGCATTtagaaacaaaacaaggaacATTTAGATTTGGAGCTGGAGCCACAGCTAAGGTTAGTGCAG
GACGAGTCAGAGCAAAGCGGAGCATCGTTTCAGCTGAAGCTAACAGTTCCTCCATAGTTGATATCGAAGG
ACCATATTCTGAGGCAGAGGCAAACATGAATGTTGAAAGAACACAAAGAAATCTAGAATTTTTAGCTAAG
GCCAGAGCTAAAGTTAGTGCAGGACGGATAAGAGTAACAGAGGATGACGTTTCAGCTGAAGCTAGTGGTC
CAAATGCTTCTGCCATGTTTGTGGCTAATAATAAAGAAGCACAGATCATGGTCGGAGCTGAGATTGCAAG
TGCCTCAGTTTCCGCAGGTCCAGTTGGAATGAAGGTTGGTCTTGGGGTTGAAACCGGTGCTCGCGTTGGT
GATGATGGCGCAGAAGTGAAAGTCCTGGGCTGTGGCATGTCAATCGGCAGAACGATGGGAGCtttttttg
gaaatgaactgaaatttaGGTTATGG


transeq HOG0046344.otherPara.fa HOG0046344.otherPara.prot ; sed -i 's/_1$//g' HOG0046344.otherPara.prot


grep "Ictalurus\|Pangas\|Ameir\|Param\|Triplo\|Albula\|Param" N5.HOG0046344.prot | sed 's/>//g' > curr_id
xargs samtools faidx N5.HOG0046344.prot < curr_id > curr_id.prot 

cat curr_id.prot HOG0046344.otherPara.prot > HOG0046344.otherPara.combined.prot


muscle5.1 -align HOG0046344.otherPara.combined.prot -output HOG0046344.otherPara.combined.aln
trimal -in HOG0046344.otherPara.combined.aln -gt 0.9 -out HOG0046344.otherPara.combined.aln.trimmed
iqtree -s HOG0046344.otherPara.combined.aln -st AA -nt 8 -bb 1000 --redo -m LG+F+G4


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



#N5.HOG0029633.clade_hgt.txt=> clade where the HGT is extracted from iTOL (from N5.HOG0029633.prot.aln.treefile)

cut -f1,2 -d "_" N5.HOG0029633.clade_hgt.txt | sort | uniq  > species_to_draw.clade1
#Choose close species without the gene .. 


#echo "Anguilla_anguilla" > non_transfer_species_1
echo "Scleropages_formosus" > non_transfer_species_1
echo "Brienomyrus_brachyistius" >> non_transfer_species_1

echo "Pangasianodon_hypophthalmus" > non_transfer_species_2
echo "Silurus_meridionalis" >> non_transfer_species_2
echo "Xyrauchen_texanus" >> non_transfer_species_2



rm -r Syn_tables_dir ; mkdir Syn_tables_dir


#First order the species names 

rm species_to_draw.clade1.ordered

echo "Paramormyrops_kingsleyae" >> species_to_draw.clade1.ordered
echo "Ictalurus_punctatus" >> species_to_draw.clade1.ordered



#If there are more than 4 species from the same order, keep only 4 


#First extract genes around 

curr_OGG=N5.HOG0029633

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


grep "NW_019712620_1" Syn_tables_dir/Paramormyrops_kingsleyae.synt.final.df  > temp ; mv temp Syn_tables_dir/Paramormyrops_kingsleyae.synt.final.df


#First add Brienomyrus_brachyistius

curr_OGG=N5.HOG0029633
curr_sp=Brienomyrus_brachyistius
ref_sp=Paramormyrops_kingsleyae


grep -A10 -B10 "N5_HOG0031098" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df

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




#Now add Scleropages_formosus

curr_OGG=N5.HOG0029633
curr_sp=Scleropages_formosus
ref_sp=Paramormyrops_kingsleyae

rm Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0031098" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df


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




#Now add Silurus_meridionalis

curr_OGG=N5.HOG0029633
curr_sp=Silurus_meridionalis
ref_sp=Ictalurus_punctatus

grep -A10 -B10 "N5_HOG0029633" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df


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





#Now add Xyrauchen_texanus

curr_OGG=N5.HOG0029633
curr_sp=Xyrauchen_texanus
ref_sp=Ictalurus_punctatus

grep -A10 -B10 "N5_HOG0029633" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df


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



#Now add Pangasianodon_hypophthalmus

curr_OGG=N5.HOG0029633
curr_sp=Pangasianodon_hypophthalmus
ref_sp=Ictalurus_punctatus

grep -A10 -B10 "N5_HOG0029633" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret > Syn_tables_dir/$curr_sp.synt.df


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
##========================================== Verify that the genes are indeed absent in close species  ======================================
##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================

samtools faidx N5.HOG0046344.prot Paramormyrops_kingsleyae_rna_XM_023822671_1 > Paramormyrops_kingsleyae_rna_XM_023822671_1.prot


#B. brachyistius

tblastn -query Ictalurus_punctatus_rna_XM_053685012_1.prot -db GCF_023856365.1_BBRACH_0.4_genomic.fna > Brienomyrus_brachyistius.Ictalurus_punctatus.HOG0029633.tblastn
tblastn -query Paramormyrops_kingsleyae---rna-XM_023822687.1.prot -db GCF_023856365.1_BBRACH_0.4_genomic.fna > Brienomyrus_brachyistius.Paramormyrops_kingsleyae.HOG0029633.tblastn
tblastn -query Paramormyrops_kingsleyae_rna_XM_023822671_1.prot -db GCF_023856365.1_BBRACH_0.4_genomic.fna > Brienomyrus_brachyistius.Paramormyrops_kingsleyae.HOG0046344.tblastn


exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs" Paramormyrops_kingsleyae---rna-XM_023822687.1.prot GCF_023856365.1_BBRACH_0.4_genomic.fna > Brienomyrus_brachyistius.Paramormyrops_kingsleyae.HOG0029633.exo
exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs" Paramormyrops_kingsleyae_rna_XM_023822671_1.prot GCF_023856365.1_BBRACH_0.4_genomic.fna > Brienomyrus_brachyistius.Paramormyrops_kingsleyae.HOG0046344.exo



#Scleropages_formosus

tblastn -query Ictalurus_punctatus_rna_XM_053685012_1.prot -db GCF_900964775.1_fSclFor1.1_genomic.fna > Scleropages_formosus.Ictalurus_punctatus.HOG0029633.tblastn
tblastn -query Paramormyrops_kingsleyae---rna-XM_023822687.1.prot -db GCF_900964775.1_fSclFor1.1_genomic.fna > Scleropages_formosus.Paramormyrops_kingsleyae.HOG0029633.tblastn
tblastn -query Paramormyrops_kingsleyae_rna_XM_023822671_1.prot -db GCF_900964775.1_fSclFor1.1_genomic.fna > Scleropages_formosus.Paramormyrops_kingsleyae.HOG0046344.tblastn


exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs" Paramormyrops_kingsleyae---rna-XM_023822687.1.prot GCF_900964775.1_fSclFor1.1_genomic.fna > Scleropages_formosus.Paramormyrops_kingsleyae.HOG0029633.exo
exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs" Paramormyrops_kingsleyae_rna_XM_023822671_1.prot GCF_900964775.1_fSclFor1.1_genomic.fna > Scleropages_formosus.Paramormyrops_kingsleyae.HOG0046344.exo


##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##========================================== Try to make a dotter plot between Paramormyrops and Ictalurus  ==========================
##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================

#Extract the region between N5_HOG0029633 and N5_HOG0046344 (genes included)
grep -A3 -B3 "N5_HOG0029633"  GFF3_N5_OGGs/Paramormyrops_kingsleyae.gff.simplified.sorted.OGG.tiret
samtools faidx GCF_002872115.1_PKINGS_0.1_genomic.fna NW_019712620.1:11592-32660 > Paramormyrops_kingsleyae-N5.HOG0029633-N5.HOG0046344-region.fa

grep -A3 -B3 "N5_HOG0029633"  GFF3_N5_OGGs/Ictalurus_punctatus.gff.simplified.sorted.OGG.tiret
samtools faidx GCF_001660625.3_Coco_2.0_genomic.fna NC_030428.2:28396767-28413211 > Ictalurus_punctatus-N5.HOG0029633-N5.HOG0046344-region.fa


revseq Ictalurus_punctatus-N5.HOG0029633-N5.HOG0046344-region.fa Ictalurus_punctatus-N5.HOG0029633-N5.HOG0046344-region.fa.rev


#Extract exons of N5_HOG0029633 and N5_HOG0046344 on these locations in Paramormyrops


grep "XM_023822687" GFF3_files_per_species/Paramormyrops_kingsleyae.gff | grep "exon" | cut -f4,5 | tr "\t" "," > gene1.exons
grep "XM_023822671" GFF3_files_per_species/Paramormyrops_kingsleyae.gff | grep "exon" | cut -f4,5 | tr "\t" "," > gene2.exons


real_start=11592

for line in `cat gene1.exons` ; do exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","` ; real_exon_start=$(( exon_start - real_start )) ; real_exon_stop=$(( exon_stop - real_start )) ; echo "$real_exon_start , $real_exon_stop" ; done
for line in `cat gene2.exons` ; do exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","` ; real_exon_start=$(( exon_start - real_start )) ; real_exon_stop=$(( exon_stop - real_start )) ; echo "$real_exon_start , $real_exon_stop" ; done


#Extract exons of N5_HOG0029633 and N5_HOG0046344 on these locations in Ictalurus punctatus


grep "XM_017482960" GFF3_files_per_species/Ictalurus_punctatus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene1.exons
grep "XM_017482962" GFF3_files_per_species/Ictalurus_punctatus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene2.exons

real_start=28396767
length=$((28413211 - real_start))

for line in `cat gene2.exons` ; do 
	exon_start=`echo "$line" | cut -f1 -d ","` 
	exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - real_start )) 
	real_exon_stop=$(( exon_stop - real_start ))

	exon_start_rev=$(( length - real_exon_start ))
	exon_stop_rev=$(( length - real_exon_stop ))

	echo "$exon_start_rev , $exon_stop_rev" 

done





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



#Extract the region between N5_HOG0029633 and N5_HOG0046344 (genes included) + 5kb upstream / downstream

#Paramormyrops
grep -A3 -B3 "N5_HOG0029633"  GFF3_N5_OGGs/Paramormyrops_kingsleyae.gff.simplified.sorted.OGG.tiret
samtools faidx GCF_002872115.1_PKINGS_0.1_genomic.fna NW_019712620.1:6592-32684 > N5.HOG0029633.Paramormyrops_kingsleyae.extended.1.fa
sed -i 's/:/-/g' N5.HOG0029633.Paramormyrops_kingsleyae.extended.1.fa
makeblastdb -in N5.HOG0029633.Paramormyrops_kingsleyae.extended.1.fa -dbtype nucl

#Ictalurus
grep -A3 -B3 "N5_HOG0029633"  GFF3_N5_OGGs/Ictalurus_punctatus.gff.simplified.sorted.OGG.tiret
samtools faidx GCF_001660625.3_Coco_2.0_genomic.fna NC_030428.2:28349121-28414880 > N5.HOG0029633.Ictalurus_punctatus.extended.1.fa
samtools faidx GCF_001660625.3_Coco_2.0_genomic.fna NC_030428.2:28670482-28714356 > N5.HOG0029633.Ictalurus_punctatus.extended.2.fa
sed -i 's/:/-/g' N5.HOG0029633.Ictalurus_punctatus.extended.1.fa
makeblastdb -in N5.HOG0029633.Ictalurus_punctatus.extended.1.fa -dbtype nucl
sed -i 's/:/-/g' N5.HOG0029633.Ictalurus_punctatus.extended.2.fa
makeblastdb -in N5.HOG0029633.Ictalurus_punctatus.extended.2.fa -dbtype nucl



#Now lets find TEs using the repbase and dfam databases

tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0029633.Paramormyrops_kingsleyae.extended.1.fa -evalue 1e-5 -outfmt 6 -out TE.Paramormyrops_kingsleyae.1.tblastn -num_threads 8
sed -i 's/#//g' TE.Paramormyrops_kingsleyae.1.tblastn

tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0029633.Ictalurus_punctatus.extended.1.fa -evalue 1e-5 -outfmt 6 -out TE.Ictalurus_punctatus.1.tblastn -num_threads 8
sed -i 's/#//g' TE.Ictalurus_punctatus.1.tblastn

tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0029633.Ictalurus_punctatus.extended.2.fa -evalue 1e-5 -outfmt 6 -out TE.Ictalurus_punctatus.2.tblastn -num_threads 8
sed -i 's/#//g' TE.Ictalurus_punctatus.2.tblastn


#merge tblastn hits and find the best TE match by doing a blastx


Rscript Rscript_merge_blast_hits.R TE.Paramormyrops_kingsleyae.1.tblastn TE.Paramormyrops_kingsleyae.1.tblastn.merged
Rscript Rscript_merge_blast_hits.R TE.Ictalurus_punctatus.1.tblastn TE.Ictalurus_punctatus.1.tblastn.merged
Rscript Rscript_merge_blast_hits.R TE.Ictalurus_punctatus.2.tblastn TE.Ictalurus_punctatus.2.tblastn.merged


xargs samtools faidx N5.HOG0029633.Paramormyrops_kingsleyae.extended.1.fa < TE.Paramormyrops_kingsleyae.1.tblastn.merged > TE.Paramormyrops_kingsleyae.1.BEST.fa
blastx -query TE.Paramormyrops_kingsleyae.1.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Paramormyrops_kingsleyae.1.BEST.blastx -max_target_seqs 1
cut -f1 TE.Paramormyrops_kingsleyae.1.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Paramormyrops_kingsleyae.1.BEST.blastx >> temp  ; done ; mv temp TE.Paramormyrops_kingsleyae.1.BEST.blastx


xargs samtools faidx N5.HOG0029633.Ictalurus_punctatus.extended.1.fa < TE.Ictalurus_punctatus.1.tblastn.merged > TE.Ictalurus_punctatus.1.BEST.fa
blastx -query TE.Ictalurus_punctatus.1.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Ictalurus_punctatus.1.BEST.blastx -max_target_seqs 1
cut -f1  TE.Ictalurus_punctatus.1.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Ictalurus_punctatus.1.BEST.blastx >> temp  ; done ; mv temp TE.Ictalurus_punctatus.1.BEST.blastx


xargs samtools faidx  N5.HOG0029633.Ictalurus_punctatus.extended.2.fa < TE.Ictalurus_punctatus.2.tblastn.merged > TE.Ictalurus_punctatus.2.BEST.fa
blastx -query TE.Ictalurus_punctatus.2.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Ictalurus_punctatus.2.BEST.blastx -max_target_seqs 1
cut -f1  TE.Ictalurus_punctatus.2.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Ictalurus_punctatus.2.BEST.blastx >> temp  ; done ; mv temp TE.Ictalurus_punctatus.2.BEST.blastx



#Now find shared elements

cut -f2 TE.Paramormyrops_kingsleyae.1.BEST.blastx | sort | uniq > TE.Paramormyrops_kingsleyae.1.uniqTE
cut -f2 TE.Ictalurus_punctatus.1.BEST.blastx | sort | uniq > TE.Ictalurus_punctatus.1.uniqTE
cut -f2 TE.Ictalurus_punctatus.2.BEST.blastx | sort | uniq > TE.Ictalurus_punctatus.2.uniqTE

comm -12 TE.Paramormyrops_kingsleyae.1.uniqTE TE.Ictalurus_punctatus.1.uniqTE
comm -12 TE.Paramormyrops_kingsleyae.1.uniqTE TE.Ictalurus_punctatus.2.uniqTE


#There is one shared element : L2-5_GA_1p:ClassI:LINE:Jockey:L2

#### Now let's find every copy of this element in the genome of both species + closely related species


samtools faidx Dfam_plus_Repbase.cdhit80.prot L2-5_GA_1p:ClassI:LINE:Jockey:L2 > L2-5.prot


#Paramormyrops_kingsleyae
tblastn -query L2-5.prot -db GCF_002872115.1_PKINGS_0.1_genomic.fna -evalue 1e-2 -outfmt 6 -out L2-5.Paramormyrops_kingsleyae.tblastn -num_threads 8
sed -i 's/#//g' L2-5.Paramormyrops_kingsleyae.tblastn
Rscript Rscript_merge_blast_hits.R L2-5.Paramormyrops_kingsleyae.tblastn L2-5.Paramormyrops_kingsleyae.tblastn.merged
xargs samtools faidx GCF_002872115.1_PKINGS_0.1_genomic.fna < L2-5.Paramormyrops_kingsleyae.tblastn.merged > L2-5.Paramormyrops_kingsleyae.BEST.fa
diamond blastx --query L2-5.Paramormyrops_kingsleyae.BEST.fa --db Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out L2-5.Paramormyrops_kingsleyae.BEST.blastx --max-target-seqs 1
grep "L2-5" L2-5.Paramormyrops_kingsleyae.BEST.blastx > L2-5.Paramormyrops_kingsleyae.list

#Ictalurus_punctatus
tblastn -query L2-5.prot -db GCF_001660625.3_Coco_2.0_genomic.fna -evalue 1e-2 -outfmt 6 -out L2-5.Ictalurus_punctatus.tblastn -num_threads 8
sed -i 's/#//g' L2-5.Ictalurus_punctatus.tblastn
Rscript Rscript_merge_blast_hits.R L2-5.Ictalurus_punctatus.tblastn L2-5.Ictalurus_punctatus.tblastn.merged
xargs samtools faidx GCF_001660625.3_Coco_2.0_genomic.fna  < L2-5.Ictalurus_punctatus.tblastn.merged > L2-5.Ictalurus_punctatus.BEST.fa
diamond blastx --query L2-5.Ictalurus_punctatus.BEST.fa --db Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out L2-5.Ictalurus_punctatus.BEST.blastx --max-target-seqs 1
grep "L2-5" L2-5.Ictalurus_punctatus.BEST.blastx > L2-5.Ictalurus_punctatus.list


#Brienomyrus_brachyistius
tblastn -query L2-5.prot -db GCF_023856365.1_BBRACH_0.4_genomic.fna -evalue 1e-2 -outfmt 6 -out L2-5.Brienomyrus_brachyistius.tblastn -num_threads 8
sed -i 's/#//g' L2-5.Brienomyrus_brachyistius.tblastn
Rscript Rscript_merge_blast_hits.R L2-5.Brienomyrus_brachyistius.tblastn L2-5.Brienomyrus_brachyistius.tblastn.merged
xargs samtools faidx GCF_023856365.1_BBRACH_0.4_genomic.fna  < L2-5.Brienomyrus_brachyistius.tblastn.merged > L2-5.Brienomyrus_brachyistius.BEST.fa
diamond blastx --query L2-5.Brienomyrus_brachyistius.BEST.fa --db Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out L2-5.Brienomyrus_brachyistius.BEST.blastx --max-target-seqs 1
grep "L2-5" L2-5.Brienomyrus_brachyistius.BEST.blastx > L2-5.Brienomyrus_brachyistius.list

#Scleropages_formosus
tblastn -query L2-5.prot -db GCF_900964775.1_fSclFor1.1_genomic.fna -evalue 1e-2 -outfmt 6 -out L2-5.Scleropages_formosus.tblastn -num_threads 8
sed -i 's/#//g' L2-5.Scleropages_formosus.tblastn
Rscript Rscript_merge_blast_hits.R L2-5.Scleropages_formosus.tblastn L2-5.Scleropages_formosus.tblastn.merged
xargs samtools faidx GCF_900964775.1_fSclFor1.1_genomic.fna  < L2-5.Scleropages_formosus.tblastn.merged > L2-5.Scleropages_formosus.BEST.fa
diamond blastx --query L2-5.Scleropages_formosus.BEST.fa --db Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out L2-5.Scleropages_formosus.BEST.blastx --max-target-seqs 1
grep "L2-5" L2-5.Scleropages_formosus.BEST.blastx > L2-5.Scleropages_formosus.list


#Pangasianodon_hypophthalmus
tblastn -query L2-5.prot -db GCF_027358585.1_fPanHyp1.pri_genomic.fna -evalue 1e-2 -outfmt 6 -out L2-5.Pangasianodon_hypophthalmus.tblastn -num_threads 8
sed -i 's/#//g' L2-5.Pangasianodon_hypophthalmus.tblastn
Rscript Rscript_merge_blast_hits.R L2-5.Pangasianodon_hypophthalmus.tblastn L2-5.Pangasianodon_hypophthalmus.tblastn.merged
xargs samtools faidx GCF_027358585.1_fPanHyp1.pri_genomic.fna  < L2-5.Pangasianodon_hypophthalmus.tblastn.merged > L2-5.Pangasianodon_hypophthalmus.BEST.fa
diamond blastx --query L2-5.Pangasianodon_hypophthalmus.BEST.fa --db Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out L2-5.Pangasianodon_hypophthalmus.BEST.blastx --max-target-seqs 1
grep "L2-5" L2-5.Pangasianodon_hypophthalmus.BEST.blastx > L2-5.Pangasianodon_hypophthalmus.list


#Silurus_meridionalis
tblastn -query L2-5.prot -db GCF_014805685.1_ASM1480568v1_genomic.fna -evalue 1e-2 -outfmt 6 -out L2-5.Silurus_meridionalis.tblastn -num_threads 8
sed -i 's/#//g' L2-5.Silurus_meridionalis.tblastn
Rscript Rscript_merge_blast_hits.R L2-5.Silurus_meridionalis.tblastn L2-5.Silurus_meridionalis.tblastn.merged
xargs samtools faidx GCF_014805685.1_ASM1480568v1_genomic.fna  < L2-5.Silurus_meridionalis.tblastn.merged > L2-5.Silurus_meridionalis.BEST.fa
diamond blastx --query L2-5.Silurus_meridionalis.BEST.fa --db Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out L2-5.Silurus_meridionalis.BEST.blastx --max-target-seqs 1
grep "L2-5" L2-5.Silurus_meridionalis.BEST.blastx > L2-5.Silurus_meridionalis.list

#Xyrauchen_texanus
tblastn -query L2-5.prot -db GCF_025860055.1_RBS_HiC_50CHRs_genomic.fna -evalue 1e-2 -outfmt 6 -out L2-5.Xyrauchen_texanus.tblastn -num_threads 8
sed -i 's/#//g' L2-5.Xyrauchen_texanus.tblastn
Rscript Rscript_merge_blast_hits.R L2-5.Xyrauchen_texanus.tblastn L2-5.Xyrauchen_texanus.tblastn.merged
xargs samtools faidx GCF_025860055.1_RBS_HiC_50CHRs_genomic.fna  < L2-5.Xyrauchen_texanus.tblastn.merged > L2-5.Xyrauchen_texanus.BEST.fa
diamond blastx --query L2-5.Xyrauchen_texanus.BEST.fa --db Dfam_plus_Repbase.cdhit80 --outfmt 6 -p 8 --out L2-5.Xyrauchen_texanus.BEST.blastx --max-target-seqs 1
grep "L2-5" L2-5.Xyrauchen_texanus.BEST.blastx > L2-5.Xyrauchen_texanus.list


wc -l L2-5.Scleropages_formosus.list
wc -l L2-5.Brienomyrus_brachyistius.list
wc -l L2-5.Paramormyrops_kingsleyae.list
wc -l L2-5.Ictalurus_punctatus.list
wc -l L2-5.Pangasianodon_hypophthalmus.list
wc -l L2-5.Silurus_meridionalis.list 
wc -l L2-5.Xyrauchen_texanus.list


### In a FINAL step we will make another synteny plot, only between Paramormyrops and Ictalurus to show a zoom on the region


## Paramormyrops

rm seq_clustered_infos_ogg.TE.txt ; rm clusters_ID_TE.txt


scaffold=`grep ">" N5.HOG0029633.Paramormyrops_kingsleyae.extended.1.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0029633.Paramormyrops_kingsleyae.extended.1.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0029633.Paramormyrops_kingsleyae.extended.1.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Paramormyrops_kingsleyae,$scaffold,$length" > clusters_ID_TE.txt

grep "N5_HOG0029633"  GFF3_N5_OGGs/Paramormyrops_kingsleyae.gff.simplified.sorted.OGG.tiret
grep "N5_HOG0046344"  GFF3_N5_OGGs/Paramormyrops_kingsleyae.gff.simplified.sorted.OGG.tiret

grep "XM_023822687" GFF3_files_per_species/Paramormyrops_kingsleyae.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons
grep "XM_023822671" GFF3_files_per_species/Paramormyrops_kingsleyae.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_2.exons


nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene1.$nbexon,N5_HOG0029633,+,Paramormyrops_kingsleyae" >> seq_clustered_infos_ogg.TE.txt
done
nbexon=0
for line in `cat gene_2.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene2.$nbexon,N5_HOG0046344,-,Paramormyrops_kingsleyae" >> seq_clustered_infos_ogg.TE.txt
done


cut -f1 TE.Paramormyrops_kingsleyae.1.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Paramormyrops_kingsleyae.1.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0029633.Paramormyrops_kingsleyae.extended.1.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0029633.Paramormyrops_kingsleyae.extended.1.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold,/g" | sed "s/$/,$strand,Paramormyrops_kingsleyae/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done




## Ictalurus punctatus


scaffold=`grep ">" N5.HOG0029633.Ictalurus_punctatus.extended.1.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0029633.Ictalurus_punctatus.extended.1.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0029633.Ictalurus_punctatus.extended.1.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Ictalurus_punctatus,$scaffold,$length" >> clusters_ID_TE.txt

grep -A3 -B3 "N5_HOG0029633"  GFF3_N5_OGGs/Ictalurus_punctatus.gff.simplified.sorted.OGG.tiret

grep "XM_053684837" GFF3_files_per_species/Ictalurus_punctatus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons
grep "XM_017482959" GFF3_files_per_species/Ictalurus_punctatus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_2.exons
grep "XM_017482962" GFF3_files_per_species/Ictalurus_punctatus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_3.exons
grep "XM_017482960" GFF3_files_per_species/Ictalurus_punctatus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_4.exons


nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene1.$nbexon,N5_HOG0029633,-,Ictalurus_punctatus" >> seq_clustered_infos_ogg.TE.txt
done

nbexon=0
for line in `cat gene_2.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene2.$nbexon,N5_HOG0029633,+,Ictalurus_punctatus" >> seq_clustered_infos_ogg.TE.txt
done
nbexon=0
for line in `cat gene_3.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene3.$nbexon,N5_HOG0046344,+,Ictalurus_punctatus" >> seq_clustered_infos_ogg.TE.txt
done
nbexon=0
for line in `cat gene_4.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene4.$nbexon,N5_HOG0029633,-,Ictalurus_punctatus" >> seq_clustered_infos_ogg.TE.txt
done


cut -f1 TE.Ictalurus_punctatus.1.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Ictalurus_punctatus.1.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0029633.Ictalurus_punctatus.extended.1.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0029633.Ictalurus_punctatus.extended.1.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold,/g" | sed "s/$/,$strand,Ictalurus_punctatus/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done



#Second region of Ictalurus


scaffold=`grep ">" N5.HOG0029633.Ictalurus_punctatus.extended.2.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0029633.Ictalurus_punctatus.extended.2.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0029633.Ictalurus_punctatus.extended.2.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Ictalurus_punctatus,$scaffold.sec,$length" >> clusters_ID_TE.txt

grep "N5_HOG0029633"  GFF3_N5_OGGs/Ictalurus_punctatus.gff.simplified.sorted.OGG.tiret
grep "XM_053685012" GFF3_files_per_species/Ictalurus_punctatus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons
grep "XM_053685197" GFF3_files_per_species/Ictalurus_punctatus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_2.exons



nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold.sec,$real_exon_start,$real_exon_stop,gene5.$nbexon,N5_HOG0029633,+,Ictalurus_punctatus" >> seq_clustered_infos_ogg.TE.txt
done
nbexon=0
for line in `cat gene_2.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold.sec,$real_exon_start,$real_exon_stop,gene6.$nbexon,N5_HOG0029633,-,Ictalurus_punctatus" >> seq_clustered_infos_ogg.TE.txt
done



cut -f1 TE.Ictalurus_punctatus.2.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Ictalurus_punctatus.2.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0029633.Ictalurus_punctatus.extended.2.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0029633.Ictalurus_punctatus.extended.2.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold.sec,/g" | sed "s/$/,$strand,Ictalurus_punctatus/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done


#Extract the conserved non exonic element observed in the dotter plot :


samtools faidx N5.HOG0029633.Paramormyrops_kingsleyae.extended.1.fa.rev NW_019712620.1-6592-32684:21093-22800 > similarity_region_1.Paramormyrops_kingsleyae.fa


header=`grep ">" N5.HOG0029633.Ictalurus_punctatus.extended.1.fa | sed 's/>//'`
samtools faidx N5.HOG0029633.Ictalurus_punctatus.extended.1.fa $header:59289-61289 > similarity_region_1.Ictalurus_punctatus.fa

cat similarity_region_1.Paramormyrops_kingsleyae.fa similarity_region_1.Ictalurus_punctatus.fa > similarity_region_1.combined.fa

needle -asequence similarity_region_1.Paramormyrops_kingsleyae.fa -bsequence similarity_region_1.Ictalurus_punctatus.fa -outfile similarity_region_1.aln -gapopen 10.0 -gapextend 0.5


muscle5.1 -align  similarity_region_1.combined.fa -output similarity_region_1.combined.aln

trimal -in similarity_region_1.combined.aln -out similarity_region_1.combined.aln.trimmed -gt 0.9


sbatch --qos=6hours -c 10 --mem=20G launch_blast_allgenomes.sh similarity_region_1.Paramormyrops_kingsleyae.fa 1e-20 region1.blastn.tsv





##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##========================================== Selective pressure analysis - HOG0029633 ============================================================
##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================



#Compute the dN/dS on every branches
sbatch --qos=1day -c 4 --mem=10G -e error.fitMG4.out -o slurm.fitMG4.out launch_fitMG4.sh N5.HOG0029633


Recipient branches : 
Paramormyrops_kingsleyae_rna_XM_023822687_1


#Test positive selection and relaxed selection on receiver branch


sbatch --qos=1day -c 4 --mem=10G -e error.absrel.cand.out -o slurm.absrel.cand.out launch_absrel_cand.sh N5.HOG0029633
#Relax not launched, as no branch under positive selection detected with aBSREL


#Extract dN/dS to table

grep "LB\":" N5.HOG0029633.cds.aln.FITTER.json | sed 's/.*://g' | sed 's/,.*//g' > curr_LB_values.txt
grep "MLE\":" N5.HOG0029633.cds.aln.FITTER.json | sed 's/.*://g' | sed 's/,.*//g' > curr_MLE_values.txt
grep "UB\":" N5.HOG0029633.cds.aln.FITTER.json | sed 's/.*://g' | sed 's/,.*//g' > curr_UB_values.txt
grep "\"dN\"" N5.HOG0029633.cds.aln.FITTER.json  | grep -v "{" | sed 's/.*://g' | sed 's/,//g' > curr_dN_values.txt
grep "\"dS\"" N5.HOG0029633.cds.aln.FITTER.json  | grep -v "{" | sed 's/.*://g' | sed 's/,//g' > curr_dS_values.txt
grep -B2 "LB\":" N5.HOG0029633.cds.aln.FITTER.json | grep -v "\-\-" | grep -v "Confidence" | grep -v "LB\":"  | sed 's/\"//g' | sed 's/:.*//g' | sed 's/^ *//g' > curr_labels

paste -d "," curr_labels curr_LB_values.txt curr_MLE_values.txt curr_UB_values.txt curr_dN_values.txt curr_dS_values.txt > N5.HOG0029633.dN_dS.csv



##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##========================================== Selective pressure analysis -- HOG0046344 ============================================================
##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================



#Compute the dN/dS on every branches
sbatch --qos=1day -c 4 --mem=10G -e error.fitMG4.HOG0046344.out -o slurm.fitMG4.HOG0046344.out launch_fitMG4.HOG0046344.sh N5.HOG0046344


Recipient branches : 
Paramormyrops_kingsleyae_rna_XM_023822671_1

#Test positive selection and relaxed selection on receiver branch


sbatch --qos=1day -c 4 --mem=10G -e error.absrel.cand.HOG0046344.out -o slurm.absrel.cand.HOG0046344.out launch_absrel_cand.HOG0046344.sh N5.HOG0046344

#Now lets launch RELAX

sbatch --qos=1week -c 4 --mem=10G -e error.relax.HOG0046344.out -o slurm.relax.HOG0046344.out --job-name=HOG0046344 launch_RELAX.HOG0046344.sh N5.HOG0046344


#Extract dN/dS to table

grep "LB\":" N5.HOG0046344.cds.aln.FITTER.json | sed 's/.*://g' | sed 's/,.*//g' > curr_LB_values.txt
grep "MLE\":" N5.HOG0046344.cds.aln.FITTER.json | sed 's/.*://g' | sed 's/,.*//g' > curr_MLE_values.txt
grep "UB\":" N5.HOG0046344.cds.aln.FITTER.json | sed 's/.*://g' | sed 's/,.*//g' > curr_UB_values.txt
grep "\"dN\"" N5.HOG0046344.cds.aln.FITTER.json  | grep -v "{" | sed 's/.*://g' | sed 's/,//g' > curr_dN_values.txt
grep "\"dS\"" N5.HOG0046344.cds.aln.FITTER.json  | grep -v "{" | sed 's/.*://g' | sed 's/,//g' > curr_dS_values.txt
grep -B2 "LB\":" N5.HOG0046344.cds.aln.FITTER.json | grep -v "\-\-" | grep -v "Confidence" | grep -v "LB\":"  | sed 's/\"//g' | sed 's/:.*//g' | sed 's/^ *//g' > curr_labels

paste -d "," curr_labels curr_LB_values.txt curr_MLE_values.txt curr_UB_values.txt curr_dN_values.txt curr_dS_values.txt > N5.HOG0046344.dN_dS.csv
















