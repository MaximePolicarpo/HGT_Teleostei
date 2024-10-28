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


samtools faidx Proteomes_BUSCO80/Paramormyrops_kingsleyae.fa Paramormyrops_kingsleyae---rna-XM_023831520.1 > Paramormyrops_kingsleyae---rna-XM_023831520.1.prot

#blastp against the fish proteome (50 best hits) + Uniprot database without any fish species (10 best hits). Here I will only take 20 catffishes otherwise there are way too much of them

blastp -query Paramormyrops_kingsleyae---rna-XM_023831520.1.prot -db concatenated_proteomes.fa -outfmt 6 -out Gene_vs_FishProteome.blastp -num_threads 8 -max_target_seqs 51
blastp -query Paramormyrops_kingsleyae---rna-XM_023831520.1.prot -db non_actino_uniprot.fa -outfmt 6 -out Gene_vs_Uniprot.blastp -num_threads 8  -max_target_seqs 10 -evalue 1e-5


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




##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##========================================== Alternative assembly search  ============================================================
##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================

#Sarch the gene on an alternative assembly of Paramormyrops kingslaye + Paramormyrops_hopkinsi 

makeblastdb -in GCA_025177845.1_MSU160_NANOPORE_genomic.fna -dbtype nucl
makeblastdb -in GCA_026780785.1_ASM2678078v1_genomic.fna -dbtype nucl
makeblastdb -in GCA_027576875.1_ASM2757687v1_genomic.fna -dbtype nucl

samtools faidx Proteomes_BUSCO80/Paramormyrops_kingsleyae.fa Paramormyrops_kingsleyae---rna-XM_023831520.1 > Paramormyrops_kingsleyae---rna-XM_023831520.1.prot

tblastn -query Paramormyrops_kingsleyae---rna-XM_023831520.1.prot -db GCA_025177845.1_MSU160_NANOPORE_genomic.fna
samtools faidx GCA_025177845.1_MSU160_NANOPORE_genomic.fna JAMQYA010004842.1:203440-283440 > Paramormyrops_kingsleyae.alternative.fa
exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs"  Paramormyrops_kingsleyae---rna-XM_023831520.1.prot Paramormyrops_kingsleyae.alternative.fa > Paramormyrops_kingsleyae.alternative.exo

tblastn -query Paramormyrops_kingsleyae---rna-XM_023831520.1.prot -db GCA_026780785.1_ASM2678078v1_genomic.fna
samtools faidx GCA_026780785.1_ASM2678078v1_genomic.fna JAODKI010002750.1:1-51396 > Paramormyrops_hopkinsi.fa
exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs"  Paramormyrops_kingsleyae---rna-XM_023831520.1.prot Paramormyrops_hopkinsi.fa > Paramormyrops_hopkinsi.exo


transeq other_paramoryrops_gene.fa other_paramoryrops_gene.prot ; sed -i 's/_1$//g' other_paramoryrops_gene.prot

grep "Pangasianodon_hypophthalmus_rna_XM_053232595_1\|Ictalurus_punctatus_rna_XM_017463853_3\|Pygocentrus_nattereri_rna_XM_017714086_2\|Electrophorus_electricus_rna_XM_027006769_2" Coding_sequences_alignments/N5.HOG0010622.prot | sed 's/>//g' > curr_seq.id
xargs samtools faidx Coding_sequences_alignments/N5.HOG0010622.prot < curr_seq.id > curr_seq.prot


cat other_paramoryrops_gene.prot Paramormyrops_kingsleyae---rna-XM_023831520.1.prot curr_seq.prot > OGG_plus_ManuallyFound.prot

muscle5.1 -align OGG_plus_ManuallyFound.prot -output OGG_plus_ManuallyFound.aln


iqtree -s OGG_plus_ManuallyFound.aln -st AA -nt 8 -m TEST -mrate G4 --redo



#other_paramoryrops_gene.fa
>Paramormyrops_kingsleyae.Alternative_assembly.GCA_025177845.1
ataGAAGAACTTACCAGCACTCCCATTGATAAATGGACAGAGTCCATGGTGAGCTCATGGCTGGCTTCAA
TAGGAATTAAAGATCAATACATTGAAACACTTCATGAAAAGGAAGTTAATGGCAAAGTCCTCcttaaaat
cacagagcaaTTTCTGAAAGAGACTGAAATGAAACCTGGTCCTGTACATCTGATAATAGAGAGCAGAAAT
GAGTTAGTTAAGACTCAAAAAGCTCAGAAGCAGCAGAATAAAGCACCCAAGAGCACTGCTGAAAACACAA
ACCAAGGAGCTTCATTCATCCCAGAGCTTCCAATATCCAGCAATTCCCAGAAAGAACAACAGATGAACCA
AGACACTATAGCAAGAAGAGATTCAAAGCCACGACCTTTTGGCAAACCAGGCATTGACCACACATATGTA
AAGCATGATGTTCTTCAGCCTGAAACAGGGGTCAGTAATCTTATAACTCCATGCCATGAGTACAAGGCAG
ATACAGCTGCGACTTTAGATGCCAGAAGAATTAAAGCCAAGCTGGCATATGAAGTCCTTAAATTTGCCAC
AGGATGCATTAACATGAGAACAAATGGCACAATACACTTTGGGGTAATGGACAACCAGGAAAAGTCTGGC
TACAAACACGGTGAAATCATTGGCATTCCAATTAAGGAACAATGTGTGTATACTGATGCACTGGATTACA
TAGAAAGGTGTTTTGATAGTGACAGAGAGCTTGTACGGCAGTGTATAAGACCAGAGTTCATTCTGGTAAC
TGAGCCCAACAGGAAAGAACAACACTATGTGGTTGAATATGATGTCGAACCTTCAGTTAGCATAGTGAGA
AATGGGGTGTTTTCTGTCAGGTTGGTGAAGTTCAGTGAAGAGTCTGGAAGAATTGAACAGGAACcaaaaa
tgtattgcagAGTTGGAGCCTCAACAAAACCAGTTGAAGGTATGAGTGGATTTTTAGAAGGAGTCAATCA
CAGAGATGCTCGAAGAGAGAAGGCAGAGAAGTCCTTTCCAGAACCCTGCCAAGATCTTGAAAGGAAACTC
ATCATGCTTATAACAAGTGGgaagaaacaaatggaaaaggAGAAATGGTACATACTTGTTACAAACAAGT
TTTCAGGGAAGATCTTAAGAATGTGGATTTTTATGAACATGAAGTTATTCTGTGTGTTTGACTTTGATCC
AGATTCCAAGGTGTCTGGATTATGCCATGAATATGCTAAGCACCATGCAGTGAACCTTCACTTCATGCAA
AACTACAAAATTCCAAGTGACAAGAACATCAGAGAATTTGAGAGTCATCTACATCTGTTTGAACAAACCA
GTTGGATATTTTGTAATGGACGAAATGATTTCAGTGGGAATGAAACTTGTGATGAGAATACATGGTGTAG
AACAAGGAGGACTTTCCTTAAAGATTGTGTGTCACTGATCTGCAAGGATATCTTACCCGGAACCTTCCTT
GTCATCTTCCTCCTCATGTCTCCTGTTGAGATACCTCTACTGAAGACATTCGATGAGTTTTTCACTGACA
TGCAAGGTCATGAAGACATTATATGCATCTCAGAATCAGAAGGAAACTTTCAGAAATGGCGAGCATCTGC
TGTTGAATTCTGTGATGGGGAAACAGTGAACAATTCGAGTATAGTAGGTTTGAAAATGAGCCACGTCCAG
TGCAACTGTCCAGCAAATCCAGGCCCTAATACTCGTGTGAATAAACTCTTACCTGTGTCTGTCAAAAAAT
GCCATCTACTGACACGTGATGAGGAAACTATGTCGTCTTTGGAGGTTTTAGGTGTGAATCAGTGTGAAGA
AATCAGTGCAGAGTTCATTGagtcaaagaaagaaaaaatagacaGAGACTTTTATCGAGGAGGAAAAGTT
AAATGGATGAACTTGTGGCTTGCAGAAGTTGGAGAAGTTGGGGAAGTGATTCAAAGAGATGCCTATCATG
AGGTGATTAATCTTCTGGATGCCTCTCTTAAATTGAGCTCAGAACAAAAGCCTGTCAAATGCATCAACAT
ATACCATTACGCTGGCAGTGGAGGGAGCACAGTGGCAAGGCAGGTTCTGTGGAAGTACAGAAAGGATCTG
AGGTGTGCAGTTGTGAAACCATCATATGCTGTTGGCACAGTTTCAAAACATGCTGTGATGCTTGAGTATG
AGGAAAAGATCCAGAGAAATGTCTCTGTTCTTTTACTCATTGACGATTGTGAGAGGGAATATTTAGAAGA
GTTAAAGCATGAGTTAGAAACAGCTATCAACACAAAGAAGATTGCAAATGGAATACCGTGCTTCATTCTC
CTCAGCTGTAGAAGATGCCACAATCCAGAGAAAATGTCGAAGGAGTCGTTACTGAGTGTCTCTGTGACTC
ACAAGCTTTCACCAACCGAAAACACATTTTGCTAGAAAGCGACAAAACTCGAGAAACAGTTTCAGCCAGA
GTTCATTCTGACATTTGTCTTAATGAGTGAAGAATTTGAACATGAATATGTGGAAAAGTTTGTAGAGCAT
TTATTACAGGACATTGTTCTCACATCTGTTGACACCCAGCTGATTCAATTTGTAGCTCTGCTTAACACCT
ATGTGGAGAACTCGTACATTTCTCAGTCACATTGTGAAGCCCTCCTTCAACTCCAGCTCACTTTCTATGG
AGAAAGATTCGACAACACGGTTGAGAATTCTCTGAGTGAGCAGGCTAAATTGGtctttatacatttaaaa
gatgAAACAACCCACATCAACTCGATCAGAATCATCCACCAAGTGGTGGCAAAAGAAATTCTCCATCAGC
TTTTgggacaaaaacaacaaagtgaGCTTGCGCTTGAAATTCTTAGAAACAACGTTCTCTTTAATCACAG
ATTTGGAACGCAGTACAGGAAGTTCCTTCGTGAGTTGTTCATAAGACGCTACAAAATCAGCAGAGGTGAC
AAATCAGACACTTTTTTCTCACCCATTGAACATGTGAGAGAAAAAGAAACGGCAGAAAATGCTATTGAGC
TTCTTCATGAGGCATACAAAAGATTCAGCGAGGATGCATTTGCTCAACAGCTAGCTCGCCTCTATTACCG
ACATGAAAGATTTGACCTTGCAGAACGTTGGGCAGAAACTGCAGCAGCTAAACTGCCTAAGAACTCCTAC
ATTCTTGACACAAAAGGTCAGGtgtacagaaaatggttcacaacaaaacagaatagaAAATGTCAGAGAA
CACCAGAGTCCATAGCAGATGCCATTGAGACAGCTCTTAAAGCCATTGAATGTTTCAAATGTCAGTTGGT
TGCTGTTGAAGAAACTGAGACCATGAACAACTCGGGATTTGTTGGAGCTATAGAAGTTGGGTGTAACTTG
TTGAGCTAATTTCTTCACTTGATGTTCTCAAATAAATATGGGGATCATTCTGAATTACAGAAGTACTTGC
TTACAGAATACATTCCCAAAGAAGTTGAAGCACCATGGGAACACTTCCACTACAAACTCAAAAGTCTCCA
GCTCACAATGCACAAGGCATTCGAGTGGATGTCAGAAGAACTGAGTTATTTCCAAAACCAAGACACTGCT
GAGGAGGAAACCTCCAAAACCTCAGAGCTGACGATACATCGCCCCAAACACTGGCTGGCTGGCAAGTCCT
GTGTTTATGGAAAGTTCTTCTGTGAGGTCTCACTCAGCAGCACACCCTGTAACTGGCAAtctcatttaaa
taacatgAGCGACTTCAGCAAACGCATGGCTATCTACCAGCTCGGGAGGAAACATTACAACAATCTTTTT
CCATCCTGGACCAGAAAAAGAGAACAAGACAAAACTTTGGAGAATATAATTTCACTGTATCCAAAAGGTG
CTAAAATGGATCAACTGGAATTTGGCAACTACACAGCATCACATTTTGCCTTAAGTGCAATCTCTCCAGG
TTCACGTTACCTGTCTGTCCTCAAACATCTGCAGAAGCTGAGCAGACCGTTTCTTCAAGACAAATCAAAA
TGTCCATCAAGTGTTCTGTTTCTGTGCACACTACTGTTTTGGCCAGAAAATTTGATGTTGTCAAGAAAAG
AggagaaatataaaacaatccTTACTGCTGTCAAATTCCTCCAACAAACATACAAGGCCAAAATGAAGGA
TATCCCTGTCAGAAGGAGGCGGATTTACACCCACTTCTACCTGGGAAAAGGATTTGGGTATGAGAAATTT
GTCCACAAGAATAAAATTGAAACGATCAGAAAATTTTCCTCTGTTTCAGAAAAGCGTCAAAAATGGCTTG
AAGGGGGAGTGTGGAAAACACCGGAAATTGCTGGAGAGCTCATGCGTGTAAGCGGCTGGACAGAAGATGG
AAATGTGTATCTTGAAGGTCCTAAAACTGAAAAGTTCATTCATCCCCTTAATGAAAGTTCCGTGCCTGCT
GGAAATGAAAACGTCACCTTCTACGTAGGCTTCACGTTCGGACCTGTTGCATGTGACATCACAATAAGAA
AA
>Paramormyrops_hopkinsi
ATAGAAGAACTTACCAGCACTCCCATTGATAAATGGACAGAGTCCATGGTGAGCTCATGGCTGGCTTCAA
TAGGAATTAAAGATCAATACATTGAAACACTTCATGAAAAGGAAGTTAATGGCAAAGTCCTTCTTAAAAT
CACAGAGCAATTTCTGAAAAAAGAGACTGGAATGAAACCTGGTCCTGTACATCTGATAATAGAGAGCAGA
AATGACTTAGTTAAGACTCAAAAAACTCAGAAGCAGCAGAATAAAGCACCCAAAAGCACTGCTGAAAACA
CAAACCAAGGAGCTTCATTCATCCCAGAGCTTCCAATATCCAGCAATTCCCAGAAAGAAGAACAGATGAA
CCAAGACACTATAGCAAGAAAAAGAGATTCAAAGCCACGACCTTTTGGCAAACCAGGCATTGACCACACA
TATGTAAAGCATGATGTTCTTCAGCCTGAAACAGGGGTCAGTAATCTTATAACTCCATGCCATGAGTATA
AGGCATTCAATACAGCTGCAACATTAGATGCCAGAAGAATTAAAGCCAAGCTGGCATATGAAGTCCTTAA
ATTTGCCACAGGATGCATTAACATGAGAACAAATGGCACAATACACTTTGGGGTAATGGACAACCAGGAG
AAAAAATCTGGCTACAAACACGGTGAAATCATTGGCATTCCCATTAAGGAACAATGTGTGTATACTGATG
CACTGGATTACATAGAAAGGTGTTTTGATAGTGATAGAGAGCTTGTACGGCAGTGTATAAGACCCCCAGA
GTTCATTCTGGTAACTGAGCCCAACAGGAAAGAACAACACTATGTGGTTGAATATGATGTCGAACCTTCA
GTTAGCATAGTGAGAAATAGGGTGTTTTCTGTCAGCTTGGTGAAGTTCAGTGAAGAGTCTGGAAGAATTG
AACAGGAACcagaaaaaatgtattgcagaGTTGGAGCCTCAACAAAACCAGTTGAAGGTATAAGTGGATT
TTTAGAAGGAGTCAATCACAGAGATGCTCGAAGAGAGAAGGCAGAGAAGTCCTTTTTGCCAGAACCCTGC
CAAGATCTTGAAAGAAAACTCATCATGCTTATTACAAGTGGGaagaaacaaatggaaaaggaGAAATGGT
ACATACTTGTTACAAACAAGTTTTCAGAGGAAGATCTTAAGAATGTggattttttaatgaacatgaaGTT
ATTCTGTGTGTTTGACTTTGATCCAGATTCCAAGGTGTCTGGATTATGCCATGAATATGCTAAGCACCAT
GCAGTGAACCTTCACTTCATGCAAAACTACAAAATTCCAAGTGACAAGAACAGCAGAGAATTTGAGCGTC
ATCTACATCTGTTTGAACAAACCAGTTGGATATTTTGTAATGGACGAAATGATTTCAGTGGGAATGAAAC
TCCCTGTGATGAGAATACATGGTGTAGAACAAGGAGGACTTTCCTTAAGGATTGTGTGTCACTGATCTGC
AAGGATATCTTACCCAAGGGAACCTTCCTTgtcatcttcctcctcatgTCTCCTGTTGAGATACCTCTAC
TGAAGACATTTGATGAGTTTTTCACTGACATGCAAGGTCATGAAGACATTATATGCATCTCAGAATCAGA
AGGAAACTTTCAGAAATGGCGAGCATCTGCTGTTGAATTCTGTGATGGGGAAACAGTGAACAATTCGAGT
ATAGTAGGTTTGAAAATGAGCCACGTCAGTGTAACTGTCCAGCAAATCCAGGGCCCTAATACTCGTGTGA
ATAAACTCTTACCTGTGTCTGTCAAAGCAAAATGCCATCTACTGACACATGATGAGGAAACAATGTCGTC
TTTGGAGGTTTTAGGTGTGAATCAGTGTGAAGAAATCAGTGCAGAGTTCATTGagtcaaagaaagaaaaa
atagacaGAGACTTTTATCGAGGAGGAAAAGTTACATGGATGAACTTGTGGCTTGCAGAAGTTGGAGAAG
TTGGGGAAGTGATTCAAAGAGATGCCTATCATGAGGTGATTAATCTTCTGGATGCCTCTCTTAAATTGAG
CTCAGAACAAAAGCCTGTCAAATGCATCAACATATACCATTACGCTGGCAGTGGAGGGAGCACAGTGGCA
AGGCAGGTTCTGTGGAAGTACAGAAAGGATCTGAGGTGTGCAGTTGTGAAACCATCATATGCTGTTGGCA
CAGTTTCAAAACATGCTGTGATGCTTCGGGAGTATGAGGAAAAAGATCCAGAGAAATGTCTCCCTGTTCT
TTTACTCATTGACGATTGTGAGAGAGAATATTTAGAAGAGTTAAAGCATGAGTTAGAAACAGCTATCAAC
ACAAAGAAGATTGCAAATGGAACACCGTGCTTCATTCTCCTCAGCTGTAGAAGATGCCACAATCCAGAGA
AAATGTTGAAGGAGTCGCCTTTACTGAGTGTCTCTGTGACTCACAAGCTTTCACCAGCcgagaaaacaca
ttttgctaGAAAGCGACAAAACCTCGAGAAACAGTTTCAGCCAGAGTTCATTCTGACATTTGTCTTAATG
AGTGAAGAATTTGAACATGAATATGTGGAAAAGTTTGTAGAACATTTATTACAGAACATTGTTCTCACAT
CTGATGTCATCCAGCTGATTCGATTTGTAGCTCTGCTTAACACCTATGTGGAGAACTCTTACATTTCTCA
GTCACATTGTGAAGCCCTCCTTCAACTTCAGCTCACTATCTATGGAGAAAGATTCAGACAACACATGTTC
AAGAATTCTCTGAGTGAGCAGGCTAACTTGGtctttatacatttaaaagatGAAACAACCCACATCAACT
CGATCAGAATCATCCACCAAGTGGTGGCAAAAGAAATTCTCCATCAGCTTTTgggacaaaaacaacaaag
tgaGCTTGCGCTTGAAATACTTGGAAACAACATTCTCTTTAATCACAGGTTTGGGGGGACGCAGTACAGG
AAGTTCCTTCGTGAGTTGTTCATAAGACGCTACAAAATCAGCAGAGGTGACAAATCAGACACTTTTTTCT
CACCCCTCATTGAACAtgtgagagaaaaagaaaagccagaAAATGCTATTGAGCTTCTTCATGAGGCATA
CAAAAGATTCAGCGAGGATGCATTTTTTGCTCAACAGCTAGCTCGCCTCTATTACCGACATGAAAGATTT
CACCTTGCAGAACGTTGGGCAGAAACTGCAGCAGCTAAACTGCCTAAGAACTCCTACATTCTTGACACCA
AAGGTCAGGTGTACAGAAATTGGttcacaacaaaaaacagtttaataGACAAACTTCAGAAAACACCAGA
GTCCATAGCAGATGCCATTGAGACAGCTCTTAAAGCCATTGAATGTTTCCAAAAATGTCAGTTGGTTGCT
ATTGAAGAAACTGAGACCATGAACAACTCGGGATTCGTTGGAGCTATAGAAGTTGGGTGTAAGTTGTTGC
AGCTAATTTCTTCACTTGATGTTTTCTCAAATAAATATGGGGATCATTCTGAATTACAGAAGTACTTACT
TACAGAATACATTCCCAAAGAAGTTGAAGCACCATGGGAACACTTCCACTACAAACTCAAAAGTCTCCAG
CTCACAATGCACAAGGCATTCGAGTGGATGTCAGAAGAACTGAGTTATTTCCAAAACCAAGACACTGCTG
AGGAGGAAACCTCCAAAACCTCAGAGCTGACGATACATCGCCCCAAACACTGGCTGGCTGCCAAGTCCTG
TGTTTATGGAAAGTTCTTCTGTGAGGTCTCACTCAGCAGCACACCTTGTAACTGGCAatctcatttaaat
aacatgAGCGACTTCAGCAAACGCATGGCTATCTACCAGCTCGGGGGAGGAAACATTACAACAATCTTTT
CCATCCTGgaccagaaaaaaagagaacaagaCAAAACTTTGGAGAATATAATTTCACTGTATCCAAAAGG
TGCTAAAATGGATCAACTGGAATTTGGCAACTACACAGCATCACATTTTGCCTTAAGTGCAATCTCTCCA
GGTTCACGTTACCTGTCTGTCCTCAAACATCTGCAGAAGCTGAGCAGACCGTTTCTTCAAGACAAATCAA
AATGTCCATCAAGTGTTCTGTTTCTGTGCACACTACTGTTTTGGCCAGAAAAATTTGATGTTGATcaaga
aaaagaggagaaatataaaacaatccTTACTGCTGTCAAATTCCTCCAACAAACATACAAGGCCAAAATG
AAGGATATCCCTGTCAGGAGGAGGCGGATTTACACCCACTTCTACCTGGGAAAAGGATCTGGGTACGAGA
AATTTGTCCACAAGAATAAAATTGAAACGATCAGAAAATTTTCCTCTGTTTCAGAAAAGCGTCAAAAATG
GCTTGAAGGGGGAGTGTGGAAAACACCGGAAATTGCTGGAGAGCTCATGCGTGTAAACGGCTGGACAGAA
GATGGAAATGTGTATCTTGAAGGTCCTAAAACTGAAAAGTTTTTCATTCATCCCCTTAATGAAAGTTCCG
TGCCTGCTGGAAATGAAAACGTCACCTTCTACGTAGGCTTCACGTTCAGGGGACCTGTTGCATGTGACAT
CACAATAAGAAAA


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



#N5.HOG0010622.clade_hgt.txt=> clade where the HGT is extracted from iTOL (from N5.HOG0010622.prot.aln.treefile)

cut -f1,2 -d "_" N5.HOG0010622.clade_hgt.txt | sort | uniq  > species_to_draw.clade1
#Choose close species without the gene .. 


echo "Scleropages_formosus" > non_transfer_species_1
echo "Brienomyrus_brachyistius" >> non_transfer_species_1

echo "Pangasianodon_hypophthalmus" > non_transfer_species_2
echo "Electrophorus_electricus" >> non_transfer_species_2
echo "Pygocentrus_nattereri" >> non_transfer_species_2



rm -r Syn_tables_dir ; mkdir Syn_tables_dir


#First order the species names 

rm species_to_draw.clade1.ordered

echo "Paramormyrops_kingsleyae" >> species_to_draw.clade1.ordered
echo "Ictalurus_punctatus" >> species_to_draw.clade1.ordered



#If there are more than 4 species from the same order, keep only 4 


#First extract genes around 

curr_OGG=N5.HOG0010622

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


grep "NW_019713398" Syn_tables_dir/Paramormyrops_kingsleyae.synt.final.df  > temp ; mv temp Syn_tables_dir/Paramormyrops_kingsleyae.synt.final.df
grep "NW_026521115\|NC_030417" Syn_tables_dir/Ictalurus_punctatus.synt.final.df  > temp ; mv temp Syn_tables_dir/Ictalurus_punctatus.synt.final.df

#First add Brienomyrus_brachyistius

curr_OGG=N5.HOG0010622
curr_sp=Brienomyrus_brachyistius
ref_sp=Paramormyrops_kingsleyae


grep -A10 -B10 "N5_HOG0016328" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret  > Syn_tables_dir/$curr_sp.synt.df

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

curr_OGG=N5.HOG0010622
curr_sp=Scleropages_formosus
ref_sp=Paramormyrops_kingsleyae


grep -A10 -B10 "N5_HOG0016328" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret  > Syn_tables_dir/$curr_sp.synt.df



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



#Now add Pangasianodon_hypophthalmus

curr_OGG=N5.HOG0010622
curr_sp=Pangasianodon_hypophthalmus
ref_sp=Ictalurus_punctatus


grep -A10 -B10 "N5_HOG0017170" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret  > Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0018014" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret | grep "NC_069712" >> Syn_tables_dir/$curr_sp.synt.df



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




#Now add Electrophorus_electricus

curr_OGG=N5.HOG0010622
curr_sp=Electrophorus_electricus
ref_sp=Ictalurus_punctatus


grep -A10 -B10 "N5_HOG0010622" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret  > Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0045507" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret  >> Syn_tables_dir/$curr_sp.synt.df



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




#Now add Pygocentrus_nattereri

curr_OGG=N5.HOG0010622
curr_sp=Pygocentrus_nattereri
ref_sp=Ictalurus_punctatus


grep -A10 -B10 "N5_HOG0010622" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret  > Syn_tables_dir/$curr_sp.synt.df
grep -A10 -B10 "N5_HOG0045507" GFF3_N5_OGGs/$curr_sp.gff.simplified.sorted.OGG.tiret  >> Syn_tables_dir/$curr_sp.synt.df



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
##========================================== Verify that the gene is absent in close species ==========================================
##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================

samtools faidx GCF_023856365.1_BBRACH_0.4_genomic.fna NC_064539.1:7524103-8209340 > Brienomyrus_brachyistius.region.fa
samtools faidx GCF_900964775.1_fSclFor1.1_genomic.fna NC_041811.1:8506043-9224772 > Scleropages_formosus.region.fa

exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs" Paramormyrops_kingsleyae---rna-XM_023831520.1.prot Brienomyrus_brachyistius.region.fa
exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs" Paramormyrops_kingsleyae---rna-XM_023831520.1.prot Scleropages_formosus.region.fa



exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs" Paramormyrops_kingsleyae---rna-XM_023831520.1.prot GCF_023856365.1_BBRACH_0.4_genomic.fna > Brienomyrus_brachyistius.whole.fa
exonerate  --showtargetgff TRUE --model protein2genome --ryo "%tcs" Paramormyrops_kingsleyae---rna-XM_023831520.1.prot GCF_900964775.1_fSclFor1.1_genomic.fna > Scleropages_formosus.whole.fa



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


#Extract the region between N5_HOG0010622 and N5_HOG0046344 (genes included) + 5kb upstream / downstream

#Paramormyrops
grep -A3 -B3 "N5_HOG0010622"  GFF3_N5_OGGs/Paramormyrops_kingsleyae.gff.simplified.sorted.OGG.tiret
samtools faidx GCF_002872115.1_PKINGS_0.1_genomic.fna NW_019713398.1:637791-652729 > N5.HOG0010622.Paramormyrops_kingsleyae.extended.1.fa
sed -i 's/:/-/g' N5.HOG0010622.Paramormyrops_kingsleyae.extended.1.fa
makeblastdb -in N5.HOG0010622.Paramormyrops_kingsleyae.extended.1.fa -dbtype nucl

#Ictalurus
grep -A3 -B3 "N5_HOG0010622"  GFF3_N5_OGGs/Ictalurus_punctatus.gff.simplified.sorted.OGG.tiret
samtools faidx GCF_001660625.3_Coco_2.0_genomic.fna NC_030417.2:23089172-23102712 > N5.HOG0010622.Ictalurus_punctatus.extended.1.fa
samtools faidx GCF_001660625.3_Coco_2.0_genomic.fna NW_026521115.1:118979-133731 > N5.HOG0010622.Ictalurus_punctatus.extended.2.fa
sed -i 's/:/-/g' N5.HOG0010622.Ictalurus_punctatus.extended.1.fa
makeblastdb -in N5.HOG0010622.Ictalurus_punctatus.extended.1.fa -dbtype nucl
sed -i 's/:/-/g' N5.HOG0010622.Ictalurus_punctatus.extended.2.fa
makeblastdb -in N5.HOG0010622.Ictalurus_punctatus.extended.2.fa -dbtype nucl



#Now lets find TEs using the repbase and dfam databases

tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0010622.Paramormyrops_kingsleyae.extended.1.fa -evalue 1e-5 -outfmt 6 -out TE.Paramormyrops_kingsleyae.1.tblastn -num_threads 8
sed -i 's/#//g' TE.Paramormyrops_kingsleyae.1.tblastn

tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0010622.Ictalurus_punctatus.extended.1.fa -evalue 1e-5 -outfmt 6 -out TE.Ictalurus_punctatus.1.tblastn -num_threads 8
sed -i 's/#//g' TE.Ictalurus_punctatus.1.tblastn

tblastn -query Dfam_plus_Repbase.cdhit80.prot -db N5.HOG0010622.Ictalurus_punctatus.extended.2.fa -evalue 1e-5 -outfmt 6 -out TE.Ictalurus_punctatus.2.tblastn -num_threads 8
sed -i 's/#//g' TE.Ictalurus_punctatus.2.tblastn


#merge tblastn hits and find the best TE match by doing a blastx

cp ../Rscript_merge_blast_hits.R ./

Rscript Rscript_merge_blast_hits.R TE.Paramormyrops_kingsleyae.1.tblastn TE.Paramormyrops_kingsleyae.1.tblastn.merged
Rscript Rscript_merge_blast_hits.R TE.Ictalurus_punctatus.1.tblastn TE.Ictalurus_punctatus.1.tblastn.merged
Rscript Rscript_merge_blast_hits.R TE.Ictalurus_punctatus.2.tblastn TE.Ictalurus_punctatus.2.tblastn.merged


xargs samtools faidx N5.HOG0010622.Paramormyrops_kingsleyae.extended.1.fa < TE.Paramormyrops_kingsleyae.1.tblastn.merged > TE.Paramormyrops_kingsleyae.1.BEST.fa
blastx -query TE.Paramormyrops_kingsleyae.1.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Paramormyrops_kingsleyae.1.BEST.blastx -max_target_seqs 1
cut -f1 TE.Paramormyrops_kingsleyae.1.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Paramormyrops_kingsleyae.1.BEST.blastx >> temp  ; done ; mv temp TE.Paramormyrops_kingsleyae.1.BEST.blastx


xargs samtools faidx N5.HOG0010622.Ictalurus_punctatus.extended.1.fa < TE.Ictalurus_punctatus.1.tblastn.merged > TE.Ictalurus_punctatus.1.BEST.fa
blastx -query TE.Ictalurus_punctatus.1.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Ictalurus_punctatus.1.BEST.blastx -max_target_seqs 1
cut -f1  TE.Ictalurus_punctatus.1.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Ictalurus_punctatus.1.BEST.blastx >> temp  ; done ; mv temp TE.Ictalurus_punctatus.1.BEST.blastx


xargs samtools faidx  N5.HOG0010622.Ictalurus_punctatus.extended.2.fa < TE.Ictalurus_punctatus.2.tblastn.merged > TE.Ictalurus_punctatus.2.BEST.fa
blastx -query TE.Ictalurus_punctatus.2.BEST.fa -db Dfam_plus_Repbase.cdhit80.prot -outfmt 6 -num_threads 8 -out TE.Ictalurus_punctatus.2.BEST.blastx -max_target_seqs 1
cut -f1  TE.Ictalurus_punctatus.2.BEST.blastx  | sort | uniq > uniq_regions.txt ; rm temp ; for region in `cat uniq_regions.txt` ; do grep -m1 "$region" TE.Ictalurus_punctatus.2.BEST.blastx >> temp  ; done ; mv temp TE.Ictalurus_punctatus.2.BEST.blastx



#Now find shared elements

cut -f2 TE.Paramormyrops_kingsleyae.1.BEST.blastx | sort | uniq > TE.Paramormyrops_kingsleyae.1.uniqTE
cut -f2 TE.Ictalurus_punctatus.1.BEST.blastx | sort | uniq > TE.Ictalurus_punctatus.1.uniqTE
cut -f2 TE.Ictalurus_punctatus.2.BEST.blastx | sort | uniq > TE.Ictalurus_punctatus.2.uniqTE

comm -12 TE.Paramormyrops_kingsleyae.1.uniqTE TE.Ictalurus_punctatus.1.uniqTE
comm -12 TE.Paramormyrops_kingsleyae.1.uniqTE TE.Ictalurus_punctatus.2.uniqTE




#NO SHARED TEs


### In a FINAL step we will make another synteny plot, only between Paramormyrops and Ictalurus to show a zoom on the region


## Paramormyrops

rm seq_clustered_infos_ogg.TE.txt ; rm clusters_ID_TE.txt


scaffold=`grep ">" N5.HOG0010622.Paramormyrops_kingsleyae.extended.1.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0010622.Paramormyrops_kingsleyae.extended.1.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0010622.Paramormyrops_kingsleyae.extended.1.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Paramormyrops_kingsleyae,$scaffold,$length" > clusters_ID_TE.txt

grep "N5_HOG0010622"  GFF3_N5_OGGs/Paramormyrops_kingsleyae.gff.simplified.sorted.OGG.tiret

grep "XM_023831520" GFF3_files_per_species/Paramormyrops_kingsleyae.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons


nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene1.$nbexon,N5_HOG0010622,-,Paramormyrops_kingsleyae" >> seq_clustered_infos_ogg.TE.txt
done


cut -f1 TE.Paramormyrops_kingsleyae.1.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Paramormyrops_kingsleyae.1.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0010622.Paramormyrops_kingsleyae.extended.1.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0010622.Paramormyrops_kingsleyae.extended.1.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold,/g" | sed "s/$/,$strand,Paramormyrops_kingsleyae/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done




## Ictalurus punctatus


scaffold=`grep ">" N5.HOG0010622.Ictalurus_punctatus.extended.1.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0010622.Ictalurus_punctatus.extended.1.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0010622.Ictalurus_punctatus.extended.1.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Ictalurus_punctatus,$scaffold,$length" >> clusters_ID_TE.txt

grep -A3 -B3 "N5_HOG0010622"  GFF3_N5_OGGs/Ictalurus_punctatus.gff.simplified.sorted.OGG.tiret

grep "XM_017491464" GFF3_files_per_species/Ictalurus_punctatus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons


nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene1.$nbexon,N5_HOG0010622,-,Ictalurus_punctatus" >> seq_clustered_infos_ogg.TE.txt
done

cut -f1 TE.Ictalurus_punctatus.1.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Ictalurus_punctatus.1.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0010622.Ictalurus_punctatus.extended.1.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0010622.Ictalurus_punctatus.extended.1.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold,/g" | sed "s/$/,$strand,Ictalurus_punctatus/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done



#Second region of Ictalurus


scaffold=`grep ">" N5.HOG0010622.Ictalurus_punctatus.extended.2.fa | sed 's/>//g' | cut -f1 -d "-"`
start_region=`grep ">" N5.HOG0010622.Ictalurus_punctatus.extended.2.fa | sed 's/>//g' | cut -f2 -d "-"`
stop_region=`grep ">" N5.HOG0010622.Ictalurus_punctatus.extended.2.fa | sed 's/>//g' | cut -f3 -d "-"`
length=$(( stop_region -  start_region + 100))

echo "Ictalurus_punctatus,$scaffold,$length" >> clusters_ID_TE.txt

grep "N5_HOG0010622"  GFF3_N5_OGGs/Ictalurus_punctatus.gff.simplified.sorted.OGG.tiret
grep "XM_017463853" GFF3_files_per_species/Ictalurus_punctatus.gff | grep "CDS" | cut -f4,5 | tr "\t" "," > gene_1.exons



nbexon=0
for line in `cat gene_1.exons` ; do 
	nbexon=$(( nbexon +1 ))
	exon_start=`echo "$line" | cut -f1 -d ","` ; exon_stop=`echo "$line" | cut -f2 -d ","`
	real_exon_start=$(( exon_start - start_region )) ; real_exon_stop=$(( exon_stop - start_region ))
	echo "$scaffold,$real_exon_start,$real_exon_stop,gene2.$nbexon,N5_HOG0010622,-,Ictalurus_punctatus" >> seq_clustered_infos_ogg.TE.txt
done


cut -f1 TE.Ictalurus_punctatus.2.BEST.blastx | sed 's/.*://g' | sed 's/-/,/g' > locations_TE.txt
cut -f2 TE.Ictalurus_punctatus.2.BEST.blastx > names_TE.txt
paste -d "," locations_TE.txt names_TE.txt names_TE.txt > locations_names.false.txt

header=`grep ">" N5.HOG0010622.Ictalurus_punctatus.extended.2.fa | sed 's/>//g'`
for curr_TE_line in `cat locations_names.false.txt` ;  do 
	start=`echo "$curr_TE_line" | cut -f1 -d ","`
	stop=`echo "$curr_TE_line" | cut -f2 -d ","`
	TE_name=`echo "$curr_TE_line" | cut -f3 -d ","`
	samtools faidx Dfam_plus_Repbase.cdhit80.prot $TE_name > curr_TE.fa
	samtools faidx N5.HOG0010622.Ictalurus_punctatus.extended.2.fa $header:$start-$stop > curr_region.fa
	makeblastdb -in curr_region.fa -dbtype nucl
	tblastn -query curr_TE.fa -db curr_region.fa -outfmt "6 sstart send" -max_target_seqs 1 -out curr_blast.txt
	TE_start=`head -1 curr_blast.txt | cut -f1`
	TE_stop=`head -1 curr_blast.txt | cut -f2`

	if [ $TE_start -ge $TE_stop ] ; then strand="-" ; else strand="+" ; fi

	echo "$curr_TE_line" | sed "s/^/$scaffold,/g" | sed "s/$/,$strand,Ictalurus_punctatus/g" | sed 's/#//g' >> seq_clustered_infos_ogg.TE.txt

done





## Look at an alignment of gene introns

header=`grep ">" N5.HOG0010622.Paramormyrops_kingsleyae.extended.1.fa | sed 's/>//g'`
samtools faidx N5.HOG0010622.Paramormyrops_kingsleyae.extended.1.fa $header:9782-9890 > Paramormyrops_kingsleyae.intron.fa

header=`grep ">" N5.HOG0010622.Ictalurus_punctatus.extended.1.fa | sed 's/>//g'`
samtools faidx N5.HOG0010622.Ictalurus_punctatus.extended.1.fa $header:7472-7561 > Ictalurus_punctatus.intron.1.fa

header=`grep ">" N5.HOG0010622.Ictalurus_punctatus.extended.2.fa | sed 's/>//g'`
samtools faidx N5.HOG0010622.Ictalurus_punctatus.extended.2.fa $header:9869-11354 > Ictalurus_punctatus.intron.2.fa



cat Paramormyrops_kingsleyae.intron.fa  Ictalurus_punctatus.intron.1.fa Ictalurus_punctatus.intron.2.fa > introns.combined.fa

muscle5.1 -align  introns.combined.fa -output introns.combined.aln




##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##========================================== Selective pressure analysis ============================================================
##====================================================================================================================================
##====================================================================================================================================
#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================


#Compute the dN/dS on every branches
sbatch --qos=1day -c 4 --mem=10G launch_fitMG4.sh N5.HOG0010622

Recipient branches : 
Paramormyrops_kingsleyae_rna_XM_023831520_1



#Test positive selection and relaxed selection on receiver branch
sbatch --qos=1week -c 4 --mem=10G -e error.absrel.cand.out -o slurm.absrel.cand.out --job-name=HOG0010622 launch_absrel_cand.sh N5.HOG0010622


#Now lets launch RELAX on the same branch
sbatch --qos=1day -c 4 --mem=10G -e error.relax.out -o slurm.relax.out --job-name=HOG0010622 launch_RELAX.sh N5.HOG0010622



#Extract dN/dS to table

HOG=N5.HOG0010622

grep "LB\":" $HOG.cds.aln.FITTER.json | sed 's/.*://g' | sed 's/,.*//g' > curr_LB_values.txt
grep "MLE\":" $HOG.cds.aln.FITTER.json | sed 's/.*://g' | sed 's/,.*//g' > curr_MLE_values.txt
grep "UB\":" $HOG.cds.aln.FITTER.json | sed 's/.*://g' | sed 's/,.*//g' > curr_UB_values.txt
grep "\"dN\"" $HOG.cds.aln.FITTER.json  | grep -v "{" | sed 's/.*://g' | sed 's/,//g' > curr_dN_values.txt
grep "\"dS\"" $HOG.cds.aln.FITTER.json  | grep -v "{" | sed 's/.*://g' | sed 's/,//g' > curr_dS_values.txt
grep -B2 "LB\":" $HOG.cds.aln.FITTER.json | grep -v "\-\-" | grep -v "Confidence" | grep -v "LB\":"  | sed 's/\"//g' | sed 's/:.*//g' | sed 's/^ *//g' > curr_labels

paste -d "," curr_labels curr_LB_values.txt curr_MLE_values.txt curr_UB_values.txt curr_dN_values.txt curr_dS_values.txt > $HOG.dN_dS.csv












