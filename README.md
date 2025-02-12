# HGT_Teleostei

All the files necessary to run the scripts below, as well as intermediate and resulting files, can be found in Zenodo : 10.5281/zenodo.14011660 + 10.5281/zenodo.14012814 + 10.5281/zenodo.14012963

## I - Ks pipeline to identify candidate HGTs

[Ks_pipeline.sh](Ks_pipeline.sh) 

This scripts allows to :
- Compute alignments for each orthogroup (HOG) and for each BUSCO group.
- Compute Ks statistics between all pairs of genes among each HOG.
- Extract Ks statistics (quantiles, mean, median) for each pair of species.
- Compute micro-syntenic scores between all pairs of genes among each HOG.
- Identify genes matching to transposable elements.

Accessory scripts : [Extract_candidate_HGT.R](Extract_candidate_HGT.R) ; [Extract_ks_stats_OGG.R](Extract_ks_stats_OGG.R) ; [Compute_aligned_CDS_stats_FAST.R](Compute_aligned_CDS_stats_FAST.R) ; [measure_Ks_stats_HGT.sh](measure_Ks_stats_HGT.sh) ; [Compute_synteny_score.sh](Compute_synteny_score.sh)


## II - BLAST pipeline to identify candidate HGTs

[Blast_pipeline.sh](Blast_pipeline.sh) 

This scripts allows to :
- Run an All-vs-all blastn
- Extract BLAST outliers with three different strategies (A: different order, B: high phylogenetic distance, C: high HGTindex)

Accessory scripts : [run_blast_proc.sh](run_blast_proc.sh) ; [Parse_blast_results.R](Parse_blast_results.R)


## III - Filter and analyze results of the Ks and BLAST pipelines

[HGT_analysis.R](HGT_analysis.R) 

This scripts allows to :
- Compute some statistics about Ks values (correlation between Ks computed on BUSCO and HOG, correlation between Ks and divergence time)
- Filter Ks pipeline results to obtain a final list of candidate HGTs and associated HOGs
- Filter BLAST pipeline results to obtain a final list of candidate HGTs and associated HOGs
- Compare candidate HGTs and HOGs obtained with the Ks and BLAST pipelines
- Draw a species tree with links between manually curated HGTs
- Draw a circos plot between H. transpacificus and C. harengus
- Compute some statistics about micro-syntenic scores and divergence times
- Analyze the results of HyPhy (fitMG94, aBSRREL and RELAX)
- Compute ENCprime values for all teleost species
- Compute RSCU values for a subset of species
- Compute GC3 distributions for a subset of species


## IV - Compute maximum likelihood phylogenies

[Compute_ML_tree.sh](Compute_ML_tree.sh)  : allows to compute ML trees for each HOG obtained through the Ks and Blast pipelines. 


## V - Further analysis of HGTs

For each HOG with a HGT manually retained, two files named "HOG.sh" and "HOG.R" (example: [N5.HOG0001647.sh](N5.HOG0001647.sh) and [N5.HOG0001647.R](N5.HOG0001647.R)) can be used to further analyse the transfer, its region and neighboring TEs. It also allows to perform another ML tree, based on BLASTP best match (using the transfered gene as query) found in the uniprot database + all teleost proteomes.

[Rscript_merge_blast_hits.R](Rscript_merge_blast_hits.R) : Rscript to parse BLAST matches of TE against genomes/genomic regions
[launch_blast_allgenomes.sh](launch_blast_allgenomes.sh) : Script to launch tblastn of a sequence against all ray-finned fishes genomes 


These script need a file called "non_actino_uniprot.fa" which correspond to the Uniprot protein database, with no representative fish sequences
This can easily be generated using the taxid of ray-finned fishes https://www.ncbi.nlm.nih.gov/taxonomy/?term=txid7898[Subtree] and removing related sequences
from uniprot (https://www.uniprot.org/uniprotkb?query=reviewed:true)

Commands to run HyPhy can be found here ( [launch_absrel_cand.sh](launch_absrel_cand.sh) ; [launch_RELAX.sh](launch_RELAX.sh) ; [launch_fitMG4.sh](launch_fitMG4.sh) ). Selection of branches to test in aBSREL and RELAX were performed on https://phylotree.hyphy.org/.
FitMG94.bf can be found here : https://github.com/veg/hyphy-analyses/tree/master/FitMG94

Genome fasta files, to be found on RefSeq and GenBank, are also needed for these analysis.  

Alignments of the genes in HGT clades (in DNA), and associated maximum likelihood phylogenies can be found in the repository "HGT_clade_DNA_phylogenies".


All the other files needed to run these scripts are available on Zenodo.

## VI - Mitochondrial genome analysis (check for contamination)

[Mitochondrial_genome_analysis.sh](Mitochondrial_genome_analysis.sh) : Script to make a maximum likelihood phylogeny of mitochondrial genes, to map raw reads of 11 species to their own mitochondrial genome and to the mitochondrial genome of C. harengus, and to compute the number of mapped reads and the mean read depth. 

Associated files are in the folder "Mitochondrial_analysis"

## VII - Long reads analysis

[Long_reads_analysis.sh](Long_reads_analysis.sh) : Script to (i) map long reads to the genome of C. harnegus and H. transpacificus, and to compute the number of long reads entirely spanning transferred genes and (ii) extract long reads spanning the genes and predict exons locations on those reads using EXONERATE (associated files in the folder "Long_reads_files"). 


## VIII - Go term analysis

[GOterm_analysis.sh](GOterm_analysis.sh) : Script to extract the longest sequence per HOG and BLAST it against the UniProt database
[GOterm_analysis.R](GOterm_analysis.R) : Script to extract GO terms from the UniProt database and run the over-representation enrichment anlysis

Associated files are in the folder "GO_term_analysis"

## IX - Map the geographical distributions of species involved in HGTs

[Geographical_distributions.R](Geographical_distributions.R) 




