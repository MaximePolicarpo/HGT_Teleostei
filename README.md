# HGT_Teleostei

All the files necessary to run the scripts below, as well as intermediate and resulting files, can be found in FigShare : XXXX

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

For each HOG with a HGT manually retained, two files named "HOG.sh" and "HOG.R" (example: [N5.HOG0001647.sh](N5.HOG0001647.sh) and [N5.HOG0001647.R](N5.HOG0001647.R)) can be used to further analyse the transfer, its region and neighboring TEs. 




