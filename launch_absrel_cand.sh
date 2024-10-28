#!/bin/bash


#SBATCH --job-name=abSREL_candidate   # Job name



#conda = miniprot ; cd-hit ; seqkit

eval "$(conda shell.bash hook)"
conda activate HyPhy_Env

HOG=$1

hyphy absrel --alignment $HOG.cds.aln --tree $HOG.prot.aln.treefile.SelecMarked --code Universal --branches test ENV=TOLERATE_NUMERICAL_ERRORS=1

mv $HOG.cds.aln.ABSREL.json $HOG.cds.aln.ABSREL.cand.json