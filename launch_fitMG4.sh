#!/bin/bash


#SBATCH --job-name=FitMG4   # Job name



#conda = miniprot ; cd-hit ; seqkit

eval "$(conda shell.bash hook)"
conda activate HyPhy_Env

HOG=$1

hyphy FitMG94.bf --alignment $HOG.cds.aln --tree $HOG.prot.aln.treefile --code Universal --type local ENV=TOLERATE_NUMERICAL_ERRORS=1