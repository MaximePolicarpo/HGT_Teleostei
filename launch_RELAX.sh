#!/bin/bash


#SBATCH --job-name=RELAX   # Job name



#conda = miniprot ; cd-hit ; seqkit

eval "$(conda shell.bash hook)"
conda activate HyPhy_Env

HOG=$1

hyphy relax --alignment $HOG.cds.aln --tree $HOG.prot.aln.treefile.SelecMarked --code Universal --test test --reference background ENV=TOLERATE_NUMERICAL_ERRORS=1
