#!/bin/bash


#SBATCH --job-name=BLAST-ALL   # Job name

module load BLAST

myquery=$1
myevalue=$2
myoutput=$3

blastn -query $myquery -db Concatenated_assemblies.fa -outfmt 6 -out $myoutput -num_threads 10 -evalue $myevalue

