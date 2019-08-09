#!/bin/bash
#SBATCH -N 1 -n 8
#SBATCH --mem=24G
export PATH=/opt/apps/labs/gtac/opt/R-3.4.1/bin:$PATH
#export R_LIBS=

mkdir /scratch/gtac/analysis/rna_seq/8998_3_magee-s4483/A_F_gene-level_eXpresso
cd /scratch/gtac/analysis/rna_seq/8998_3_magee-s4483/A_F_gene-level_eXpresso
Rscript /scratch/gtac/analysis/rna_seq/8998_3_magee-s4483/code/eXpresso.A.R -k /scratch/gtac/analysis/rna_seq/8998_3_magee-s4483/A_H.pair_file -m fit -b FALSE -f gene -s mouse -p rnaseq -c voom
