#!/usr/bin/env bash
#SBATCH --export=ALL
#SBATCH -N 1
#SBATCH --mem=6G
##SBATCH -e samtools_bam2fq.sh.e
##SBATCH -o samtools_bam2fq.sh.o
# Time: 2021-01-17 20:41:21

module load samtools
input_bam=$1
input="${input_bam%.*}"
sorted_bam=${input}.sorted.bam
samtools sort -n $input_bam -o $sorted_bam

# save fastq reads in separate R1 and R2 files
samtools fastq -@ 8 $sorted_bam \
    -1 ${input}.R1.fastq.gz \
    -2 ${input}.R2.fastq.gz \
    -0 /dev/null -s /dev/null -n

rm $sorted_bam

module load fastqc
fastqc ${input}.R*fastq.gz
