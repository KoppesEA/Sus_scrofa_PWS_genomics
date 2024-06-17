#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH -J Bedtoolz
#SBATCH --output=Bedtoolz-%A_%a.out
#SBATCH --cpus-per-task=1 #
#SBATCH --mem=32g # Memory pool for all cores (see also --mem-per-cpu)

module load gcc/8.2.0
module load bedtools/2.29.0

for VCF_file in *hard-filtered.vcf
do
	VCF_base=`basename $VCF_file .vcf`

	bedtools intersect \
	-a $VCF_file \
	-b gRNA5_1_OT.bed > ${VCF_base}_5_1_OT.vcf

	bedtools intersect \
	-a $VCF_file \
	-b gRNA6_2_OT.bed > ${VCF_base}_6_2_OT.vcf
done