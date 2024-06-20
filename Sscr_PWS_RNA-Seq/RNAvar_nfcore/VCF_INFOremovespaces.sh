#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH -J VCFedit
#SBATCH --output=htseq_PWSSscrofav98-%A_%a.out
#SBATCH --cpus-per-task=1
#SBATCH --mem=16g #

module load gcc/8.2.0
module load
module load bcftools/1.15.1

bcftools view sus_scrofa.vcf.gz | awk -v FS='\t' -v OFS='\t' '{gsub(/ /, "_", $8); print}' | bcftools view -Oz -o sus_scrofa_forGATK.vcf.gz
bcftools index -t sus_scrofa_forGATK.vcf.gz