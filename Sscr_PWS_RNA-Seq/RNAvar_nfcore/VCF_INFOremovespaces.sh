#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH -J VCFedit
#SBATCH --output=htseq_PWSSscrofav98-%A_%a.out
#SBATCH --cpus-per-task=1
#SBATCH --mem=16g #

module load gcc/8.2.0
module load htslib/1.9

bgzip -d -c | awk -v FS='\t' -v OFS='\t' '{gsub(/ /, "_", $8); print}' | bgzip -o sus_scrofa_forGATK.vcf.gz
tabix -p sus_scrofa_forGATK.vcf.gz