#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH -J VCFtoolz
#SBATCH --output=VCFtoolz-%A_%a.out
#SBATCH --cpus-per-task=1 #
#SBATCH --mem=32g # Memory pool for all cores (see also --mem-per-cpu)

echo "loading VCFtools 0.1.16"
module load gcc/8.2.0
module load vcftools/0.1.16

WKDIR=/bgfs/rnicholls/PWS_Sus_2022/30-764669944/samples
OUTDIR=$WKDIR/VCF_sgRNA

for VCF_gz in $WKDIR/*/*.vcf.gz
do
	echo $VCF_gz
	VCF_base=`basename $VCF_gz .vcf.gz`
	vcftools --gzvcf $VCF_gz --chr AEMK02000602.1 --recode --out $OUTDIR/${VCF_base}_AEMK
done

