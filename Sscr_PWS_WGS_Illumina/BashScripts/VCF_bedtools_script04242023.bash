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

WKDIR=/bgfs/rnicholls/PWS_Sus_2022/30-764669944/samples
OUTDIR=$WKDIR/VCF_sgRNA

rm -r $OUTDIR
mkdir $OUTDIR

for VCF_gz in $WKDIR/*/*.vcf.gz
do
	VCF_file=`basename $VCF_gz .gz`
	VCF_base=`basename $VCF_file .vcf`
	echo "extracting " $VCF_gz "to " $VCF_file ".vcf"
	gunzip -c $VCF_gz > $WKDIR/$VCF_file
	bedtools intersect \
	-a $WKDIR/$VCF_file \
	-b $WKDIR/SusPWS_sgRNA_1.bed > $OUTDIR/${VCF_base}_sgRNA1_OT.vcf

	bedtools intersect \
	-a $WKDIR/$VCF_file \
	-b $WKDIR/SusPWS_sgRNA_5.bed > $OUTDIR/${VCF_base}_sgRNA5_OT.vcf
done