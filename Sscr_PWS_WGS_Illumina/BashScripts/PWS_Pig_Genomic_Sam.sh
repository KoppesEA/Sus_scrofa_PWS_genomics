#!/bin/bash
#
#SBATCH -N 1
#SBATCH -t 1-00:00
#SBATCH -J SAM_RIP
#SBATCH --output=PWS_SamWigPrep.txt
#SBATCH --cpus-per-task=4


##Script to extract PWS mapped .bam files and then convert to wiggle files

WKDIR=/bgfs/rnicholls/PWS_Sus_2022/30-764669944/samples
OUTDIR=$WKDIR/wiggle

module load gcc/8.2.0
module load samtools/1.14

for bamfile in ${WKDIR}/*/*.aln.bam
do
echo $bamfile

filebase=`basename $bamfile .aln.bam`
echo $filebase

#index sorted BAM files (Genewiz Enesemble Sscrofa 11.1 Alignment)
## Chr1 contains IPW/SNORD115/UBE3A and the proximal  NDN and unannotated MKRN3 and MAGEL2
## AEMK02000602.1 SNURF, SNRPN and 16 annotated SNORD116 copies on a contig 
#extract all chr1 alignments from sorted BAM file and redirect output to IGV-chr1 folder
samtools view -@ 4 -b -o $OUTDIR/${filebase}_Chr1.bam $bamfile "1"
samtools view -@ 4 -b -o $OUTDIR/${filebase}_Chr1PWS.bam $bamfile "1:141500000-143000000"
samtools view -@ 4 -b -o $OUTDIR/${filebase}_AEMK02000602.bam $bamfile "AEMK02000602.1"

#index extract alignment files
samtools index -@ 4 $OUTDIR/${filebase}_Chr1.bam
samtools index -@ 4 $OUTDIR/${filebase}_Chr1PWS.bam
samtools index -@ 4 $OUTDIR/${filebase}_AEMK02000602.bam

done
