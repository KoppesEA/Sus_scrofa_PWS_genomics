#!/bin/bash
#
#SBATCH -N 1
#SBATCH -t 1-00:00
#SBATCH -J SAM_WIG_RIP
#SBATCH --output=PWS_wiggleprep.txt
#SBATCH --cpus-per-task=4


##Script to extract PWS mapped .bam files and then convert to wiggle files

WKDIR=/bgfs/rnicholls/PWS_Sus_2022/30-764669944/samples
OUTDIR=$WKDIR/wiggle

module load deeptools/3.3.0

for bamfile in ${WKDIR}/*/*.aln.bam
do
echo $bamfile

filebase=`basename $bamfile .aln.bam`
echo $filebase

#index sorted BAM files (Genewiz Enesemble Sscrofa 11.1 Alignment)
## Chr1 contains IPW/SNORD115/UBE3A and the proximal  NDN and unannotated MKRN3 and MAGEL2
## AEMK02000602.1 SNURF, SNRPN and 16 annotated SNORD116 copies on a contig 

## Convert pre-sorted and indexed .bam files to bigwig
bamCoverage -p 4 -b $OUTDIR/${filebase}_Chr1.bam -o $OUTDIR/${filebase}_Chr1.bw
bamCoverage -p 4 -b $OUTDIR/${filebase}_Chr1PWS.bam -o $OUTDIR/${filebase}_Chr1PWS.bw
bamCoverage -p 4 -b $OUTDIR/${filebase}_AEMK02000602.bam -o $OUTDIR/${filebase}_AEMK02000602.bw

done
