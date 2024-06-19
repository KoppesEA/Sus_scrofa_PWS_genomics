#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 3-00:00 # Runtime in D-HH:MM
#SBATCH -J STARSscrofa
#SBATCH --output=STARSscrofa-%A_%a.out
#SBATCH --cpus-per-task=4 # Request that ncpus be allocated per process.
#SBATCH --mem=48g # Memory pool for all cores (see also --mem-per-cpu)
##array should start from zero
#SBATCH --array=0-2

module load star/2.7.0e
 
names=($(cat PWS_Sscrofa_ExpList))

WKDIR=/bgfs/rnicholls/PWS_Sus_2017/STAR_Susv98
READ_BASE=${names[${SLURM_ARRAY_TASK_ID}]}
INPUT_FASTQ1_trim=$WKDIR/$READ_BASE/${READ_BASE}_R1_001_val_1.fq.gz
INPUT_FASTQ2_trim=$WKDIR/$READ_BASE/${READ_BASE}_R2_001_val_2.fq.gz
ANNOTDIR=/bgfs/rnicholls/REFGenomes/SscrofaSus_v11.1.98/STAR_v98


echo $READ_BASE
echo $INPUT_FASTQ1_trim
echo $INPUT_FASTQ2_trim

# Aligning by STAR
rm -rf $WKDIR/${READ_BASE}tmp/*
rmdir $WKDIR/${READ_BASE}tmp
mkdir $WKDIR/${READ_BASE}

STAR \
--outTmpDir $WKDIR/${READ_BASE}tmp \
--runThreadN 4 \
--outFilterIntronMotifs RemoveNoncanonicalUnannotated \
--outSAMprimaryFlag AllBestScore \
--chimSegmentMin 15 \
--chimOutType Junctions \
--twopassMode Basic \
--readFilesIn $INPUT_FASTQ1_trim $INPUT_FASTQ2_trim \
--outFileNamePrefix $WKDIR/${READ_BASE}/${READ_BASE} \
--quantMode GeneCounts \
--outStd Log \
--genomeDir $ANNOTDIR \
--readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate