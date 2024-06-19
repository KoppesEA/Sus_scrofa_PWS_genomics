#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH -J SusSTAR
#SBATCH --output=genSTARGenome-%A_%a.out
#SBATCH --cpus-per-task=16 # Request that ncpus be allocated per process.
#SBATCH --mem=256g # Memory pool for all cores (see also --mem-per-cpu)

module load star/2.7.9a

WKDIR=/bgfs/rnicholls/REFGenomes/Sscrofa_v11.1_v104
STARDIR=$WKDIR/Sscr11.1v104_STAR
FASTA=$WKDIR/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa
GTF=$WKDIR/Sus_scrofa.Sscrofa11.1.104.gtf

mkdir -p $STARDIR

STAR \
--runThreadN 16 \
--runMode genomeGenerate \
--genomeDir $STARDIR \
--genomeFastaFiles $FASTA \
--sjdbGTFfile $GTF \
--sjdbOverhang 74