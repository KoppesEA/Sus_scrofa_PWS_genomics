#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH -J htseq-PWSSscrofav98
#SBATCH --output=htseq_PWSSscrofav98-%A_%a.out
#SBATCH --cpus-per-task=1 # HTseq does not support multithread
#SBATCH --mem=64g # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --array=0-1

names=($(cat PWS_Sscrofa_ExpList))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

WKDIR=/bgfs/rnicholls/PWS_Sus_2017/STAR_Susv98
READ_BASE=${names[${SLURM_ARRAY_TASK_ID}]}
INPUT_BAM=/$WKDIR/$READ_BASE/${READ_BASE}Aligned.sortedByCoord.out.bam
ANNOT=/bgfs/rnicholls/REFGenomes/SscrofaSus_v11.1.98/Sus_scrofa.Sscrofa11.1.98.gtf

echo $READ_BASE
echo $INPUT_BAM
echo $ANNOT

module load htseq/0.11.2
#htseq with union overlap and non-unique all options
htseq-count \
--format bam \
--order pos \
--mode union \
--nonunique all \
--minaqual 1 \
--stranded reverse \
--type exon \
--idattr gene_id \
--additional-attr gene_name \
$INPUT_BAM $ANNOT > $WKDIR/htseq_unionall/${READ_BASE}.tsv
