#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH -J TrimO
#SBATCH --output=SusTrim-%A_%a.out
#SBATCH --cpus-per-task=3 # Request that ncpus be allocated per process.
#SBATCH --mem=32g # Memory pool for all cores (see also --mem-per-cpu)
##array should start from zero
#SBATCH --array=0-1

module load trimgalore/0.5.0

names=($(cat PWS_Sscrofa_ExpList))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

#refer to --cpus-per-task
FASTQDIR=/bgfs/rnicholls/PWS_Sus_2017/raw_fastq
WKDIR=/bgfs/rnicholls/PWS_Sus_2017/STAR_Susv98
INPUT_FASTQ1=$FASTQDIR/${names[${SLURM_ARRAY_TASK_ID}]}_R1_001.fastq.gz
INPUT_FASTQ2=$FASTQDIR/${names[${SLURM_ARRAY_TASK_ID}]}_R2_001.fastq.gz
OUTPUT=$WKDIR/${names[${SLURM_ARRAY_TASK_ID}]}


echo $INPUT_FASTQ1
echo $INPUT_FASTQ2
echo $OUTPUT


#Trim_galore: cutadapt
mkdir $OUTPUT
trim_galore \
--paired \
--retain_unpaired \
--fastqc \
--gzip \
--output $OUTPUT \
$INPUT_FASTQ1 $INPUT_FASTQ2

