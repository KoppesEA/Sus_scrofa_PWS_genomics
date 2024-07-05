#!/bin/bash
#
#SBATCH -N 1 # 1 node
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH -J Sscr_VEP104
#SBATCH --output=Sscr_VEP104-%A_%a.txt
#SBATCH --cpus-per-task=8 # 8-vep forks
#SBATCH --mem=32g # (see also --mem-per-cpu)
#SBATCH --array=0-1 # Two PWS pig fibroblast RNA-Seq samples

#set path RNAvar VCF directories and output directories
DIRvcf=/bgfs/rnicholls/PWS_Sus_2017/RNAvar-nfcore/variant_calling
DIRvepcache=/bgfs/rnicholls/REFGenomes/Sscrofa_v11.1_v104/sus_scrofa/104_Sscrofa11.1

#make output directory
mkdir $DIRdedup

#Load VEP v95
module load vep/95

#Extract fq sample names from list text file
names=($(tail -n +2 samplesheet.csv | cut -d, -f1))
PWS_FIBR_SAMPLE=${names[${SLURM_ARRAY_TASK_ID}]}
echo $GSE
VCFin=$DIRvcf/$PWS_FIBR_SAMPLE/${PWS_FIBR_SAMPLE}haplotypecaller.filtered.vcf.gz
VCFout=$DIRvcf/$PWS_FIBR_SAMPLE/${PWS_FIBR_SAMPLE}haplotypecaller.filtered.VEP104.vcf.gz
TBI=$DIRvcf/$PWS_FIBR_SAMPLE/${PWS_FIBR_SAMPLE}haplotypecaller.filtered.vcf.gz.tbi
STATSout=$DIRvcf/$PWS_FIBR_SAMPLE/${PWS_FIBR_SAMPLE}haplotypecaller.filtered.VEP104.summary.html
echo $VCF
VEP_CACHE=
echo $VEP_CACHE

echo "Performing VEP annotation"

vep \
      -i $VCFin \
      -o $VCFout \
      --everything \
      --species pig \
      --fork 8 \
      --vcf \
      --stats_file $STATSout