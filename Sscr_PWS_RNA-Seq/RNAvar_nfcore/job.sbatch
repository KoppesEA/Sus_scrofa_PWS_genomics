#!/bin/bash
#
#SBATCH --job-name=nextflow_head
#SBATCH -c 1
#SBATCH --mem=8g
#SBATCH -t 3-00:00 # Runtime in D-HH:MM
#SBATCH --output=nextflow_head.out
#SBATCH --mail-user=eak37@pitt.edu
#SBATCH --mail-type=ALL

unset TMPDIR

module purge
module load nextflow/23.10.1
module load singularity/3.9.6
export NXF_SINGULARITY_CACHEDIR=/ihome/crc/install/genomics_nextflow/nf-core-rnavar-1.0.0/singularity-images
export SINGULARITY_CACHEDIR=/ihome/crc/install/genomics_nextflow/nf-core-rnavar-1.0.0/singularity-images


nextflow run /ihome/crc/install/genomics_nextflow/nf-core-rnavar-1.0.0/workflow \
	-resume \
	-profile htc \
	-name EAK_Sscr_RNAvar_attempt17 \
	-work-dir /bgfs/rnicholls/PWS_Sus_2017/RNAvar-nfcore/work \
	-params-file nf-params.json

