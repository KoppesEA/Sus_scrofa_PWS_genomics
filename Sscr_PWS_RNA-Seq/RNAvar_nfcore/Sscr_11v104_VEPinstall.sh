#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH -J VEPcacheinstall
#SBATCH --output=VEPcacheinstall-%A_%a.out
#SBATCH --cpus-per-task=1
#SBATCH --mem=16g #

#load VEP v95
module load vep/95

DIRvepcache=/bgfs/rnicholls/PWS_Sus_2017/RNAvar-nfcore/variant_calling/vep
mkdir $DIRvepcache

vep_install \
	-a cf \
	-s sus_scrofa \
	-y 11
	-v 104 \
	-c $DIRvepcache \
	--assembly 11
