#!/bin/bash
#
for VCF_file in *.gz
do
	VCF_base=`basename $VCF_file .gz`
	gunzip -c $VCF_file > $VCF_base
done