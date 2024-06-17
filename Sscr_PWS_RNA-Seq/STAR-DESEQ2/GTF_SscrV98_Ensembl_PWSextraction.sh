#BASH script to print PWS [chr1 and chr] Ensembl_ID and gene_names from Sus scrofa v98
#!bin/bash

cat Sus_scrofa.Sscrofa11.1.98.gtf  | awk '{ if ($1==1) print $_; }' | \
awk '{ if ($3== "gene") print $_; }' | \
awk '{ if ($4 >= 141500000) print $_; }' | \
awk '{ if ($4 <= 143000000) print $_; }' | \
cut -f9 | \
cut -d ";" -f1,3 | \
cut -d \" -f2,4 | \
sed 's/\"/\t/g' > Sscr_PWS_Ensembl_chr1PWS.tsv

cat Sus_scrofa.Sscrofa11.1.98.gtf  | awk '{ if ($1== "AEMK02000602.1") print $_; }' | \
awk '{ if ($3== "gene") print $_; }' | \
cut -f9 | \
cut -d ";" -f1,3 | \
cut -d \" -f2,4 | \
sed 's/\"/\t/g' > Sscr_PWS_Ensembl_AEMK02000602PWS.tsv