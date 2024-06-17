#BASH script to print ALL Ensembl_ID and gene_names from Sus scrofa v98
#!bin/bash
cat Sus_scrofa.Sscrofa11.1.98.gtf  |
awk '{ if ($3== "gene") print $_; }' | \
cut -f9 | \
cut -d ";" -f1,3 | \
cut -d \" -f2,4 | \
sed 's/\"/\t/g' > Sscr_GeneID_Ensembl.tsv
