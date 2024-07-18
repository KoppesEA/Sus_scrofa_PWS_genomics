#BASH ZFS script to Ensembl_ID and gene_names from Mmus GRCm39 Gencode
#!bin/bash
gzcat Sus_scrofa.Sscrofa10.2.87.gtf.gz |
grep -v "^##" |
awk -F'\t' '$1 == "MT" { print $0 }' |
awk -F'\t' '$3 == "gene" { print $9 }' |
cut -d ";" -f1,3 |
cut -d \" -f2,4 |
sed 's/\"/\t/g' |
cut -f2 | sort -u > Sus_scrofa.mito.Sscrofa10.2.87.tsv