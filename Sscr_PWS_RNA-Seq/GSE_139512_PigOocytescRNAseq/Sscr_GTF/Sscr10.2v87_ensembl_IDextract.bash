#BASH ZFS script to Ensembl_ID and gene_names from Mmus GRCm39 Gencode
#!bin/bash
zcat Sus_scrofa.Sscrofa10.2.87.gtf.gz |
grep -v "^##" |
awk -F'\t' '$3 == "gene" { print $9 }' |
cut -d ";" -f1,3 |
cut -d \" -f2,4 |
sed 's/\"/\t/g' > Sus_scrofa.Sscrofa10.2.87.tsv

cat Sus_scrofa.Sscrofa10.2.87.tsv | cut -f2 | sort | uniq -d > Sus_scrofa.Sscrofa10.2.87_nonuniq.tsv