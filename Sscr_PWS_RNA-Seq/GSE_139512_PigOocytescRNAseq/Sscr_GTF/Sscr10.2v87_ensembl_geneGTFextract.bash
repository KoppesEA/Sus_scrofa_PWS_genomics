#BASH ZFS script to Ensembl_ID and gene_names from Mmus GRCm39 Gencode
#!bin/bash
gzcat Sus_scrofa.Sscrofa10.2.87.gtf.gz |
grep -v "^##" |
awk -F'\t' '$3 == "gene" { print $9 }' |
cut -d ";" -f1,3 |
cut -d \" -f2,4 |
sed 's/\"/\t/g' > Sus_scrofa.Sscrofa10.2.87.tsv

awk -F'\t' 'BEGIN {OFS="\t"} {
    if ($2 == "ensembl" || $2 == "ensembl_havana" || $2 == "havana" || $2 == "insdc") {
        print $1, $1
    } else {
        print $1, $2
    }
}' Sus_scrofa.Sscrofa10.2.87.tsv > Sus_scrofa.Sscrofa10.2.87.geneIDconv.tsv

##need to make file of nonunique ensembl tracking IDs