#BASH script to select PWS syntenic domains in Sus sscrofa 10.87 chromosome 7 PWS region from mouse annotation GRCm39.104
#!bin/bash
gzcat Sus_scrofa.Sscrofa10.2.87.gtf.gz |
grep -v "^##" |
awk '{ if ($1==1) print $0 }' | \
awk '{ if ($3== "gene") print $0 }' | \
awk '{ if ($4 >= 157783500) print $0 }' | \
awk '{ if ($4 <= 158900000) print $0 }' | \
cut -f9 | \
cut -d ";" -f1,3 | \
cut -d \" -f2,4 | \
sed 's/\"/\t/g' > Ssus_chrom1_PWS_Ensembl.tsv

gzcat Sus_scrofa.Sscrofa10.2.87.gtf.gz |
grep -v "^##" |
awk '{ if ($1=="GL896379.1") print $0 }' | \
awk '{ if ($3== "gene") print $0 }' | \
cut -f9 | \
cut -d ";" -f1,3 | \
cut -d \" -f2,4 | \
sed 's/\"/\t/g' > Ssus_GL896379_PWS_Ensembl.tsv

cat Ssus_chrom1_PWS_Ensembl.tsv Ssus_GL896379_PWS_Ensembl.tsv > Ssus_PWS_Ensembl.tsv

awk -F'\t' '{
    if ($2 == "ensembl" || $2 == "SNORD115" || $2 == "SNORD116") {
        print $1
    } else {
        print $2
    }
}' Ssus_PWS_Ensembl.tsv > Ssus_PWS_Ensembl_geneIDs.tsv


# 1:157783500-158900000 UBE3A + SNORD115 + SNORD115 + MAGEL2
# Chromosome GL896379: 1-22,821 SNRPN + SNORD107