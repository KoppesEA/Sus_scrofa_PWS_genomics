## GSE139512 Pig (Sus scrofa) scRNA-Seq analysis of oocyte and embryos
## Erik Koppes, PhD; University of Pittsburgh
## 6/25/2024
#rm(list=ls())

## load dependencies
library(readr)
library(tibble)
library(stringr)
library(dplyr)
library(Seurat)
#install.packages("Seurat")
library(patchwork)
## ggplot2 to export graphs
library(ggplot2)

GSE_139512_SraRunTbl <- read_csv("GSE_139512_SraRunTable.txt")

# Import GSE57249 matrice
GSE139512_df <- read_tsv("./GSE139512%5Ffpkm.txt") #note row1 is gene names, cols are sample GSE

## Ensembl GTF for Sscrofa 10.2.87
# set Ensembl Gene IDs to rownames rather than column [see awk/sort commands]
Sscr10.2v87_IDs <- read_tsv("./Sscr_GTF/Sus_scrofa.Sscrofa10.2.87.tsv",
                          col_names = c("Ensembl_ID", "Gene_Name"))
# set nonuniq Gene_Names [see sort/uniq commands]
Sscr10.2v87_nonuniq <- read_tsv("./Sscr_GTF/Sus_scrofa.Sscrofa10.2.87_nonuniq.tsv",
                              col_names = c("Gene_Name"))
# get corresponding Ensembl Ensemble_IDs (4013) that match nonuniq Gene_Names (293); mostly s7k, snords, U6 etc.
Sscr10.2v87_nonuniq_IDs <- Sscr10.2v87_IDs %>% filter(Gene_Name %in% Sscr10.2v87_nonuniq$Gene_Name)


# filter out nonuniq ensembl_IDs from count matrix then join to get Gene_Names
GSE139512_df <- GSE139512_df %>% filter(!(ID %in% Sscr10.2v87_nonuniq_IDs$Ensembl_ID)) %>%
  left_join(Sscr10.2v87_IDs, by = join_by(ID == Ensembl_ID)) %>%
  select(-ID)

## and set column Gene_Name to rownames
GSE139512_df <- column_to_rownames(GSE139512_df, var = "Gene_Name")

# Extract metadata from CPB metadata
source("GSE139512_metaScript.R")

# Extract metadata from matrix
series_full <- CPB_meta %>% column_to_rownames(var = "sample")

# Bindcols and Transpose metadata; Cell names as rows and  Description as columns
series_full_t <- data.frame(t(series_full))
colnames(series_full_t) <- c("Celltype", "Stage", "Meta_ID")
rownames(series_full_t) == colnames(GSE57249_df)  ##all TRUE--goodtogo

# Create Seurat Object
GSE57249_seurat <- CreateSeuratObject(counts = GSE57249_df, meta.data = series_full_t)