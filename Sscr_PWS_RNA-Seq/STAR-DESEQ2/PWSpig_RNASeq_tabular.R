##  Tabulation of RNA-Seq of PWS 120-1 clonal fibroblast lines
##  Erik Koppes, University of Pittsburgh
##  October 24th, 2023
##  Study Design::  #120_1_1 is PWS patDel; #120_1_2 is PWS matinv

## rm(list = ls())
## Load R Library Dependencies
library(dplyr)
library(readr)
library(openxlsx)

## Import ht-seq union-all .tsv count tables (Sus scrofa v98)
PWSpig120_1_1_df <- read_tsv("./120-1_subline_1_S1.tsv",col_names = c("Sscr_Ens98_ID", "GeneName", "HTSeqCount" ))
PWSpig120_2_2_df <- read_tsv("./120-1_subline_2_S2.tsv",col_names = c("Sscr_Ens98_ID", "GeneName", "HTSeqCount" ))

## Import and Join count tables
PWSpig120_tbl_all <- left_join(
  read_tsv("./120-1_subline_1_S1.tsv",col_names = c("Sscr_Ens98_ID", "GeneName", "HTSeqCount" )),
  read_tsv("./120-1_subline_2_S2.tsv",col_names = c("Sscr_Ens98_ID", "GeneName", "HTSeqCount" )),
  by = c("Sscr_Ens98_ID", "GeneName"), suffix = c("_120_1_1", "_120_1_2")) 

##  Calculate fold change and log2FC (#120_1_1 is PWS patDel; #120_1_2 is PWS matinv)
PWSpig120_tbl_all <- PWSpig120_tbl_all %>% mutate (FoldChange = HTSeqCount_120_1_1/HTSeqCount_120_1_2) %>%
  mutate(Log2FC = log2(FoldChange)) %>%
  arrange(FoldChange)

## Filter >= 10 counts
PWSpig120_tbl_10c <- PWSpig120_tbl_all %>% filter(HTSeqCount_120_1_2 >= 10)

## Table for PWS genes see GTF_SscrV98_Ensembl_PWSextraction.sh for GTF to .tsv
Sscr_PWS_genes <- rbind(read_tsv("./Sscr_PWS_Ensembl_chr1PWS.tsv", col_names = c("Sscr_Ens98_ID", "GeneName")), # if GeneName=Ensembl no geneName was specified in v98
                      read_tsv("./Sscr_PWS_Ensembl_AEMK02000602PWS.tsv", col_names = c("Sscr_Ens98_ID", "GeneName")))
PWSpig120_tbl_all_PWSgenes <- PWSpig120_tbl_all %>% filter(Sscr_Ens98_ID %in% unique(Sscr_PWS_genes$Sscr_Ens98_ID)) %>% arrange(desc(HTSeqCount_120_1_2))
PWSpig120_tbl_10c_PWSgenes <- PWSpig120_tbl_10c %>% filter(Sscr_Ens98_ID %in% unique(Sscr_PWS_genes$Sscr_Ens98_ID)) %>% arrange(desc(HTSeqCount_120_1_2))

##  Save as Excel
#write.xlsx()
wb <- createWorkbook()
addWorksheet(wb, sheetName = "PWS120_1_allgenes")
addWorksheet(wb, sheetName = "PWS120_1_genes10c")
addWorksheet(wb, sheetName = "PWSgenes_43")
addWorksheet(wb, sheetName = "PWSgenes_7")
writeDataTable(wb, sheet = 1, x = PWSpig120_tbl_all, colNames = TRUE, rowNames = FALSE, tableStyle = "TableStyleLight12")
writeDataTable(wb, sheet = 2, x = PWSpig120_tbl_10c, colNames = TRUE, rowNames = FALSE, tableStyle = "TableStyleLight12")
writeDataTable(wb, sheet = 3, x = PWSpig120_tbl_all_PWSgenes, colNames = TRUE, rowNames = FALSE, tableStyle = "TableStyleLight12")
writeDataTable(wb, sheet = 4, x = PWSpig120_tbl_10c_PWSgenes, colNames = TRUE, rowNames = FALSE, tableStyle = "TableStyleLight12")
saveWorkbook(wb, "PWSpig120_1_Sscrv98RNAseq.xlsx", overwrite = TRUE)  ## save to working directory


#options(openxlsx.borderColour = "#4F80BD")
#options(openxlsx.borderStyle = "thin")
#modifyBaseFont(wb, fontSize = 10, fontName = "Arial Narrow")

