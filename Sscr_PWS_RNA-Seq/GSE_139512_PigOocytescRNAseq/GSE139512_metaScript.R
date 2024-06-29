## script to import SRA Run meta data from SRARunTable.txt

library(readr)
library(dplyr)

GSE_139512_meta <- read_csv("GSE_139512_SraRunTable.txt") %>%
  select(Run, BioProject, BioSample, cell_type, Developmental_stage,  `Sample Name`)
