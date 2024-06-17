##script to import .gtf and then extract Snord115 and Snor116 sequences
##Erik Koppes 3/9/20

library(tidyverse)

##Extracting GTF
S_scrofa_gtf <- read_tsv("Sus_scrofa.Sscrofa11.98.sorted.gtf",
                        comment = "#",
                        col_names = c("seqname", "source", "feature", "start", "end",
                                      "score", "strand", "frame", "attribute"))
S_scrofa_gtf_1 <- S_scrofa_gtf %>% separate(attribute, letters[1:15], sep = ";")

##Extracting SNORD115 gene
S_scrofa_gtf_SNORD115 <- S_scrofa_gtf_1[str_detect(S_scrofa_gtf_1$c, pattern = fixed("SNORD115")),] %>%
  select(seqname, start, end, strand, a , c) %>%
  mutate(to_print = str_c(seqname, ":", start, "-", end, sep =""))

write_csv(S_scrofa_gtf_SNORD115, "S_scrofa_gtf_SNORD115.csv")

S_scrofa_gtf_SNORD115 %>% select(to_print) %>%
  write_tsv("S_scrofa_SNORD115.tsv", col_names = F)

##use below cmd line w/ samtools to extract, actually on negative strand
#samtools faidx --reverse-complement -o F_catus_6.0_SNORD115.fa Felis_catus.Felis_catus_9.0.dna.toplevel.fa F_catus_SNORD115.tsv


                          
##Extracting SNORD116
S_scrofa_gtf_SNORD116 <- S_scrofa_gtf_1[str_detect(S_scrofa_gtf_1$c, pattern = fixed("SNORD116")),] %>%
  select(seqname, start, end, strand, a, c) %>%
  mutate(seqname = "AEMK02000602.1") %>% #need to put col_types
  mutate(to_print = str_c(seqname, ":", start, "-", end, sep =""))
                          
write_csv(S_scrofa_gtf_SNORD116, "S_scrofa_gtf_SNORD116.csv")
                          
S_scrofa_gtf_SNORD116 %>% select(to_print) %>%
  write_tsv("S_scrofa_SNORD116.tsv", col_names = F)
