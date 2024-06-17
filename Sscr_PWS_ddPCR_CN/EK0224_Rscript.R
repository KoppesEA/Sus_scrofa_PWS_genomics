##  EK0224 PWSICdel Pig ddPCR CNV and Methylation assay
##  Script for pig PWSIC CNV and methylation ddPCR Graphing Copy of PWSIC/UBE3A or SNORD107/UBE3A using RE digests for EcoRI, EcoRI+HpaII, McrBC, EcoRI + BstUI
##  Plate 1 EcoRI and EcoRI/HpaII (PWSIC/UBE3A) and EcoRI (SNORD107/UBE3A)
##  Plate 2 McrBC and EcoRI/BstUI (PWSIC/UBE3A)
##  Dec 9th 2021; revused 2/27/2023

##  library dependencies
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(forcats)
library(stringr)


## Vec list for designating RE digests
Vec_01_to_04 <- c(paste0("A", 0, 1:4), paste0("B", 0, 1:4), paste0("C", 0, 1:4), paste0("D", 0, 1:4),
                  paste0("E", 0, 1:4), paste0("F", 0, 1:4), paste0("G", 0, 1:4), paste0("H", 0, 1:4))
Vec_05_to_08 <- c(paste0("A", 0, 5:8), paste0("B", 0, 5:8), paste0("C", 0, 5:8), paste0("D", 0, 5:8),
                  paste0("E", 0, 5:8), paste0("F", 0, 5:8), paste0("G", 0, 5:8), paste0("H", 0, 5:8))
Vec_09_to_12 <- c(c(paste0("A", "09"), paste0("A", 10:12)), c(paste0("B", "09"), paste0("B", 10:12)),
                  c(paste0("C", "09"), paste0("C", 10:12)), c(paste0("D", "09"), paste0("D", 10:12)),
                  c(paste0("E", "09"), paste0("E", 10:12)), c(paste0("F", "09"), paste0("F", 10:12)),
                  c(paste0("G", "09"), paste0("G", 10:12)), c(paste0("H", "09"), paste0("H", 10:12)))

##  Import and clean data
EK0224_plate_1 <- read_csv("EK0224_PWSpig_probe_manual.csv") %>%
  select(Well, Sample, Target, ReferenceUsed, "Conc(copies/µL)", CNV, PoissonCNVMax, PoissonCNVMin, Ratio, PoissonRatioMax, PoissonRatioMin) %>% 
  filter(Target %in% c("PWSIC", "SNORD107")) %>%
  filter(!(Well %in% c("A01", "B01", "A05", "B05", "A09", "B09"))) %>% ## remove buffer-control and NTC reactions
#  mutate(Target = factor(Target, levels =  c("ACADVL_ex10", "ACADVL_in11", "ACADVL_ex15", "ACADVL_ex20"))) %>%
  mutate(Conc = as.numeric(`Conc(copies/µL)`)) %>%
  mutate(RE_dig = if_else(Well %in% Vec_01_to_04, "EcoRI",
                         if_else(Well %in% Vec_05_to_08, "EcoRI-HpaII",
                                 if_else(Well %in% Vec_09_to_12, "EcoRI", "error")))) %>%
  mutate(plate = "plate_1")

EK0224_plate_2 <- read_csv("EK0224_PWSpig_probe_plate2_manual.csv") %>%
  select(Well, Sample, Target, ReferenceUsed, "Conc(copies/µL)", CNV, PoissonCNVMax, PoissonCNVMin, Ratio, PoissonRatioMax, PoissonRatioMin) %>% 
  filter(Target == "PWSIC") %>%
  filter(!(Well %in% c("A01","B01", "A05", "B05"))) %>% ## remove buffer-control and NTC reactions
  #  mutate(Target = factor(Target, levels =  c("ACADVL_ex10", "ACADVL_in11", "ACADVL_ex15", "ACADVL_ex20"))) %>%
  mutate(Conc = as.numeric(`Conc(copies/µL)`)) %>%
  mutate(RE_dig = if_else(Well %in% Vec_01_to_04, "McrBC",
                          if_else(Well %in% Vec_05_to_08, "EcoRI-BstUI", "error"))) %>%
  mutate(plate = "plate_1")
  
##combine plates 1 and 2
EK0224_full <- rbind(EK0224_plate_1, EK0224_plate_2) %>%
  mutate(Targ_RE_combo = factor(str_c(Target, RE_dig, sep = "-")))




##adjust names for Sample
Sample_Conv <- data.frame(
  Sample = c(paste0("EK", 0, 2:9), paste0("EK", 10:31)),
  Sample_ID = factor(c("P039_Testis", "P039_Intestine", "P039_Pancreas", "P039_Liver", "P039_Kidney", "P039_Brain", "P039_Pituitary", "P039_Lung",
                       "P039_Heart", "P039_Muscle", "P040_Pancreas", "P040_Liver", "P041_Pancreas", "P041_Ovary", "P041_Liver",
                       "119-1", "120-1", "120-2", "120-3", "P039", "P040", "P041",
                       "CHO", "SCH-8", "SCH-16", "SCH-18", "SCH-19", "Yuc_Boar", "Yuc_GFP", "Duroc_Fetal")))


EK0224_full_rename <- EK0224_full %>% left_join(Sample_Conv) %>%
  select(Sample, Sample_ID, RE_dig, Target, Targ_RE_combo, ReferenceUsed, CNV, PoissonCNVMax, PoissonCNVMin, Ratio, PoissonRatioMax, PoissonRatioMin) %>%
  arrange(Sample) 



##plot of PWS_pigs by Target-RE_dig ACADVL exon
EK0224_full_rename %>%
  ggplot(aes(x = Sample_ID, y = CNV)) +
  geom_col(aes(fill = Targ_RE_combo)) +
  facet_wrap(.~Targ_RE_combo, ncol = 1) +
  theme(axis.text.x = element_text(angle =45, hjust =1)) +
  geom_errorbar(aes(ymin = PoissonCNVMin, ymax = PoissonCNVMax), width = 0.25, size =0.125) +
  scale_y_continuous(
    name = "Genomic Copy Number",
    limits = c(0, 2.5),
    breaks = c(0, 0.5, 1.0, 1.5, 2, 2.5)
  ) +
  #  scale_x_discrete(NULL, labels = NULL) +
#  labs(fill = "Cell Line") +
  scale_fill_brewer(palette = "Set3")


ggsave("EK0224_PWSpigddPCR_1.pdf",
       device = "pdf",
       units = "cm",
       width =  20,
       height = 25,
       dpi = "retina",
       useDingbats = FALSE)

## plot to just Adult pig DNA
EK0224_full_rename %>%
  mutate(Sample_ID = str_replace_all(!!sym("Sample_ID"),"P040", "120-2")) %>%  ##chatGPT code :)
  filter(Sample %in% c("EK17", "EK18", "EK22", "EK20")) %>%
  mutate(Sample_ID = as_factor(Sample_ID)) %>%
  mutate(Sample_ID = fct_relevel(Sample_ID, c("119-1", "120-1", "120-2", "120-3"))) %>%
  mutate(Targ_RE_combo = fct_relevel(Targ_RE_combo,
        c("SNORD107-EcoRI", "PWSIC-EcoRI", "PWSIC-EcoRI-BstUI", "PWSIC-EcoRI-HpaII", "PWSIC-McrBC"))) %>%
  ggplot(aes(x = Sample_ID, y = CNV, fill = Targ_RE_combo)) +
  geom_col(color = "black",
           position = position_dodge(0.8), width = 0.675) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust =1)) +
  geom_errorbar(aes(ymin = CNV, ymax = PoissonCNVMax),
                position = position_dodge(width=0.8),
                width = 0.25, size =0.125) +
  scale_y_continuous(
    name = "Genomic Copy Number",
    limits = c(0, 2.25),
    breaks = c(0, 0.5, 1.0, 1.5, 2),
    expand = c(0,0)
  ) +
  scale_fill_brewer(palette = "Set3")

ggsave("EK0224_PWSpigddPCR_1.pdf",
       device = "pdf",
       units = "cm",
       width =  9,
       height = 6,
       dpi = "retina",
       useDingbats = FALSE)

## plot to just CHO and SCH lines
EK0224_full_rename %>%
#  mutate(Sample_ID = str_replace_all(!!sym("Sample_ID"),"P040", "120-2")) %>%  ##chatGPT code :)
  filter(Sample %in% paste0("EK", 24:28)) %>%
#  mutate(Sample_ID = as_factor(Sample_ID)) %>%
#  mutate(Sample_ID = fct_relevel(Sample_ID, c("119-1", "120-1", "120-2", "120-3"))) %>%
  mutate(Targ_RE_combo = fct_relevel(Targ_RE_combo,
                                     c("SNORD107-EcoRI", "PWSIC-EcoRI", "PWSIC-EcoRI-BstUI", "PWSIC-EcoRI-HpaII", "PWSIC-McrBC"))) %>%
  ggplot(aes(x = Sample_ID, y = CNV, fill = Targ_RE_combo)) +
  geom_col(color = "black",
           position = position_dodge(0.8), width = 0.675) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust =1)) +
  geom_errorbar(aes(ymin = CNV, ymax = PoissonCNVMax),
                position = position_dodge(width=0.8),
                width = 0.25, size =0.125) +
  scale_y_continuous(
    name = "Genomic Copy Number",
    limits = c(0, 1.0),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    expand = c(0,0)
  ) +
  scale_fill_brewer(palette = "Set3")

ggsave("EK0224_PWSpigddPCR_2.pdf",
       device = "pdf",
       units = "cm",
       width =  9,
       height = 6,
       dpi = "retina",
       useDingbats = FALSE)

##plot all pig tissues
EK0224_full_rename %>%
  #  mutate(Sample_ID = str_replace_all(!!sym("Sample_ID"),"P040", "120-2")) %>%  ##chatGPT code :)
  filter(!(Sample %in% c(paste0("EK", 24:28), c("EK17", "EK18","EK19", "EK20")))) %>%
  mutate(Sample_ID = as_factor(Sample_ID)) %>%
  mutate(Sample_ID = fct_relevel(Sample_ID, c("Duroc_Fetal", "Yuc_GFP", "Yuc_Boar"))) %>% ##yucboar should be 1n
  mutate(Targ_RE_combo = fct_relevel(Targ_RE_combo,
                                     c("SNORD107-EcoRI", "PWSIC-EcoRI", "PWSIC-EcoRI-BstUI", "PWSIC-EcoRI-HpaII", "PWSIC-McrBC"))) %>%
  ggplot(aes(x = Sample_ID, y = CNV, fill = Targ_RE_combo)) +
  geom_col(color = "black",
           position = position_dodge(0.8), width = 0.675) +
  theme_bw() +
  theme(axis.text.x = element_text(angle =45, hjust =1)) +
  geom_errorbar(aes(ymin = CNV, ymax = PoissonCNVMax),
                position = position_dodge(width=0.8),
                width = 0.25, size =0.125) +
  scale_y_continuous(
    name = "Genomic Copy Number",
    limits = c(0, 2.40),
    breaks = c(0, 0.5, 1.0, 1.5, 2),
    expand = c(0,0)
  ) +
  scale_fill_brewer(palette = "Set3")

ggsave("EK0224_PWSpigddPCR_3.pdf",
       device = "pdf",
       units = "cm",
       width =  18,
       height = 9,
       dpi = "retina",
       useDingbats = FALSE)

