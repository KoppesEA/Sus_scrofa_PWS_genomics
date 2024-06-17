##Erik Koppes  Feb 23th 2020
##Script to examin single nucleotide variants in PKU swine

library(readr)
library(dplyr)
library(tidyr)
library(forcats)
library(openxlsx)

##import vcf files from PKU pigs hard filtered
MJ01_SNV_filt <- read_tsv("MJ01-116-1.hard-filtered.vcf", comment = "#", col_names = c("CHROM", "POS","ID", "REF", "ALT", "QUAL", "FILTER", "INFO"))
MJ02_SNV_filt <- read_tsv("MJ02-116-2.hard-filtered.vcf", comment = "#", col_names = c("CHROM", "POS","ID", "REF", "ALT", "QUAL", "FILTER", "INFO"))
MJ03_SNV_filt <- read_tsv("MJ03-YucGFP.hard-filtered.vcf", comment = "#", col_names = c("CHROM", "POS","ID", "REF", "ALT", "QUAL", "FILTER", "INFO"))

##subset for PAH region
MJ01_SNV_PAH <- MJ01_SNV_filt %>% filter(CHROM == 5 & (81236096< POS  & POS < 81463452)) %>% mutate(Sample = "116-1")
MJ01_SNV_PAH_ex6 <- MJ01_SNV_filt %>% filter(CHROM == 5 & (81435000< POS  & POS < 81445000)) %>% mutate(Sample = "116-1")
MJ02_SNV_PAH <- MJ02_SNV_filt %>% filter(CHROM == 5 & (81236096< POS  & POS < 81463452)) %>% mutate(Sample = "116-2")
MJ02_SNV_PAH_ex6 <- MJ02_SNV_filt %>% filter(CHROM == 5 & (81435000< POS  & POS < 81445000)) %>% mutate(Sample = "116-2")
MJ03_SNV_PAH <- MJ03_SNV_filt %>% filter(CHROM == 5 & (81236096< POS  & POS < 81463452)) %>% mutate(Sample = "Yuc:GFP")
MJ03_SNV_PAH_ex6 <- MJ03_SNV_filt %>% filter(CHROM == 5 & (81435000< POS  & POS < 81445000)) %>% mutate(Sample = "Yuc:GFP")
##add Sample column as removed _filt object (very big)
MJ01_SNV_PAH <- MJ01_SNV_PAH %>% mutate(Sample = "116-1")
MJ01_SNV_PAH_ex6 <- MJ01_SNV_PAH_ex6 %>% mutate(Sample = "116-1")
MJ02_SNV_PAH <- MJ02_SNV_PAH %>% mutate(Sample = "116-2")
MJ02_SNV_PAH_ex6 <- MJ02_SNV_PAH_ex6 %>% mutate(Sample = "116-2")
MJ03_SNV_PAH <- MJ03_SNV_PAH %>% mutate(Sample = "Yuc:GFP")
MJ03_SNV_PAH_ex6 <- MJ03_SNV_PAH_ex6 %>% mutate(Sample = "Yuc:GFP")


##combine samples to one tab
SNV_PAH_combo <- bind_rows(MJ01_SNV_PAH, MJ02_SNV_PAH, MJ03_SNV_PAH) %>% arrange(CHROM, POS) %>% select(Sample, everything())
SNV_PAH_ex6_combo <- bind_rows(MJ01_SNV_PAH_ex6, MJ02_SNV_PAH_ex6, MJ03_SNV_PAH_ex6) %>% arrange(CHROM, POS) %>% select(Sample, everything())


#Create Excel Table of SNVs
SNV_Pig_PAH_summary <- createWorkbook("SNV_Pig_PAH_summary")
addWorksheet(SNV_Pig_PAH_summary, "MJ01_116-1_PAH")
writeData(SNV_Pig_PAH_summary, "MJ01_116-1_PAH", MJ01_SNV_PAH, rowNames = F, colNames = T)
addWorksheet(SNV_Pig_PAH_summary, "MJ02_116-2_PAH")
writeData(SNV_Pig_PAH_summary, "MJ02_116-2_PAH", MJ02_SNV_PAH, rowNames = F, colNames = T)
addWorksheet(SNV_Pig_PAH_summary, "MJ03_YucGFP_PAH")
writeData(SNV_Pig_PAH_summary, "MJ03_YucGFP_PAH", MJ03_SNV_PAH, rowNames = F, colNames = T)
addWorksheet(SNV_Pig_PAH_summary, "All3_PAH")
writeData(SNV_Pig_PAH_summary, "All3_PAH", SNV_PAH_combo, rowNames = F, colNames = T)
saveWorkbook(SNV_Pig_PAH_summary, "SNV_Pig_PAH_summary.xlsx", overwrite = T)

SNV_Pig_PAHex6_summary <- createWorkbook("SNV_Pig_PAHex6_summary")
addWorksheet(SNV_Pig_PAHex6_summary, "MJ01_116-1_PAHex6")
writeData(SNV_Pig_PAHex6_summary, "MJ01_116-1_PAHex6", MJ01_SNV_PAH_ex6, rowNames = F, colNames = T)
addWorksheet(SNV_Pig_PAHex6_summary, "MJ02_116-2_PAHex6")
writeData(SNV_Pig_PAHex6_summary, "MJ02_116-2_PAHex6", MJ02_SNV_PAH_ex6, rowNames = F, colNames = T)
addWorksheet(SNV_Pig_PAHex6_summary, "MJ03_YucGFP_PAHex6")
writeData(SNV_Pig_PAHex6_summary, "MJ03_YucGFP_PAHex6", MJ03_SNV_PAH_ex6, rowNames = F, colNames = T)
addWorksheet(SNV_Pig_PAHex6_summary, "All3_PAHex6")
writeData(SNV_Pig_PAHex6_summary, "All3_PAHex6", SNV_PAH_ex6_combo, rowNames = F, colNames = T)
saveWorkbook(SNV_Pig_PAHex6_summary, "SNV_Pig_PAHex6_summary.xlsx", overwrite = T)



##function to extract gRNA5_1 off Targets
gRNA5_1_filter <- function(df) {
  filter(df,
         CHROM == 6 & (153946924 < POS & POS < 153946946) | CHROM == 2 & (36380787 < POS & POS < 36380809) |
           CHROM == 2 & (118033140 < POS & POS < 118033162) | CHROM == 8 & (129310246 < POS & POS < 129310268) |
           CHROM == 7 & (108651595 < POS & POS < 108651617) | CHROM == 2 & (28813230 < POS & POS < 28813252) |
           CHROM == 13 & (173341952 < POS & POS < 173341974) | CHROM == 6 & (154108660 < POS & POS < 154108682) |
           CHROM == 11 & (66518022 < POS & POS < 66518044) | CHROM == 8 & (117499059 < POS & POS < 117499081) |
           CHROM == 2 & (76478135 < POS & POS < 76478157) | CHROM == 1 & (242428459 < POS & POS < 242428481) |
           CHROM == 15 & (93947762 < POS & POS < 93947784) | CHROM == 6 & (105755606 < POS & POS < 105755628) |
           CHROM == 3 & (123499712 < POS & POS < 123499734) | CHROM == 5 & (48137882 < POS & POS < 48137904) |
           CHROM == 13 & (162836805 < POS & POS < 162836827) | CHROM == 1 & (115677643 < POS & POS < 115677665) |
           CHROM == 15 & (112846739 < POS & POS < 112846761) | CHROM == 7 & (65486048 < POS & POS < 65486070) |
           CHROM == 14 & (66673359 < POS & POS < 66673381) | CHROM == 2 & (76135945 < POS & POS < 76135967) |
           CHROM == "X" & (73572208 < POS & POS < 73572230) | CHROM == 6 & (122462984 < POS & POS < 122463006) |
           CHROM == 2 & (98178155 < POS & POS < 98178177) | CHROM == 3 & (26932699 < POS & POS < 26932721) |
           CHROM == "X" & (52011179 < POS & POS < 52011201) | CHROM == 7 & (117164797 < POS & POS < 117164819) |
           CHROM == 14 & (32396361 < POS & POS < 32396383) | CHROM == 9 & (59180736 < POS & POS < 59180758) |
           CHROM == 17 & (2120479 < POS & POS < 2120501) | CHROM == 3 & (60156435 < POS & POS < 60156457) |
           CHROM == 14 & (23708758 < POS & POS < 23708780) | CHROM == 12 & (59079545 < POS & POS < 59079567) |
           CHROM == 1 & (69601949 < POS & POS < 69601971) | CHROM == 4 & (56962314 < POS & POS < 56962336) |
           CHROM == 14 & (27083895 < POS & POS < 27083917) | CHROM == 10 & (33882830 < POS & POS < 33882852) |
           CHROM == 9 & (42683328 < POS & POS < 42683350) | CHROM == 15 & (119050505 < POS & POS < 119050527) |
           CHROM == 12 & (46613008 < POS & POS < 46613030) | CHROM == 1 & (182583347 < POS & POS < 182583369) |
           CHROM == 6 & (85663420 < POS & POS < 85663442) | CHROM == 15 & (126697576 < POS & POS < 126697598) |
           CHROM == 11 & (75412590 < POS & POS < 75412612) | CHROM == 1 & (87817240 < POS & POS < 87817262) |
           CHROM == 17 & (26858166 < POS & POS < 26858188) | CHROM == 13 & (134877233 < POS & POS < 134877255) |
           CHROM == "X" & (71315069 < POS & POS < 71315091) | CHROM == 13 & (84399412 < POS & POS < 84399434) |
           CHROM == 9 & (76710332 < POS & POS < 76710354) | CHROM == 4 & (113915288 < POS & POS < 113915310) |
           CHROM == "X" & (43642835 < POS & POS < 43642857) | CHROM == 3 & (22588384 < POS & POS < 22588406) |
           CHROM == 2 & (130864112 < POS & POS < 130864134) | CHROM == 1 & (244153529 < POS & POS < 244153551) |
           CHROM == "X" & (5297567 < POS & POS < 5297589) | CHROM == "Y" & (3944293 < POS & POS < 3944315) |
           CHROM == 15 & (101178624 < POS & POS < 101178646) | CHROM == 1 & (75147516 < POS & POS < 75147538) |
           CHROM == 2 & (42721700 < POS & POS < 42721722) | CHROM == 10 & (51477707 < POS & POS < 51477729) |
           CHROM == 2 & (118769205 < POS & POS < 118769227) | CHROM == 13 & (174649711 < POS & POS < 174649733) |
           CHROM == "X" & (94457408 < POS & POS < 94457430) | CHROM == 14 & (20180373 < POS & POS < 20180395) |
           CHROM == 1 & (52370837 < POS & POS < 52370859) | CHROM == 1 & (202941333 < POS & POS < 202941355) |
           CHROM == 5 & (41928234 < POS & POS < 41928256) | CHROM == 4 & (62333732 < POS & POS < 62333754) |
           CHROM == 8 & (137796325 < POS & POS < 137796347) | CHROM == 13 & (197326569 < POS & POS < 197326591) |
           CHROM == 9 & (62387634 < POS & POS < 62387656) | CHROM == 14 & (96236359 < POS & POS < 96236381) |
           CHROM == 4 & (41101468 < POS & POS < 41101490) | CHROM == 18 & (41621861 < POS & POS < 41621883) |
           CHROM == 9 & (107731483 < POS & POS < 107731505) | CHROM == 14 & (73006306 < POS & POS < 73006328) |
           CHROM == 3 & (46876326 < POS & POS < 46876348) | CHROM == 1 & (62991757 < POS & POS < 62991779) |
           CHROM == 13 & (82292789 < POS & POS < 82292811) | CHROM == 3 & (17468641 < POS & POS < 17468663) |
           CHROM == "X" & (72392768 < POS & POS < 72392790) | CHROM == 9 & (77426981 < POS & POS < 77427003) |
           CHROM == 1 & (93900530 < POS & POS < 93900552) | CHROM == 15 & (29503983 < POS & POS < 29504005) |
           CHROM == 6 & (139299264 < POS & POS < 139299286) | CHROM == 18 & (34071686 < POS & POS < 34071708) |
           CHROM == "X" & (63949081 < POS & POS < 63949103) | CHROM == 6 & (40438289 < POS & POS < 40438311) |
           CHROM == 1 & (52203201 < POS & POS < 52203223) | CHROM == 5 & (100222474 < POS & POS < 100222496) |
           CHROM == 13 & (145382790 < POS & POS < 145382812) | CHROM == 18 & (53511452 < POS & POS < 53511474) |
           CHROM == 8 & (67415497 < POS & POS < 67415519) | CHROM == "AEMK02000261.1" & (394016 < POS & POS < 394038) |
           CHROM == 9 & (109554902 < POS & POS < 109554924) | CHROM == 6 & (33685982 < POS & POS < 33686004) |
           CHROM == 9 & (20532618 < POS & POS < 20532640) | CHROM == 12 & (25271649 < POS & POS < 25271671) |
           CHROM == 12 & (25408492 < POS & POS < 25408514) | CHROM == 2 & (32600758 < POS & POS < 32600780) |
           CHROM == 1 & (164175676 < POS & POS < 164175698) | CHROM == 9 & (60189972 < POS & POS < 60189994) |
           CHROM == 16 & (58492001 < POS & POS < 58492023) | CHROM == 13 & (157006860 < POS & POS < 157006882) |
           CHROM == 11 & (34858585 < POS & POS < 34858607) | CHROM == 15 & (40099672 < POS & POS < 40099694) |
           CHROM == 3 & (117449662 < POS & POS < 117449684) | CHROM == 15 & (119713084 < POS & POS < 119713106) |
           CHROM == 8 & (55116877 < POS & POS < 55116899) | CHROM == 14 & (104214797 < POS & POS < 104214819) |
           CHROM == 4 & (121535956 < POS & POS < 121535978) | CHROM == 3 & (123926378 < POS & POS < 123926400) |
           CHROM == 15 & (30672400 < POS & POS < 30672422) | CHROM == 7 & (8650065 < POS & POS < 8650087) |
           CHROM == 13 & (170167843 < POS & POS < 170167865) | CHROM == 13 & (20808662 < POS & POS < 20808684) |
           CHROM == 9 & (13568870 < POS & POS < 13568892) | CHROM == 16 & (58034432 < POS & POS < 58034454) |
           CHROM == 6 & (117480327 < POS & POS < 117480349) | CHROM == 8 & (105215554 < POS & POS < 105215576) |
           CHROM == 2 & (76566458 < POS & POS < 76566480) | CHROM == 9 & (11406805 < POS & POS < 11406827) |
           CHROM == 9 & (11759438 < POS & POS < 11759460) | CHROM == 14 & (13006317 < POS & POS < 13006339) |
           CHROM == 1 & (147750500 < POS & POS < 147750522) | CHROM == 4 & (88805579 < POS & POS < 88805601) |
           CHROM == 4 & (9776721 < POS & POS < 9776743))
}

##function to extract gRNA6-2 offtargets
gRNA6_2_filter <- function(df) {
  filter(df,
         CHROM == 13 & (150485290 < POS & POS < 150485312) | CHROM == 16 & (20443576 < POS & POS < 20443598) |
           CHROM == 11 & (28997578 < POS & POS < 28997600) | CHROM == 11 & (78854769 < POS & POS < 78854791) |
           CHROM == 9 & (6314984 < POS & POS < 6315006) | CHROM == 12 & (22343582 < POS & POS < 22343604) |
           CHROM == 14 & (107862953 < POS & POS < 107862975) | CHROM == 13 & (75612002 < POS & POS < 75612024) |
           CHROM == 18 & (52343454 < POS & POS < 52343476) | CHROM == 4 & (107113143 < POS & POS < 107113165) |
           CHROM == 1 & (200728324 < POS & POS < 200728346) | CHROM == 9 & (7171086 < POS & POS < 7171108) |
           CHROM == 7 & (70568800 < POS & POS < 70568822) | CHROM == "AEMK02000598.1" & (1713433 < POS & POS < 1713455) |
           CHROM == 7 & (6141500 < POS & POS < 6141522) | CHROM == 10 & (56478114 < POS & POS < 56478136) |
           CHROM == 1 & (263663985 < POS & POS < 263664007) | CHROM == 1 & (262810902 < POS & POS < 262810924) |
           CHROM == "AEMK02000598.1" & (894275 < POS & POS < 894297) | CHROM == "AEMK02000598.1" & (181817 < POS & POS < 181839) |
           CHROM == 1 & (263022478 < POS & POS < 263022500) | CHROM == 3 & (40192690 < POS & POS < 40192712) |
           CHROM == 2 & (38569500 < POS & POS < 38569522) | CHROM == "AEMK02000510.1" & (882378 < POS & POS < 882400) |
           CHROM == "X" & (14295561 < POS & POS < 14295583) | CHROM == 16 & (50889716 < POS & POS < 50889738) |
           CHROM == 6 & (55686958 < POS & POS < 55686980) | CHROM == 1 & (140064990 < POS & POS < 140065012) |
           CHROM == 6 & (10032004 < POS & POS < 10032026) | CHROM == 17 & (18161353 < POS & POS < 18161375) |
           CHROM == 8 & (115803639 < POS & POS < 115803661) | CHROM == 11 & (43240872 < POS & POS < 43240894) |
           CHROM == 9 & (43098985 < POS & POS < 43099007) | CHROM == 1 & (62701420 < POS & POS < 62701442) |
           CHROM == 14 & (23336185 < POS & POS < 23336207) | CHROM == "X" & (40933252 < POS & POS < 40933274) |
           CHROM == 13 & (46715841 < POS & POS < 46715863) | CHROM == 1 & (109825362 < POS & POS < 109825384) |
           CHROM == 1 & (263136102 < POS & POS < 263136124) | CHROM == 8 & (68070827 < POS & POS < 68070849) |
           CHROM == "AEMK02000598.1" & (1061634 < POS & POS < 1061656) | CHROM == 15 & (41392812 < POS & POS < 41392834) |
           CHROM == 4 & (42619832 < POS & POS < 42619854) | CHROM == 15 & (102980449 < POS & POS < 102980471) |
           CHROM == 10 & (41550575 < POS & POS < 41550597) | CHROM == 15 & (29120786 < POS & POS < 29120808) |
           CHROM == 2 & (17821842 < POS & POS < 17821864) | CHROM == 3 & (51500596 < POS & POS < 51500618) |
           CHROM == 14 & (89357146 < POS & POS < 89357168) | CHROM == 4 & (8100050 < POS & POS < 8100072) |
           CHROM == 3 & (20951693 < POS & POS < 20951715) | CHROM == 5 & (30700210 < POS & POS < 30700232) |
           CHROM == 14 & (86614833 < POS & POS < 86614855) | CHROM == "AEMK02000598.1" & (851480 < POS & POS < 851502) |
           CHROM == 16 & (37469948 < POS & POS < 37469970) | CHROM == 9 & (16703939 < POS & POS < 16703961) |
           CHROM == 7 & (3617273 < POS & POS < 3617295) | CHROM == 6 & (66055353 < POS & POS < 66055375) |
           CHROM == 13 & (162352541 < POS & POS < 162352563) | CHROM == 9 & (65860317 < POS & POS < 65860339) |
           CHROM == 11 & (77812300 < POS & POS < 77812322) | CHROM == 1 & (153773562 < POS & POS < 153773584) |
           CHROM == 13 & (46340639 < POS & POS < 46340661))
}

##Extract SNVs from gRNA5-1 OTs
MJ01_5_1_OT <- gRNA5_1_filter(df= MJ01_SNV_filt)
MJ02_5_1_OT <- gRNA5_1_filter(df= MJ02_SNV_filt)
MJ03_5_1_OT <- gRNA5_1_filter(df= MJ03_SNV_filt)
##add sample column
MJ01_5_1_OT <- MJ01_5_1_OT %>% mutate(Sample = "116-1")
MJ02_5_1_OT <- MJ02_5_1_OT %>% mutate(Sample = "116-2")
MJ03_5_1_OT <- MJ03_5_1_OT %>% mutate(Sample = "Yuc:GFP")

##Extract SNVs from gRNA6-2 OTs
MJ01_6_2_OT <- gRNA6_2_filter(df= MJ01_SNV_filt)
MJ02_6_2_OT <- gRNA6_2_filter(df= MJ02_SNV_filt)
MJ03_6_2_OT <- gRNA6_2_filter(df= MJ03_SNV_filt)
##add sample column
MJ01_6_2_OT <- MJ01_6_2_OT %>% mutate(Sample = "116-1")
MJ02_6_2_OT <- MJ02_6_2_OT %>% mutate(Sample = "116-2")
MJ03_6_2_OT <- MJ03_6_2_OT %>% mutate(Sample = "Yuc:GFP")


#combine SNV OTs
SNV_5_1_combo <- bind_rows(MJ01_5_1_OT, MJ02_5_1_OT, MJ03_5_1_OT) %>% arrange(CHROM, POS) %>% select(Sample, everything())
SNV_6_2_combo <- bind_rows(MJ01_6_2_OT, MJ02_6_2_OT, MJ03_6_2_OT) %>% arrange(CHROM, POS) %>% select(Sample, everything())


#Create Excel Table of SNVs in OT 5-1
SNV_5_1_OT_summary <- createWorkbook("SNV_5_1_OT_summary")
addWorksheet(SNV_5_1_OT_summary, "MJ01_116-1_5_1_OT")
writeData(SNV_5_1_OT_summary, "MJ01_116-1_5_1_OT", MJ01_5_1_OT, rowNames = F, colNames = T)
addWorksheet(SNV_5_1_OT_summary, "MJ02_116-2_5_1_OT")
writeData(SNV_5_1_OT_summary, "MJ02_116-2_5_1_OT", MJ02_5_1_OT, rowNames = F, colNames = T)
addWorksheet(SNV_5_1_OT_summary, "MJ03_YucGFP_5_1_OT")
writeData(SNV_5_1_OT_summary, "MJ03_YucGFP_5_1_OT", MJ03_5_1_OT, rowNames = F, colNames = T)
addWorksheet(SNV_5_1_OT_summary, "All3_5_1_OT")
writeData(SNV_5_1_OT_summary, "All3_5_1_OT", SNV_5_1_combo, rowNames = F, colNames = T)
saveWorkbook(SNV_5_1_OT_summary, "SNV_5_1_OT_summary.xlsx", overwrite = T)

#Create Excel Table of SNVs in OT 6-2
SNV_6_2_OT_summary <- createWorkbook("SNV_6_2_OT_summary")
addWorksheet(SNV_6_2_OT_summary, "MJ01_116-1_6_2_OT")
writeData(SNV_6_2_OT_summary, "MJ01_116-1_6_2_OT", MJ01_6_2_OT, rowNames = F, colNames = T)
addWorksheet(SNV_6_2_OT_summary, "MJ02_116-2_6_2_OT")
writeData(SNV_6_2_OT_summary, "MJ02_116-2_6_2_OT", MJ02_6_2_OT, rowNames = F, colNames = T)
addWorksheet(SNV_6_2_OT_summary, "MJ03_YucGFP_6_2_OT")
writeData(SNV_6_2_OT_summary, "MJ03_YucGFP_6_2_OT", MJ03_6_2_OT, rowNames = F, colNames = T)
addWorksheet(SNV_6_2_OT_summary, "All3_6_2_OT")
writeData(SNV_6_2_OT_summary, "All3_6_2_OT", SNV_6_2_combo, rowNames = F, colNames = T)
saveWorkbook(SNV_6_2_OT_summary, "SNV_6_2_OT_summary.xlsx", overwrite = T)


