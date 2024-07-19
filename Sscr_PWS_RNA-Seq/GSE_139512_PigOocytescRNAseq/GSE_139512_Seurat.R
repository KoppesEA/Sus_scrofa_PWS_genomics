## GSE139512 Pig (Sus scrofa) scRNA-Seq analysis of oocyte and embryos (Kong et al., FASEB, 2019)
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
GSE139512_df <- read_tsv("./GSE139512%5Ffpkm.txt") #note col1 is gene names, next 91 cols are GSE_139512 sample ; 25414 features

## Ensembl GTF for Sscrofa 10.2.87
# set Ensembl Gene IDs to rownames rather than column [see awk/sort commands] (25322 IDs)
Sscr10.2v87_IDs <- read_tsv("./Sscr_GTF/Sus_scrofa.Sscrofa10.2.87.geneIDconv.tsv",
                          col_names = c("Ensembl_ID", "Gene_Name"))

# Find Ensembl TrackingID elements that have nonunique Gene_Name: 1928 total tracking IDs: mostly s7k, snords, U6 etc.
Sscr10.2v87_nonuniq <- Sscr10.2v87_IDs$Ensembl_ID[duplicated(Sscr10.2v87_IDs$Gene_Name, fromLast = TRUE)]
Sscr10.2v87_nonuniq_uniq <- unique(Sscr10.2v87_IDs$Gene_Name[duplicated(Sscr10.2v87_IDs$Gene_Name, fromLast = TRUE)]) #includes ZFP57
Sscr10.2v87_IDs <- Sscr10.2v87_IDs %>% mutate(Gene_Name_adj = if_else(Ensembl_ID %in% Sscr10.2v87_nonuniq, Ensembl_ID, Gene_Name))

## Remove 7SK, U snRNAs,   
remove_genes <- c(paste0("U", c(1:8)), "U6atac", "snoU6-53", "7SK", "Metazoa_SRP", "5_8S_rRNA", "5S_rRNA") #list of 14 Gene_Names
remove_tracking <- filter(Sscr10.2v87_IDs, Gene_Name %in% remove_genes) # list of 193 tracking IDs

# Join gene_names_adj to GSE139512_df
GSE139512_df <- GSE139512_df %>% left_join(Sscr10.2v87_IDs, by = join_by(tracking_id == Ensembl_ID)) %>% select(tracking_id, Gene_Name, Gene_Name_adj, everything())
GSE139512_df %>% filter(Gene_Name == "7SK") #check that repetitive gene features retain tracking_IDs (generally low counts anyway for 7SK)

# 92 additional tracking IDs starting with ERCC-; set Gene_Name_adj to tracking_ID
GSE139512_df %>% filter(is.na(Gene_Name_adj))
GSE139512_df <- GSE139512_df %>%
  mutate(Gene_Name_adj = if_else(str_detect(tracking_id, "^ERCC-\\d{5}$"), tracking_id, Gene_Name_adj))

# filter out remove_tracking genes (S7K, rRNA etc) nonuniq ensembl_IDs from count matrix then join to get Gene_Names 
GSE139512_df <- GSE139512_df %>% filter(!(tracking_id %in% remove_tracking$Ensembl_ID)) # down to 23966 observations (genes)

# remove extra columns
GSE139512_df <- GSE139512_df %>% select(-tracking_id, -Gene_Name)

## No more nonunique Gene_Name_adj found
# Find nonunique character class values in the 'values' column
nonunique_values <- GSE139512_df %>%
  group_by(Gene_Name_adj) %>%
  filter(n() > 1) %>%
  pull(Gene_Name_adj) %>%
  unique()
GSE139512_df %>% filter(is.na(Gene_Name_adj))

## and set column Gene_Name_adj to rownames
GSE139512_df <- column_to_rownames(GSE139512_df, var = "Gene_Name_adj")

# Extract metadata from CPB metadata
source("GSE139512_metaScript.R")

## However it appears that the the Development Stage variable in the meta-data is not the sample name, also there is a 2,4,8 and Morula -cell pull sample (Pools?)
## Make Metadata manuallty
GSE_139512_meta <- bind_cols(
  Sample_Name = colnames(GSE139512_df),
  Cell_Type = c(rep("Oocyte", 2), rep("1-Cell", 3), rep("2-Cell", 6), "2-Cell-Pool", rep("4-Cell", 13), "4-Cell-Pool", rep("8-Cell", 25), "8-Cell-Pool",
                rep("Morula", 30), "Morula-Pool", rep("Blast-1cell", 2), rep("Blast-ICM", 3), rep("Blast-TE",3)),
  Experiment = c(GSE_139512_SraRunTbl$Experiment[79:80], GSE_139512_SraRunTbl$Experiment[86:88], GSE_139512_SraRunTbl$Experiment[1:7], # 7 2-cell includes pooled
                 GSE_139512_SraRunTbl$Experiment[8:21], GSE_139512_SraRunTbl$Experiment[22:45], rep(NA, 2), GSE_139512_SraRunTbl$Experiment[49:78], NA, # 14 4-cell includes pooled, 24 8-cell is short an 8-cell and pool, 30 morula no pooled?
                 GSE_139512_SraRunTbl$Experiment[81:82], GSE_139512_SraRunTbl$Experiment[46:48], GSE_139512_SraRunTbl$Experiment[83:85]),
  Run = c(GSE_139512_SraRunTbl$Run[79:80], GSE_139512_SraRunTbl$Run[86:88], GSE_139512_SraRunTbl$Run[1:7], # 7 2-cell includes pooled
               GSE_139512_SraRunTbl$Run[8:21], GSE_139512_SraRunTbl$Run[22:45], rep(NA, 2), GSE_139512_SraRunTbl$Run[49:78], NA, # 14 4-cell includes pooled, 24 8-cell is short an 8-cell and pool, 30 morula no pooled?
               GSE_139512_SraRunTbl$Run[81:82], GSE_139512_SraRunTbl$Run[46:48], GSE_139512_SraRunTbl$Run[83:85]))
  
# Set series_full a meta_df with Sample_Name as rowname
series_full <- GSE_139512_meta %>% column_to_rownames(var = "Sample_Name")

all(rownames(series_full) == colnames(GSE139512_df))  ##all TRUE--goodtogo

# Create Seurat Object
GSE139512_seurat <- CreateSeuratObject(counts = GSE139512_df, meta.data = series_full)
# Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
# Warning: Data is of class data.frame. Coercing to dgCMatrix.





# Does not appear to be Mitochondrial Features or Read counts in this gencode annotated trancriptomic data
GSE139512_seurat[["percent.mt"]] <- PercentageFeatureSet(GSE139512_seurat, pattern = "^M") ## Ensembl uses "^MT-" to denote mitochondrial genes usually but is left of in this .gtf
rownames(GSE139512_seurat)[rownames(GSE139512_df) %>% str_detect("ND1")]
##this shows ND1 is present, does not have MT- or M- marking as mitochondrial so will need manual mitogene list extracted from gtf [see awk bash script]
Sscr_mito_genes <- read_tsv("./Sscr_GTF/Sus_scrofa.mito.Sscrofa10.2.87.uniq.tsv", col_names = "Gene_Name")
Sscr_mito_regex <- paste0("\\b(", paste(Sscr_mito_genes$Gene_Name, collapse = "|"), ")\\b")
# Sscr_mito_regex # "\\b(ATP6|ATP8|COX1|COX2|COX3|CYTB|ND1|ND2|ND3|ND4|ND4L|ND5|ND6|...)\\b" # the \\b: first backslash to escape second, and \b for word boundary

# Use regex to set pattern
GSE139512_seurat[["percent.mt"]] <- PercentageFeatureSet(GSE139512_seurat, pattern = Sscr_mito_regex)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.5041  8.3736 11.8782 19.0077 17.8335 98.7157 
head(GSE139512_seurat@meta.data, 10) # check metadata updated
summary(GSE139512_seurat@meta.data$percent.mt)
## Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## 0.5041  8.3736 11.8782 19.0077 17.8335 98.7157 
## mostly OK but some very high mito % should be removed
head(GSE139512_seurat@meta.data %>% arrange(desc(percent.mt)), 25) # top 25 cells with highest percent.mt, top6 have >90% mt reads; top 15 are > 20% mt Reads to remove

# Visualize nFeature, RNA count and percent mt QC metrics as a violin plot; need to remove high percent mt
VlnPlot(GSE139512_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(GSE139512_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = 'Cell_Type')

# ICM-2, 8-cell_1_5 and M-1_1 should be removed due to mt and M_4_15 removed due to low nfeature_RNA and 2-cell_1_2 with high nCountRNA
# Also want to remove and "pull" samples; these are pooled
Cells_Remove <- c("2-cell_1_2", "4-cell-3_2", "1-cell_3", "1-cell_1", "S_Bla-1", "4-cell-2_3", "ICM-2", "8-cell_1_5", "8-cell_2_3", "M-1_1", "4-cell-1_3", "TE-2")
Cell_Type_Remove <- c("2-Cell-Pool", "4-Cell-Pool", "8-Cell-Pool", "Blast-1cell", "Morula-Pool")
GSE139512_seurat_filter <- subset(GSE139512_seurat, cells = Cells_Remove, invert = TRUE)
GSE139512_seurat_filter <- subset(GSE139512_seurat_filter, subset = Cell_Type %in% Cell_Type_Remove, invert = TRUE)

# Visualize nFeature, RNA count and percent mt QC metrics as a violin plot; need to remove high percent mt
VlnPlot(GSE139512_seurat_filter, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(GSE139512_seurat_filter, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = 'Cell_Type')

#######################-----------#################

# FeatureScatter plot to see correlation of nFeature and nCount
FeatureScatter(GSE139512_seurat_filter, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
## One slight outlier with high Ncount (4-cell-4)

# Normalization
GSE139512_seurat_filter <- NormalizeData(GSE139512_seurat_filter, normalization.method = "LogNormalize", scale.factor = 10000) #default norm param

# Find variable features
GSE139512_seurat_filter <- FindVariableFeatures(GSE139512_seurat_filter, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(GSE139512_seurat_filter), 10)
top100 <- head(VariableFeatures(GSE139512_seurat_filter), 100)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(GSE139512_seurat_filter)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Scale the data
all.genes <- rownames(GSE139512_seurat_filter)
GSE139512_seurat_filter <- ScaleData(GSE139512_seurat_filter, features = all.genes)

# Run PCA
GSE139512_seurat_filter <- RunPCA(GSE139512_seurat_filter, features = VariableFeatures(object = GSE139512_seurat_filter))

# Examine and visualize PCA results a few different ways
print(GSE139512_seurat_filter[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(GSE139512_seurat_filter, dims = 1:2, reduction = "pca")
DimPlot(GSE139512_seurat_filter, reduction = "pca", group.by = 'Cell_Type')

# Heatmaps
DimHeatmap(GSE139512_seurat_filter, dims = 1, cells = 74, balanced = TRUE) # for PC1 only
DimHeatmap(GSE139512_seurat_filter, dims = 1:4, cells = 74, balanced = TRUE)

#Jackstrawplot
GSE139512_seurat_filter <- JackStraw(GSE139512_seurat_filter, num.replicate = 100)
GSE139512_seurat_filter <- ScoreJackStraw(GSE139512_seurat_filter, dims = 1:20)
JackStrawPlot(GSE139512_seurat_filter, dims = 1:20) ## PC1-8 significant, and we have 8 known Cell_Type groups
ElbowPlot(GSE139512_seurat_filter) ##pivot is on PC7, we know there are 8 cell types: Oocyte, 1-cell, 2-cell, 4-cell, 8-cell, Blast-ICM, Blast-TE, Morula, Oocyte

# Cluster the Cells
# FindNeighbors; dims=8
GSE139512_seurat_filter <- FindNeighbors(GSE139512_seurat_filter, dims = 1:8) #set to 8 dims (see comments above); elbow plot

# FindClusters; using avoce dims=1:8 and resolution=1.625 correctly predicts the 5 cell types
GSE139512_seurat_filter <- FindClusters(GSE139512_seurat_filter, resolution = 3.65) # 4 groups: Oo, 1-cell and 2-cell together; 4 + 8-vell togetherl Morula diverse ICM/TE together; not very good
head(Idents(GSE139512_seurat_filter), 100)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
GSE139512_seurat_filter <- RunUMAP(GSE139512_seurat_filter, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(GSE139512_seurat_filter, reduction = "umap")

## To Edit below for tsne with goup.by=Cell_type
# # note that you can set `label = TRUE` or use the LabelClusters function to help label
# # UMAP cluster visualization with cluster-defined, known-group and litter labels
# DimPlot(CPB_seurat_filter, reduction = "umap", group.by = 'ident') + ggtitle("UMAP Plot with cluster-defined em-stage(ZNAT/ZSO/2CNAT/2CSO) labels")
# DimPlot1 = last_plot()
# DimPlot(CPB_seurat_filter, reduction = "umap", group.by = 'em_mat') + ggtitle("UMAP Plot with group-known em-stage(ZNAT/ZSO/2CNAT/2CSO) labels")
# DimPlot2 = last_plot()
# DimPlot(CPB_seurat_filter, reduction = "umap", group.by = 'litter') + ggtitle("UMAP Plot with group-known mouse litter labels")
# DimPlot3 = last_plot()
# 
# # TSNE cluster visualization with cluster-defined, known-group and litter labels
# DimPlot(CPB_seurat_filter, reduction = "tsne",  group.by = 'ident') + ggtitle("tSNE Plot with cluster-defined em-stage(ZNAT/ZSO/2CNAT/2CSO) labels")
# DimPlot4 = last_plot()
# DimPlot(CPB_seurat_filter, reduction = "tsne",  group.by = 'em_mat') + ggtitle("tSNE Plot with group-known em-stage(ZNAT/ZSO/2CNAT/2CSO) labels")
# DimPlot5 = last_plot()
# DimPlot(CPB_seurat_filter, reduction = "tsne", group.by = 'litter') + ggtitle("tSNE Plot with mouse litter labels")
# DimPlot6 = last_plot()
# 
# # Loop through the plots and save each UMAP/TSNE Plot
# DimPlots <- paste0("DimPlot", c(1:6))
# DimPath="./CPB_cluster_UMAP_TSNE"
# DimPlotUMAP_names <- c("CPB_UMAP_clusterident.jpeg", "CPB_UMAP_knowngroup.jpeg", "CPB_UMAP_knownlitter.jpeg",
#                        "CPB_TSNE_clusterident.jpeg", "CPB_TSNE_knowngroup.jpeg", "CPB_TSNE_knownlitter.jpeg")
# for (i in 1:length(DimPlots)) {
#   plot_name = DimPlots[i]
#   plot = get(plot_name)
#   ggsave(filename = DimPlotUMAP_names[i], plot = plot, path = DimPath,
#          device = "jpeg", units = "cm", scale = 1, dpi = 300, quality = 50)
# }

#############################################
#############################################
## instead of clusters will use known groups

# Set Idents to metadata embroy-mating 4-level category
GSE139512_seurat_filter_origID <- SetIdent(GSE139512_seurat_filter, value = GSE139512_seurat_filter@meta.data$Cell_Type)
# > unique(GSE139512_seurat_filter@meta.data$Cell_Type)
# [1] "Oocyte"    "1-Cell"    "2-Cell"    "4-Cell"    "8-Cell"    "Morula"    "Blast-ICM" "Blast-TE"  

# Find All markers for every cluster compared to all remaining cells, report only the positive
# Runs for ~2min on M1 macbook
GSE139512_seurat_group_markers <- FindAllMarkers(GSE139512_seurat_filter_origID, only.pos = TRUE)

## Set a df with top 250 markers from each cluster
top250_group_df <- GSE139512_seurat_group_markers  %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 250) %>%
  ungroup()

# Filter by cluster; No markers fouhd for Oocyte, One-Cell, Blast-ICM and Blast-TE (too few samples?)
Oocyte_250 <- top250_group_df %>% dplyr::filter(cluster == "Oocyte") %>% column_to_rownames(var = "gene")
One_Cell_250 <- top250_group_df %>% dplyr::filter(cluster == "1-Cell") %>% column_to_rownames(var = "gene")
Two_Cell_250 <- top250_group_df %>% dplyr::filter(cluster == "2-Cell") %>% column_to_rownames(var = "gene")
Four_Cell_250 <- top250_group_df %>% dplyr::filter(cluster == "4-Cell") %>% column_to_rownames(var = "gene")
Eight_Cell_250 <- top250_group_df %>% dplyr::filter(cluster == "8-Cell") %>% column_to_rownames(var = "gene")
Morula_250 <- top250_group_df %>% dplyr::filter(cluster == "Morula") %>% column_to_rownames(var = "gene")
Blast_ICM_250 <- top250_group_df %>% dplyr::filter(cluster == "Blast-ICM") %>% column_to_rownames(var = "gene")
Blast_TE_250 <- top250_group_df %>% dplyr::filter(cluster == "Blast-TE") %>% column_to_rownames(var = "gene")


#Create Excel Table of Early mammalian cell division lineage markers
library(openxlsx)
GSE139512_scRNAseq_markergenes <- createWorkbook("GSE139512_scRNAseq_markergenes")
addWorksheet(GSE139512_scRNAseq_markergenes, "Oocyte")
writeData(GSE139512_scRNAseq_markergenes, "Oocyte", Oocyte_250, rowNames = T, colNames = T)
addWorksheet(GSE139512_scRNAseq_markergenes, "One_Cell")
writeData(GSE139512_scRNAseq_markergenes, "One_Cell", One_Cell_250, rowNames = T, colNames = T)
addWorksheet(GSE139512_scRNAseq_markergenes, "Two_Cell")
writeData(GSE139512_scRNAseq_markergenes, "Two_Cell", Two_Cell_250, rowNames = T, colNames = T)
addWorksheet(GSE139512_scRNAseq_markergenes, "Four_Cell")
writeData(GSE139512_scRNAseq_markergenes, "Four_Cell", Four_Cell_250, rowNames = T, colNames = T)
addWorksheet(GSE139512_scRNAseq_markergenes, "Eight_Cell")
writeData(GSE139512_scRNAseq_markergenes, "Eight_Cell", Eight_Cell_250, rowNames = T, colNames = T)
addWorksheet(GSE139512_scRNAseq_markergenes, "Morula")
writeData(GSE139512_scRNAseq_markergenes, "Morula", Morula_250, rowNames = T, colNames = T)
addWorksheet(GSE139512_scRNAseq_markergenes, "Blast_ICM")
writeData(GSE139512_scRNAseq_markergenes, "Blast_ICM", Blast_ICM_250, rowNames = T, colNames = T)
addWorksheet(GSE139512_scRNAseq_markergenes, "Blast_TE")
writeData(GSE139512_scRNAseq_markergenes, "Blast_TE", Blast_TE_250, rowNames = T, colNames = T)
saveWorkbook(GSE139512_scRNAseq_markergenes, "GSE139512_scRNAseq_markergenes.xlsx", overwrite = T)


#######################################
##Adding graphs for known groups

## Known Maternal Effect genes::
matEffect <- c("DNMT1", "ZFP57", "PADI6")
VlnPlot(GSE139512_seurat_filter_origID, features = matEffect)
ggsave("CPB_matEffect_violin.jpeg", plot = last_plot(), device = "jpeg", path = "./CPB_group_violinPlots", units = "cm",
       scale=1, dpi=300, quality=50)
FeaturePlot(GSE139512_seurat_filter_origID,features  = matEffect)
ggsave("CPB_matEffect_feature.jpeg", plot = last_plot(), device = "jpeg", path = "./CPB_group_featurePlots", units = "cm",
       scale=1, dpi=300, quality=50)

## Ref :: KHDC3L, NLRP2, NLRP4F, NLRP9A, NLRP9B, NLRP9C, NLRP7 not found
OoCortex<- c("KHDC3L", "OOEP", "PADI6", "TLE6", "ZBED3", "NLRP2", "NLRP4F",
             "NLRP5", "NLRP9A", "NLRP9B", "NLRP9C", "NLRP7")
VlnPlot(GSE139512_seurat_filter_origID, features = OoCortex)
ggsave("CPB_matEffect_OoCortex_violin.jpeg", plot = last_plot(), device = "jpeg", path = "./CPB_group_violinPlots", units = "cm",
       scale=1, dpi=300, quality=50)


## Ref   [Nlxt2, Fb, Pol2f, Uch3, Uch4 not found mm39]
Biase_bivalent <- c("BBS4", "NLXT2",
                    "PDLIM7", "SCAMP4",
                    "TCF3","BHMT22",
                    "TRP53","ZFP951",
                    "ASPM", "PARP12",
                    "LGR4",  "FB",
                    "CDK1", "POL2F",
                    "UCH3", "UCH4")
VlnPlot(GSE139512_seurat_filter_origID, features = Biase_bivalent)
ggsave("CPB_bivalent_violin.jpeg", plot = last_plot(), device = "jpeg", path = "./CPB_group_violinPlots", units = "cm",
       scale=1, dpi=300, quality=50)

## Nup genes; Elys=Ahctf1; ; NUP37, NUP50, NUP96, NUP107, NUP160, SEC13 not found
NupGenes <- c("NUP37", "NUP43", "NUP50", "NUP62", "NUP85", "NUP93", "NUP96", "NUP98", "NUP107", 
              "NUP133", "NUP153", "NUP160", "NUP214", "SEC13", "SEH1L", "AHCTF1")
VlnPlot(GSE139512_seurat_filter_origID, features = NupGenes)
ggsave("CPB_NupGenes_violin.jpeg", plot = last_plot(), device = "jpeg", path = "./CPB_group_violinPlots", units = "cm",
       scale=1, dpi=300, quality=50)
FeaturePlot(GSE139512_seurat_filter_origID,features  = NupGenes)
ggsave("CPB_NupGenes_features.jpeg", plot = last_plot(), device = "jpeg", path = "./CPB_group_featurePlots", units = "cm",
       scale=1, dpi=300, quality=50)

## Imprint TFs; Only ZSCAN4 found
imprintTF <- c("YY1", "NRF1", "ZSCAN4", "ZFP445", "ZFP274", "ZSCAN4C", "ZSCAN4D", "ZSCAN4F")
VlnPlot(GSE139512_seurat_filter_origID, features = imprintTF)
ggsave("CPB_imprintTFs_violin.jpeg", plot = last_plot(), device = "jpeg", path = "./CPB_group_violinPlots", units = "cm",
       scale=1, dpi=300, quality=50)

## Zscan4 isoforms
Zscan4_isoforms <- GRCm39gencode_IDs %>% filter(grepl("Zscan4", Gene_Name)) %>% pull(Gene_Name)
VlnPlot(GSE139512_seurat_filter_origID, features = Zscan4_isoforms)
ggsave("CPB_Zscan4_violin.jpeg", plot = last_plot(), device = "jpeg", path = "./CPB_group_violinPlots", units = "cm",
       scale=1, dpi=300, quality=50)

## DNA Methylation; DNMT3L, TET1, UHRF2, TRIM28
DNAmethylation <- c("DNMT1", "DNMT3A", "DNMT3B", "DNMT3L",
                    "TET1", "TET2", "TET3", "SETDB1",
                    "UHRF1", "UHRF2", "CDCA7", "HELLS",
                    "ZFP57", "TRIM28", "DMAP1", "ZBTB24")
VlnPlot(GSE139512_seurat_filter_origID, features  = DNAmethylation)
ggsave("CPB_DNAme_violin.jpeg", plot = last_plot(), device = "jpeg", path = "./CPB_group_violinPlots", units = "cm",
       scale=1, dpi=300, quality=50)

#Ehmt2=G9a  #Ls1 not found
HistoneModifiers <- c("SETDB1", "SIN3A", "LSD1", "SUV39H1", "EHMT2", "HDAC4", "HDAC6")
VlnPlot(GSE139512_seurat_filter_origID, features  = HistoneModifiers)
ggsave("CPB_Histone_violin.jpeg", plot = last_plot(), device = "jpeg", path = "./CPB_group_violinPlots", units = "cm",
       scale=1, dpi=300, quality=50)

MoreofInterest <- c("ZFP688", "GADD45A", "GADD45B",
                    "STAT3", "SIRT1", "SIRT3")
VlnPlot(GSE139512_seurat_filter_origID, features  = MoreofInterest)
ggsave("CPB_moreofInterest_violin.jpeg", plot = last_plot(), device = "jpeg", path = "./CPB_group_violinPlots", units = "cm",
       scale=1, dpi=300, quality=50)

# Lineage Marker Definitions:: 
ICM_marker <- c("POU5F1", "NANOG", "CARM1", "KLF4", "OTX2", "STAT3", "NR5A2", "TDGF1")
TE_marker <- c("CDX2", "EOMES", "SOX2", "GATA3", "PAR3", "PARD6B", "ELF5", "KRT7",
               "KRT8", "COLA4", "CDH1", "ECAD", "GH", "TEAD4", "TFEB", "ITGB5")
XEN_marker <- c("GATA6", "GATA4", "INS", "IGF2R", "MEST")
Zygote <- c("ZP1", "ZP2", "ZP3", "PADI6", "OOEP", "DNMT1", "DNMT3A")
Oocyte <- c("ZP1", "ZP2", "ZP3", "PADI6", "OOEP", "DNMT1", "DNMT3A")

## Lineage Marker ViolinPlots
VlnPlot(GSE139512_seurat_filter_origID,features  = ICM_marker)
ggsave("CPB_ICM_violin.jpeg", plot = last_plot(), device = "jpeg", path = "./CPB_group_violinPlots", units = "cm",
       scale=1, dpi=300, quality=50)
VlnPlot(GSE139512_seurat_filter_origID,features  = TE_marker) #Par3, Cola4, Ecad, GH not found
ggsave("CPB_TE_violin.jpeg", plot = last_plot(), device = "jpeg", path = "./CPB_group_violinPlots", units = "cm",
       scale=1, dpi=300, quality=50)
VlnPlot(GSE139512_seurat_filter_origID,features  = XEN_marker)
ggsave("CPB_XEN_violin.jpeg", plot = last_plot(), device = "jpeg", path = "./CPB_group_violinPlots", units = "cm",
       scale=1, dpi=300, quality=50)
VlnPlot(GSE139512_seurat_filter_origID,features  = Zygote) #Zp1 not found
VlnPlot(GSE139512_seurat_filter_origID,features  = Oocyte) #Zp1 not found; same as Zygote

## Mellissa List from R01
MelList1 <- c("ZFP57", "PADI6", "BIRC3", "EIF4E1B",
              "RSPO2", "GJA4", "BMP15", "MLLT3",
              "DAGLB", "HTRA4", "PHC1", "NPM2")
MelList2 <- c("ESR2", "GALM", "UHRF1", "B4GALT4",
              "OAS1H", "PCGF1", "RALBP1", "MEIS2",
              "NUP214", "TXNIP", "AKAP17B", "NLRP5")

VlnPlot(GSE139512_seurat_filter_origID,features  = MelList1)
ggsave("CPB_Mann1_violin.jpeg", plot = last_plot(), device = "jpeg", path = "./CPB_group_violinPlots", units = "cm",
       scale=1, dpi=300, quality=50)

VlnPlot(GSE139512_seurat_filter_origID,features  = MelList2)
ggsave("CPB_Mann2_violin.jpeg", plot = last_plot(), device = "jpeg", path = "./CPB_group_violinPlots", units = "cm",
       scale=1, dpi=300, quality=50)

#################### PWS genes ####################
PWS_genes <- read_tsv("./Sscr_GTF/Ssus_PWS_Ensembl_geneIDs.tsv", col_names = "PWS_Gene")
VlnPlot(GSE139512_seurat_filter_origID,features  = PWS_genes$PWS_Gene) #ENSSSCG00000018361, ENSSSCG00000027296 not found; expression detected only in UBE3A, MKRN3, SNRPN, ENSSSCG00000004834, ENSSSCG00000026997
# ENSSSCG00000004834 + ENSSSCG00000026997 == Ndn or pseuedogenes near Mkrn3/Magel2 chr1
#