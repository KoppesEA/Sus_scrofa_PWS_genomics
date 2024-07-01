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
# set Ensembl Gene IDs to rownames rather than column [see awk/sort commands] (25322 IDs)
Sscr10.2v87_IDs <- read_tsv("./Sscr_GTF/Sus_scrofa.Sscrofa10.2.87.tsv",
                          col_names = c("Ensembl_ID", "Gene_Name"))

# set nonuniq Gene_Names [see sort/uniq commands] (197 genes)
Sscr10.2v87_nonuniq <- read_tsv("./Sscr_GTF/Sus_scrofa.Sscrofa10.2.87_nonuniq.tsv",
                              col_names = c("Gene_Name"))
# remove nonuniq Gene_names that are "Ensembl", "ensembl_havana", "havana" or "insdc" as gene_name not given in v87 GTF (194 genes)
Sscr10.2v87_nonuniq <- Sscr10.2v87_nonuniq %>% filter(!(Gene_Name %in% c("ensembl", "ensembl_havana", "havana", "insdc")))

# get corresponding Ensembl Ensemble_IDs (4013) that match nonuniq Gene_Names (293); mostly s7k, snords, U6 etc. (2121 IDs)
Sscr10.2v87_nonuniq_IDs <- Sscr10.2v87_IDs %>% filter(Gene_Name %in% Sscr10.2v87_nonuniq$Gene_Name)

# filter out nonuniq ensembl_IDs from count matrix then join to get Gene_Names 
GSE139512_df <- GSE139512_df %>% filter(!(tracking_id %in% Sscr10.2v87_nonuniq_IDs$Ensembl_ID)) %>%
  left_join(Sscr10.2v87_IDs, by = join_by(tracking_id == Ensembl_ID))

# make new Gene_annot column to fill in gene_IDs that are "ensembl" etc. with the gene_ID instead
GSE139512_df <- GSE139512_df %>%
  mutate(Gene_annot = if_else(Gene_Name %in% c("ensembl", "ensembl_havana", "havana", "insdc") | is.na(Gene_Name), tracking_id, Gene_Name)) %>%
  select(-tracking_id, -Gene_Name, Gene_annot)

## No more nonunique Gene_annot found
# Find nonunique character class values in the 'values' column
nonunique_values <- GSE139512_df %>%
  group_by(Gene_annot) %>%
  filter(n() > 1) %>%
  pull(Gene_annot) %>%
  unique()

## and set column Gene_annot to rownames
GSE139512_df <- column_to_rownames(GSE139512_df, var = "Gene_annot")

# Extract metadata from CPB metadata
source("GSE139512_metaScript.R")

## However it appears that the the Devolpment Stage variable in the meta-data is not the sample name, also there is a 2,4,8 and Morula -cell pull sample (Pools?)
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
GSE57249_seurat <- CreateSeuratObject(counts = GSE139512_df, meta.data = series_full)

# Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
# Warning: Data is of class data.frame. Coercing to dgCMatrix.

# Does not appear to be Mitochondrial Features or Read counts in this gencode annotated trancriptomic data
GSE57249_seurat[["percent.mt"]] <- PercentageFeatureSet(GSE57249_seurat, pattern = "^MT-1") ## Gencode uses "^MT-" to denote mitochondrial genes
rownames(GSE57249_seurat)[rownames(GSE139512_df) %>% str_detect("MT-1")] #or possibly just not annotated by gene_name in Sscr 10.2 v87
# character(0)
head(GSE57249_seurat@meta.data, 10) # unfiltered but <10% in most cases
summary(GSE57249_seurat@meta.data$percent.mt)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0       0       0       0       0       0

#######################-----------#################

# Visualize QC metrics as a violin plot; at least one zygote should be removed due to low feature and ncounts
VlnPlot(GSE57249_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(GSE57249_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = 'Cell_Type')

# FeatureScatter plot to see correlation of nFeature and nCount
FeatureScatter(GSE57249_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
## perhaps need to remove pull/pool samples or combine into cells
##m Extreme outlier in 2-cell_1_2 nCount_RNA

# # Filter bad zygote --to revise filtering here for bad 2-cell
# GSE57249_seurat@meta.data %>% filter(Cell_Type == "zygote") %>% arrange(desc(nCount_RNA)) #GSM1377865 sample Z7 has abberrantly low nCount and nFeature
GSE57249_seurat_filter <- subset(GSE57249_seurat, subset = Experiment != "SRX7066509")

# Idents(GSE57249_seurat)  %>% grepl(pattern = "2-cell")
# identities <- Idents(GSE57249_seurat)
# identities <- as.character(identities)
# print(any(identities == "2-cell_1_2"))
##easiest just to addback a full label as orig.ident to metadata....

# FeatureScatter plot of filtered to see correlation of nFeature and nCount
FeatureScatter(GSE57249_seurat_filter, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
subset(GSE57249_seurat_filter, subset, Idents %in% c("1-cell","2-Cell", "4-Cell")) %>%  FeatureScatter(feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Normalization
CPB_seurat <- NormalizeData(CPB_seurat, normalization.method = "LogNormalize", scale.factor = 10000) #default norm param

# Find variable features
CPB_seurat <- FindVariableFeatures(CPB_seurat, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(CPB_seurat), 10)
top100 <- head(VariableFeatures(CPB_seurat), 100)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(CPB_seurat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Scale the data
all.genes <- rownames(CPB_seurat)
CPB_seurat <- ScaleData(CPB_seurat, features = all.genes)

# Run PCA
CPB_seurat <- RunPCA(CPB_seurat, features = VariableFeatures(object = CPB_seurat))

# Examine and visualize PCA results a few different ways
print(CPB_seurat[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(CPB_seurat, dims = 1:2, reduction = "pca")
DimPlot(CPB_seurat, reduction = "pca", group.by = 'emstage')
DimPlot(CPB_seurat, reduction = "pca", group.by = 'mating')

# Heatmaps
DimHeatmap(CPB_seurat, dims = 1, cells = 55, balanced = TRUE) # for PC1 only
DimHeatmap(CPB_seurat, dims = 1:4, cells = 55, balanced = TRUE)

#Jackstrawplot
CPB_seurat <- JackStraw(CPB_seurat, num.replicate = 100)
CPB_seurat <- ScoreJackStraw(CPB_seurat, dims = 1:20)
JackStrawPlot(CPB_seurat, dims = 1:20) ## PC1-4 significant, PC5 less so, PC6 more so,
ElbowPlot(CPB_seurat) ##pivot is on PC4, we know there are 4 cell types: 2C-NAT, 2C-SO, Z-NAT, Z-SO

# Cluster the Cells
# FindNeighbors; dims=10
CPB_seurat <- FindNeighbors(CPB_seurat, dims = 1:10) #set to 10 dims (see comments above); elbo plot

# FindClusters; using avoce dims=1:10 and resolution=1.625 correctly predicts the 5 cell types
CPB_seurat <- FindClusters(CPB_seurat, resolution = 1.375) #1.3 start to see Z_NAT vs Z_SO diff
head(Idents(CPB_seurat), 100)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
CPB_seurat <- RunUMAP(CPB_seurat, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(CPB_seurat, reduction = "umap")

## quick save, although with <100 cells was quick to run
saveRDS(CPB_seurat, file = "./CPB_seurat_seurat_v1.rds")

## rename and re-level Clusters
new_cluster_IDs <- c("2CSO", "ZNAT", "ZSO", "2CNAT")
names(new_cluster_IDs) <- levels(CPB_seurat)
CPB_seurat <- RenameIdents(CPB_seurat, new_cluster_IDs)

# Relevel order for violin plots and FindMarkers
levels(CPB_seurat) <-  c("ZNAT", "ZSO", "2CNAT", "2CSO")

## Dim plot w/new names  
DimPlot(CPB_seurat, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

## Find markers
# find  markers of SO vs NAT Zygote
cluster_ZSO_markers <- FindMarkers(CPB_seurat, ident.1 = "ZSO", ident.2 = "ZNAT", min.pct = 0.25)
cluster_2CSO_markers <- FindMarkers(CPB_seurat, ident.1 = "2CSO", ident.2 = "2CNAT", min.pct = 0.25)

# find markers of ZNAT (cmpared to 2CNAT) and 2CNAT (compared to ZNAT) for cell lineaege markers; should be reciprocal
cluster_ZNAT_markers <- FindMarkers(CPB_seurat, ident.1 = "ZNAT", ident.2 = "2CNAT", min.pct = 0.25)
cluster_2CNAT_markers <- FindMarkers(CPB_seurat, ident.1 = "2CNAT", ident.2 = "ZNAT", min.pct = 0.25)

## Filter each for p_val_adj < 0.05
cluster_ZSO_markers_padj0.05 <- cluster_ZSO_markers[cluster_ZSO_markers$p_val_adj < .05,]
cluster_2CSO_markers_padj0.05 <- cluster_2CSO_markers[cluster_2CSO_markers$p_val_adj < .05,]
cluster_ZNAT_markers_padj0.05 <- cluster_ZNAT_markers[cluster_ZNAT_markers$p_val_adj < .05,]
cluster_2CNAT_markers_padj0.05 <- cluster_2CNAT_markers[cluster_2CNAT_markers$p_val_adj < .05,]

# find markers for every cluster compared to all remaining cells, report only the positive
# Runs for ~2min on M1 macbook
CPB_seurat_markers <- FindAllMarkers(CPB_seurat, only.pos = TRUE)

## Set a df with top 250 markers from each cluster
top250_df <- CPB_seurat_markers  %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 250) %>%
  ungroup()

# Filter by cluster
ZNAT_250 <- top250_df %>% dplyr::filter(cluster == "ZNAT")
ZSO_250 <- top250_df %>% dplyr::filter(cluster == "ZSO")
TwoCNAT_250 <- top250_df %>% dplyr::filter(cluster == "2CNAT")
TwoCSO_250 <- top250_df %>% dplyr::filter(cluster == "2CSO")

#Create Excel Table of Early mammalian cell division lineage markers
library(openxlsx)
CPB_scRNAseq_markergenes <- createWorkbook("CPB_scRNAseq_markergenes")
addWorksheet(CPB_scRNAseq_markergenes, "ZSOvZNAT")
writeData(CPB_scRNAseq_markergenes, "ZSOvZNAT", cluster_ZSO_markers_padj0.05, rowNames = F, colNames = T)
addWorksheet(CPB_scRNAseq_markergenes, "2CSOv2CNAT")
writeData(CPB_scRNAseq_markergenes, "2CSOv2CNAT", cluster_2CSO_markers_padj0.05, rowNames = F, colNames = T)
addWorksheet(CPB_scRNAseq_markergenes, "ZNATv2CNAT")
writeData(CPB_scRNAseq_markergenes, "ZNATv2CNAT", cluster_ZNAT_markers_padj0.05, rowNames = F, colNames = T)
addWorksheet(CPB_scRNAseq_markergenes, "2CNATvZNAT")
writeData(CPB_scRNAseq_markergenes, "2CNATvZNAT", cluster_2CNAT_markers_padj0.05, rowNames = F, colNames = T)
## to add cluster vs all top 250 lists
addWorksheet(CPB_scRNAseq_markergenes, "TE")
writeData(CPB_scRNAseq_markergenes, "TE", TE_df, rowNames = F, colNames = T)
saveWorkbook(CPB_scRNAseq_markergenes, "CPB_scRNAseq_markergeness.xlsx", overwrite = T)


## Known Maternal Effect genes::
matEffect <- c("Dnmt1", "Zfp57", "Padi6")
VlnPlot(GSE57249_seurat, features = matEffect)
ggsave("Biase_matEffect_violin.jpeg", plot = last_plot(), device = "jpeg", path = "./Biase_scRNAseq_violinPlots", units = "cm",
       scale=1, dpi=300, quality=50)
FeaturePlot(GSE57249_seurat,features  = matEffect)
ggsave("Biase_matEffect_feature.jpeg", plot = last_plot(), device = "jpeg", path = "./Biase_scRNAseq_featurePlots", units = "cm",
       scale=1, dpi=300, quality=50)


## Ref :: Khdc3l and Nlrp7 missing
OoCortex<- c("Khdc3l", "Ooep", "Padi6", "Tle6", "Zbed3", "Nlrp2", "Nlrp4f",
             "Nlrp5", "Nlrp9a", "Nlrp9b", "Nlrp9c", "Nlrp7")
VlnPlot(GSE57249_seurat, features = OoCortex)
ggsave("Biase_OoCortex_violin.jpeg", plot = last_plot(), device = "jpeg", path = "./Biase_scRNAseq_violinPlots", units = "cm",
       scale=1, dpi=300, quality=50)


## Ref   [Nlxt2, Fb, Pol2f, Uch3, Uch4 not found mm9]
Biase_bivalent <- c("Bbs4", "Nlxt2",
                    "Pdlim7", "Scamp4",
                    "Tcf3","Bhmt2",
                    "Trp53","Zfp951",
                    "Aspm", "Parp12",
                    "Lgr4",  "Fb",
                    "Cdk1", "Pol2f",
                    "Uch3", "Uch4")
VlnPlot(GSE57249_seurat, features = Biase_bivalent)
ggsave("Biase_bivalent_violin.jpeg", plot = last_plot(), device = "jpeg", path = "./Biase_scRNAseq_violinPlots", units = "cm",
       scale=1, dpi=300, quality=50)

## Nup genes; Elys=Ahctf1; 
NupGenes <- c("Nup37", "Nup43","Nup50", "Nup62", "Nup85",  "Nup93", "Nup96",
              "Nup98", "Nup107", "Nup133", "Nup153", "Nup160", "Nup214", "Sec13",
              "Seh1l", "Ahctf1")
VlnPlot(GSE57249_seurat, features = NupGenes)
ggsave("Biase_NupGenes_violin.jpeg", plot = last_plot(), device = "jpeg", path = "./Biase_scRNAseq_violinPlots", units = "cm",
       scale=1, dpi=300, quality=50)
FeaturePlot(GSE57249_seurat,features  = NupGenes)
ggsave("Biase_NupGenes_features.jpeg", plot = last_plot(), device = "jpeg", path = "./Biase_scRNAseq_violinPlots", units = "cm",
       scale=1, dpi=300, quality=50)

## Imprint TFs #Zfp274 not found
imprintTF <- c("Yy1", "Nrf1", "Zfp445", "Zfp274", "Zscan4c", "Zscan4d", "Zscan4f")
VlnPlot(GSE57249_seurat, features = imprintTF)
ggsave("Biase_imprintTFs_violin.jpeg", plot = last_plot(), device = "jpeg", path = "./Biase_scRNAseq_violinPlots", units = "cm",
       scale=1, dpi=300, quality=50)

#Tcl1a=Tcl1
## The following requested variables were not found: Dab2, Cdc73 [##check with reference build and GTF gene name]
Camilo_of_interest <- c("Pdia3", "Tcl1", "Ooep", "Rras2",
                        "Snrpf", "Dab2", "Cops5", "Sbds",
                        "Pelo", "Cdc73", "Emg1", "Ctnna1")
VlnPlot(GSE57249_seurat,features  = Camilo_of_interest)
ggsave("Biase_Camillo_violin.jpeg", plot = last_plot(), device = "jpeg", path = "./Biase_scRNAseq_violinPlots", units = "cm",
       scale=1, dpi=300, quality=50)

## DNA Methylation
DNAmethylation <- c("Dnmt1", "Dnmt3a", "Dnmt3b", "Dnmt3l",
                    "Tet1", "Tet2", "Tet3", "Setdb1",
                    "Uhrf1", "Uhrf2", "Cdca7", "Hells",
                    "Zfp57", "Trim28", "Dmap1", "Zbtb24")
VlnPlot(GSE57249_seurat,features  = DNAmethylation)
ggsave("Biase_DNAme_violin.jpeg", plot = last_plot(), device = "jpeg", path = "./Biase_scRNAseq_violinPlots", units = "cm",
       scale=1, dpi=300, quality=50)

#Ehmt2=G9a 
HistoneModifiers <- c("Setdb1", "Sin3a", "Lsd1", "Suv39h1", "Ehmt2", "Hdac4", "Hdac6") #tooadd ask MRW
VlnPlot(GSE57249_seurat,features  = HistoneModifiers)
ggsave("Biase_Histone_violin.jpeg", plot = last_plot(), device = "jpeg", path = "./Biase_scRNAseq_violinPlots", units = "cm",
       scale=1, dpi=300, quality=50)


MoreofInterest <- c("Zfp688", "Gadd45a", "Gadd45b",
                    "Stat3", "Sirt1", "Sirt3")
VlnPlot(GSE57249_seurat,features  = MoreofInterest)
ggsave("Biase_moreofInterest_violin.jpeg", plot = last_plot(), device = "jpeg", path = "./Biase_scRNAseq_violinPlots", units = "cm",
       scale=1, dpi=300, quality=50)


# Lineage Marker Definitions:: 
ICM_marker <- c("Pou5f1", "Nanog", "Carm1", "Klf4", "Otx2", "Stat3", "Nr5a2", "Tdgf1") #Oct4=Pou5f1 ZGA <-
TE_marker <- c("Cdx2", "Eomes", "Sox2", "Gata3", "Par3", "Pard6b", "Elf5", "Krt7", #Par3, Cola4, E-cadherin, cGH, Dab2 not found; Igc2 removed
               "Krt8", "Cola4", "Cdh1", "Ecad", "GH", "Tead4", "Tfeb", "Itgb5")
XEN_marker <- c("Gata6", "Gata4", "Ins1", "Ins2", "Igf2r", "Mest")
Zygote <- c("Zp1", "Zp2", "Zp3", "Padi6", "Ooep", "Dnmt1", "Dnmt3a")
# Oocyte <- c("Zp1", "Zp2", "Zp3", "Padi6", "Ooep" "Dnmt1", "Dnmt3a")

## Lineage Marker ViolinPlots
VlnPlot(GSE57249_seurat,features  = ICM_marker)
ggsave("Biase_ICM_violin.jpeg", plot = last_plot(), device = "jpeg", path = "./Biase_scRNAseq_violinPlots", units = "cm",
       scale=1, dpi=300, quality=50)
VlnPlot(GSE57249_seurat,features  = TE_marker) #Par3, Cola4, Ecad, GH not found
ggsave("Biase_TE_violin.jpeg", plot = last_plot(), device = "jpeg", path = "./Biase_scRNAseq_violinPlots", units = "cm",
       scale=1, dpi=300, quality=50)
VlnPlot(GSE57249_seurat,features  = XEN_marker)
ggsave("Biase_XEN_violin.jpeg", plot = last_plot(), device = "jpeg", path = "./Biase_scRNAseq_violinPlots", units = "cm",
       scale=1, dpi=300, quality=50)

## Mellissa List from R01
MelList1 <- c("Zfp57", "Padi6", "Birc3", "Eif4e1b",
              "Rspo2", "Gja4", "Bmp15", "Mllt3",
              "Daglb", "Htra4", "Phc1", "Npm2")
MelList2 <- c("Esr2", "Galm", "Uhrf1", "B4galt4",
              "Oas1h", "Pcgf1", "Ralbp1", "Meis2",
              "Nup214", "Txnip", "Akap17b", "Nlrp5")
VlnPlot(GSE57249_seurat,features  = MelList1)
ggsave("Biase_Mann1_violin.jpeg", plot = last_plot(), device = "jpeg", path = "./Biase_scRNAseq_violinPlots", units = "cm",
       scale=1, dpi=300, quality=50)

VlnPlot(GSE57249_seurat,features  = MelList2)
ggsave("Biase_Mann2_violin.jpeg", plot = last_plot(), device = "jpeg", path = "./Biase_scRNAseq_violinPlots", units = "cm",
       scale=1, dpi=300, quality=50)

# Jiang 2012 Cell Research; Table S1
JiZscan4 <- c( "Zscan4f", "Npm1", "Npm2", "Dyrk2", "Kdm1b", "Rad51", "Zar1", "Dppa3", "Ezh2", "Pms2")
VlnPlot(GSE57249_seurat,features  = JiZscan4)
ggsave("Biase_Zscan4Ji_violin.jpeg", plot = last_plot(), device = "jpeg", path = "./Biase_scRNAseq_violinPlots", units = "cm",
       scale=1, dpi=300, quality=50)
