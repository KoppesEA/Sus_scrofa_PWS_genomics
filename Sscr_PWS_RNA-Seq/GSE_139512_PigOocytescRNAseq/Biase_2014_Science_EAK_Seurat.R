### Seurat-based (re)-analysis of Biase et al (2014, Genome Research)
### Series GSE57249 .fpkm count matrix
### https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
### 
### Erik Koppes, July 2023
### M. RW Mann Lab, MWRI

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

# Import GSE57249 matrice
GSE57249_df <- read_tsv("./GSE57249_fpkm.txt") #note row1 is gene names, cols are sample GSE
# set Ensembl Gene IDs to rownames rather than column [see awk/sort commands]
GRCm38v70_IDs <- read_tsv("./Mmus_GRCm38v70_GeneID_Ensembl.tsv",
                          col_names = c("Ensembl_ID", "Gene_Name"))
# set nonuniq Gene_Names [see sort/uniq commands]
GRCm38v70_nonuniq <- read_tsv("./Mmus_GRCm38v70_GeneID_Ensembl_nonuniq.tsv",
                          col_names = c("Gene_Name"))
# get corresponding Ensembl Ensemble_IDs (4013) that match nonuniq Gene_Names (293); mostly s7k, snords, U6 etc.
GRCm38v70_nonuniq_IDs <- GRCm38v70_IDs %>% filter(Gene_Name %in% GRCm38v70_nonuniq$Gene_Name)

# filter out nonuniq ensembl_IDs from FPKM matrix then join to get Gene_Names
GSE57249_df <- GSE57249_df %>% filter(!(ID %in% GRCm38v70_nonuniq_IDs$Ensembl_ID)) %>%
  left_join(GRCm38v70_IDs, by = join_by(ID == Ensembl_ID)) %>%
  select(-ID)

## and set column Gene_Name to rownames
GSE57249_df <- column_to_rownames(GSE57249_df, var = "Gene_Name")


# Extract metadata from two matrices
series_1 <- read_tsv("./GSE57249-GPL13112_series_matrix.txt", skip=35)[c(6,9,16),] %>% select(-1)
series_2 <- read_tsv("./GSE57249-GPL17021_series_matrix.txt", skip=35)[c(6,9,16),] %>% select(-1)
series_full <- bind_cols(series_1, series_2) %>% select(colnames(GSE57249_df)) ##reorder using Sample ID from FPKM matrix

# Bindcols and Transpose metadata; Cell names as rows and  Description as columns
series_full_t <- data.frame(t(series_full))
colnames(series_full_t) <- c("Celltype", "Stage", "Meta_ID")
rownames(series_full_t) == colnames(GSE57249_df)  ##all TRUE--goodtogo

# Create Seurat Object
GSE57249_seurat <- CreateSeuratObject(counts = GSE57249_df, meta.data = series_full_t)

# Add mito % column to data (although this is not 10x, rather bulk from physically separated cells)
GSE57249_seurat[["percent.mt"]] <- PercentageFeatureSet(GSE57249_seurat, pattern = "^MT-")
head(GSE57249_seurat@meta.data, 10) # all 0, mt transcripts were filtered

# Visualize QC metrics as a violin plot; at least one zygote should be removed due to low feature and ncounts
VlnPlot(GSE57249_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(GSE57249_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = 'Celltype')

# FeatureScatter plot to see correlation of nFeature and nCount
FeatureScatter(GSE57249_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Filter bad zygote
GSE57249_seurat@meta.data %>% filter(Celltype == "zygote") %>% arrange(desc(nCount_RNA)) #GSM1377865 sample Z7 has abberrantly low nCount and nFeature
GSE57249_seurat <- subset(GSE57249_seurat, subset = Meta_ID != "Z7")

# Normalization
GSE57249_seurat <- NormalizeData(GSE57249_seurat, normalization.method = "LogNormalize", scale.factor = 10000) #default norm param

# Find variable features
GSE57249_seurat <- FindVariableFeatures(GSE57249_seurat, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(GSE57249_seurat), 10)
top100 <- head(VariableFeatures(GSE57249_seurat), 100)



# plot variable features with and without labels
plot1 <- VariableFeaturePlot(GSE57249_seurat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Scale the data
all.genes <- rownames(GSE57249_seurat)
GSE57249_seurat <- ScaleData(GSE57249_seurat, features = all.genes)

# Run PCA
GSE57249_seurat <- RunPCA(GSE57249_seurat, features = VariableFeatures(object = GSE57249_seurat))

# Examine and visualize PCA results a few different ways
print(GSE57249_seurat[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(GSE57249_seurat, dims = 1:2, reduction = "pca")
DimPlot(GSE57249_seurat, reduction = "pca", group.by = 'Celltype')
?DimPlot()

# Heatmaps
DimHeatmap(GSE57249_seurat, dims = 1, cells = 55, balanced = TRUE) # for PC1 only
DimHeatmap(GSE57249_seurat, dims = 1:4, cells = 55, balanced = TRUE)

#Jackstrawplot
GSE57249_seurat <- JackStraw(GSE57249_seurat, num.replicate = 100)
GSE57249_seurat <- ScoreJackStraw(GSE57249_seurat, dims = 1:20)
JackStrawPlot(GSE57249_seurat, dims = 1:20) ## dropoff between 4 (p<=0.000293) distinct trajectories as compared to 5 (p=0.0411)+allclose.
ElbowPlot(GSE57249_seurat) ##pivot is on PC5, we know there are 4 cell types: 2E, 4E, Z, TE; but use >5dim

# Cluster the Cells
# FindNeighbors; dims=10
GSE57249_seurat <- FindNeighbors(GSE57249_seurat, dims = 1:10) #set to 10 dims (see comments above); elbo plot

# FindClusters; using avoce dims=1:10 and resolution=1.625 correctly predicts the 5 cell types
GSE57249_seurat <- FindClusters(GSE57249_seurat, resolution = 1.625) #see tutorial, study of blastomeres have unusually low cell count as 0.4-1.2 range optimized for 3k cells? 
head(Idents(GSE57249_seurat), 100)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
GSE57249_seurat <- RunUMAP(GSE57249_seurat, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(GSE57249_seurat, reduction = "umap")

## quick save, although with <100 cells was quick to run
saveRDS(GSE57249_seurat, file = "./output/GSE57249_seurat_v1.rds")

## rename and re-level Clusters
new_cluster_IDs <- c("2-Cell", "4-Cell", "Zygote", "ICM", "TE")
names(new_cluster_IDs) <- levels(GSE57249_seurat)
GSE57249_seurat <- RenameIdents(GSE57249_seurat, new_cluster_IDs)

# Relevel order for violin plots
levels(GSE57249_seurat) <-  c( "Zygote", "2-Cell", "4-Cell", "ICM", "TE")

## Dim plot w/new names  
DimPlot(GSE57249_seurat, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

## Find markers
# find all markers of cluster 1-4
# cluster1_markers <- FindMarkers(GSE57249_seurat, ident.1 = 1, min.pct = 0.25)
# cluster2_markers <- FindMarkers(GSE57249_seurat, ident.2 = 2, min.pct = 0.25)
# cluster3_markers <- FindMarkers(GSE57249_seurat, ident.1 = 3, min.pct = 0.25)
# cluster4_markers <- FindMarkers(GSE57249_seurat, ident.1 = 4, min.pct = 0.25)


# find markers for every cluster compared to all remaining cells, report only the positive
# Runs for ~2min on M1 macbook
GSE57249_markers <- FindAllMarkers(GSE57249_seurat, only.pos = TRUE)

## Set a df with top 250 markers from each cluster
top250_df <- GSE57249_markers  %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 250) %>%
  ungroup()

# Filter by cluster
Zygotic_df <- top250_df %>% dplyr::filter(cluster == "Zygote")
TwoCell_df <- top250_df %>% dplyr::filter(cluster == "2-Cell")
FourCell_df <- top250_df %>% dplyr::filter(cluster == "4-Cell")
ICM_df <- top250_df %>% dplyr::filter(cluster == "ICM")
TE_df <-  top250_df %>% dplyr::filter(cluster == "TE")

                                                                  
                                                                  
#Create Excel Table of Early mammalian cell division lineage markers
library(openxlsx)
Biase_scRNAseq_markergenes <- createWorkbook("Biase_scRNAseq_markergenes")
addWorksheet(Biase_scRNAseq_markergenes, "Zygote")
writeData(Biase_scRNAseq_markergenes, "Zygote", Zygotic_df, rowNames = F, colNames = T)
addWorksheet(Biase_scRNAseq_markergenes, "2-Cell")
writeData(Biase_scRNAseq_markergenes, "2-Cell", TwoCell_df, rowNames = F, colNames = T)
addWorksheet(Biase_scRNAseq_markergenes, "4-Cell")
writeData(Biase_scRNAseq_markergenes, "4-Cell", FourCell_df, rowNames = F, colNames = T)
addWorksheet(Biase_scRNAseq_markergenes, "ICM")
writeData(Biase_scRNAseq_markergenes, "ICM", ICM_df, rowNames = F, colNames = T)
addWorksheet(Biase_scRNAseq_markergenes, "TE")
writeData(Biase_scRNAseq_markergenes, "TE", TE_df, rowNames = F, colNames = T)
saveWorkbook(Biase_scRNAseq_markergenes, "Biase_scRNAseq_markergenes.xlsx", overwrite = T)


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

## References 