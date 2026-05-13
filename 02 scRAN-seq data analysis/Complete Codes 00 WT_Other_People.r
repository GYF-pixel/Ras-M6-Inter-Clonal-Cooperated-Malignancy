#Step 3 整合不同Genotype的数据
library(devtools)
library(DropletUtils)
library(Seurat)
library(scran)
library(scDblFinder)
library(scater)
library(DoubletFinder)
library(scDblFinder)
library(ggplot2)
library(ggsci)
library(cowplot)
library(tidyverse)
library(ggunchull)
library(SCENIC)
library(scales)
library(scCustomize) # 需要Seurat版本4.3.0
library(viridis)
library(RColorBrewer)
library(gridExtra)
set.seed(1234)

getwd()
setwd("I:/M6")
setwd("H:/M6")
####################################################Bravo González-Blas et al., 2019;
#The count matrix can be downloaded from GEO (GSE141589).
library(data.table)
count.matrix <- suppressWarnings(data.frame(fread('00WT_Bravo_EAD//GSM4209275_10X_SCRNASEQ_WT_EA_AGGR.tsv', sep='\t', verbose=F), row.names=1))

#Dimensionality reduction and cluster annotation - Figure 1b
library(Seurat)
# Initialize object
EAdisc <- CreateSeuratObject(counts = count.matrix, min.cells = 1, min.features = 1, project = "EAdisc")
# Normalize data
EAdisc <- NormalizeData(object = EAdisc, normalization.method = "LogNormalize", scale.factor = 1e4)
# Find variable features
EAdisc <- FindVariableFeatures(object = EAdisc, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))
# Scale data
EAdisc <- ScaleData(object = EAdisc, features = rownames(x = EAdisc), vars.to.regress = c("nCount_RNA"))
# Run PCA
EAdisc <- RunPCA(object = EAdisc, features = VariableFeatures(object = EAdisc), verbose = FALSE, npcs=150)

# Select number of PCs
source("00WT_Bravo_EAD//Seurat_Utils.R")
data.use <- PrepDR(object = EAdisc, genes.use = VariableFeatures(object = EAdisc), use.imputed = F, assay.type = "RNA")
nPC <- PCA_estimate_nPC(data.use, whereto="00WT_Bravo_EAD//nPC_selection.Rds", to.nPC = 150)

# Find clusters
EAdisc <- FindNeighbors(object = EAdisc, dims = 1:nPC, k.param=30)
EAdisc <- FindClusters(object = EAdisc, resolution = 1.2)
# Generate t-SNE
EAdisc <- RunTSNE(object = EAdisc, dims = 1:nPC)
# Rename IDs
new.cluster.ids <- c('AMF_prog_prec', 'PM_medial', 'Hemocytes', 'Antenna_A3_Arista', 'PMF_PR_Late/CC', 'Head_vertex', 'PMF_PR_Early', 'PMF_Interommatidial', 'Antenna_A1', 'PM_lateral', 'Glia', 'Antenna_A2', 'Brain_A', 'MF_Morphogenetic_Furrow', 'Brain_B', 'Brain_C', 'twi_cells')
names(x = new.cluster.ids) <- levels(x = EAdisc)
EAdisc <- RenameIdents(object = EAdisc, new.cluster.ids)
# Set colors
colors <- c("#FFC29C", "#01E6B3", "#01FAF7", "#F98E19", "#02B4E2", "#6CA7FF", "#01D3F8", "#67EAFF", "#85BC05", "#BB90FE", "#E5D006", "#59F889", "#FD79A3", "#9EF75F", "#5CB8F7", "#B5DEFF", "#F56DE0")
#col5 <- colorRampPalette((pal_npg(palette = c("nrc"))(7)))(17)
names(colors) <- levels(x = EAdisc)
# Save object
saveRDS(EAdisc, file='00WT_Bravo_EAD//Bravo_10X_SeuratObject.Rds')

#Load
EAdisc <- readRDS('00WT_Bravo_EAD//Bravo_10X_SeuratObject.Rds')
#Plot
library(Seurat)
DimPlot(object = EAdisc, reduction = 'tsne', cols=colors, label=TRUE, label.size = 3.5, repel=TRUE) + NoLegend() + NoAxes()

##提取其中的Hemocytes
WT_Bravo_Hemocytes <- subset(x = EAdisc, ident = "Hemocytes")
saveRDS(WT_Bravo_Hemocytes, file='00WT_Bravo_EAD//WT_Bravo_Hemocytes.Rds')

###############################################
# Get eye disc data and filter
library(SCopeLoomR)
loom <- open_loom("H:\\M6\\00WT_Ariss_EAD\\EAD_Ariss_WT_Seurat_SCENIC.loom")
count.matrix <- get_dgem(loom)
cell.data <- get_cellAnnotation(loom)[,2,drop=FALSE]
EyeDisc_Ariss_integration <- CreateSeuratObject(counts = count.matrix, min.cells = 3, min.features = 1, project = "Frolov - WT", meta.data=cell.data)
EyeDisc_Ariss_integration <- subset(x = EyeDisc_Ariss_integration, nCount_RNA > 1000 & nFeature_RNA > 500)
rm(loom)
# Remove brain cells for integration
EAdisc_integration <- EAdisc
EAdisc_integration <- subset(x = EAdisc, idents = c("Brain_A","Brain_B","Brain_C"),invert = TRUE)

# Normalize and find variable features in each set
EyeDisc_Ariss_integration <-NormalizeData(object = EyeDisc_Ariss_integration, verbose = FALSE)
EyeDisc_Ariss_integration <- FindVariableFeatures(object = EyeDisc_Ariss_integration, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
EAdisc_integration <-NormalizeData(object = EAdisc_integration, verbose = FALSE)
EAdisc_integration <- FindVariableFeatures(object = EAdisc_integration, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
# Determine anchors and transfer labels
EAdisc.anchors <- FindTransferAnchors(reference = EyeDisc_Ariss_integration, query = EAdisc_integration, 
    dims = 1:30)
Transfer_labels <- TransferData(anchorset = EAdisc.anchors, refdata = EyeDisc_Ariss_integration$Ariss.labels, 
    dims = 1:30)
# Add transferred labels as metadata
Ariss_labels <- as.vector(unlist(EAdisc@active.ident))
names(Ariss_labels) <- names(EAdisc@active.ident)
Ariss_labels[rownames(Transfer_labels)] <- Transfer_labels[,1]
Ariss_labels[grep('Brain', Ariss_labels)] <- 'Brain (Not transferred)'
EAdisc <- AddMetaData(object = EAdisc, metadata = as.data.frame(Ariss_labels))
#Set colors
colors <- c("#59F889", "#FD79A3", "#BB90FE", "#85BC05", "#01D3F8", "#01FAF7", "#67EAFF", "#02B4E2", "#9EF75F", "#6CA7FF", "#000000", "#996515", "#FFC29C", "#A020F0", "#01E6B3", "#E5D006")
names(colors) <- levels(x = EAdisc@meta.data$Ariss_labels)
saveRDS(colors, file='Figure_1/Processed_data/Seurat/ArissTo10X_ColVars.Rds')
# Save object
saveRDS(EAdisc, file='Figure_1/Processed_data/Seurat/10X_SeuratObject.Rds')

#Remove other cells for visualization
cells <- rownames(EAdisc@meta.data[which(EAdisc@meta.data$Ariss_labels != 'other'),])
#Load
colors <- readRDS('Figure_1/Processed_data/Seurat/ArissTo10X_ColVars.Rds')
EAdisc <- readRDS('Figure_1/Processed_data/Seurat/10X_SeuratObject.Rds')
#Plot
DimPlot(object = EAdisc, reduction = 'tsne', group.by='Ariss_labels', cells=cells, cols=colors, label=TRUE, label.size = 3.5) + NoLegend() + NoAxes()