
library(optparse)
library(tidyverse)
library(future.apply)
library(Seurat)
library(RColorBrewer)
library(reshape2)
library(network)
library(igraph)
setwd("I:\\FlyPhoneDB-master\\FlyPhoneDB-master")
setwd("I:/M6")
setwd("I:\\M6\\00 Step6 Ligand&Receptor Analysis\\GFP_Hemocytes")

# heatmap
exprMat <- read.csv("expMatrix.csv", row.names = 1, check.names = FALSE)
# exprMat <- read.csv(opt$matrix, row.names = 1, check.names = FALSE)
cellInfo <- read.csv("meta.csv", row.names = 1)
# cellInfo <- read.csv(opt$metadata, row.names = 1)
cellInfo$celltype <- as.character(cellInfo$celltype)
# str(cellInfo)

exprMat <- exprMat[ , row.names(cellInfo)]

# create seuratObj
seuratObj <- CreateSeuratObject(counts = exprMat)
seuratObj <- NormalizeData(seuratObj)
seuratObj <- FindVariableFeatures(seuratObj, selection.method = "vst", nfeatures = 2000)
all_genes <- rownames(seuratObj)
seuratObj <- ScaleData(seuratObj, features = all_genes)
seuratObj$celltype <- as.factor(cellInfo$celltype)
Idents(seuratObj) <- "celltype"
# clusterMetadataTable <- table(seuratObj@meta.data[ , "celltype"]) %>% as.data.frame()
# colnames(clusterMetadataTable) <- c("celltype", "count")
# print(clusterMetadataTable$celltype)

exprMat <- sweep(exprMat, 2, Matrix::colSums(exprMat), FUN = "/") * 10000



avgexp <- AverageExpression(seuratObj, assay = "RNA", return.seurat = TRUE)

Pathway_core_components <- read.table("./annotation/Pathway_core_components_2021vs1_clean.txt", sep = "\t", header = TRUE)

df <- Pathway_core_components[Pathway_core_components$pathway == "PVR RTK signaling pathway",]
genes <- df$gene

p <- DoHeatmap(avgexp, features = genes, label = TRUE ,draw.lines = FALSE, raster = FALSE, angle = 90) +
    scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu"))) # & NoLegend()
print(p)
ggsave(p, file = "Heatmap.pdf",   # The directory you want to save the file in
         width = 8, # The width of the plot in inches
         height = 12)


mat <- GetAssayData(avgexp,slot = 'scale.data')
mat <- as.data.frame(mat)
input <- mat[genes,]
write.csv(input, file="input.csv")

#Heatmap 画热图
library(pheatmap)
library(RColorBrewer)

colors <- colorRampPalette(c("navy", "white", "firebrick3"))(100)
heatmap=pheatmap(input,color = colors,
                 main="",
                 fontsize = 15,
                 scale="row",
                 border_color = "black",
                 na_col = "grey",
                 cluster_rows = T,cluster_cols = T,
                 show_rownames = T,show_colnames = T,
                 treeheight_row = 20,treeheight_col = 10,
                 cellheight = 12,cellwidth = 20,
                 cutree_row=2,cutree_col=2,
                 display_numbers = F,legend = T,
)
heatmap

ggsave("Part3_Heatmap_genes.pdf", plot = heatmap, width = 8, height = 8) 

#pbmc <- pbmc_Epithelium
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
#library(scCustomize) # 需要Seurat版本4.3.0
library(viridis)
library(RColorBrewer)
library(gridExtra)
set.seed(1234)
setwd("I:/M6")
pbmc <- readRDS("CCA_WTRM6_UniqueClusters_reclustered_Epithelium_GFP_and_Hemocytes.Rds")
DefaultAssay(pbmc) <- "RNA"

a <- FeaturePlot(pbmc,features = c("Pvf1", "Pvf2","Pvf3","Pvr"),cols = c('lightgrey','red'),ncol = 4)
ggsave("Part4_Pvf_genes_pbmcUMAP1_Hemocytes_plot.pdf", plot = plot_grid(a), width = 16, height = 3)


pbmc <- readRDS("CCA_WTRasM6_UniqueClusters_reclustered_Hemocytes_reclustered.Rds")
a <- FeaturePlot(pbmc,features = c("Pvr"),cols = c('lightgray','red'),pt.size = 1.6)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Pvr")
ggsave("Part4_Pvr_genes_pbmcUMAP1_Hemocytes_plot.pdf", plot = plot_grid(a), width = 5, height = 4)
pbmc <- readRDS("CCA_WT_RasM6_Bravo_XGY_Hemocytes_Reculstered.Rds")
a <- FeaturePlot(pbmc,features = c("Pvr"),cols = c('lightgray','red'),pt.size = 1.6)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Pvr")
ggsave("Part4_Pvr_genes_pbmcUMAP1_All_Hemocytes_plot.pdf", plot = plot_grid(a), width = 5, height = 4)
