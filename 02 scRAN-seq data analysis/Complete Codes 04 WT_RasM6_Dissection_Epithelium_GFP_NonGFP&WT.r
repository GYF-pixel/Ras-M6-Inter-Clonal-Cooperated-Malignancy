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

getwd()
setwd("F:\\05Cooperated_Project\\02M6\\scRNAseq")
setwd("H:/M6")

pbmc_GFP <- readRDS("CCA_WTRM6_UniqueClusters_reclustered_Epithelium_GFP.Rds")

pbmc_NonGFP <- readRDS("CCA_WTRM6_UniqueClusters_reclustered_Epithelium_NonGFP.Rds")

pbmc_WT <- readRDS("CCA_WT_Epithelium.Rds")

###WT&GFP
DefaultAssay(pbmc_WT) <- "RNA"
DefaultAssay(pbmc_GFP) <- "RNA"

WTGFP.anchors <- FindIntegrationAnchors(object.list = list(pbmc_WT, pbmc_GFP), anchor.features = 2000, dims = 1:50)
WTGFP.combined <- IntegrateData(anchorset = WTGFP.anchors, dims = 1:50)
pbmc <-  WTGFP.combined

#pbmc <- NormalizeData(pbmc, verbose = FALSE, normalization.method = "LogNormalize", scale.factor = 1e4)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2500)
pbmc <- ScaleData(pbmc, vars.to.regress = c("nCount_RNA"), verbose = TRUE)
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc), npcs = 40, nfeature.print = 10, ndims.print = 1:5, verbose = T)
pc.num=1:40
pbmc <- RunUMAP(pbmc, dims=pc.num)
pbmc <- FindNeighbors(pbmc, dims = pc.num)
pbmc = FindClusters(pbmc,resolution = 1.2)

saveRDS(pbmc, file="CCA_WTRM6_UniqueClusters_reclustered_Epithelium_GFP&WT_reclustered.Rds") 


##检查default assay
DefaultAssay(pbmc)

##查看有多少个cluster
levels(pbmc)

##设置对应cluster的配色数目
col5 <- colorRampPalette((pal_npg(palette = c("nrc"))(7)))(19)
show_col(col5)

##正常的UMAP降维数据结果
UMPplot_label <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "seurat_clusters", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("")
ggsave("P1_UMPplot_label.pdf", plot = plot_grid(UMPplot_label), width = 10, height = 8)
UMPplot_unlabel <- DimPlot(pbmc, reduction = "umap", label = FALSE, group.by = "seurat_clusters", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("")
ggsave("P1_UMPplot_unlabel.pdf", plot = plot_grid(UMPplot_unlabel), width = 10, height = 8)


##相同Genotype的分组结果
table(pbmc@meta.data[["group"]])
##按group区别
UMPplot_label_group <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "group", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("")
ggsave("P2_UMPplot_label_group.pdf", plot = plot_grid(UMPplot_label_group), width = 10, height = 8)
UMPplot_unlabel_group <- DimPlot(pbmc, reduction = "umap", label = FALSE, group.by = "group", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("")
ggsave("P2_UMPplot_unlabel_group.pdf", plot = plot_grid(UMPplot_unlabel_group), width = 10, height = 8)

##相同Genotype的分组结果
table(pbmc@meta.data[["sec_group"]])
##按sec_group区别
UMPplot_label_sec_group <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "sec_group", pt.size = 1.0, label.size = 7, cols = c("#4DBBD5","#E64B35","#00A087","#3C5488"), raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("")
ggsave("P2_UMPplot_label_sec_group.pdf", plot = plot_grid(UMPplot_label_sec_group), width = 10, height = 8)
UMPplot_unlabel_sec_group <- DimPlot(pbmc, reduction = "umap", label = FALSE, group.by = "sec_group", pt.size = 1.0, label.size = 7, cols = c("#4DBBD5","#E64B35","#00A087","#3C5488"), raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("")
ggsave("P2_UMPplot_unlabel_sec_group.pdf", plot = plot_grid(UMPplot_unlabel_sec_group), width = 10, height = 8)

##各cell cluster的细胞占比
table(pbmc@meta.data$group, pbmc@meta.data$integrated_snn_res.1.2)
table(pbmc@meta.data$sec_group, pbmc@meta.data$integrated_snn_res.1.2)

##生成文件画出堆积柱状图
Cells_in_different_Clusters_at_different_group <- table(pbmc@meta.data$sec_group, pbmc@meta.data$integrated_snn_res.1.2)
write.csv(Cells_in_different_Clusters_at_different_group,file = "P3_Cells_in_different_Clusters_at_different_group.csv",row.names = T)
##计算单个cluster中两种cell type的占比，设置Unique Cluster筛选条件

#改为RNA
DefaultAssay(pbmc) <- "RNA"

#Old Codes
a <- FeaturePlot(pbmc,features = c("GFP"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("GFP")
ggsave("P4_pbmcUMAP1_GFPplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(pbmc,features = c("Ras85D"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Ras85D")
ggsave("P4_pbmcUMAP1_Ras85Dplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(pbmc,features = c("M6"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("M6")
ggsave("P4_pbmcUMAP1_M6plot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(pbmc,features = c("GAL80"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("GAL80")
ggsave("P4_pbmcUMAP1_GAL80plot2.pdf", plot = plot_grid(a), width = 6, height = 5)

#Select Unique Cluster >120 cells & >80%
pbmc@meta.data$UniqueCluster = NA
Cells.sub <- subset(pbmc@meta.data, integrated_snn_res.1.2 %in% c("5","11","12","13","14"))
table(pbmc@meta.data$integrated_snn_res.1.2)
table(Cells.sub$integrated_snn_res.1.2)

for (i in row.names(Cells.sub)) {
    pbmc@meta.data[[i, "UniqueCluster"]] <- as.character(pbmc@meta.data[[i,"integrated_snn_res.1.2"]])
}
table(pbmc@meta.data$UniqueCluster)
length(unique(pbmc@meta.data$UniqueCluster))
##as.numeric来实现cluster排列顺序的调整
pbmc@meta.data$UniqueCluster <- as.numeric(pbmc@meta.data$UniqueCluster)
table(pbmc@meta.data$UniqueCluster)
length(unique(pbmc@meta.data$UniqueCluster))

##设置对应cluster的配色数目
col5 <- colorRampPalette((pal_npg(palette = c("nrc"))(7)))(6)
show_col(col5)


##单独显示Unique Cluster，其他是NA
UMP1plot_UniqeClusters <- DimPlot(pbmc, reduction = "umap", pt.size = 1.2, cols = col5, group.by ="UniqueCluster", raster=FALSE, label = TRUE)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Uniqe Clusters")
ggsave("P5_UMP1plot_UniqeClusters.pdf", plot = plot_grid(UMP1plot_UniqeClusters), width = 7, height = 6)
UMP1plot_UniqeClusters_unlabel <- DimPlot(pbmc, reduction = "umap", pt.size = 1.2, cols = col5, group.by ="UniqueCluster", raster=FALSE, label = FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Uniqe Clusters")
ggsave("P5_UMP1plot_UniqeClusters_unlabel.pdf", plot = plot_grid(UMP1plot_UniqeClusters_unlabel), width = 7, height = 6)


# GFP Cluster 7 16 vs WT_Epithelial
clusterGFP.markers <- FindMarkers(pbmc, ident.1 = c("5","11","12","13","14"), ident.2 = c("0","1","2","4","6","7","8","9","10","15","16","17","18"), min.pct = 0.25)
head(clusterGFP.markers, n = 5)
write.csv(clusterGFP.markers, "P5_Unique_GFP.markers_vs_WT.csv", row.names = T)

clusterGFP5.markers <- FindMarkers(pbmc, ident.1 = c("5"), ident.2 = c("0","1","2","4","6","7","8","9","10","15","16","17","18"), min.pct = 0.25)
head(clusterGFP5.markers, n = 5)
write.csv(clusterGFP5.markers, "P5_Unique_GFP5.markers_vs_WT.csv", row.names = T)

clusterGFP11.markers <- FindMarkers(pbmc, ident.1 = c("11"), ident.2 = c("0","1","2","4","6","7","8","9","10","15","16","17","18"), min.pct = 0.25)
head(clusterGFP11.markers, n = 5)
write.csv(clusterGFP11.markers, "P5_Unique_GFP11.markers_vs_WT.csv", row.names = T)

clusterGFP12.markers <- FindMarkers(pbmc, ident.1 = c("12"), ident.2 = c("0","1","2","4","6","7","8","9","10","15","16","17","18"), min.pct = 0.25)
head(clusterGFP12.markers, n = 5)
write.csv(clusterGFP12.markers, "P5_Unique_GFP12.markers_vs_WT.csv", row.names = T)

clusterGFP13.markers <- FindMarkers(pbmc, ident.1 = c("13"), ident.2 = c("0","1","2","4","6","7","8","9","10","15","16","17","18"), min.pct = 0.25)
head(clusterGFP13.markers, n = 5)
write.csv(clusterGFP13.markers, "P5_Unique_GFP13.markers_vs_WT.csv", row.names = T)

clusterGFP14.markers <- FindMarkers(pbmc, ident.1 = c("14"), ident.2 = c("0","1","2","4","6","7","8","9","10","15","16","17","18"), min.pct = 0.25)
head(clusterGFP14.markers, n = 5)
write.csv(clusterGFP14.markers, "P5_Unique_GFP14.markers_vs_WT.csv", row.names = T)

#Genes of interest
a <- FeaturePlot(pbmc,features = c("Dif", "dl","Myd88"),cols = c('lightgrey','red'),ncol = 1)
ggsave("P7_pbmcUMAP1_Toll_Feature_plot.pdf", plot = plot_grid(a), width = 3, height = 7)
a <- FeaturePlot(pbmc,features = c("kay","Jra","puc","Mmp1","bsk","msn"),cols = c('lightgrey','red'),ncol=2)
ggsave("P7_pbmcUMAP1_MAPK_Feature_plot.pdf", plot = plot_grid(a), width = 6, height = 7)
a <- FeaturePlot(pbmc,features = c("fj", "Diap1","sd"),cols = c('lightgrey','red'),ncol = 1)
ggsave("P7_pbmcUMAP1_Hippo_Feature_plot.pdf", plot = plot_grid(a), width = 3, height = 7)

a <- FeaturePlot(scRNAsub,features = c("kay"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("kay")
ggsave("P6_scRNAsubUMAP1_kayplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(scRNAsub,features = c("Jra"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Jra")
ggsave("P6_scRNAsubUMAP1_Jraplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(scRNAsub,features = c("Mmp1"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Mmp1")
ggsave("P6_scRNAsubUMAP1_Mmp1plot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(scRNAsub,features = c("nmo"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("nmo")
ggsave("P6_scRNAsubUMAP1_nmoplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(scRNAsub,features = c("puc"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("puc")
ggsave("P6_scRNAsubUMAP1_pucplot2.pdf", plot = plot_grid(a), width = 6, height = 5)




###WT&NonGFP
DefaultAssay(pbmc_WT) <- "RNA"
DefaultAssay(pbmc_NonGFP) <- "RNA"

WTNonGFP.anchors <- FindIntegrationAnchors(object.list = list(pbmc_WT, pbmc_NonGFP), anchor.features = 2000, dims = 1:50)
WTNonGFP.combined <- IntegrateData(anchorset = WTNonGFP.anchors, dims = 1:50)
pbmc <-  WTNonGFP.combined

#pbmc <- NormalizeData(pbmc, verbose = FALSE, normalization.method = "LogNormalize", scale.factor = 1e4)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2500)
pbmc <- ScaleData(pbmc, vars.to.regress = c("nCount_RNA"), verbose = TRUE)
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc), npcs = 40, nfeature.print = 10, ndims.print = 1:5, verbose = T)
pc.num=1:40
pbmc <- RunUMAP(pbmc, dims=pc.num)
pbmc <- FindNeighbors(pbmc, dims = pc.num)
pbmc = FindClusters(pbmc,resolution = 1.0)

saveRDS(pbmc, file="CCA_WTRM6_UniqueClusters_reclustered_Epithelium_NonGFP&WT_reclustered.Rds") 

#pbmc <- readRDS("CCA_WTRM6_UniqueClusters_reclustered_Epithelium_NonGFP&WT_reclustered.Rds")
##检查default assay
DefaultAssay(pbmc)

##查看有多少个cluster
levels(pbmc)

##设置对应cluster的配色数目
col5 <- colorRampPalette((pal_npg(palette = c("nrc"))(7)))(21)
show_col(col5)

##正常的UMAP降维数据结果
UMPplot_label <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "seurat_clusters", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("")
ggsave("P1_UMPplot_label.pdf", plot = plot_grid(UMPplot_label), width = 10, height = 8)
UMPplot_unlabel <- DimPlot(pbmc, reduction = "umap", label = FALSE, group.by = "seurat_clusters", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("")
ggsave("P1_UMPplot_unlabel.pdf", plot = plot_grid(UMPplot_unlabel), width = 10, height = 8)


##相同Genotype的分组结果
table(pbmc@meta.data[["group"]])
##按group区别
UMPplot_label_group <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "group", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("")
ggsave("P2_UMPplot_label_group.pdf", plot = plot_grid(UMPplot_label_group), width = 10, height = 8)
UMPplot_unlabel_group <- DimPlot(pbmc, reduction = "umap", label = FALSE, group.by = "group", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("")
ggsave("P2_UMPplot_unlabel_group.pdf", plot = plot_grid(UMPplot_unlabel_group), width = 10, height = 8)

##相同Genotype的分组结果
table(pbmc@meta.data[["sec_group"]])
##按sec_group区别
UMPplot_label_sec_group <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "sec_group", pt.size = 1.0, label.size = 7, cols = c("#4DBBD5","#E64B35","#00A087","#3C5488"), raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("")
ggsave("P2_UMPplot_label_sec_group.pdf", plot = plot_grid(UMPplot_label_sec_group), width = 10, height = 8)
UMPplot_unlabel_sec_group <- DimPlot(pbmc, reduction = "umap", label = FALSE, group.by = "sec_group", pt.size = 1.0, label.size = 7, cols = c("#4DBBD5","#E64B35","#00A087","#3C5488"), raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("")
ggsave("P2_UMPplot_unlabel_sec_group.pdf", plot = plot_grid(UMPplot_unlabel_sec_group), width = 10, height = 8)

##各cell cluster的细胞占比
table(pbmc@meta.data$group, pbmc@meta.data$integrated_snn_res.1)
table(pbmc@meta.data$sec_group, pbmc@meta.data$integrated_snn_res.1)

##生成文件画出堆积柱状图
Cells_in_different_Clusters_at_different_group <- table(pbmc@meta.data$sec_group, pbmc@meta.data$integrated_snn_res.1)
write.csv(Cells_in_different_Clusters_at_different_group,file = "P3_Cells_in_different_Clusters_at_different_group.csv",row.names = T)
##计算单个cluster中两种cell type的占比，设置Unique Cluster筛选条件


#改为RNA
DefaultAssay(pbmc) <- "RNA"

#Old Codes
a <- FeaturePlot(pbmc,features = c("GFP"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("GFP")
ggsave("P4_pbmcUMAP1_GFPplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(pbmc,features = c("Ras85D"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Ras85D")
ggsave("P4_pbmcUMAP1_Ras85Dplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(pbmc,features = c("M6"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("M6")
ggsave("P4_pbmcUMAP1_M6plot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(pbmc,features = c("GAL80"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("GAL80")
ggsave("P4_pbmcUMAP1_GAL80plot2.pdf", plot = plot_grid(a), width = 6, height = 5)


#Select Unique Cluster >120 cells & >80%
pbmc@meta.data$UniqueCluster = NA
Cells.sub <- subset(pbmc@meta.data, integrated_snn_res.1 %in% c("2","7","14","18"))
table(pbmc@meta.data$integrated_snn_res.1)
table(Cells.sub$integrated_snn_res.1)

for (i in row.names(Cells.sub)) {
    pbmc@meta.data[[i, "UniqueCluster"]] <- as.character(pbmc@meta.data[[i,"integrated_snn_res.1"]])
}
table(pbmc@meta.data$UniqueCluster)
length(unique(pbmc@meta.data$UniqueCluster))
##as.numeric来实现cluster排列顺序的调整
pbmc@meta.data$UniqueCluster <- as.numeric(pbmc@meta.data$UniqueCluster)
table(pbmc@meta.data$UniqueCluster)
length(unique(pbmc@meta.data$UniqueCluster))

##设置对应cluster的配色数目
col5 <- colorRampPalette((pal_npg(palette = c("nrc"))(7)))(5)
show_col(col5)

##单独显示Unique Cluster，其他是NA
UMP1plot_UniqeClusters <- DimPlot(pbmc, reduction = "umap", pt.size = 1.2, cols = col5, group.by ="UniqueCluster", raster=FALSE, label = TRUE)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Uniqe Clusters")
ggsave("P5_UMP1plot_UniqeClusters.pdf", plot = plot_grid(UMP1plot_UniqeClusters), width = 7, height = 6)
UMP1plot_UniqeClusters_unlabel <- DimPlot(pbmc, reduction = "umap", pt.size = 1.2, cols = col5, group.by ="UniqueCluster", raster=FALSE, label = FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Uniqe Clusters")
ggsave("P5_UMP1plot_UniqeClusters_unlabel.pdf", plot = plot_grid(UMP1plot_UniqeClusters_unlabel), width = 7, height = 6)

# NonGFP vs WT_Epithelial
clusterNonGFP.markers <- FindMarkers(pbmc, ident.1 = c("2","7","14","18"), ident.2 = c("0","1","3","4","5","6","8","9","10","11","12","13","15","16","17","19","20"), min.pct = 0.25)
head(clusterNonGFP.markers, n = 5)
write.csv(clusterNonGFP.markers, "P5_Unique_NonGFP.markers_vs_WT.csv", row.names = T)

clusterGFP2.markers <- FindMarkers(pbmc, ident.1 = c("2"), ident.2 = c("0","1","3","4","5","6","8","9","10","11","12","13","15","16","17","19","20"), min.pct = 0.25)
head(clusterGFP2.markers, n = 5)
write.csv(clusterGFP2.markers, "P5_Unique_NonGFP2.markers_vs_WT.csv", row.names = T)

clusterGFP7.markers <- FindMarkers(pbmc, ident.1 = c("7"), ident.2 = c("0","1","3","4","5","6","8","9","10","11","12","13","15","16","17","19","20"), min.pct = 0.25)
head(clusterGFP7.markers, n = 5)
write.csv(clusterGFP7.markers, "P5_Unique_NonGFP7.markers_vs_WT.csv", row.names = T)

clusterGFP14.markers <- FindMarkers(pbmc, ident.1 = c("14"), ident.2 = c("0","1","3","4","5","6","8","9","10","11","12","13","15","16","17","19","20"), min.pct = 0.25)
head(clusterGFP14.markers, n = 5)
write.csv(clusterGFP14.markers, "P5_Unique_NonGFP14.markers_vs_WT.csv", row.names = T)

clusterGFP18.markers <- FindMarkers(pbmc, ident.1 = c("18"), ident.2 = c("0","1","3","4","5","6","8","9","10","11","12","13","15","16","17","19","20"), min.pct = 0.25)
head(clusterGFP18.markers, n = 5)
write.csv(clusterGFP18.markers, "P5_Unique_NonGFP18.markers_vs_WT.csv", row.names = T)


#Genes of interest
a <- FeaturePlot(pbmc,features = c("kay","Jra","puc","Egfr","sty","kek1","nmo","pnt","aop"),cols = c('lightgray','red'),ncol=3)
ggsave("P7_pbmcUMAP1_Apoptosis&MAPK_Feature_plot.pdf", plot = plot_grid(a), width = 9, height = 7)




gene_sets <- read.table("P7_pbmcUMAP1_MAPK_signaling.txt")
gene_sets2 <- list(gene_sets[,1])
gene_sets2

DefaultAssay(pbmc) <- "RNA"

sce_T <- AddModuleScore(pbmc, features = gene_sets, nbin = 24, ctrl = 100,name = "MAPK_signaling_score")
#Toll_Imd_signaling_score AddModuleScore之后就是Toll_Imd_signaling_score1
FeaturePlot(sce_T,'MAPK_signaling_score1',cols=rev(brewer.pal(10, name = "RdBu")))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("MAPK signaling score")
ggsave("P7_pbmcUMAP1_MAPK_signaling_score.pdf")

table(sce_T@meta.data[["sec_group"]])


Cells.sub <- subset(sce_T@meta.data, sec_group %in% c("RM6"))
scRNAsub <- subset(sce_T, cells=row.names(Cells.sub))

FeaturePlot(scRNAsub,'MAPK_signaling_score1',order=TRUE,cols=rev(brewer.pal(10, name = "RdBu")))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("MAPK signaling score")
ggsave("P7_pbmcUMAP1_MAPK_signaling_score_NonGFP.pdf")


Cells.sub <- subset(sce_T@meta.data, integrated_snn_res.1 %in% c("2","7","14","18"))
scRNAsub <- subset(sce_T, cells=row.names(Cells.sub))

FeaturePlot(scRNAsub,'MAPK_signaling_score1',order=TRUE,cols=rev(brewer.pal(10, name = "RdBu")))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("MAPK signaling score")
ggsave("P7_pbmcUMAP1_MAPK_signaling_score_NonGFP_2_7_14_18.pdf")



