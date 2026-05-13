
#scRNAsub_Epithelium重新cluster
#scRNAsub <- scRNAsub_Epithelium
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
scRNAsub <- readRDS("CCA_WTRM6_UniqueClusters_reclustered_Epithelium.Rds")

DefaultAssay(scRNAsub) <- "integrated"

scRNAsub <- FindVariableFeatures(scRNAsub, selection.method = "vst", nfeatures = 2500)
scRNAsub <- ScaleData(scRNAsub, vars.to.regress = c("nCount_RNA"), verbose = TRUE)
scRNAsub <- RunPCA(scRNAsub, features = VariableFeatures(scRNAsub), npcs = 60, nfeature.print = 10, ndims.print = 1:5, verbose = T)
pc.num=1:60
scRNAsub <- RunUMAP(scRNAsub, dims=pc.num)
scRNAsub <- FindNeighbors(scRNAsub, dims = pc.num)
scRNAsub = FindClusters(scRNAsub,resolution = 1.3)

##查看有多少个cluster
levels(scRNAsub)

##设置对应cluster的配色数目
col5 <- colorRampPalette((pal_npg(palette = c("nrc"))(7)))(19)
show_col(col5)

##正常的UMAP降维数据结果
UMPplot_label <- DimPlot(scRNAsub, reduction = "umap", label = TRUE, group.by = "seurat_clusters", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("")
ggsave("P1_UMPplot_label.pdf", plot = plot_grid(UMPplot_label), width = 10, height = 8)
UMPplot_unlabel <- DimPlot(scRNAsub, reduction = "umap", label = FALSE, group.by = "seurat_clusters", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("")
ggsave("P1_UMPplot_unlabel.pdf", plot = plot_grid(UMPplot_unlabel), width = 10, height = 8)

##相同Genotype的分组结果
table(scRNAsub@meta.data[["group"]])
##按group区别
UMPplot_label_group <- DimPlot(scRNAsub, reduction = "umap", label = TRUE, group.by = "group", pt.size = 1.0, label.size = 7, cols = c("#4DBBD5","#E64B35","#00A087","#3C5488"), raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("")
ggsave("P2_UMPplot_label_group.pdf", plot = plot_grid(UMPplot_label_group), width = 10, height = 8)
UMPplot_unlabel_group <- DimPlot(scRNAsub, reduction = "umap", label = FALSE, group.by = "group", pt.size = 1.0, label.size = 7, cols = c("#4DBBD5","#E64B35","#00A087","#3C5488"), raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("")
ggsave("P2_UMPplot_unlabel_group.pdf", plot = plot_grid(UMPplot_unlabel_group), width = 10, height = 8)

##这里需要将默认格式输出调整为“RNA”
DefaultAssay(scRNAsub) <- "RNA"
scRNAsub <- ScaleData(scRNAsub, vars.to.regress = c("nCount_RNA"), verbose = TRUE)

## report only the positive ones
scRNAsub.markers <- FindAllMarkers(scRNAsub,
                               only.pos = FALSE,
                               min.pct = 0.4,
                               logfc.threshold = 0.25)

Allmarkers_scRNAsub = scRNAsub.markers %>% select(gene, everything()) %>% subset(p_val<0.05)
write.csv(Allmarkers_scRNAsub, "P3_Allmarkers__scRNAsub_wilcox.csv", row.names = T)

## 使用Top 5基因画热图
top5scRNAsub.markers <- scRNAsub.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)
##做Heatmap的热图
Heatmap_markers_scRNAsub <- DoHeatmap(scRNAsub,features = top5scRNAsub.markers$gene,
          group.colors = col5) + theme(text = element_text (size = 18)) + 
		  scale_fill_gradient2(low = '#0099CC',mid = 'white',high = '#CC0033',name = 'Z-score')
ggplot2::ggsave(filename = 'P3_Heatmap_markers_scRNAsub.pdf',width =16,height = 19)

##做jjVolcano的火山图
library(scRNAtoolVis)
jjVolcano_scRNAsub.markers <- jjVolcano(diffData = scRNAsub.markers,
          log2FC.cutoff = 0.25, 
          size  = 3.5, #设置点的大小
          fontface = 'italic', #设置字体形式
          aesCol = c('#00468B','#ED0000'), #设置点的颜色
          tile.col = col5, #设置cluster的颜色
          #col.type = "adjustP", #设置矫正方式
          topGeneN = 5 #设置展示topN的基因
         )
jjVolcano_scRNAsub.markers
ggsave(filename = 'P3_jjVolcano_scRNAsub.markers.pdf', plot = plot_grid(jjVolcano_scRNAsub.markers), width = 15, height = 6)


##Old codes
a <- FeaturePlot(scRNAsub,features = c("GFP"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("GFP")
ggsave("P3_scRNAsubUMAP1_GFPplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(scRNAsub,features = c("Ras85D"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Ras85D")
ggsave("P3_scRNAsubUMAP1_Ras85Dplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(scRNAsub,features = c("M6"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("M6")
ggsave("P3_scRNAsubUMAP1_M6plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(scRNAsub,features = c("GAL80"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("GAL80")
ggsave("P3_scRNAsubUMAP1_GAL80plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(scRNAsub,features = c("puc"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("puc")
ggsave("P3_scRNAsubUMAP1_pucplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(scRNAsub,features = c("Jra"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Jra")
ggsave("P3_scRNAsubUMAP1_Jraplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(scRNAsub,features = c("kay"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("kay")
ggsave("P3_scRNAsubUMAP1_kayplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)

a <- FeaturePlot(scRNAsub,features = c("CycE"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("CycE")
ggsave("P3_scRNAsubUMAP1_CycEplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(scRNAsub,features = c("Diap1"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Diap1")
ggsave("P3_scRNAsubUMAP1_Diap1plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(scRNAsub,features = c("Myc"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Myc")
ggsave("P3_scRNAsubUMAP1_Mycplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(scRNAsub,features = c("ex"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("ex")
ggsave("P3_scRNAsubUMAP1_explot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(scRNAsub,features = c("dl"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("dl")
ggsave("P3_scRNAsubUMAP1_dlplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(scRNAsub,features = c("Rel"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Rel")
ggsave("P3_scRNAsubUMAP1_Relplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(scRNAsub,features = c("Dif"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Dif")
ggsave("P3_scRNAsubUMAP1_Difplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(scRNAsub,features = c("wg"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("wg")
ggsave("P3_scRNAsubUMAP1_wgplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
##GFP细胞在各个Cluster中的占比

#计算GFP 细胞占比
table_GFP <- as.data.frame(scRNAsub@assays$RNA@counts["GFP",])
table(table_GFP[,1] != 0)
a_GFP <- subset(table_GFP, table_GFP[,1]!=0)

scRNAsub_GFP <- subset(scRNAsub, cells=row.names(a_GFP))
table(scRNAsub_GFP@meta.data$group)

GFP_seurat_clusters <- table(scRNAsub_GFP@meta.data$sec_group, scRNAsub_GFP@meta.data$seurat_clusters)
GFP_seurat_clusters 
write.csv(GFP_seurat_clusters, file = "P4_GFP_seurat_clusters_WTRasM6.csv",row.names = T)

#各个Cluster中细胞的占比
Cells_in_different_Clusters_at_different_group <- table(scRNAsub@meta.data$sec_group, scRNAsub@meta.data$seurat_clusters)
write.csv(Cells_in_different_Clusters_at_different_group,file = "P4_Cells_in_different_Clusters_at_different_group.csv",row.names = T)

saveRDS(scRNAsub, file="CCA_WTRM6_UniqueClusters_reclustered_Epithelium_Reclustered.Rds") 






#画出占比图

#提取GFP细胞和NonGFP细胞

levels(scRNAsub)
Cells.sub <- subset(scRNAsub@meta.data, integrated_snn_res.1.2 %in% c("2","4","7","8","11","12","16","17"))

scRNAsub_Epithelium_GFP <- subset(scRNAsub, cells=rownames(Cells.sub))

saveRDS(scRNAsub_Epithelium_GFP, file="CCA_WTRM6_UniqueClusters_reclustered_Epithelium_GFP.Rds") 

genes = c("Tl",				#Toll-1
				"18w",			#Toll-2
				#"MstProx", 	#Toll-3		
				"Toll-4",
				#"Tehao",		#Toll-5
				#"Toll-6",
				"Toll-7",
				"Tollo",			#Toll-8	
				#"Toll-9",
				"spz",	
				"NT1",			#spz2
				"spz3",	
				"spz4",	
				"spz5",	
				"spz6"
				)

a <- FeaturePlot(pbmc,features = genes,cols = c('lightgrey','red'),order = T, ncol = 3)
ggsave("P3_pbmcUMAP1_Epithelium_GFP_Toll_Feature_plot.pdf", plot = plot_grid(a), width = 11, height = 12)








levels(scRNAsub)
Cells.sub <- subset(scRNAsub@meta.data, integrated_snn_res.1.2 %in% c("0","1","3","5","6","9","10","13","15","18"))

scRNAsub_Epithelium_NonGFP <- subset(scRNAsub, cells=rownames(Cells.sub))

saveRDS(scRNAsub_Epithelium_NonGFP, file="CCA_WTRM6_UniqueClusters_reclustered_Epithelium_NonGFP.Rds") 
