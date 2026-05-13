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

################
pbmc <- readRDS("CCA_WTR79.Rds")

##检查default assay
DefaultAssay(pbmc)

##查看有多少个cluster
levels(pbmc)

##设置对应cluster的配色数目
col5 <- colorRampPalette((pal_npg(palette = c("nrc"))(7)))(26)
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

#Select Unique Cluster >120 cells & >80%
#Cluster Unique in R79 c("6","10","13","18","20","21","24","25")
pbmc@meta.data$UniqueCluster = NA

Cells.sub <- subset(pbmc@meta.data, integrated_snn_res.1 %in% c("6","10","13","18","20","21","24","25"))
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
col5 <- colorRampPalette((pal_npg(palette = c("nrc"))(7)))(8)
show_col(col5)

##单独显示Unique Cluster，其他是NA
UMP1plot_UniqeClusters <- DimPlot(pbmc, reduction = "umap", pt.size = 1.2, cols = col5, group.by ="UniqueCluster", raster=FALSE, label = TRUE)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Uniqe Clusters")
ggsave("P6_UMP1plot_UniqeClusters.pdf", plot = plot_grid(UMP1plot_UniqeClusters), width = 7, height = 6)
UMP1plot_UniqeClusters_unlabel <- DimPlot(pbmc, reduction = "umap", pt.size = 1.2, cols = col5, group.by ="UniqueCluster", raster=FALSE, label = FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Uniqe Clusters")
ggsave("P6_UMP1plot_UniqeClusters_unlabel.pdf", plot = plot_grid(UMP1plot_UniqeClusters_unlabel), width = 7, height = 6)

##保存含有Unique Cluster信息的pbmc
saveRDS(pbmc, file = "CCA_WTR79.Rds")

##################
#提取Unique Cluster的全部细胞
scRNAsub <- subset(pbmc, cells=row.names(Cells.sub))

##检查default assay
DefaultAssay(scRNAsub)

##剔除WT细胞
table(scRNAsub@meta.data$sec_group)

Cells.sub2 <- subset(scRNAsub@meta.data, sec_group %in% c("R79"))
scRNAsub <- subset(scRNAsub, cells=row.names(Cells.sub2))

table(scRNAsub@meta.data$sec_group)

##重新聚类分析
DefaultAssay(scRNAsub) <- "integrated"

scRNAsub <- FindVariableFeatures(scRNAsub, selection.method = "vst", nfeatures = 2500)
scRNAsub <- ScaleData(scRNAsub, vars.to.regress = c("nCount_RNA"), verbose = TRUE)
scRNAsub <- RunPCA(scRNAsub, features = VariableFeatures(scRNAsub), npcs = 60, nfeature.print = 10, ndims.print = 1:5, verbose = T)
pc.num=1:60
scRNAsub <- RunUMAP(scRNAsub, dims=pc.num)
scRNAsub <- FindNeighbors(scRNAsub, dims = pc.num)
scRNAsub = FindClusters(scRNAsub,resolution = 1.2)

##查看有多少个cluster
levels(scRNAsub)

##设置对应cluster的配色数目
col5 <- colorRampPalette((pal_npg(palette = c("nrc"))(7)))(19)
show_col(col5)

##正常的UMAP降维数据结果
UMPplot_label <- DimPlot(scRNAsub, reduction = "umap", label = TRUE, group.by = "seurat_clusters", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("")
ggsave("P7_UMPplot_label.pdf", plot = plot_grid(UMPplot_label), width = 10, height = 8)
UMPplot_unlabel <- DimPlot(scRNAsub, reduction = "umap", label = FALSE, group.by = "seurat_clusters", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("")
ggsave("P7_UMPplot_unlabel.pdf", plot = plot_grid(UMPplot_unlabel), width = 10, height = 8)

##相同Genotype的分组结果
table(scRNAsub@meta.data[["group"]])
##按group区别
UMPplot_label_group <- DimPlot(scRNAsub, reduction = "umap", label = TRUE, group.by = "group", pt.size = 1.0, label.size = 7, cols = c("#4DBBD5","#E64B35","#00A087","#3C5488"), raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("")
ggsave("P8_UMPplot_label_group.pdf", plot = plot_grid(UMPplot_label_group), width = 10, height = 8)
UMPplot_unlabel_group <- DimPlot(scRNAsub, reduction = "umap", label = FALSE, group.by = "group", pt.size = 1.0, label.size = 7, cols = c("#4DBBD5","#E64B35","#00A087","#3C5488"), raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("")
ggsave("P8_UMPplot_unlabel_group.pdf", plot = plot_grid(UMPplot_unlabel_group), width = 10, height = 8)

##各cell cluster--Find markers for every cluster compared to all remaining cells
##这里需要将默认格式输出调整为“RNA”
DefaultAssay(scRNAsub) <- "RNA"
scRNAsub <- ScaleData(scRNAsub, vars.to.regress = c("nCount_RNA"), verbose = TRUE)

## report only the positive ones
scRNAsub.markers <- FindAllMarkers(scRNAsub,
                               only.pos = FALSE,
                               min.pct = 0.4,
                               logfc.threshold = 0.25)

Allmarkers_scRNAsub = scRNAsub.markers %>% select(gene, everything()) %>% subset(p_val<0.05)
write.csv(Allmarkers_scRNAsub, "P9_Allmarkers__scRNAsub_wilcox.csv", row.names = T)

## 使用Top 5基因画热图
top5scRNAsub.markers <- scRNAsub.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)
##做Heatmap的热图
Heatmap_markers_scRNAsub <- DoHeatmap(scRNAsub,features = top5scRNAsub.markers$gene,
          group.colors = col5) + theme(text = element_text (size = 18)) + 
		  scale_fill_gradient2(low = '#0099CC',mid = 'white',high = '#CC0033',name = 'Z-score')
ggplot2::ggsave(filename = 'P9_Heatmap_markers_scRNAsub.pdf',width =16,height = 19)

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
ggsave(filename = 'P9_jjVolcano_scRNAsub.markers.pdf', plot = plot_grid(jjVolcano_scRNAsub.markers), width = 15, height = 6)

#Old Codes
a <- FeaturePlot(scRNAsub,features = c("GFP"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("GFP")
ggsave("P10_scRNAsubUMAP1_GFPplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(scRNAsub,features = c("Ras85D"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Ras85D")
ggsave("P10_scRNAsubUMAP1_Ras85Dplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(scRNAsub,features = c("M6"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("M6")
ggsave("P10_scRNAsubUMAP1_M6plot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(scRNAsub,features = c("GAL80"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("GAL80")
ggsave("P10_scRNAsubUMAP1_GAL80plot2.pdf", plot = plot_grid(a), width = 6, height = 5)

a <- FeaturePlot(scRNAsub,features = c("pnr", "Hml", "pigs","Ppn","kuz","He","srp","NimC1","atilla","lz","Mmp1","IM18","CecA2","CecC","Mtk","DptB","Drs","Prx2540-1","Prx2540-2","CG12896","Abl","Snoo","Ubx","CG15550","CG6023","mthl7","Cys","CG8860","COX8"),cols = c('lightgrey','red'))
ggsave("P10_scRNAsubUMAP1_Hemocytes_Feature_plot.pdf", plot = plot_grid(a), width = 18, height = 16)
a <- FeaturePlot(scRNAsub,features = c("repo", "CG3168","Gli", "nrv2", "sty","moody", "NK7.1", "CG9336"),cols = c('lightgrey','red'))
ggsave("P10_scRNAsubUMAP1_Nerve_Feature_plot.pdf", plot = plot_grid(a), width = 12, height = 10)

saveRDS(scRNAsub, file="CCA_WT_R79_Unique_Cluster.Rds")   # "RNA" scaled

###剔除Hemocytes和Glia细胞
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

################
pbmc <- readRDS("CCA_WT_R79_Unique_Cluster.Rds")

#剔除Hemocytes和Glia细胞
levels(pbmc)
Cells.sub <- subset(pbmc@meta.data, integrated_snn_res.1.2 %in% c("0","1","2","3","4","5","7","8","9","11","12","13","14","16"))

scRNAsub_Epithelium <- subset(pbmc, cells=rownames(Cells.sub))

saveRDS(scRNAsub_Epithelium, file="CCA_WTR79_UniqueClusters_reclustered_Epithelium.Rds") 



#提取Hemocytes
Cells.sub <- subset(pbmc@meta.data, integrated_snn_res.1.2 %in% c("19","20","21"))

scRNAsub_Hemocytes <- subset(pbmc, cells=rownames(Cells.sub))

saveRDS(scRNAsub_Hemocytes, file="CCA_WTRasM6_UniqueClusters_reclustered_Hemocytes.Rds") 

#提取Nerve cells
Cells.sub <- subset(pbmc@meta.data, integrated_snn_res.1.2 %in% c("13"))

scRNAsub_Nerve <- subset(pbmc, cells=rownames(Cells.sub))

saveRDS(scRNAsub_Nerve, file="CCA_WTRasM6_UniqueClusters_reclustered_Nerve.Rds") 


