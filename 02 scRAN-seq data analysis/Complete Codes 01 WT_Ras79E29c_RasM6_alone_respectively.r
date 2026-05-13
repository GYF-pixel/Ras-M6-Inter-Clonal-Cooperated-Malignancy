#Step 2 WT_Ras79E29c_RasM6_alone_respectively
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

pbmc_WT <- readRDS("CCA_WT.Rds")
pbmc_R79 <- readRDS("CCA_Ras79E29c.Rds")
pbmc_RM6 <- readRDS("CCA_RasM6.Rds")

#WT Group
pbmc <- pbmc_WT

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
UMPplot_label_group <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "group", pt.size = 1.0, label.size = 7, cols = c("#4DBBD5","#E64B35","#00A087","#3C5488"), raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("")
ggsave("P2_UMPplot_label_group.pdf", plot = plot_grid(UMPplot_label_group), width = 10, height = 8)
UMPplot_unlabel_group <- DimPlot(pbmc, reduction = "umap", label = FALSE, group.by = "group", pt.size = 1.0, label.size = 7, cols = c("#4DBBD5","#E64B35","#00A087","#3C5488"), raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("")
ggsave("P2_UMPplot_unlabel_group.pdf", plot = plot_grid(UMPplot_unlabel_group), width = 10, height = 8)

##各cell cluster的细胞占比
table(pbmc@meta.data$group, pbmc@meta.data$integrated_snn_res.0.7)
##生成文件画出堆积柱状图
Cells_in_different_Clusters_at_different_group <- table(pbmc@meta.data$group, pbmc@meta.data$integrated_snn_res.0.7)
write.csv(Cells_in_different_Clusters_at_different_group,file = "P3_Cells_in_different_Clusters_at_different_group.csv",row.names = T)

##各cell cluster--Find markers for every cluster compared to all remaining cells
##这里需要将默认格式输出调整为“RNA”
DefaultAssay(pbmc) <- "RNA"
pbmc <- ScaleData(pbmc, vars.to.regress = c("nCount_RNA"), verbose = TRUE)

## report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc,
                               only.pos = FALSE,
                               min.pct = 0.4,
                               logfc.threshold = 0.25)

Allmarkers_pbmc = pbmc.markers %>% select(gene, everything()) %>% subset(p_val<0.05)
write.csv(Allmarkers_pbmc, "P4_Allmarkers__pbmc_wilcox.csv", row.names = T)

## 使用Top 5基因画热图
top5pbmc.markers <- pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)
##做Heatmap的热图
Heatmap_markers_pbmc <- DoHeatmap(pbmc,features = top5pbmc.markers$gene,
          group.colors = col5) + theme(text = element_text (size = 18)) + 
		  scale_fill_gradient2(low = '#0099CC',mid = 'white',high = '#CC0033',name = 'Z-score')
ggplot2::ggsave(filename = 'P4_Heatmap_markers_pbmc.pdf',width =16,height = 19)

##做jjVolcano的火山图
library(scRNAtoolVis)
jjVolcano_scRNAsub.markers <- jjVolcano(diffData = pbmc.markers,
          log2FC.cutoff = 0.25, 
          size  = 3.5, #设置点的大小
          fontface = 'italic', #设置字体形式
          aesCol = c('#00468B','#ED0000'), #设置点的颜色
          tile.col = col5, #设置cluster的颜色
          #col.type = "adjustP", #设置矫正方式
          topGeneN = 5 #设置展示topN的基因
         )
jjVolcano_scRNAsub.markers
ggsave(filename = 'P4_jjVolcano_scRNAsub.markers.pdf', plot = plot_grid(jjVolcano_scRNAsub.markers), width = 15, height = 6)

##做特定基因的FeaturePlot
library(Seurat)
library(tidyverse)
library(scCustomize) # 需要Seurat版本4.3.0
library(viridis)
library(RColorBrewer)
library(gridExtra)

#scCustomize的使用方法
https://mp.weixin.qq.com/s/Y8sPLN05SKXHtjFQtFTBrA
https://mp.weixin.qq.com/s/RMUTJ3xRKPoJZJ787o32TQ
https://mp.weixin.qq.com/s/wkhmYO7kjC7dp7usxu4HPA
https://mp.weixin.qq.com/s/g9Zu7o9bZ-lx_QRdBlidkw
https://mp.weixin.qq.com/s/pdNlIV7FIv-mFUhsG30QmA

PalettePlot(viridis_plasma_dark_high)
PalettePlot(viridis_light_high)
#默认na_cutoff = 1e-09, na_color = "lightgray"
a <- FeaturePlot_scCustom(pbmc,features = c("GFP"),colors_use = viridis_light_high,na_color = "lightgray",na_cutoff = 0.5)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("GFP")
ggsave("P5_pbmcUMAP1_GFPplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot_scCustom(pbmc,features = c("Ras85D"),colors_use = viridis_light_high,na_color = "lightgray",na_cutoff = 0.5)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Ras85D")
ggsave("P5_pbmcUMAP1_Ras85Dplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot_scCustom(pbmc,features = c("M6"),colors_use = viridis_light_high,na_color = "lightgray",na_cutoff = 0.5)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("M6")
ggsave("P5_pbmcUMAP1_M6plot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot_scCustom(pbmc,features = c("GAL80"),colors_use = viridis_light_high,na_color = "lightgray",na_cutoff = 0.5)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("GAL80")
ggsave("P5_pbmcUMAP1_GAL80plot2.pdf", plot = plot_grid(a), width = 6, height = 5)

#Old Codes
a <- FeaturePlot(pbmc,features = c("GFP"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("GFP")
ggsave("P6_pbmcUMAP1_GFPplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(pbmc,features = c("Ras85D"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Ras85D")
ggsave("P6_pbmcUMAP1_Ras85Dplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(pbmc,features = c("M6"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("M6")
ggsave("P6_pbmcUMAP1_M6plot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(pbmc,features = c("GAL80"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("GAL80")
ggsave("P6_pbmcUMAP1_GAL80plot2.pdf", plot = plot_grid(a), width = 6, height = 5)

a <- FeaturePlot(pbmc,features = c("pnr", "Hml", "pigs","Ppn","kuz","He","srp","NimC1","atilla","lz","Mmp1","IM18","CecA2","CecC","Mtk","DptB","Drs","Prx2540-1","Prx2540-2","CG12896","Abl","Snoo","Ubx","CG15550","CG6023","mthl7","Cys","CG8860","COX8"),cols = c('lightgrey','red'))
ggsave("P6_pbmcUMAP1_Hemocytes_Feature_plot.pdf", plot = plot_grid(a), width = 20, height = 14)
a <- FeaturePlot(pbmc,features = c("repo", "CG3168","Gli", "nrv2", "sty","moody", "NK7.1", "CG9336"),cols = c('lightgrey','red'))
ggsave("P6_pbmcUMAP1_Nerve_Feature_plot.pdf", plot = plot_grid(a), width = 12, height = 10)

saveRDS(pbmc, file="CCA_WT.Rds")   # "RNA" scaled

#R79 Group
pbmc <- pbmc_R79

##检查default assay
DefaultAssay(pbmc)

##查看有多少个cluster
levels(pbmc)

##设置对应cluster的配色数目
col5 <- colorRampPalette((pal_npg(palette = c("nrc"))(7)))(28)
show_col(col5)

##正常的UMAP降维数据结果
UMPplot_label <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "seurat_clusters", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("")
ggsave("P1_UMPplot_label.pdf", plot = plot_grid(UMPplot_label), width = 10, height = 8)
UMPplot_unlabel <- DimPlot(pbmc, reduction = "umap", label = FALSE, group.by = "seurat_clusters", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("")
ggsave("P1_UMPplot_unlabel.pdf", plot = plot_grid(UMPplot_unlabel), width = 10, height = 8)

##相同Genotype的分组结果
table(pbmc@meta.data[["group"]])
##按group区别
UMPplot_label_group <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "group", pt.size = 1.0, label.size = 7, cols = c("#4DBBD5","#E64B35","#00A087","#3C5488"), raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("")
ggsave("P2_UMPplot_label_group.pdf", plot = plot_grid(UMPplot_label_group), width = 10, height = 8)
UMPplot_unlabel_group <- DimPlot(pbmc, reduction = "umap", label = FALSE, group.by = "group", pt.size = 1.0, label.size = 7, cols = c("#4DBBD5","#E64B35","#00A087","#3C5488"), raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("")
ggsave("P2_UMPplot_unlabel_group.pdf", plot = plot_grid(UMPplot_unlabel_group), width = 10, height = 8)

##各cell cluster的细胞占比
table(pbmc@meta.data$group, pbmc@meta.data$integrated_snn_res.0.7)
##生成文件画出堆积柱状图
Cells_in_different_Clusters_at_different_group <- table(pbmc@meta.data$group, pbmc@meta.data$integrated_snn_res.0.7)
write.csv(Cells_in_different_Clusters_at_different_group,file = "P3_Cells_in_different_Clusters_at_different_group.csv",row.names = T)

##各cell cluster--Find markers for every cluster compared to all remaining cells
##这里需要将默认格式输出调整为“RNA”
DefaultAssay(pbmc) <- "RNA"
pbmc <- ScaleData(pbmc, vars.to.regress = c("nCount_RNA"), verbose = TRUE)

## report only the positive ones--FALSE
pbmc.markers <- FindAllMarkers(pbmc,
                               only.pos = FALSE,
                               min.pct = 0.4,
                               logfc.threshold = 0.25)

Allmarkers_pbmc = pbmc.markers %>% select(gene, everything()) %>% subset(p_val<0.05)
write.csv(Allmarkers_pbmc, "P4_Allmarkers__pbmc_wilcox.csv", row.names = T)

## 使用Top 5基因画热图
top5pbmc.markers <- pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)
##做Heatmap的热图
Heatmap_markers_pbmc <- DoHeatmap(pbmc,features = top5pbmc.markers$gene,
          group.colors = col5) + theme(text = element_text (size = 18)) + 
		  scale_fill_gradient2(low = '#0099CC',mid = 'white',high = '#CC0033',name = 'Z-score')
ggplot2::ggsave(filename = 'P4_Heatmap_markers_pbmc.pdf',width =16,height = 19)

##做jjVolcano的火山图
library(scRNAtoolVis)
jjVolcano_scRNAsub.markers <- jjVolcano(diffData = pbmc.markers,
          log2FC.cutoff = 0.25, 
          size  = 3.5, #设置点的大小
          fontface = 'italic', #设置字体形式
          aesCol = c('#00468B','#ED0000'), #设置点的颜色
          tile.col = col5, #设置cluster的颜色
          #col.type = "adjustP", #设置矫正方式
          topGeneN = 5 #设置展示topN的基因
         )
jjVolcano_scRNAsub.markers
ggsave(filename = 'P4_jjVolcano_scRNAsub.markers.pdf', plot = plot_grid(jjVolcano_scRNAsub.markers), width = 15, height = 6)

##做特定基因的FeaturePlot
a <- FeaturePlot_scCustom(pbmc,features = c("GFP"),colors_use = viridis_light_high,na_color = "lightgray",na_cutoff = 0.5)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("GFP")
ggsave("P5_pbmcUMAP1_GFPplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot_scCustom(pbmc,features = c("Ras85D"),colors_use = viridis_light_high,na_color = "lightgray",na_cutoff = 0.5)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Ras85D")
ggsave("P5_pbmcUMAP1_Ras85Dplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot_scCustom(pbmc,features = c("M6"),colors_use = viridis_light_high,na_color = "lightgray",na_cutoff = 0.5)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("M6")
ggsave("P5_pbmcUMAP1_M6plot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot_scCustom(pbmc,features = c("GAL80"),colors_use = viridis_light_high,na_color = "lightgray",na_cutoff = 0.5)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("GAL80")
ggsave("P5_pbmcUMAP1_GAL80plot2.pdf", plot = plot_grid(a), width = 6, height = 5)

#Old Codes
a <- FeaturePlot(pbmc,features = c("GFP"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("GFP")
ggsave("P6_pbmcUMAP1_GFPplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(pbmc,features = c("Ras85D"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Ras85D")
ggsave("P6_pbmcUMAP1_Ras85Dplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(pbmc,features = c("M6"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("M6")
ggsave("P6_pbmcUMAP1_M6plot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(pbmc,features = c("GAL80"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("GAL80")
ggsave("P6_pbmcUMAP1_GAL80plot2.pdf", plot = plot_grid(a), width = 6, height = 5)

a <- FeaturePlot(pbmc,features = c("pnr", "Hml", "pigs","Ppn","kuz","He","srp","NimC1","atilla","lz","Mmp1","IM18","CecA2","CecC","Mtk","DptB","Drs","Prx2540-1","Prx2540-2","CG12896","Abl","Snoo","Ubx","CG15550","CG6023","mthl7","Cys","CG8860","COX8"),cols = c('lightgrey','red'))
ggsave("P6_pbmcUMAP1_Hemocytes_Feature_plot.pdf", plot = plot_grid(a), width = 18, height = 18)
a <- FeaturePlot(pbmc,features = c("repo", "CG3168","Gli", "nrv2", "sty","moody", "NK7.1", "CG9336"),cols = c('lightgrey','red'))
ggsave("P6_pbmcUMAP1_Nerve_Feature_plot.pdf", plot = plot_grid(a), width = 12, height = 10)

saveRDS(pbmc, file="CCA_Ras79E29c.Rds")   # "RNA" scaled

#RM6 Group
pbmc <- pbmc_RM6

##检查default assay
DefaultAssay(pbmc)

##查看有多少个cluster
levels(pbmc)

##设置对应cluster的配色数目
col5 <- colorRampPalette((pal_npg(palette = c("nrc"))(7)))(25)
show_col(col5)

##正常的UMAP降维数据结果
UMPplot_label <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "seurat_clusters", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("")
ggsave("P1_UMPplot_label.pdf", plot = plot_grid(UMPplot_label), width = 10, height = 8)
UMPplot_unlabel <- DimPlot(pbmc, reduction = "umap", label = FALSE, group.by = "seurat_clusters", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("")
ggsave("P1_UMPplot_unlabel.pdf", plot = plot_grid(UMPplot_unlabel), width = 10, height = 8)

##相同Genotype的分组结果
table(pbmc@meta.data[["group"]])
##按group区别
UMPplot_label_group <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "group", pt.size = 1.0, label.size = 7, cols = c("#4DBBD5","#E64B35","#00A087","#3C5488"), raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("")
ggsave("P2_UMPplot_label_group.pdf", plot = plot_grid(UMPplot_label_group), width = 10, height = 8)
UMPplot_unlabel_group <- DimPlot(pbmc, reduction = "umap", label = FALSE, group.by = "group", pt.size = 1.0, label.size = 7, cols = c("#4DBBD5","#E64B35","#00A087","#3C5488"), raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("")
ggsave("P2_UMPplot_unlabel_group.pdf", plot = plot_grid(UMPplot_unlabel_group), width = 10, height = 8)

##各cell cluster的细胞占比
table(pbmc@meta.data$group, pbmc@meta.data$integrated_snn_res.0.7)
##生成文件画出堆积柱状图
Cells_in_different_Clusters_at_different_group <- table(pbmc@meta.data$group, pbmc@meta.data$integrated_snn_res.0.7)
write.csv(Cells_in_different_Clusters_at_different_group,file = "P3_Cells_in_different_Clusters_at_different_group.csv",row.names = T)

##各cell cluster--Find markers for every cluster compared to all remaining cells
##这里需要将默认格式输出调整为“RNA”
DefaultAssay(pbmc) <- "RNA"
pbmc <- ScaleData(pbmc, vars.to.regress = c("nCount_RNA"), verbose = TRUE)

## report only the positive ones--FALSE
pbmc.markers <- FindAllMarkers(pbmc,
                               only.pos = FALSE,
                               min.pct = 0.4,
                               logfc.threshold = 0.25)

Allmarkers_pbmc = pbmc.markers %>% select(gene, everything()) %>% subset(p_val<0.05)
write.csv(Allmarkers_pbmc, "P4_Allmarkers__pbmc_wilcox.csv", row.names = T)

## 使用Top 5基因画热图
top5pbmc.markers <- pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)
##做Heatmap的热图
Heatmap_markers_pbmc <- DoHeatmap(pbmc,features = top5pbmc.markers$gene,
          group.colors = col5) + theme(text = element_text (size = 18)) + 
		  scale_fill_gradient2(low = '#0099CC',mid = 'white',high = '#CC0033',name = 'Z-score')
ggplot2::ggsave(filename = 'P4_Heatmap_markers_pbmc.pdf',width =16,height = 19)

##做jjVolcano的火山图
library(scRNAtoolVis)
jjVolcano_scRNAsub.markers <- jjVolcano(diffData = pbmc.markers,
          log2FC.cutoff = 0.25, 
          size  = 3.5, #设置点的大小
          fontface = 'italic', #设置字体形式
          aesCol = c('#00468B','#ED0000'), #设置点的颜色
          tile.col = col5, #设置cluster的颜色
          #col.type = "adjustP", #设置矫正方式
          topGeneN = 5 #设置展示topN的基因
         )
jjVolcano_scRNAsub.markers
ggsave(filename = 'P4_jjVolcano_scRNAsub.markers.pdf', plot = plot_grid(jjVolcano_scRNAsub.markers), width = 15, height = 6)

##做特定基因的FeaturePlot
a <- FeaturePlot_scCustom(pbmc,features = c("GFP"),colors_use = viridis_light_high,na_color = "lightgray",na_cutoff = 0.5)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("GFP")
ggsave("P5_pbmcUMAP1_GFPplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot_scCustom(pbmc,features = c("Ras85D"),colors_use = viridis_light_high,na_color = "lightgray",na_cutoff = 0.5)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Ras85D")
ggsave("P5_pbmcUMAP1_Ras85Dplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot_scCustom(pbmc,features = c("M6"),colors_use = viridis_light_high,na_color = "lightgray",na_cutoff = 0.5)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("M6")
ggsave("P5_pbmcUMAP1_M6plot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot_scCustom(pbmc,features = c("GAL80"),colors_use = viridis_light_high,na_color = "lightgray",na_cutoff = 0.5)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("GAL80")
ggsave("P5_pbmcUMAP1_GAL80plot2.pdf", plot = plot_grid(a), width = 6, height = 5)

a <- FeaturePlot(pbmc,features = c("GFP"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("GFP")
ggsave("P6_pbmcUMAP1_GFPplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(pbmc,features = c("Ras85D"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Ras85D")
ggsave("P6_pbmcUMAP1_Ras85Dplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(pbmc,features = c("M6"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("M6")
ggsave("P6_pbmcUMAP1_M6plot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(pbmc,features = c("GAL80"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("GAL80")
ggsave("P6_pbmcUMAP1_GAL80plot2.pdf", plot = plot_grid(a), width = 6, height = 5)

a <- FeaturePlot(pbmc,features = c("pnr", "Hml", "pigs","Ppn","kuz","He","srp","NimC1","atilla","lz","Mmp1","IM18","CecA2","CecC","Mtk","DptB","Drs","Prx2540-1","Prx2540-2","CG12896","Abl","Snoo","Ubx","CG15550","CG6023","mthl7","Cys","CG8860","COX8"),cols = c('lightgrey','red'))
ggsave("P6_pbmcUMAP1_Hemocytes_Feature_plot.pdf", plot = plot_grid(a), width = 20, height = 18)
a <- FeaturePlot(pbmc,features = c("repo", "CG3168","Gli", "nrv2", "sty","moody", "NK7.1", "CG9336"),cols = c('lightgrey','red'))
ggsave("P6_pbmcUMAP1_Nerve_Feature_plot.pdf", plot = plot_grid(a), width = 12, height = 10)

saveRDS(pbmc, file="CCA_RasM6.Rds")  # "RNA" scaled
##单细胞基因集打分AddModuleScore
https://mp.weixin.qq.com/s/S_BnQ6qZsqJEkpL89xufGQ
https://www.jianshu.com/p/cef5663888ff
Hemocytes_anno <- AddModuleScore(object = pbmc, features = c("pnr", "Hml", "pigs","Ppn","kuz","He","srp","NimC1","atilla","lz","Mmp1","IM18","CecA2","CecC","Mtk","DptB","Drs","Prx2540-1","Prx2540-2","CG12896","Abl","Snoo","Ubx","CG15550","CG6023","mthl7","Cys","CG8860","COX8"), 
                           name =c("Hemocytes_Markers"), 
                           nbin = 12)
