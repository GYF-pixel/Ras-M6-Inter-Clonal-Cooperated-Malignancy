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
setwd("J:\\05Cooperated_Project\\02M6\\scRNAseq")


################WTRM6_Hemocytes alone
pbmc_WTRM6_Hemocytes <- readRDS("CCA_WTRasM6_UniqueClusters_reclustered_Hemocytes.Rds")
pbmc_WTRM6_Hemocytes <- FindVariableFeatures(pbmc_WTRM6_Hemocytes, selection.method = "vst", nfeatures = 2500)
pbmc_WTRM6_Hemocytes <- ScaleData(pbmc_WTRM6_Hemocytes, vars.to.regress = c("nCount_RNA"), verbose = TRUE)
pbmc_WTRM6_Hemocytes <- RunPCA(pbmc_WTRM6_Hemocytes, features = VariableFeatures(pbmc_WTRM6_Hemocytes), npcs = 40, nfeature.print = 10, ndims.print = 1:5, verbose = T)
pc.num=1:40
pbmc_WTRM6_Hemocytes <- RunUMAP(pbmc_WTRM6_Hemocytes, dims=pc.num)
pbmc_WTRM6_Hemocytes <- FindNeighbors(pbmc_WTRM6_Hemocytes, dims = pc.num)
pbmc_WTRM6_Hemocytes = FindClusters(pbmc_WTRM6_Hemocytes,resolution = 1)

pbmc <- pbmc_WTRM6_Hemocytes
##检查default assay
DefaultAssay(pbmc)

##查看有多少个cluster
levels(pbmc)

##设置对应cluster的配色数目
col5 <- colorRampPalette((pal_npg(palette = c("nrc"))(7)))(19)
show_col(col5)

##正常的UMAP降维数据结果
UMPplot_label <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "seurat_clusters", pt.size = 1.4, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Hemocytes")
ggsave("P1_UMPplot_label.pdf", plot = plot_grid(UMPplot_label), width = 5.5, height = 5)

##正常的UMAP降维数据结果 revision
UMPplot_label <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "seurat_clusters", pt.size = 1.4, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Hemocytes")
ggsave("P1_UMPplot_label.pdf", plot = plot_grid(UMPplot_label), width = 4.5, height = 4)

UMPplot_unlabel <- DimPlot(pbmc, reduction = "umap", label = FALSE, group.by = "seurat_clusters", pt.size = 1.4, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Hemocytes")
ggsave("P1_UMPplot_unlabel.pdf", plot = plot_grid(UMPplot_unlabel), width = 5.5, height =5)

##相同Genotype的分组结果
table(pbmc@meta.data[["group"]])
##按group区别
UMPplot_label_group <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "group", pt.size = 1.4, label.size = 7, cols = c("#4DBBD5","#E64B35","#00A087","#3C5488"), raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("")
ggsave("P2_UMPplot_label_group.pdf", plot = plot_grid(UMPplot_label_group), width = 10, height = 8)
UMPplot_unlabel_group <- DimPlot(pbmc, reduction = "umap", label = FALSE, group.by = "group", pt.size = 1.4, label.size = 7, cols = c("#4DBBD5","#E64B35","#00A087","#3C5488"), raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("")
ggsave("P2_UMPplot_unlabel_group.pdf", plot = plot_grid(UMPplot_unlabel_group), width = 10, height = 8)


##检查default assay
DefaultAssay(pbmc) <- "RNA"

pbmc.markers <- FindAllMarkers(pbmc,
                               only.pos = TRUE,
                               min.pct = 0.4,
                               logfc.threshold = 0.25)

Allmarkers_pbmc = pbmc.markers %>% select(gene, everything()) %>% subset(p_val<0.05)
write.csv(Allmarkers_pbmc, "P2_Allmarkers__pbmc_wilcox_WTRM6_Hemocytes.csv", row.names = T)
# get top 10 genes
top5pbmc.markers <- pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)
# plot
Heatmap_markers_pbmc <- DoHeatmap(pbmc,features = top5pbmc.markers$gene,
          group.colors = col5) + theme(text = element_text (size = 18)) + 
		  scale_fill_gradient2(low = '#0099CC',mid = 'white',high = '#CC0033',name = 'Z-score')
ggplot2::ggsave(filename = 'P2_Heatmap_markers_WTRM6_Hemocytes.pdf',width =16,height = 19)

a <- FeaturePlot(pbmc,features = c("Tl", "spz","Toll-4","egr","Toll-6","Toll-7","Toll-9"),cols = c('lightgrey','red'),ncol = 2)
ggsave("P3_pbmcUMAP1_Hemocytes_Feature_plot.pdf", plot = plot_grid(a), width = 6, height = 7)

saveRDS(pbmc, file="CCA_WTRasM6_UniqueClusters_reclustered_Hemocytes_reclustered.Rds") 

###hemocytes classification

pbmc_RM6 <- readRDS("CCA_WTRasM6_UniqueClusters_reclustered_Hemocytes_reclustered.Rds")
pbmc <- pbmc_RM6

#各个cluster的细胞占比
table(pbmc@meta.data$RNA_snn_res.1)
##生成文件画出堆积柱状图
Cells_in_different_Clusters_at_different_group <- table(pbmc@meta.data$RNA_snn_res.1)
write.csv(Cells_in_different_Clusters_at_different_group,file = "P3_Cells_in_different_Clusters_at_different_group.csv",row.names = T)

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
				#"spz4",	
				#"spz5",	
				"spz6"
				)


a <- FeaturePlot(pbmc,features = genes,cols = c('lightgrey','red'),order = T, ncol = 3)
ggsave("P3_pbmcUMAP1_Hemocytes_Feature_plot.pdf", plot = plot_grid(a), width = 11, height = 9)

a <- FeaturePlot(pbmc,features = c("egr"),cols = c('lightgray','red'),order = T)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("egr")
ggsave("P3_pbmcUMAP1_egrplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)

genes = c("spz",	"egr")	#Toll-1

a <- FeaturePlot(pbmc,features = genes,cols = c('lightgrey','red'),order = T, ncol = 2)
ggsave("P3_pbmcUMAP1_spz_egrplot2.pdf", plot = plot_grid(a), width = 8, height = 3.5)



#pan-hemocytes
a <- FeaturePlot(pbmc,features = c("srp"),cols = c('lightgray','red'),order = T)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("srp")
ggsave("P3_pbmcUMAP1_srpplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("He"),cols = c('lightgray','red'),order = T)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("He")
ggsave("P3_pbmcUMAP1_Heplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)

"He","srp",


#pan-hemocytes
a <- FeaturePlot(pbmc,features = c("Col4a1"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Col4a1")
ggsave("P3_pbmcUMAP1_Col4a1plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("Col4a1"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Col4a1")
ggsave("P3_pbmcUMAP1_Col4a1plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)


#Prohemocytes/thanacytes
a <- FeaturePlot(pbmc,features = c("Tep4"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Tep4")
ggsave("P3_pbmcUMAP1_Tep4plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("Ance"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Ance")
ggsave("P3_pbmcUMAP1_Anceplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)


a <- FeaturePlot(pbmc,features = c("CG15506"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("CG15506")
ggsave("P3_pbmcUMAP1_CG15506plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("CG1648"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("CG1648")
ggsave("P3_pbmcUMAP1_CG1648plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("Ance-2"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Ance-2")
ggsave("P3_pbmcUMAP1_Ance-2plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)

"Tep4","Ance","CG1648","CG15506","Ance-2"

a <- FeaturePlot(pbmc,features = c("IM18"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Tep4")
ggsave("P3_pbmcUMAP1_Tep4plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("CecA2"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Ance")
ggsave("P3_pbmcUMAP1_Anceplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("NimB3"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("NimB3")
ggsave("P3_pbmcUMAP1_NimB3plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("Ppa"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Ppa")
ggsave("P3_pbmcUMAP1_Ppaplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("shg"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("shg")
ggsave("P3_pbmcUMAP1_shgplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("E(spl)m3-HLH"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("E(spl)m3-HLH")
ggsave("P3_pbmcUMAP1_E(spl)m3-HLHplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("E(spl)m4-BFM"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("E(spl)m4-BFM")
ggsave("P3_pbmcUMAP1_E(spl)m4-BFMplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("Dl"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Dl")
ggsave("P3_pbmcUMAP1_Dlplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("Ilp6"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Ilp6")
ggsave("P3_pbmcUMAP1_Ilp6plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)


#Cluster1
a <- FeaturePlot(pbmc,features = c("GstD1"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("GstD1")
ggsave("P3_pbmcUMAP1_GstD1plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("Lst"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Lst")
ggsave("P3_pbmcUMAP1_Lstplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("et"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("et")
ggsave("P3_pbmcUMAP1_etplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)

a <- FeaturePlot(pbmc,features = c("Idgf6"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Idgf6")
ggsave("P3_pbmcUMAP1_Idgf6plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("Eip55E"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Eip55E")
ggsave("P3_pbmcUMAP1_Eip55Eplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("HmgZ"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("HmgZ")
ggsave("P3_pbmcUMAP1_HmgZplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("Inos"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Inos")
ggsave("P3_pbmcUMAP1_Inosplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("CG4793"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("CG4793")
ggsave("P3_pbmcUMAP1_CG4793plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("Nplp2"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Nplp2")
ggsave("P3_pbmcUMAP1_Nplp2plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("stg"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("stg")
ggsave("P3_pbmcUMAP1_stgplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("Incenp"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Incenp")
ggsave("P3_pbmcUMAP1_Incenpplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("Chd64"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Chd64")
ggsave("P3_pbmcUMAP1_Chd64plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("Pcd"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Pcd")
ggsave("P3_pbmcUMAP1_Pcdplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)

a <- FeaturePlot(pbmc,features = c("CG15347"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("CG15347")
ggsave("P3_pbmcUMAP1_CG15347plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("CG6770"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("CG6770")
ggsave("P3_pbmcUMAP1_CG6770plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("CG42394"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("CG42394")
ggsave("P3_pbmcUMAP1_CG42394plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("CG1572"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("CG1572")
ggsave("P3_pbmcUMAP1_CG1572plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("Gdh"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Gdh")
ggsave("P3_pbmcUMAP1_Gdhplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("cathD"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("cathD")
ggsave("P3_pbmcUMAP1_cathDplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("sesB"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("sesB")
ggsave("P3_pbmcUMAP1_sesBplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)

#Plasmatocytes
#"NimC1","Cg","Pxn","Hml","eater","crq"
a <- FeaturePlot(pbmc,features = c("NimC1"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("NimC1")
ggsave("P3_pbmcUMAP1_NimC1plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)  # Not included
a <- FeaturePlot(pbmc,features = c("Cg"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Cg")
ggsave("P3_pbmcUMAP1_Cgplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)  # Not included
a <- FeaturePlot(pbmc,features = c("Pxn"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Pxn")
ggsave("P3_pbmcUMAP1_Pxnplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("Hml"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Hml")
ggsave("P3_pbmcUMAP1_Hmlplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("eater"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("eater")
ggsave("P3_pbmcUMAP1_eaterplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("crq"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("crq")
ggsave("P3_pbmcUMAP1_crqplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("Ppn"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Ppn")
ggsave("P3_pbmcUMAP1_Ppnplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)

a <- FeaturePlot(pbmc,features = c("CAH7"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("CAH7")
ggsave("P3_pbmcUMAP1_CAH7plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("Lsp2"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Lsp2")
ggsave("P3_pbmcUMAP1_Lsp2plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("CG2444"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("CG2444")
ggsave("P3_pbmcUMAP1_CG2444plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)

#Crystal cells
#Early markers: lozenge (lz) and pebbled (peb); later markers:  PPO1 and PPO2

a <- FeaturePlot(pbmc,features = c("PPO2"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("PPO2")
ggsave("P3_pbmcUMAP1_PPO2plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)  # Not included
a <- FeaturePlot(pbmc,features = c("PPO1"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("PPO1")
ggsave("P3_pbmcUMAP1_PPO1plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)  # Not included
a <- FeaturePlot(pbmc,features = c("lz"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("lz")
ggsave("P3_pbmcUMAP1_lzplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("peb"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("peb")
ggsave("P3_pbmcUMAP1_pebplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("CG9119"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("CG9119")
ggsave("P3_pbmcUMAP1_CG9119plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4) 		# Not included
a <- FeaturePlot(pbmc,features = c("fok"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("fok")
ggsave("P3_pbmcUMAP1_fokplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)

#Lamellocytes
#mys is also expressed by hemocyte progenitors and plasmatocytes
a <- FeaturePlot(pbmc,features = c("atilla"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("atilla")
ggsave("P3_pbmcUMAP1_atillaplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("PPO3"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("PPO3")
ggsave("P3_pbmcUMAP1_PPO3plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("ItgaPS4"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("ItgaPS4")
ggsave("P3_pbmcUMAP1_ItgaPS4plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)  # Not included
a <- FeaturePlot(pbmc,features = c("msn"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("msn")
ggsave("P3_pbmcUMAP1_msnplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("mys"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("mys")
ggsave("P3_pbmcUMAP1_mysplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("CG33225"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("CG33225")
ggsave("P3_pbmcUMAP1_CG33225plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("alphaTub85E"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("alphaTub85E")
ggsave("P3_pbmcUMAP1_alphaTub85Eplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("betaTub60D"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("betaTub60D")
ggsave("P3_pbmcUMAP1_betaTub60Dplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)

#JGG novel markers for Lamellocytes
a <- FeaturePlot(pbmc,features = c("mthl4"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("mthl4")
ggsave("P3_pbmcUMAP1_mthl4plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("CAP"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("CAP")
ggsave("P3_pbmcUMAP1_CAPplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("shot"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("shot")
ggsave("P3_pbmcUMAP1_shotplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)

#thanacyte	Tep4


#PSC/Primocytes/PM11/PL−ImpL2
# the early hemocyte marker srp, no detectable levels of more mature hemocyte markers like He, Pxn and crq

a <- FeaturePlot(pbmc,features = c("Antp"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Antp")
ggsave("P3_pbmcUMAP1_Antpplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("kn"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("kn")
ggsave("P3_pbmcUMAP1_knplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("tau"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("tau")
ggsave("P3_pbmcUMAP1_tauplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("ham"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("ham")
ggsave("P3_pbmcUMAP1_hamplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("Mad"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Mad")
ggsave("P3_pbmcUMAP1_Madplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("CG15550"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("CG15550")
ggsave("P3_pbmcUMAP1_CG15550plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("CG6023"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("CG6023")
ggsave("P3_pbmcUMAP1_CG6023plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("beat-IIb"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("beat-IIb")
ggsave("P3_pbmcUMAP1_beat-IIbplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("Fas1"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Fas1")
ggsave("P3_pbmcUMAP1_Fas1plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("ImpL2"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("ImpL2")
ggsave("P3_pbmcUMAP1_ImpL2plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)






"Antp","crq","eater","LpR2","Lsd-2","Sirup"
a <- FeaturePlot(pbmc,features = c("Klp61F","Klp67A","Ncd","Ncd80","Nuf2","CycE","CycB"),cols = c('lightgrey','red'))
ggsave("P10_scRNAsubUMAP1_Hemocytes_PL_Inos&PL_Proliff_Feature_plot.pdf", plot = plot_grid(a), width = 10, height = 8)

a <- FeaturePlot(pbmc,features = c("pnr", "Hml", "pigs","Ppn","kuz","He","srp","NimC1","atilla","lz","Mmp1","IM18","CecA2","CecC","Mtk","DptB","Drs","Prx2540-1","Prx2540-2","CG12896","Abl","Snoo","Ubx","CG15550","CG6023","mthl7","Cys","CG8860","COX8"),cols = c('lightgrey','red'))
ggsave("P10_scRNAsubUMAP1_Hemocytes_Feature_plot.pdf", plot = plot_grid(a), width = 18, height = 16)
a <- FeaturePlot(pbmc,features = c("Gal","CG10038","tau","CG15550","Obp99b","Lsp1beta","Lsp1alpha","CG34166","Fbp2","Drs","vir-1","CG42807","CG43236","Arc1",
"Mmp1","Rel","CG34296","ImpL2","Obp99a","CecC","CecA2","CecA1","CG5399","Ten-a","Eip93F","E(spl)mbeta-HLH","Nplp2","Hml","NimB4","Ance-5","NLaz","Snp","Eip71CD","TwdlE",
"Sur-8","CG4793","Inos","HmgZ","Eip55E","Idgf6","robo2","NimC1","fat-spondin","NimC2","Kn","CG34437","Stg","Incenp","Ppn","Chd64","CG31431","Ance","Pcd","CG9119","Men","MtnA","PPO1","PPO2"),cols = c('lightgrey','red'))
ggsave("P10_scRNAsubUMAP1_Hemocytes_Feature_plot.pdf", plot = plot_grid(a), width = 18, height = 25)

a <- FeaturePlot(pbmc,features = c("Timp","CecC","sesB","cathD","Gdh","CG31777","MtnA","CG1572","ATP8B","CG42394","CG34330","CG34253","Tsp42Ed","NLaz","Drat",
"CG6770","AANATL2","CG9336","CG34423","CG14072","Drip","shot","RapGAP1","ItgaPS4","Treh","rhea","CG1208","CG15347","alphaTub85E"),cols = c('lightgrey','red'))
ggsave("P10_scRNAsubUMAP1_Hemocytes_Lamanocytes_Feature_plot.pdf", plot = plot_grid(a), width = 18, height = 17)



#Plasmatocytes PM3
a <- FeaturePlot(pbmc,features = c("Hml"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Hml")
ggsave("P3_pbmcUMAP1_Hmlplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("CG5399"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("CG5399")
ggsave("P3_pbmcUMAP1_CG5399plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("miple2"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("miple2")
ggsave("P3_pbmcUMAP1_miple2plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)



#PM11
a <- FeaturePlot(pbmc,features = c("CG15550"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("CG15550")
ggsave("P3_pbmcUMAP1_CG15550plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("CG8860"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("CG8860")
ggsave("P3_pbmcUMAP1_CG8860plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("CG6023"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("CG6023")
ggsave("P3_pbmcUMAP1_CG6023plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("mthl12"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("mthl12")
ggsave("P3_pbmcUMAP1_mthl12plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)

#Editor required figure
a <- FeaturePlot(pbmc,features = c("Nplp2"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("nplp2")
ggsave("P3_pbmcUMAP1_nplp2plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("Tep4"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Tep4")
ggsave("P3_pbmcUMAP1_Tep4plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("Idgf6"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Idgf6")
ggsave("P3_pbmcUMAP1_Idgf6plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("et"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("et")
ggsave("P3_pbmcUMAP1_etplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)

#PL-ImpL2
a <- FeaturePlot(pbmc,features = c("Nlaz"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("nplp2")
ggsave("P3_pbmcUMAP1_nplp2plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("Eip93F"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Idgf6")
ggsave("P3_pbmcUMAP1_Eip93Fplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)


a <- FeaturePlot(pbmc,features = c("kn "),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("kn")
ggsave("P3_pbmcUMAP1_knplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("CG15550 "),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("CG15550")
ggsave("P3_pbmcUMAP1_CG15550plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("CG10038 "),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("CG10038")
ggsave("P3_pbmcUMAP1_CG10038plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)



########################################RM6 Hemocyte 细胞注释 DotPlot
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
setwd("I:/M6")

pbmc_RM6 <- readRDS("CCA_WTRasM6_UniqueClusters_reclustered_Hemocytes_reclustered.Rds")
pbmc <- pbmc_RM6

#每种celltype的marker基因/或者需要展示的基因
DefaultAssay(pbmc) <- "RNA"

Pan_hemocytes_genes <- c("srp","He","Cg")
Plasmatocytes_genes <- c("NimC1","Pxn","Hml","eater","crq","Ppn","mys")
Crystal_cell_genes <- c("lz","peb","fok","PPO1","PPO2","CG9119")
Lamellocytes_genes <- c("atilla","PPO3","alphaTub85E","ItgaPS4","msn","CG33225","shot","mthl4","CAP")   	#"betaTub60D", excluded
PSC_Primocytes_genes <- c("Antp","kn","tau","CG15550","CG6023","CG4793","CG6770","CG10038","ImpL2","HmgZ","Fas1","beat-IIb")
PL_prolif_genes <- c("Ance","CG42394","Tep4","GstD1","Idgf6","et","Eip55E","Lst","CG1648","Nplp2")
PL_Inos_genes <- c("Inos","Chd64","CG1572","Gdh","Incenp","stg")
PL_0_PL_robo2_genes <- c("Eip71CD","fat-spondin","Nlaz","Gal","Sur-8","TwdlE","E(spl)mbeta-HLH","Drs","CG5399","Snp")

#"cathD" excluded
#

##设置celltype
levels(pbmc)
levels(pbmc@active.ident)

pbmc$celltype <- pbmc@active.ident

new.cluster.ids <- c("0"="PSC_1",
								"3"="PSC_2",
								"1"="PL_prolif-like cells",
								"2"="Plasmatocytes")  #(PL_0&PL_robo2-like cells)

pbmc <- RenameIdents(pbmc, new.cluster.ids)    

levels(pbmc) 
#pbmc$celltype <- pbmc@active.ident

#基因储存为一个list，这个list是DotPlot的输入，才会实现分面的效果
features <- list(#"Pan_hemocytes" = Pan_hemocytes_genes,
                 "PSC_Primocytes" = PSC_Primocytes_genes,
				 "PL_Inos_genes" = PL_Inos_genes,
				 "PL_prolif" = PL_prolif_genes,
				 "Plasmatocytes" = Plasmatocytes_genes,
				 "PL_0_PL_robo2_genes" = PL_0_PL_robo2_genes,	 
				 "Lamellocytes" = Lamellocytes_genes,
                 "Crystal_cell" = Crystal_cell_genes
)

#设置celltype展示顺序
levels(pbmc)  <- c("PSC_1","PSC_2","PL_prolif-like cells","Plasmatocytes")


#作图修饰
Adotplot <- DotPlot(object = pbmc, features=features)&
  theme_bw()& #设置主题
  geom_point(shape=21, aes(size=pct.exp),stroke=1)& #给点添加边框
  theme(axis.title = element_blank(),
        axis.text.x = element_text(color = 'black', size = 10, angle = 90, hjust = 1, vjust = 0.5, face = "bold"),
        axis.text.y = element_text(color = 'black', size = 12, face = "bold"),
        # panel.grid.major = element_blank(), 
        # panel.grid.minor = element_blank(),
        strip.background = element_blank(), #去除分面图背景
        strip.text = element_blank(),#去除分面图文字
        plot.margin=unit(c(1, 1, 1, 1),'cm'),#缩小作图范围
        panel.border = element_rect(color="black",size = 1.2, linetype="solid"),#修改边框大小
        panel.spacing = unit(0.12, "cm"),#修改分面图间距
        # legend.frame = element_rect(colour = "black"),#修改legend边框
        # legend.ticks = element_line(colour = "black", linewidth  = 0),
        legend.key.width = unit(0.3, "cm"),#修改legend宽度
        legend.key.height = unit(0.5, "cm"),#修改legend高度
        legend.title = element_text(color = 'black', face = "bold", size=9))& #修改legend标题字体大小
  scale_color_gradientn(colours = colorRampPalette(c("navy","white","firebrick3"))(100))& #修改连续legend颜色
  labs(tag = "Cell type annotation")& #添加标签
  theme(plot.tag.position = c(0.3, 1.05), #设置标签位置，自行调整
        plot.tag = element_text(size = 12,face = "bold"))& #设置标签文字大小
  guides(size=guide_legend(title="Proportion of\nexpressing cells"), #设置legend size的标题
         colour=guide_colorbar(title="Average\nexpression"))#设置legend color bar的标题

Adotplot

#45 marker genes 14 cm
ggsave("P3_Hemocytes_marker_gene_Adotplot.pdf",plot=Adotplot,width = 14, height = 4.5)





















################
pbmc_WT_Bravo <- readRDS("00WT_Bravo_EAD//WT_Bravo_Hemocytes.Rds")
pbmc_WT_Bravo@meta.data$group = "WT_Bravo"
pbmc_WT_Bravo@meta.data$sec_group = "WT"

pbmc_WT_XGY <- readRDS("CCA_WT_Hemocytes.Rds")

pbmc_WT <- merge(pbmc_WT_XGY, y=c(pbmc_WT_Bravo), 
             add.cell.ids = c("WT_XGY","WT_Bravo"),
             merge.data = TRUE)
			 
pbmc_RM6 <- readRDS("CCA_WTRasM6_UniqueClusters_reclustered_Hemocytes_reclustered.Rds")

###WT and RM6
DefaultAssay(pbmc_WT) <- "RNA"
DefaultAssay(pbmc_RM6) <- "RNA"

WTRM6.anchors <- FindIntegrationAnchors(object.list = list(pbmc_WT, pbmc_RM6), anchor.features = 2500, dims = 1:50)
WTRM6.combined <- IntegrateData(anchorset = WTRM6.anchors, dims = 1:50)
pbmc <-  WTRM6.combined

#pbmc <- NormalizeData(pbmc, verbose = FALSE, normalization.method = "LogNormalize", scale.factor = 1e4)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2500)
pbmc <- ScaleData(pbmc, vars.to.regress = c("nCount_RNA"), verbose = TRUE)
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc), npcs = 40, nfeature.print = 10, ndims.print = 1:5, verbose = T)
pc.num=1:40
pbmc <- RunUMAP(pbmc, dims=pc.num)
pbmc <- FindNeighbors(pbmc, dims = pc.num)
pbmc = FindClusters(pbmc,resolution = 1.0)

saveRDS(pbmc, file="CCA_WT_RasM6_Bravo_XGY_Hemocytes_Reculstered.Rds") 

#pbmc <- readRDS("CCA_WT_RasM6_Bravo_XGY_Hemocytes_Reculstered.Rds")

##检查default assay
DefaultAssay(pbmc)
##查看有多少个cluster
levels(pbmc)
##设置对应cluster的配色数目
col5 <- colorRampPalette((pal_npg(palette = c("nrc"))(7)))(8)
show_col(col5)

##正常的UMAP降维数据结果
UMPplot_label <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "seurat_clusters", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("WT & RasM6 Integrated Hemocytes")
ggsave("P1_UMPplot_label.pdf", plot = plot_grid(UMPplot_label), width = 6, height = 5)
UMPplot_unlabel <- DimPlot(pbmc, reduction = "umap", label = FALSE, group.by = "seurat_clusters", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("WT & RasM6 Integrated Hemocytes")
ggsave("P1_UMPplot_unlabel.pdf", plot = plot_grid(UMPplot_unlabel), width = 6, height = 5)
table(pbmc@meta.data[["seurat_clusters"]])

##相同Genotype的分组结果
table(pbmc@meta.data[["group"]])
##按group区别
UMPplot_label_group <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "group", pt.size = 1.0, label.size = 7,  cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("WT & RasM6 Integrated Hemocytes")
ggsave("P2_UMPplot_label_group.pdf", plot = plot_grid(UMPplot_label_group), width = 7, height = 5)
UMPplot_unlabel_group <- DimPlot(pbmc, reduction = "umap", label = FALSE, group.by = "group", pt.size = 1.0, label.size = 7,  cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("WT & RasM6 Integrated Hemocytes")
ggsave("P2_UMPplot_unlabel_group.pdf", plot = plot_grid(UMPplot_unlabel_group), width = 7, height = 5)

##相同Genotype的分组结果
table(pbmc@meta.data[["sec_group"]])
##按sec_group区别
UMPplot_label_sec_group <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "sec_group", pt.size = 1.0, label.size = 7, cols = c("#4DBBD5","#E64B35","#00A087","#3C5488"), raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("WT & RasM6 Integrated Hemocytes")
ggsave("P2_UMPplot_label_sec_group.pdf", plot = plot_grid(UMPplot_label_sec_group), width = 6, height = 5)
UMPplot_unlabel_sec_group <- DimPlot(pbmc, reduction = "umap", label = FALSE, group.by = "sec_group", pt.size = 1.0, label.size = 7, cols = c("#4DBBD5","#E64B35","#00A087","#3C5488"), raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("WT & RasM6 Integrated Hemocytes")
ggsave("P2_UMPplot_unlabel_sec_group.pdf", plot = plot_grid(UMPplot_unlabel_sec_group), width = 6, height = 5)

##这里需要将默认格式输出调整为“RNA”
DefaultAssay(pbmc) <- "RNA"
pbmc <- ScaleData(pbmc, vars.to.regress = c("nCount_RNA"), verbose = TRUE)

## report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc,
                               only.pos = FALSE,
                               min.pct = 0.4,
                               logfc.threshold = 0.25)

Allmarkers_pbmc = pbmc.markers %>% select(gene, everything()) %>% subset(p_val<0.05)
write.csv(Allmarkers_pbmc, "P3_Allmarkers__pbmc_wilcox.csv", row.names = T)

## 使用Top 5基因画热图
top5pbmc.markers <- pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 6, wt = avg_log2FC)
##做Heatmap的热图
Heatmap_markers_pbmc <- DoHeatmap(pbmc,features = top5pbmc.markers$gene,
          group.colors = col5) + theme(text = element_text (size = 20)) + 
		  scale_fill_gradient2(low = '#00468B',mid = 'white',high = '#CC0033',name = 'Z-score')
ggplot2::ggsave(filename = 'P3_Heatmap_markers_pbmc.pdf',width =14,height = 10)

##做jjVolcano的火山图
library(scRNAtoolVis)
jjVolcano_pbmc.markers <- jjVolcano(diffData = pbmc.markers,
          log2FC.cutoff = 0.25, 
          size  = 3.5, #设置点的大小
          fontface = 'italic', #设置字体形式
          aesCol = c('#00468B','#ED0000'), #设置点的颜色
          tile.col = col5, #设置cluster的颜色
		  legend.position = c(0.7,0), #图注位置
          #col.type = "adjustP", #设置矫正方式
          topGeneN = 5 #设置展示topN的基因
         )
jjVolcano_pbmc.markers
ggsave(filename = 'P3_jjVolcano_pbmc.markers.pdf', plot = plot_grid(jjVolcano_pbmc.markers), width = 6, height = 6)


#RM6 Hemocytes Cluster 0 2 5 6 7 vs WT_Epithelial
clusterHem.markers <- FindMarkers(pbmc, ident.1 = c("0","2","5","6","7"), ident.2 = c("1","3","4"), min.pct = 0.25)
head(clusterHem.markers, n = 5)
write.csv(clusterHem.markers, "P5_Unique_Hem.markers_vs_WT.csv", row.names = T)

clusterHem0.markers <- FindMarkers(pbmc, ident.1 = c("0"), ident.2 = c("1","3","4"), min.pct = 0.25)
head(clusterHem0.markers, n = 5)
write.csv(clusterHem0.markers, "P5_Unique_Hem0.markers_vs_WT.csv", row.names = T)

clusterHem2.markers <- FindMarkers(pbmc, ident.1 = c("2"), ident.2 = c("1","3","4"), min.pct = 0.25)
head(clusterHem2.markers, n = 5)
write.csv(clusterHem2.markers, "P5_Unique_Hem2.markers_vs_WT.csv", row.names = T)

clusterHem5.markers <- FindMarkers(pbmc, ident.1 = c("5"), ident.2 = c("1","3","4"), min.pct = 0.25)
head(clusterHem5.markers, n = 5)
write.csv(clusterHem5.markers, "P5_Unique_Hem5.markers_vs_WT.csv", row.names = T)

clusterHem6.markers <- FindMarkers(pbmc, ident.1 = c("6"), ident.2 = c("1","3","4"), min.pct = 0.25)
head(clusterHem6.markers, n = 5)
write.csv(clusterHem6.markers, "P5_Unique_Hem6.markers_vs_WT.csv", row.names = T)

clusterHem7.markers <- FindMarkers(pbmc, ident.1 = c("7"), ident.2 = c("1","3","4"), min.pct = 0.25)
head(clusterHem7.markers, n = 5)
write.csv(clusterHem7.markers, "P5_Unique_Hem7.markers_vs_WT.csv", row.names = T)

#Genes of interest
a <- FeaturePlot(pbmc,features = c("grass","pirk","imd","Jra","PGRP-LF","Ank","key","ben","tub","pll","Rel"),cols = c('lightgray','red'),ncol=3)
ggsave("P6_pbmcUMAP1_Toll&Imd_Enriched_Feature_plot.pdf", plot = plot_grid(a), width = 10, height = 10)

a <- FeaturePlot(pbmc,features = c("Tl", "spz","Toll-4","egr","Toll-6","Toll-7","Toll-9"),cols = c('lightgrey','red'),ncol = 2)
ggsave("P6_pbmcUMAP1_Hemocytes_Toll_Feature_plot.pdf", plot = plot_grid(a), width = 6, height = 7)

#40 genes
a <- FeaturePlot(pbmc,features = c("cic", "kek1", "raw", "alph", "Rac2", "Btk29A", "Jra", "Rac1", "RhoL", 
"edl", "puc", "ttk", "msn", "drk", "hid", "hep", "gro", "vap", "slpr", "sty", "sina", "pnt", "Cka", "Shc", "msk", "aop", 
"rl", "chic", "eff", "Mtl", "Ras85D", "Sod2", "bsk", "RasGAP1", "kay", "lic", "mts", "p38a", "tkv", "Mkk4"),cols = c('lightgrey','red'),ncol = 7)
ggsave("P6_pbmcUMAP1_Hemocytes_MAPK signaling pathway_plot.pdf", plot = plot_grid(a), width = 26, height = 18)

#29 genes
a <- FeaturePlot(pbmc,features = c("cact", "dl", "Myd88", "ben", "imd", "Jra", "Ank", "PGRP-SA", "grass", 
"PGRP-LF", "Tl", "pirk", "key", "hep", "Diap2", "tub", "pll", "Rel", "Dif", "spz", "slmb", "eff", "bsk", "kay", 
"lic", "Uev1A", "p38a", "PGRP-LC", "GNBP3"),cols = c('lightgrey','red'),ncol = 5)
ggsave("P6_pbmcUMAP1_Hemocytes_Toll and Imd signaling pathway_plot.pdf", plot = plot_grid(a), width = 20, height = 18)

#####AddModuleScore
library(RColorBrewer)

gene_sets <- read.table("P7_pbmcUMAP1_Toll&Imd_signaling.txt")
gene_sets2 <- list(gene_sets[,1])
gene_sets2


gene_sets <- c("cact", "dl", "Myd88", "ben", "imd", "Jra", "Ank", "PGRP-SA", "grass", "PGRP-LF", "Tl", "pirk", "key", "hep", 
"Diap2", "tub", "pll", "Rel", "Dif", "spz", "slmb", "eff", "bsk", "kay", "lic", "Uev1A", "p38a", "PGRP-LC", "GNBP3")

sce_T <- AddModuleScore(pbmc, features = gene_sets, nbin = 24, ctrl = 100,name = "Toll_Imd_signaling_score")
#Toll_Imd_signaling_score AddModuleScore之后就是Toll_Imd_signaling_score1
FeaturePlot(sce_T,'Toll_Imd_signaling_score1',cols=rev(brewer.pal(10, name = "RdBu")))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Toll & Imd signaling score")
ggsave("P7_pbmcUMAP1_Toll_Imd_signaling_score.pdf")

gene_sets <- read.table("P7_pbmcUMAP1_MAPK_signaling.txt")
gene_sets2 <- list(gene_sets[,1])
gene_sets2

sce_T <- AddModuleScore(pbmc, features = gene_sets, nbin = 24, ctrl = 100,name = "MAPK_signaling_score")
#Toll_Imd_signaling_score AddModuleScore之后就是Toll_Imd_signaling_score1
FeaturePlot(sce_T,'MAPK_signaling_score1',cols=rev(brewer.pal(10, name = "RdBu")))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("MAPK signaling score")
ggsave("P7_pbmcUMAP1_MAPK_signaling_score.pdf")

###pbmc AddModuleScore 计算相关性 
#使用pbmc@assays$RNA@data调取基因表达scale之后的数据，而不是pbmc@assays$RNA@counts调取原始数据
gene_sets <- read.table("P7_pbmcUMAP1_Toll&Imd_signaling.txt")
gene_sets2 <- list(gene_sets[,1])
gene_sets2
pbmc <- AddModuleScore(pbmc, features = gene_sets, nbin = 24, ctrl = 100,name = "Toll_Imd_signaling_score")

gene_sets <- read.table("P7_pbmcUMAP1_MAPK_signaling.txt")
gene_sets2 <- list(gene_sets[,1])
gene_sets2
pbmc <- AddModuleScore(pbmc, features = gene_sets, nbin = 24, ctrl = 100,name = "MAPK_signaling_score")

head(pbmc@meta.data)

saveRDS(pbmc, file="CCA_WT_RasM6_Bravo_XGY_Hemocytes_Reculstered.Rds") 


Cell.sub <- subset(pbmc@meta.data, sec_group %in% c("RM6"))

Cell.sub <- subset(pbmc@meta.data, sec_group %in% c("WT","RM6"))

a <- as.data.frame(pbmc@meta.data[rownames(Cell.sub),c("sec_group","Toll_Imd_signaling_score1","MAPK_signaling_score1")])

write.csv(a,file="P8_Toll_Imd_signaling_score_cor_MAPK_signaling_score.csv")

Cell.sub <- subset(pbmc@meta.data, seurat_clusters %in% c("0","2","5","6","7"))

a <- as.data.frame(pbmc@meta.data[rownames(Cell.sub),c("seurat_clusters","Toll_Imd_signaling_score1","MAPK_signaling_score1")])

write.csv(a,file="P8_Toll_Imd_signaling_score_cor_MAPK_signaling_score_All_RM6_cluster.csv")

exprSet=as.data.frame(pbmc@assays$RNA@data)
#初始得到的基因表达矩阵 colnames是细胞id rownames是基因
b <- t(exprSet[c("spz"),])   #转置
a$gene = b[rownames(a),]   #增加spz的表达

write.csv(a,file="P8_Toll_Imd_signaling_score_cor_MAPK_signaling_score_All_RM6_cluster.csv")



#提取spz 细胞
table_spz <- as.data.frame(pbmc@assays$RNA@data["spz",])
table(table_spz[,1] != 0)
a_spz <- subset(table_spz, table_spz[,1]!=0)

table_puc <- as.data.frame(pbmc@assays$RNA@data["imd",])
table(table_puc[,1] != 0)
a_puc <- subset(table_puc, table_puc[,1]!=0)

scRNAsub <- subset(pbmc, cells=intersect(row.names(a_spz),row.names(a_puc)))
table(scRNAsub@meta.data$group)

Cell.sub <- subset(scRNAsub@meta.data, sec_group %in% c("WT","RM6"))
a <- as.data.frame(scRNAsub@meta.data[rownames(Cell.sub),c("sec_group","Toll_Imd_signaling_score1","MAPK_signaling_score1")])
a
exprSet=as.data.frame(scRNAsub@assays$RNA@data)
#初始得到的基因表达矩阵 colnames是细胞id rownames是基因
b <- t(exprSet[c("spz","imd"),])   #转置
a$gene = b[rownames(a),]   #增加spz的表达

write.csv(a,file="P8_spz_cor_MAPK_signaling_score.csv")




scRNAsub <- subset(pbmc, sec_group %in% c("RM6"))

a <- FeaturePlot(scRNAsub,features = c("grass","pirk","imd","Jra","PGRP-LF","Ank","key","ben","tub","pll","Rel"),cols = c('lightgray','red'),ncol=3)
ggsave("P7_scRNAsubUMAP1_Toll&Imd_Enriched_Feature_plot.pdf", plot = plot_grid(a), width = 10, height = 10)

a <- FeaturePlot(scRNAsub,features = c("Tl", "spz","Toll-4","egr","Toll-6","Toll-7","Toll-9"),cols = c('lightgrey','red'),ncol = 2)
ggsave("P7_scRNAsubUMAP1_Hemocytes_Toll_Feature_plot.pdf", plot = plot_grid(a), width = 6, height = 7)




################################################################################################WT alone IntegrateData for dotplot

pbmc_WT_Bravo <- readRDS("00WT_Bravo_EAD//WT_Bravo_Hemocytes.Rds")
pbmc_WT_Bravo@meta.data$group = "WT_Bravo"
pbmc_WT_Bravo@meta.data$sec_group = "WT"

pbmc_WT_XGY <- readRDS("CCA_WT_Hemocytes.Rds")

###WT and RM6
DefaultAssay(pbmc_WT_Bravo) <- "RNA"
DefaultAssay(pbmc_WT_XGY) <- "RNA"

WTRM6.anchors <- FindIntegrationAnchors(object.list = list(pbmc_WT_Bravo, pbmc_WT_XGY), anchor.features = 2500, dims = 1:30)
WTRM6.combined <- IntegrateData(anchorset = WTRM6.anchors, dims = 1:30, k.weight = 40)
#IntegrateData在整合对象时，如果有对象细胞数小于100就会报错。k.weight默认为100Q，需要调整k.weight为最小细胞数，假如为40，
pbmc <-  WTRM6.combined

#pbmc <- NormalizeData(pbmc, verbose = FALSE, normalization.method = "LogNormalize", scale.factor = 1e4)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2500)
pbmc <- ScaleData(pbmc, vars.to.regress = c("nCount_RNA"), verbose = TRUE)
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc), npcs = 30, nfeature.print = 10, ndims.print = 1:5, verbose = T)
pc.num=1:30
pbmc <- RunUMAP(pbmc, dims=pc.num)
pbmc <- FindNeighbors(pbmc, dims = pc.num)
pbmc = FindClusters(pbmc,resolution = 0.6)

saveRDS(pbmc, file="CCA_WT_Bravo_WT_XGY_Hemocytes_Reculstered.Rds") 

#pbmc <- readRDS("CCA_WT_Bravo_WT_XGY_Hemocytes_Reculstered.Rds")

##检查default assay
DefaultAssay(pbmc)
##查看有多少个cluster
levels(pbmc)
##设置对应cluster的配色数目
col5 <- colorRampPalette((pal_npg(palette = c("nrc"))(7)))(8)
show_col(col5)

##正常的UMAP降维数据结果
UMPplot_label <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "seurat_clusters", pt.size = 1.8, label.size = 10, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("WT Integrated Hemocytes")
ggsave("P1_UMPplot_label.pdf", plot = plot_grid(UMPplot_label), width = 6, height = 5)
UMPplot_unlabel <- DimPlot(pbmc, reduction = "umap", label = FALSE, group.by = "seurat_clusters", pt.size = 1.8, label.size = 10, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("WT Integrated Hemocytes")
ggsave("P1_UMPplot_unlabel.pdf", plot = plot_grid(UMPplot_unlabel), width = 6, height = 5)
table(pbmc@meta.data[["seurat_clusters"]])

##相同Genotype的分组结果
table(pbmc@meta.data[["group"]])
##按group区别
UMPplot_label_group <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "group", pt.size = 1.8, label.size = 10,  cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("WT Integrated Hemocytes")
ggsave("P2_UMPplot_label_group.pdf", plot = plot_grid(UMPplot_label_group), width = 7, height = 5)
UMPplot_unlabel_group <- DimPlot(pbmc, reduction = "umap", label = FALSE, group.by = "group", pt.size = 1.8, label.size = 10,  cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("WT Integrated Hemocytes")
ggsave("P2_UMPplot_unlabel_group.pdf", plot = plot_grid(UMPplot_unlabel_group), width = 7, height = 5)

##相同Genotype的分组结果
table(pbmc@meta.data[["sec_group"]])
##按sec_group区别
UMPplot_label_sec_group <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "sec_group", pt.size = 1.8, label.size = 10, cols = c("#4DBBD5","#E64B35","#00A087","#3C5488"), raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("WT Integrated Hemocytes")
ggsave("P2_UMPplot_label_sec_group.pdf", plot = plot_grid(UMPplot_label_sec_group), width = 6, height = 5)
UMPplot_unlabel_sec_group <- DimPlot(pbmc, reduction = "umap", label = FALSE, group.by = "sec_group", pt.size = 1.8, label.size = 10, cols = c("#4DBBD5","#E64B35","#00A087","#3C5488"), raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("WT Integrated Hemocytes")
ggsave("P2_UMPplot_unlabel_sec_group.pdf", plot = plot_grid(UMPplot_unlabel_sec_group), width = 6, height = 5)

##这里需要将默认格式输出调整为“RNA”
DefaultAssay(pbmc) <- "RNA"
pbmc <- ScaleData(pbmc, vars.to.regress = c("nCount_RNA"), verbose = TRUE)

## report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc,
                               only.pos = FALSE,
                               min.pct = 0.4,
                               logfc.threshold = 0.25)

Allmarkers_pbmc = pbmc.markers %>% select(gene, everything()) %>% subset(p_val<0.05)
write.csv(Allmarkers_pbmc, "P3_Allmarkers__pbmc_wilcox.csv", row.names = T)

## 使用Top 5基因画热图
top5pbmc.markers <- pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 6, wt = avg_log2FC)
##做Heatmap的热图
Heatmap_markers_pbmc <- DoHeatmap(pbmc,features = top5pbmc.markers$gene,
          group.colors = col5) + theme(text = element_text (size = 20)) + 
		  scale_fill_gradient2(low = '#00468B',mid = 'white',high = '#CC0033',name = 'Z-score')
ggplot2::ggsave(filename = 'P3_Heatmap_markers_pbmc.pdf',width =14,height = 10)

#各个cluster的细胞占比
table(pbmc@meta.data$integrated_snn_res.0.6)
##生成文件画出堆积柱状图
Cells_in_different_Clusters_at_different_group <- table(pbmc@meta.data$integrated_snn_res.0.6)
write.csv(Cells_in_different_Clusters_at_different_group,file = "P3_Cells_in_different_Clusters_at_different_group.csv",row.names = T)

########################################DotPlot

library(Seurat)
library(ggplot2)


#pbmc_RM6 <- readRDS("CCA_WTRasM6_UniqueClusters_reclustered_Hemocytes_reclustered.Rds")
pbmc <- readRDS("CCA_WT_Bravo_WT_XGY_Hemocytes_Reculstered.Rds")

#每种celltype的marker基因/或者需要展示的基因
DefaultAssay(pbmc) <- "RNA"

Pan_hemocytes_genes <- c("srp","He")
Plasmatocytes_genes <- c("NimC1","Cg","Pxn","Hml","eater","crq","Ppn","mys")
PL_0_genes <- c("Eip71CD","fat-spondin","Nlaz","Gal","Sur-8")
PL_Inos_genes <- c("Inos","Chd64","Incenp","stg","CG1572","Gdh")
PL_robo2_genes <- c("TwdlE","E(spl)mbeta-HLH","Drs","CG5399","Snp")
Crystal_cell_genes <- c("lz","peb","PPO1","PPO2","CG9119")
Lamellocytes_genes <- c("atilla","PPO3","alphaTub85E","mthl4")
PSC_Primocytes_genes <- c("Antp","kn","tau","CG15550","CG6023","CG4793","CG6770","CG10038","ImpL2","HmgZ")

PL_Prolif_genes <- c("Ance","CG42394","Tep4","GstD1","Idgf6","et","Eip55E","Lst","CG15347","Nplp2")

#"cathD" excluded

##设置celltype
levels(pbmc)
levels(pbmc@active.ident)

pbmc$celltype <- pbmc@active.ident

new.cluster.ids <- c("0"="Plasmatocytes_1 (PL_0)",
							"1"="Plasmatocytes_2 (PL_Inos)",
							"2"="Plasmatocytes_3 (PL_robo2)",
							"3"="PSC (PL_ImpL2)")

pbmc <- RenameIdents(pbmc, new.cluster.ids)    

levels(pbmc) 


#基因储存为一个list，这个list是DotPlot的输入，才会实现分面的效果
features <- list("Pan_hemocytes" = Pan_hemocytes_genes,
				 "Plasmatocytes" = Plasmatocytes_genes,	 
				 "PL_0_genes" = PL_0_genes,
				 "PL_robo2_genes" = PL_robo2_genes,
				 "PL_Inos_genes" = PL_Inos_genes,
				 "PL_Prolif_genes" = PL_Prolif_genes,
                 "PSC_Primocytes" = PSC_Primocytes_genes,
				 "Lamellocytes" = Lamellocytes_genes,
                 "Crystal_cell" = Crystal_cell_genes
)

#设置celltype展示顺序
levels(pbmc)  <- c("Plasmatocytes_2 (PL_Inos)","Plasmatocytes_1 (PL_0)","Plasmatocytes_3 (PL_robo2)","PSC (PL_ImpL2)")


#作图修饰
Adotplot <- DotPlot(object = pbmc, features=features)&
  theme_bw()& #设置主题
  geom_point(shape=21, aes(size=pct.exp),stroke=1)& #给点添加边框
  theme(axis.title = element_blank(),
        axis.text.x = element_text(color = 'black', size = 10, angle = 90, hjust = 1, vjust = 0.5, face = "bold"),
        axis.text.y = element_text(color = 'black', size = 12, face = "bold"),
        # panel.grid.major = element_blank(), 
        # panel.grid.minor = element_blank(),
        strip.background = element_blank(), #去除分面图背景
        strip.text = element_blank(),#去除分面图文字
        plot.margin=unit(c(1, 1, 1, 1),'cm'),#缩小作图范围
        panel.border = element_rect(color="black",size = 1.2, linetype="solid"),#修改边框大小
        panel.spacing = unit(0.12, "cm"),#修改分面图间距
        # legend.frame = element_rect(colour = "black"),#修改legend边框
        # legend.ticks = element_line(colour = "black", linewidth  = 0),
        legend.key.width = unit(0.3, "cm"),#修改legend宽度
        legend.key.height = unit(0.5, "cm"),#修改legend高度
        legend.title = element_text(color = 'black', face = "bold", size=9))& #修改legend标题字体大小
#  scale_color_gradientn(colours = colorRampPalette(c("navy","white","firebrick3"))(100))& #修改连续legend颜色
scale_color_gradientn(colours = colorRampPalette(c("white","firebrick3"))(100))& #修改连续legend颜色
  labs(tag = "Cell type annotation")& #添加标签
  theme(plot.tag.position = c(0.3, 1.05), #设置标签位置，自行调整
        plot.tag = element_text(size = 12,face = "bold"))& #设置标签文字大小
  guides(size=guide_legend(title="Proportion of\nexpressing cells"), #设置legend size的标题
         colour=guide_colorbar(title="Average\nexpression"))#设置legend color bar的标题

Adotplot

#45 marker genes 14 cm
ggsave("P3_Hemocytes_marker_gene_Adotplot.pdf",plot=Adotplot,width = 14, height = 4.5)


a <- FeaturePlot(pbmc,features = c("Nlaz"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("nplp2")
ggsave("P3_pbmcUMAP1_nplp2plot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("Eip93F"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Idgf6")
ggsave("P3_pbmcUMAP1_Eip93Fplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)
a <- FeaturePlot(pbmc,features = c("spz"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("spz")
ggsave("P3_pbmcUMAP1_spzplot2.pdf", plot = plot_grid(a), width = 4.5, height = 4)


















