#Complete Codes 02 WT_Isolating Hemocytes_Glia_Epithelium

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

DefaultAssay(pbmc_WT)

pbmc <- pbmc_WT

a <- FeaturePlot(pbmc,features = c("GFP"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("GFP")
ggsave("P4_pbmcUMAP1_GFPplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(pbmc,features = c("Ras85D"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Ras85D")
ggsave("P4_pbmcUMAP1_Ras85Dplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(pbmc,features = c("M6"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("M6")
ggsave("P4_pbmcUMAP1_M6plot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(pbmc,features = c("GAL80"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("GAL80")
ggsave("P4_pbmcUMAP1_GAL80plot2.pdf", plot = plot_grid(a), width = 6, height = 5)


a <- FeaturePlot(pbmc,features = c("pnr", "Hml", "pigs","Ppn","kuz","He","srp","NimC1","atilla","lz","Mmp1","IM18","CecA2","CecC","Mtk","DptB","Drs","Prx2540-1","Prx2540-2","CG12896","Abl","Snoo","Ubx","CG15550","CG6023","mthl7","Cys","CG8860","COX8"),cols = c('lightgrey','red'))
ggsave("P10_pbmcUMAP1_Hemocytes_Feature_plot.pdf", plot = plot_grid(a), width = 18, height = 16)
a <- FeaturePlot(pbmc,features = c("repo", "CG3168","Gli", "nrv2", "sty","moody", "NK7.1", "CG9336"),cols = c('lightgrey','red'))
ggsave("P10_pbmcUMAP1_Nerve_Feature_plot.pdf", plot = plot_grid(a), width = 12, height = 10)


Cells.sub <- subset(pbmc@meta.data, integrated_snn_res.0.7 %in% c("0","1","2","3","4","5","6","7","8","9","11","12","13","14","15","17"))
scRNAsub_Epithelium <- subset(pbmc, cells=rownames(Cells.sub))

saveRDS(scRNAsub_Epithelium, file="CCA_WT_Epithelium.Rds") 



