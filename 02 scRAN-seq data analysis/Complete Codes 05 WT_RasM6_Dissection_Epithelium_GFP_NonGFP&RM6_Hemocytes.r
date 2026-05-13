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
library(viridis)
library(RColorBrewer)
library(gridExtra)
#library(scCustomize) # 需要Seurat版本5.0

set.seed(1234)

getwd()
setwd("I:/M6")
setwd("H:/M6")

#读入GFP NonGFP 和 Hemocytes
pbmc_GFP <- readRDS("CCA_WTRM6_UniqueClusters_reclustered_Epithelium_GFP.Rds")

pbmc_NonGFP <- readRDS("CCA_WTRM6_UniqueClusters_reclustered_Epithelium_NonGFP.Rds")

pbmc_Hemocytes <- readRDS("CCA_WTRasM6_UniqueClusters_reclustered_Hemocytes_reclustered.Rds")

table(pbmc_GFP@meta.data[["sec_group"]])
table(pbmc_NonGFP@meta.data[["sec_group"]])
table(pbmc_Hemocytes@meta.data[["sec_group"]])

pbmc_GFP@meta.data$Type = "GFP"
pbmc_NonGFP@meta.data$Type = "NonGFP"
pbmc_Hemocytes@meta.data$Type = "Hemocytes"

DefaultAssay(pbmc_GFP) <- "RNA"
DefaultAssay(pbmc_NonGFP) <- "RNA"
DefaultAssay(pbmc_Hemocytes) <- "RNA"

###整合pbmc_GFP pbmc_Hemocytes
WTGFP.anchors <- FindIntegrationAnchors(object.list = list(pbmc_GFP, pbmc_Hemocytes), anchor.features = 2500, dims = 1:50)
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

saveRDS(pbmc, file="CCA_WTRM6_UniqueClusters_reclustered_Epithelium_GFP_and_Hemocytes.Rds") 

##出图

#pbmc <- readRDS("CCA_WTRM6_UniqueClusters_reclustered_Epithelium_GFP_and_Hemocytes.Rds")

##检查default assay
DefaultAssay(pbmc)

##查看有多少个cluster
levels(pbmc)

##设置对应cluster的配色数目
col5 <- colorRampPalette((pal_npg(palette = c("nrc"))(7)))(17)
show_col(col5)

##正常的UMAP降维数据结果
UMPplot_label <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "seurat_clusters", pt.size = 1.0, label.size = 8, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("Integrated scRNA-seq datasets")
ggsave("P1_UMPplot_label.pdf", plot = plot_grid(UMPplot_label), width = 9, height = 7)
UMPplot_unlabel <- DimPlot(pbmc, reduction = "umap", label = FALSE, group.by = "seurat_clusters", pt.size = 1.0, label.size = 8, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("Integrated scRNA-seq datasets")
ggsave("P1_UMPplot_unlabel.pdf", plot = plot_grid(UMPplot_unlabel), width = 9, height = 7)

##相同Genotype的分组结果
table(pbmc@meta.data[["Type"]])
##按group区别
UMPplot_label_group <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "Type", pt.size = 1.0, label.size = 7, cols = c("#4DBBD5","#E64B35","#00A087","#3C5488"), raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("Integrated scRNA-seq datasets")
ggsave("P2_UMPplot_label_group.pdf", plot = plot_grid(UMPplot_label_group), width = 9, height = 7)
UMPplot_unlabel_group <- DimPlot(pbmc, reduction = "umap", label = FALSE, group.by = "Type", pt.size = 1.0, label.size = 7, cols = c("#4DBBD5","#E64B35","#00A087","#3C5488"), raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("Integrated scRNA-seq datasets")
ggsave("P2_UMPplot_unlabel_group.pdf", plot = plot_grid(UMPplot_unlabel_group), width = 9, height = 7)

##各cell cluster的细胞占比
table(pbmc@meta.data$Type, pbmc@meta.data$integrated_snn_res.1.2)

##生成文件画出堆积柱状图
Cells_in_different_Clusters_at_different_group <- table(pbmc@meta.data$Type, pbmc@meta.data$integrated_snn_res.1.2)
write.csv(Cells_in_different_Clusters_at_different_group,file = "P3_Cells_in_different_Clusters_at_different_group.csv",row.names = T)
##计算单个cluster中两种cell type的占比，设置Unique Cluster筛选条件

##FlyphoneDB
#####expMatrix
expMatrix <- pbmc@assays$RNA@counts
write.csv(expMatrix, file="expMatrix2.csv")

#####meta
pbmc@meta.data$Type2 = "GFP"
meta <- pbmc@meta.data
meta$integrated_snn_res.1.2 == "11"
meta_Hemocytes = meta[meta$integrated_snn_res.1.2 == "11",]
meta[rownames(meta_Hemocytes),"Type2"] = "Hemocytes"
table(meta$Type2)
pbmc@meta.data = meta
meta <- pbmc@meta.data[,c("integrated_snn_res.1.2","Type2")]
#到Excel里去合并Meta中的data，变成GFP3 GFP4等等
write.csv(meta, file="meta2.csv")


#####Linux里去跑

#1.导入R环境
module load R/3.6.3          #导入R3.6.3
module load R/4.0.5          #导入R4.0.5
module load R/4.1.2          #导入R4.1.2
module load R/4.2.1          #导入R4.2.1 #通常用
module load R/4.3.1          #导入R4.3.1

#2.安装软件包
#用户导入R 环境后，默认会将包装在~/R/x86_64-pc-linux-gnu-library/ 目录下对应的R版本目录下，也可以设置R_LIBS_USER变量重新定义安装目录。

echo "export R_LIBS_USER=$HOME/R_LIBS" >> ~/.bash_profile
source ~/.bash_profile
module load R/4.2.1
R
install.packages("optparse")
install.packages("tidyverse")
install.packages("future.apply")
install.packages('Seurat', repos = c('https://satijalab.r-universe.dev')) # 检查版本
# packageVersion("Seurat")  # [1] ‘4.4.0’
install.packages("RColorBrewer")
install.packages("reshape2")
install.packages("network")
install.packages("igraph")

list.of.packages <- c(‘fitdistrplus’, ‘ggridges’, ‘httr’, ‘ica’, ‘igraph’, 
‘irlba’, ‘leiden’, ‘lmtest’, ‘matrixStats’, ‘miniUI’, ‘patchwork’, ‘pbapply’, 
‘plotly’, ‘png’, ‘RANN’, ‘RcppAnnoy’, ‘reticulate’, ‘ROCR’, ‘Rtsne’, 
‘scattermore’, ‘sctransform’, ‘shiny’, ‘spatstat.explore’, ‘spatstat.geom’, 
‘uwot’, ‘RcppProgress’)

 library(optparse)
 library(tidyverse)
 library(future.apply)
 library(Seurat)
 library(RColorBrewer)
 library(reshape2)
 library(network)
 library(igraph)
#运行以上命令后，在个人配置文件中,R_LIBS_USER变量被定义为用户家目录下的R_LIBS文件夹。
#导入R/3.6.3环境后，运行R命令打开R的交互式界面，这时用install.packages可以把包安装到用户自己目录下。

#3. 更换下载源
#用户可以在R命令行中输入以下命令，更换安装包过程中的下载源，以下是西湖大学本地R源示例：
options(repos=c(CRAN="https://mirrors.westlake.edu.cn/CRAN"))

#4. 修改shell脚本
#!/bin/bash

Rscript --vanilla FlyPhone_parallel_batch.R \
		--matrix ./test_dataset_60M/2021-03-25_matrix.csv \
		--metadata ./test_dataset_60M/2021-03-25_metadata.csv \
		--lrpair ./annotation/Ligand_receptor_pair_high_confident_2021vs1_clean.txt \
		--corecomponents ./annotation/Pathway_core_components_2021vs1_clean.txt \
		--cores 8 \
		--output ./output_test_dataset_60M/


Rscript --vanilla /storage/maxianjueLab/guoyifan/Software/FlyPhoneDB/FlyPhoneDB-master/FlyPhone_parallel_batch.R \
		--matrix  /storage/maxianjueLab/guoyifan/Software/FlyPhoneDB/NonGFP_Hemocytes/expMatrix.csv \
		--metadata /storage/maxianjueLab/guoyifan/Software/FlyPhoneDB/NonGFP_Hemocytes/meta.csv \
		--lrpair /storage/maxianjueLab/guoyifan/Software/FlyPhoneDB/FlyPhoneDB-master/annotation/Ligand_receptor_pair_high_confident_2021vs1_clean.txt \
		--corecomponents /storage/maxianjueLab/guoyifan/Software/FlyPhoneDB/FlyPhoneDB-master/annotation/Pathway_core_components_2021vs1_clean.txt \
		--cores 8 \
		--output /storage/maxianjueLab/guoyifan/Software/FlyPhoneDB/GFP_Hemocytes/output_dataset/


#5. 修改bash脚本

#!/bin/bash
#SBATCH -J FlyDB_NGFP
#SBATCH -p amd-ep2,intel-sc3,amd-ep2-short
#SBATCH -q normal
#SBATCH --mem=100G
#SBATCH -c 8
sh /storage/maxianjueLab/guoyifan/Software/FlyPhoneDB/GFP_Hemocytes/FlyPhone_sbatch_NonGFP_Hemocytes.sh

#6. 提交Bash文件
module load R/4.2.1

sbatch FlyPhone_sbatch_NonGFP_Hemocytes.bash


###整合pbmc_NonGFP pbmc_Hemocytes
WTGFP.anchors <- FindIntegrationAnchors(object.list = list(pbmc_NonGFP, pbmc_Hemocytes), anchor.features = 2500, dims = 1:50)
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

saveRDS(pbmc, file="CCA_WTRM6_UniqueClusters_reclustered_Epithelium_NonGFP_and_Hemocytes.Rds") 

##出图

#pbmc <- readRDS("CCA_WTRM6_UniqueClusters_reclustered_Epithelium_NonGFP_and_Hemocytes.Rds")

##检查default assay
DefaultAssay(pbmc)

##查看有多少个cluster
levels(pbmc)

##设置对应cluster的配色数目
col5 <- colorRampPalette((pal_npg(palette = c("nrc"))(7)))(16)
show_col(col5)

##正常的UMAP降维数据结果
UMPplot_label <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "seurat_clusters", pt.size = 1.0, label.size = 8, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("Integrated scRNA-seq datasets")
ggsave("P1_UMPplot_label.pdf", plot = plot_grid(UMPplot_label), width = 9, height = 7)
UMPplot_unlabel <- DimPlot(pbmc, reduction = "umap", label = FALSE, group.by = "seurat_clusters", pt.size = 1.0, label.size = 8, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("Integrated scRNA-seq datasets")
ggsave("P1_UMPplot_unlabel.pdf", plot = plot_grid(UMPplot_unlabel), width = 9, height = 7)


##相同Genotype的分组结果
table(pbmc@meta.data[["Type"]])
##按group区别
UMPplot_label_group <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "Type", pt.size = 1.0, label.size = 7, cols = c("#4DBBD5","#E64B35","#00A087","#3C5488"), raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("Integrated scRNA-seq datasets")
ggsave("P2_UMPplot_label_group.pdf", plot = plot_grid(UMPplot_label_group), width = 9, height = 7)
UMPplot_unlabel_group <- DimPlot(pbmc, reduction = "umap", label = FALSE, group.by = "Type", pt.size = 1.0, label.size = 7, cols = c("#4DBBD5","#E64B35","#00A087","#3C5488"), raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("Integrated scRNA-seq datasets")
ggsave("P2_UMPplot_unlabel_group.pdf", plot = plot_grid(UMPplot_unlabel_group), width = 9, height = 7)


##各cell cluster的细胞占比
table(pbmc@meta.data$Type, pbmc@meta.data$integrated_snn_res.1.2)

##生成文件画出堆积柱状图
Cells_in_different_Clusters_at_different_group <- table(pbmc@meta.data$Type, pbmc@meta.data$integrated_snn_res.1.2)
write.csv(Cells_in_different_Clusters_at_different_group,file = "P3_Cells_in_different_Clusters_at_different_group.csv",row.names = T)
##计算单个cluster中两种cell type的占比，设置Unique Cluster筛选条件

##FlyphoneDB
#####expMatrix
expMatrix <- pbmc@assays$RNA@counts
write.csv(expMatrix, file="expMatrix.csv")

#####meta
pbmc@meta.data$Type2 = "NonGFP"
meta <- pbmc@meta.data
meta$integrated_snn_res.1.2 == "14"
meta_Hemocytes = meta[meta$integrated_snn_res.1.2 == "14",]
meta[rownames(meta_Hemocytes),"Type2"] = "Hemocytes"
table(meta$Type2)
pbmc@meta.data = meta
meta <- pbmc@meta.data[,c("integrated_snn_res.1.2","Type2")]
#到Excel里去合并Meta中的data，变成NonGFP3 NonGFP4等等
write.csv(meta, file="meta.csv")













