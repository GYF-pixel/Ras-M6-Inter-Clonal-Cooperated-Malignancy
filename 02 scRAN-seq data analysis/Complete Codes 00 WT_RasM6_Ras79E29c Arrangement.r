#Step 1 整合不同Genotype的数据
library(devtools)
library(DropletUtils)
library(Seurat)
library(scran)
library(scDblFinder)
library(scater)
library(DoubletFinder)
library(ggplot2)
library(ggsci)
library(cowplot)
library(tidyverse)
library(ggunchull)
library(SCENIC)
library(scales)
set.seed(1234)
getwd()
setwd("I:/M6")

#WT
WT_82Ba <- Read10X(data.dir = "E:/R studio default working directory/WT_82Ba_filtered_feature_bc_matrix")
pbmc_WT_82Ba <- CreateSeuratObject(counts = WT_82Ba,
                            min.cells = 50, 
                            min.features = 200)
WT_82Bb <- Read10X(data.dir = "E:/R studio default working directory/WT_82Bb_filtered_feature_bc_matrix")
pbmc_WT_82Bb <- CreateSeuratObject(counts = WT_82Bb,
                            min.cells = 50, 
                            min.features = 200)
WT_82Bc <- Read10X(data.dir = "E:/R studio default working directory/WT_82Bc_filtered_feature_bc_matrix")
pbmc_WT_82Bc <- CreateSeuratObject(counts = WT_82Bc,
                            min.cells = 50, 
                            min.features = 200)
							
#根据每个pbmc中有多少个细胞进行标记
pbmc_WT_82Ba@meta.data$group <- rep("Wildtype_82Ba",21149)
pbmc_WT_82Bb@meta.data$group <- rep("Wildtype_82Bb",17362)
pbmc_WT_82Bc@meta.data$group <- rep("Wildtype_82Bc",20812)
pbmc = merge(pbmc_WT_82Ba, y=c(pbmc_WT_82Bb, pbmc_WT_82Bc), 
             add.cell.ids = c("Wildtype_82Ba","Wildtype_82Bb", "Wildtype_82Bc"),
             merge.data = TRUE)
head(pbmc@meta.data)
pbmc@meta.data$sec_group = c("WT")
table(pbmc@meta.data$sec_group)
table(pbmc@meta.data$group)

##mitochondrial gene 占比  Rp gene 占比 
scRNAsub = pbmc						
scRNAsub[['percent.mt']] <- PercentageFeatureSet(scRNAsub,pattern = "^mt:")
scRNAsub[["percent.rp"]] <- PercentageFeatureSet(scRNAsub, pattern = "^Rp")			
 
violin <- VlnPlot(scRNAsub,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), 
                  group.by = "group",
                  pt.size = 0, #不需要显示点，可以设置pt.size = 0
                  ncol = 4) + 
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave("vlnplot_before_qc.pdf", plot = violin, width = 12, height = 6)


scRNAsub <- subset(scRNAsub, subset = nCount_RNA > 500 & nCount_RNA < 8000 & nFeature_RNA < 1500 & nFeature_RNA > 250 & percent.mt < 30 & percent.rp < 30)


violin <- VlnPlot(scRNAsub,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), 
                  group.by = "group",
                  pt.size = 0, #不需要显示点，可以设置pt.size = 0
                  ncol = 4) + 
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave("vlnplot_after_qc.pdf", plot = violin, width = 12, height = 6)


##CellCycleScoring
cc.genes <- read.table(file = "cell_cycle_genes.txt")
g2m.genes <- cc.genes[1:68,]
s.genes <- cc.genes[69:124,]
scRNAsub <- CellCycleScoring(scRNAsub, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
scRNAsub$CC.Difference <- scRNAsub$S.Score - scRNAsub$G2M.Score

##Split for doublets analysis
### 拆分为seurat子对象(分开pca1和pca2)
table(scRNAsub@meta.data$group)
table(scRNAsub@meta.data$sec_group)
scRNAsub.list <- SplitObject(scRNAsub, split.by = "group")
scRNAsub.list

for (i in 1:length(scRNAsub.list)) {
	scRNAsub.list[[i]] <- NormalizeData(scRNAsub.list[[i]], verbose = FALSE, normalization.method = "LogNormalize", scale.factor = 1e4)
	scRNAsub.list[[i]] <- FindVariableFeatures(scRNAsub.list[[i]], selection.method = "vst", nfeatures = 2500)
	scRNAsub.list[[i]] <- ScaleData(scRNAsub.list[[i]], vars.to.regress = c("CC.Difference", "percent.mt", "nCount_RNA"), verbose = TRUE)
	scRNAsub.list[[i]] <- RunPCA(scRNAsub.list[[i]], features = VariableFeatures(scRNAsub.list[[i]]), npcs = 40, nfeature.print = 10, ndims.print = 1:5, verbose = T)
	pc.num=1:40
	scRNAsub.list[[i]] <- RunUMAP(scRNAsub.list[[i]], dims=pc.num)
	scRNAsub.list[[i]] <- FindNeighbors(scRNAsub.list[[i]], dims = pc.num)
	scRNAsub.list[[i]] = FindClusters(scRNAsub.list[[i]],resolution = 0.3)
###### 检测doublets 
#Doublets被定义为在相同细胞barcode下测序的两个细胞（例如被捕获在同一液滴中）
### define the expected number of doublet cellscells.
nExp <- round(ncol(scRNAsub.list[[i]]) * 0.05)  ### expect 20% doublets
### remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)
#这个DoubletFinder包的输入是经过预处理（包括归一化、降维，但不一定要聚类）的 Seurat 对象
scRNAsub.list[[i]] <- doubletFinder_v3(scRNAsub.list[[i]], pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:40)	
###### 找Doublets  
### DF的名字是不固定的，因此从scRNAsub@meta.data列名中提取比较保险
DF.name = colnames(scRNAsub.list[[i]]@meta.data)[grepl("DF.classification", colnames(scRNAsub.list[[i]]@meta.data))]
###### 过滤doublet
scRNAsub.list[[i]]=scRNAsub.list[[i]][, scRNAsub.list[[i]]@meta.data[, DF.name] == "Singlet"]
### Filter Mitocondrial 线粒体基因
scRNAsub.list[[i]] <- scRNAsub.list[[i]][!grepl("^mt:", rownames(scRNAsub.list[[i]]),ignore.case = T), ]
### Filter lncRNA
scRNAsub.list[[i]] <- scRNAsub.list[[i]][!grepl("^lncRNA", rownames(scRNAsub.list[[i]]),ignore.case = T), ]
### Filter FBti
scRNAsub.list[[i]] <- scRNAsub.list[[i]][!grepl("^FBti", rownames(scRNAsub.list[[i]]),ignore.case = T), ]
#过滤到此结束
}

scRNAsub.list
reference.list <- scRNAsub.list[c("Wildtype_82Ba", "Wildtype_82Bb", "Wildtype_82Bc")]
scRNAsub.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:60)
scRNAsub.integrated <- IntegrateData(anchorset = scRNAsub.anchors, dims = 1:60)

# switch to integrated assay. The variable features of this assay are automatically set during
# IntegrateData
DefaultAssay(scRNAsub.integrated) <- "integrated"

# 运行标准流程并进行可视化
scRNAsub.integrated <- ScaleData(scRNAsub.integrated, verbose = FALSE)
scRNAsub.integrated <- RunPCA(scRNAsub.integrated, npcs = 60, verbose = FALSE)
scRNAsub.integrated <- RunUMAP(scRNAsub.integrated, reduction = "pca", dims = 1:60)
scRNAsub.integrated <- FindNeighbors(scRNAsub.integrated, dims = 1:60)
scRNAsub.integrated = FindClusters(scRNAsub.integrated,resolution = 0.7)

DF.name = colnames(scRNAsub.integrated@meta.data)[!grepl("pANN_0.25", colnames(scRNAsub.integrated@meta.data))]
DF.name2 = colnames(scRNAsub.integrated@meta.data)[!grepl("DF.classifications", colnames(scRNAsub.integrated@meta.data))]
DF.name3 = colnames(scRNAsub.integrated@meta.data)[!grepl("RNA_snn", colnames(scRNAsub.integrated@meta.data))]
DF.name4 = intersect(intersect(DF.name,DF.name2),DF.name3)

scRNAsub.integrated@meta.data <- scRNAsub.integrated@meta.data[,DF.name4]

saveRDS(scRNAsub.integrated, file="CCA_WT.Rds") 

#R79E29c
R79a <- Read10X(data.dir = "I:/M6/R79A_filtered_feature_bc_matrix")
pbmc_R79a <- CreateSeuratObject(counts = R79a,
                                min.cells = 50, 
                                min.features = 200)
R79b <- Read10X(data.dir = "I:/M6/R79B_filtered_feature_bc_matrix")
pbmc_R79b <- CreateSeuratObject(counts = R79b,
                                min.cells = 50, 
                                min.features = 200)
R79c <- Read10X(data.dir = "I:/M6/R79C_filtered_feature_bc_matrix")
pbmc_R79c <- CreateSeuratObject(counts = R79c,
                                min.cells = 50, 
                                min.features = 200)
											
#根据每个pbmc中有多少个细胞进行标记
pbmc_R79a@meta.data$group <- rep("R79a",10267)
pbmc_R79b@meta.data$group <- rep("R79b",7763)
pbmc_R79c@meta.data$group <- rep("R79c",10728)
pbmc = merge(pbmc_R79a, y=c(pbmc_R79b, pbmc_R79c), 
             add.cell.ids = c("R79a","R79b", "R79c"),
             merge.data = TRUE)
head(pbmc@meta.data)
pbmc@meta.data$sec_group = c("R79")
table(pbmc@meta.data$sec_group)
table(pbmc@meta.data$group)

##mitochondrial gene 占比  Rp gene 占比 
scRNAsub = pbmc						
scRNAsub[['percent.mt']] <- PercentageFeatureSet(scRNAsub,pattern = "^mt:")
scRNAsub[["percent.rp"]] <- PercentageFeatureSet(scRNAsub, pattern = "^Rp")			
 
violin <- VlnPlot(scRNAsub,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), 
                  group.by = "group",
                  pt.size = 0, #不需要显示点，可以设置pt.size = 0
                  ncol = 4) + 
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave("vlnplot_before_qc.pdf", plot = violin, width = 12, height = 6)


scRNAsub <- subset(scRNAsub, subset = nCount_RNA > 200 & nCount_RNA < 30000 & nFeature_RNA < 4500 & nFeature_RNA > 250 & percent.mt < 35 & percent.rp < 40)


violin <- VlnPlot(scRNAsub,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), 
                  group.by = "group",
                  pt.size = 0, #不需要显示点，可以设置pt.size = 0
                  ncol = 4) + 
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave("vlnplot_after_qc.pdf", plot = violin, width = 12, height = 6)

##CellCycleScoring
cc.genes <- read.table(file = "cell_cycle_genes.txt")
g2m.genes <- cc.genes[1:68,]
s.genes <- cc.genes[69:124,]
scRNAsub <- CellCycleScoring(scRNAsub, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
scRNAsub$CC.Difference <- scRNAsub$S.Score - scRNAsub$G2M.Score

##Split for doublets analysis
### 拆分为seurat子对象(分开pca1和pca2)
table(scRNAsub@meta.data$group)
table(scRNAsub@meta.data$sec_group)
scRNAsub.list <- SplitObject(scRNAsub, split.by = "group")
scRNAsub.list

for (i in 1:length(scRNAsub.list)) {
	scRNAsub.list[[i]] <- NormalizeData(scRNAsub.list[[i]], verbose = FALSE, normalization.method = "LogNormalize", scale.factor = 1e4)
	scRNAsub.list[[i]] <- FindVariableFeatures(scRNAsub.list[[i]], selection.method = "vst", nfeatures = 2500)
	scRNAsub.list[[i]] <- ScaleData(scRNAsub.list[[i]], vars.to.regress = c("CC.Difference", "percent.mt", "nCount_RNA"), verbose = TRUE)
	scRNAsub.list[[i]] <- RunPCA(scRNAsub.list[[i]], features = VariableFeatures(scRNAsub.list[[i]]), npcs = 40, nfeature.print = 10, ndims.print = 1:5, verbose = T)
	pc.num=1:40
	scRNAsub.list[[i]] <- RunUMAP(scRNAsub.list[[i]], dims=pc.num)
	scRNAsub.list[[i]] <- FindNeighbors(scRNAsub.list[[i]], dims = pc.num)
	scRNAsub.list[[i]] = FindClusters(scRNAsub.list[[i]],resolution = 0.3)
###### 检测doublets 
#Doublets被定义为在相同细胞barcode下测序的两个细胞（例如被捕获在同一液滴中）
### define the expected number of doublet cellscells.
nExp <- round(ncol(scRNAsub.list[[i]]) * 0.05)  ### expect 20% doublets
### remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)
#这个DoubletFinder包的输入是经过预处理（包括归一化、降维，但不一定要聚类）的 Seurat 对象
scRNAsub.list[[i]] <- doubletFinder_v3(scRNAsub.list[[i]], pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:40)	
###### 找Doublets  
### DF的名字是不固定的，因此从scRNAsub@meta.data列名中提取比较保险
DF.name = colnames(scRNAsub.list[[i]]@meta.data)[grepl("DF.classification", colnames(scRNAsub.list[[i]]@meta.data))]
###### 过滤doublet
scRNAsub.list[[i]]=scRNAsub.list[[i]][, scRNAsub.list[[i]]@meta.data[, DF.name] == "Singlet"]
### Filter Mitocondrial 线粒体基因
scRNAsub.list[[i]] <- scRNAsub.list[[i]][!grepl("^mt:", rownames(scRNAsub.list[[i]]),ignore.case = T), ]
### Filter lncRNA
scRNAsub.list[[i]] <- scRNAsub.list[[i]][!grepl("^lncRNA", rownames(scRNAsub.list[[i]]),ignore.case = T), ]
### Filter FBti
scRNAsub.list[[i]] <- scRNAsub.list[[i]][!grepl("^FBti", rownames(scRNAsub.list[[i]]),ignore.case = T), ]
#过滤到此结束
}

scRNAsub.list
reference.list <- scRNAsub.list[c("R79a","R79b", "R79c")]
scRNAsub.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:60)
scRNAsub.integrated <- IntegrateData(anchorset = scRNAsub.anchors, dims = 1:60)

# switch to integrated assay. The variable features of this assay are automatically set during
# IntegrateData
DefaultAssay(scRNAsub.integrated) <- "integrated"

# 运行标准流程并进行可视化
scRNAsub.integrated <- ScaleData(scRNAsub.integrated, verbose = FALSE)
scRNAsub.integrated <- RunPCA(scRNAsub.integrated, npcs = 60, verbose = FALSE)
scRNAsub.integrated <- RunUMAP(scRNAsub.integrated, reduction = "pca", dims = 1:60)
scRNAsub.integrated <- FindNeighbors(scRNAsub.integrated, dims = 1:60)
scRNAsub.integrated = FindClusters(scRNAsub.integrated,resolution = 0.7)

DF.name = colnames(scRNAsub.integrated@meta.data)[!grepl("pANN_0.25", colnames(scRNAsub.integrated@meta.data))]
DF.name2 = colnames(scRNAsub.integrated@meta.data)[!grepl("DF.classifications", colnames(scRNAsub.integrated@meta.data))]
DF.name3 = colnames(scRNAsub.integrated@meta.data)[!grepl("RNA_snn", colnames(scRNAsub.integrated@meta.data))]
DF.name4 = intersect(intersect(DF.name,DF.name2),DF.name3)

scRNAsub.integrated@meta.data <- scRNAsub.integrated@meta.data[,DF.name4]


saveRDS(scRNAsub.integrated, file="CCA_Ras79E29c.Rds") 
							
#RasM6
RM6a <- Read10X(data.dir = "I:/M6/RM6A_filtered_feature_bc_matrix")
pbmc_RM6a <- CreateSeuratObject(counts = RM6a,
                                min.cells = 50, 
                                min.features = 200)
RM6b <- Read10X(data.dir = "I:/M6/RM6B_filtered_feature_bc_matrix")
pbmc_RM6b <- CreateSeuratObject(counts = RM6b,
                                min.cells = 50, 
                                min.features = 200)
RM6c <- Read10X(data.dir = "I:/M6/RM6C_filtered_feature_bc_matrix")
pbmc_RM6c <- CreateSeuratObject(counts = RM6c,
                                min.cells = 50, 
                                min.features = 200)
											
#根据每个pbmc中有多少个细胞进行标记
pbmc_RM6a@meta.data$group <- rep("RM6a",20133)
pbmc_RM6b@meta.data$group <- rep("RM6b",17976)
pbmc_RM6c@meta.data$group <- rep("RM6c",16193)
pbmc = merge(pbmc_RM6a, y=c(pbmc_RM6b, pbmc_RM6c), 
             add.cell.ids = c("RM6a","RM6b", "RM6c"),
             merge.data = TRUE)
head(pbmc@meta.data)
pbmc@meta.data$sec_group = c("RM6")
table(pbmc@meta.data$sec_group)
table(pbmc@meta.data$group)

##mitochondrial gene 占比  Rp gene 占比 
scRNAsub = pbmc						
scRNAsub[['percent.mt']] <- PercentageFeatureSet(scRNAsub,pattern = "^mt:")
scRNAsub[["percent.rp"]] <- PercentageFeatureSet(scRNAsub, pattern = "^Rp")			
 
violin <- VlnPlot(scRNAsub,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), 
                  group.by = "group",
                  pt.size = 0, #不需要显示点，可以设置pt.size = 0
                  ncol = 4) + 
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave("vlnplot_before_qc.pdf", plot = violin, width = 12, height = 6)


scRNAsub <- subset(scRNAsub, subset = nCount_RNA > 200 & nCount_RNA < 12000 & nFeature_RNA < 2500 & nFeature_RNA > 250 & percent.mt < 35 & percent.rp < 40)


violin <- VlnPlot(scRNAsub,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), 
                  group.by = "group",
                  pt.size = 0, #不需要显示点，可以设置pt.size = 0
                  ncol = 4) + 
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave("vlnplot_after_qc.pdf", plot = violin, width = 12, height = 6)

##CellCycleScoring
cc.genes <- read.table(file = "cell_cycle_genes.txt")
g2m.genes <- cc.genes[1:68,]
s.genes <- cc.genes[69:124,]
scRNAsub <- CellCycleScoring(scRNAsub, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
scRNAsub$CC.Difference <- scRNAsub$S.Score - scRNAsub$G2M.Score

##Split for doublets analysis
### 拆分为seurat子对象(分开pca1和pca2)
table(scRNAsub@meta.data$group)
table(scRNAsub@meta.data$sec_group)
scRNAsub.list <- SplitObject(scRNAsub, split.by = "group")
scRNAsub.list

for (i in 1:length(scRNAsub.list)) {
	scRNAsub.list[[i]] <- NormalizeData(scRNAsub.list[[i]], verbose = FALSE, normalization.method = "LogNormalize", scale.factor = 1e4)
	scRNAsub.list[[i]] <- FindVariableFeatures(scRNAsub.list[[i]], selection.method = "vst", nfeatures = 2500)
	scRNAsub.list[[i]] <- ScaleData(scRNAsub.list[[i]], vars.to.regress = c("CC.Difference", "percent.mt", "nCount_RNA"), verbose = TRUE)
	scRNAsub.list[[i]] <- RunPCA(scRNAsub.list[[i]], features = VariableFeatures(scRNAsub.list[[i]]), npcs = 40, nfeature.print = 10, ndims.print = 1:5, verbose = T)
	pc.num=1:40
	scRNAsub.list[[i]] <- RunUMAP(scRNAsub.list[[i]], dims=pc.num)
	scRNAsub.list[[i]] <- FindNeighbors(scRNAsub.list[[i]], dims = pc.num)
	scRNAsub.list[[i]] = FindClusters(scRNAsub.list[[i]],resolution = 0.3)
###### 检测doublets 
#Doublets被定义为在相同细胞barcode下测序的两个细胞（例如被捕获在同一液滴中）
### define the expected number of doublet cellscells.
nExp <- round(ncol(scRNAsub.list[[i]]) * 0.05)  ### expect 20% doublets
### remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)
#这个DoubletFinder包的输入是经过预处理（包括归一化、降维，但不一定要聚类）的 Seurat 对象
scRNAsub.list[[i]] <- doubletFinder_v3(scRNAsub.list[[i]], pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:40)	
###### 找Doublets  
### DF的名字是不固定的，因此从scRNAsub@meta.data列名中提取比较保险
DF.name = colnames(scRNAsub.list[[i]]@meta.data)[grepl("DF.classification", colnames(scRNAsub.list[[i]]@meta.data))]
###### 过滤doublet
scRNAsub.list[[i]]=scRNAsub.list[[i]][, scRNAsub.list[[i]]@meta.data[, DF.name] == "Singlet"]
### Filter Mitocondrial 线粒体基因
scRNAsub.list[[i]] <- scRNAsub.list[[i]][!grepl("^mt:", rownames(scRNAsub.list[[i]]),ignore.case = T), ]
### Filter lncRNA
scRNAsub.list[[i]] <- scRNAsub.list[[i]][!grepl("^lncRNA", rownames(scRNAsub.list[[i]]),ignore.case = T), ]
### Filter FBti
scRNAsub.list[[i]] <- scRNAsub.list[[i]][!grepl("^FBti", rownames(scRNAsub.list[[i]]),ignore.case = T), ]
#过滤到此结束
}

scRNAsub.list
reference.list <- scRNAsub.list[c("RM6a","RM6b", "RM6c")]
scRNAsub.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:60)
scRNAsub.integrated <- IntegrateData(anchorset = scRNAsub.anchors, dims = 1:60)

# switch to integrated assay. The variable features of this assay are automatically set during
# IntegrateData
DefaultAssay(scRNAsub.integrated) <- "integrated"

# 运行标准流程并进行可视化
scRNAsub.integrated <- ScaleData(scRNAsub.integrated, verbose = FALSE)
scRNAsub.integrated <- RunPCA(scRNAsub.integrated, npcs = 60, verbose = FALSE)
scRNAsub.integrated <- RunUMAP(scRNAsub.integrated, reduction = "pca", dims = 1:60)
scRNAsub.integrated <- FindNeighbors(scRNAsub.integrated, dims = 1:60)
scRNAsub.integrated = FindClusters(scRNAsub.integrated,resolution = 0.7)

DF.name = colnames(scRNAsub.integrated@meta.data)[!grepl("pANN_0.25", colnames(scRNAsub.integrated@meta.data))]
DF.name2 = colnames(scRNAsub.integrated@meta.data)[!grepl("DF.classifications", colnames(scRNAsub.integrated@meta.data))]
DF.name3 = colnames(scRNAsub.integrated@meta.data)[!grepl("RNA_snn", colnames(scRNAsub.integrated@meta.data))]
DF.name4 = intersect(intersect(DF.name,DF.name2),DF.name3)

scRNAsub.integrated@meta.data <- scRNAsub.integrated@meta.data[,DF.name4]

saveRDS(scRNAsub.integrated, file="CCA_RasM6.Rds") 
