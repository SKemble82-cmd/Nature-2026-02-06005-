###############lines 0-325 samples processing and integration. lines 326 onwards analysis for figures

library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(reticulate)
library(RColorBrewer)
library(sctransform)
library(clustree)
library(harmony)
library(BiocManager)
options(bitmapType='cairo')
library(remotes)
library(gsfisher)
library(EnhancedVolcano)
library(RColorBrewer)
library(ggsci)
library(pals)
library(wesanderson)
library(CellChat)
library(ggfortify)
library(DESeq2)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(data.table)
library(gsfisher)
library(cowplot)
library(igraph)
library(Nebulosa)
library(scProportionTest)
library(tidyr)
library(miloR)
library(scater)
library(patchwork)
library(SingleCellExperiment)
library(monocle)
library(DDRTree)
library(SeuratWrappers)
library(devtools)
library(ggthemes)
library(tvthemes)
library(pheatmap)

##########resting

CIA_r1p <- Read10X(data.dir = "/rds/projects/c/croftap-ktr-cia-sc-01/CIA_CD45_pos/Con2_CD45pos/outs/filtered_feature_bc_matrix")
CIA_r1p <- CreateSeuratObject(counts = CIA_r1p, min.cells=3, min.features=100, project="CIA_r1p")
CIA_r1p <- RenameCells(CIA_r1p, add.cell.id = "CIA_r1p")
CIA_r1p$sample_id <- "pResting1"
CIA_r1p$CD45 <- "pos"
CIA_r1p$State <- "Resting"
head(CIA_r1p)

CIA_r2p <- Read10X(data.dir = "/rds/projects/c/croftap-ktr-cia-sc-01/CIA_CD45_pos/Con3_CD45pos/outs/filtered_feature_bc_matrix")
CIA_r2p <- CreateSeuratObject(counts = CIA_r2p, min.cells=3, min.features=100, project="CIA_r2p")
CIA_r2p <- RenameCells(CIA_r2p, add.cell.id = "CIA_r2p")
CIA_r2p$sample_id <- "pResting2"
CIA_r2p$CD45 <- "pos"
CIA_r2p$State <- "Resting"
head(CIA_r2p)

CIA_r3p <- Read10X(data.dir = "/rds/projects/c/croftap-ktr-cia-sc-01/CIA_CD45_pos/Con1_CD45pos/outs/filtered_feature_bc_matrix")
CIA_r3p <- CreateSeuratObject(counts = CIA_r3p, min.cells=3, min.features=100, project="CIA_r3p")
CIA_r3p <- RenameCells(CIA_r3p, add.cell.id = "CIA_r3p")
CIA_r3p$sample_id <- "pResting3"
CIA_r3p$CD45 <- "pos"
CIA_r3p$State <- "Resting"
head(CIA_r3p)

#################non-treated
CIA_i1p <- Read10X(data.dir = "/rds/projects/c/croftap-sgcar/X204SC24074958-Z01-F001/count/CON1POS/outs/per_sample_outs/CON1POS/count/sample_filtered_feature_bc_matrix")
CIA_i1p <- CreateSeuratObject(counts = CIA_i1p, min.cells=3, min.features=100, project="CIA_i1p")
CIA_i1p <- RenameCells(CIA_i1p, add.cell.id = "CIA_i1p")
CIA_i1p$sample_id <- "pInflamed1"
CIA_i1p$CD45 <- "pos"
CIA_i1p$State <- "Inflamed"
head(CIA_i1p)

CIA_i2p <- Read10X(data.dir = "/rds/projects/c/croftap-sgcar/X204SC24074958-Z01-F001/count/CON2POS/outs/per_sample_outs/CON2POS/count/sample_filtered_feature_bc_matrix")
CIA_i2p <- CreateSeuratObject(counts = CIA_i2p, min.cells=3, min.features=100, project="CIA_i2p")
CIA_i2p <- RenameCells(CIA_i2p, add.cell.id = "CIA_i2p")
CIA_i2p$sample_id <- "pInflamed2"
CIA_i2p$CD45 <- "pos"
CIA_i2p$State <- "Inflamed"
head(CIA_i2p)

CIA_i3p <- Read10X(data.dir = "/rds/projects/c/croftap-sgcar/X204SC24074958-Z01-F001/count/CON3POS/outs/per_sample_outs/CON3POS/count/sample_filtered_feature_bc_matrix")
CIA_i3p <- CreateSeuratObject(counts = CIA_i3p, min.cells=3, min.features=100, project="CIA_i3p")
CIA_i3p <- RenameCells(CIA_i3p, add.cell.id = "CIA_i3p")
CIA_i3p$sample_id <- "pInflamed3"
CIA_i3p$CD45 <- "pos"
CIA_i3p$State <- "Inflamed"
head(CIA_i3p)

############CAR-T treated
CIA_c1p <- Read10X(data.dir = "/rds/projects/c/croftap-sgcar/X204SC24074958-Z01-F001/count/CAR1POS/outs/per_sample_outs/CAR1POS/count/sample_filtered_feature_bc_matrix")
CIA_c1p <- CreateSeuratObject(counts = CIA_c1p, min.cells=3, min.features=100, project="CIA_c1p")
CIA_c1p <- RenameCells(CIA_c1p, add.cell.id = "CIA_c1p")
CIA_c1p$sample_id <- "pCAR1"
CIA_c1p$CD45 <- "pos"
CIA_c1p$State <- "CAR"
head(CIA_c1p)

CIA_c2p <- Read10X(data.dir = "/rds/projects/c/croftap-sgcar/X204SC24074958-Z01-F001/count/CAR2POS/outs/per_sample_outs/CAR2POS/count/sample_filtered_feature_bc_matrix")
CIA_c2p <- CreateSeuratObject(counts = CIA_c2p, min.cells=3, min.features=100, project="CIA_c2p")
CIA_c2p <- RenameCells(CIA_c2p, add.cell.id = "CIA_c2p")
CIA_c2p$sample_id <- "pCAR2"
CIA_c2p$CD45 <- "pos"
CIA_c2p$State <- "CAR"
head(CIA_c2p)

CIA_c3p <- Read10X(data.dir = "/rds/projects/c/croftap-sgcar/X204SC24074958-Z01-F001/count/CAR3POS/outs/per_sample_outs/CAR3POS/count/sample_filtered_feature_bc_matrix")
CIA_c3p <- CreateSeuratObject(counts = CIA_c3p, min.cells=3, min.features=100, project="CIA_c3p")
CIA_c3p <- RenameCells(CIA_c3p, add.cell.id = "CIA_c3p")
CIA_c3p$sample_id <- "pCAR3"
CIA_c3p$CD45 <- "pos"
CIA_c3p$State <- "CAR"
head(CIA_c3p)


mito.featuresCIA_r1p <- grep(pattern="^mt-", x=rownames(x=CIA_r1p), value=T)
percent.mitoCIA_r1p <- Matrix::colSums(x = GetAssayData(object = CIA_r1p, slot = "counts")[mito.featuresCIA_r1p,]) / Matrix::colSums(x = GetAssayData(object = CIA_r1p, slot = "counts"))
CIA_r1p[["percent.mito"]] <- percent.mitoCIA_r1p
VlnPlot(object = CIA_r1p, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0)

CIA_r1p <- subset(x = CIA_r1p, subset = nFeature_RNA > 500 & nFeature_RNA <3000 & nCount_RNA > 200 & nCount_RNA < 10000 & percent.mitoCIA_r1p < 0.1)
VlnPlot(object = CIA_r1p, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
CIA_r1p <- NormalizeData(object = CIA_r1p, verbose = F)
CIA_r1p <- FindVariableFeatures(object = CIA_r1p, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesCIA_r1p <- rownames(CIA_r1p)
CIA_r1p <- ScaleData(CIA_r1p, features = all.genesCIA_r1p)
VlnPlot(object = CIA_r1p, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)


mito.featuresCIA_r2p <- grep(pattern="^mt-", x=rownames(x=CIA_r2p), value=T)
percent.mitoCIA_r2p <- Matrix::colSums(x = GetAssayData(object = CIA_r2p, slot = "counts")[mito.featuresCIA_r2p,]) / Matrix::colSums(x = GetAssayData(object = CIA_r2p, slot = "counts"))
CIA_r2p[["percent.mito"]] <- percent.mitoCIA_r2p
VlnPlot(object = CIA_r2p, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0)

CIA_r2p <- subset(x = CIA_r2p, subset = nFeature_RNA > 500 & nFeature_RNA <2000 & nCount_RNA > 500 & nCount_RNA < 6000 & percent.mitoCIA_r2p < 0.1)
VlnPlot(object = CIA_r2p, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
CIA_r2p <- NormalizeData(object = CIA_r2p, verbose = F)
CIA_r2p <- FindVariableFeatures(object = CIA_r2p, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesCIA_r2p <- rownames(CIA_r2p)
CIA_r2p <- ScaleData(CIA_r2p, features = all.genesCIA_r2p)
VlnPlot(object = CIA_r2p, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)


mito.featuresCIA_r3p <- grep(pattern="^mt-", x=rownames(x=CIA_r3p), value=T)
percent.mitoCIA_r3p <- Matrix::colSums(x = GetAssayData(object = CIA_r3p, slot = "counts")[mito.featuresCIA_r3p,]) / Matrix::colSums(x = GetAssayData(object = CIA_r3p, slot = "counts"))
CIA_r3p[["percent.mito"]] <- percent.mitoCIA_r3p
VlnPlot(object = CIA_r3p, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

CIA_r3p <- subset(x = CIA_r3p, subset = nFeature_RNA > 500 & nFeature_RNA <6000 & nCount_RNA > 500 & nCount_RNA < 30000 & percent.mitoCIA_r3p < 0.1)
VlnPlot(object = CIA_r3p, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
CIA_r3p <- NormalizeData(object = CIA_r3p, verbose = F)
CIA_r3p <- FindVariableFeatures(object = CIA_r3p, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesCIA_r3p <- rownames(CIA_r3p)
CIA_r3p <- ScaleData(CIA_r3p, features = all.genesCIA_r3p)
VlnPlot(object = CIA_r3p, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)


mito.featuresCIA_i1p <- grep(pattern="^mt-", x=rownames(x=CIA_i1p), value=T)
percent.mitoCIA_i1p <- Matrix::colSums(x = GetAssayData(object = CIA_i1p, slot = "counts")[mito.featuresCIA_i1p,]) / Matrix::colSums(x = GetAssayData(object = CIA_i1p, slot = "counts"))
CIA_i1p[["percent.mito"]] <- percent.mitoCIA_i1p
VlnPlot(object = CIA_i1p, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0)

CIA_i1p <- subset(x = CIA_i1p, subset = nFeature_RNA > 1000 & nFeature_RNA <7500 & nCount_RNA > 500 & nCount_RNA < 50000 & percent.mitoCIA_i1p < 0.1)
VlnPlot(object = CIA_i1p, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
CIA_i1p <- NormalizeData(object = CIA_i1p, verbose = F)
CIA_i1p <- FindVariableFeatures(object = CIA_i1p, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesCIA_i1p <- rownames(CIA_i1p)
CIA_i1p <- ScaleData(CIA_i1p, features = all.genesCIA_i1p)
VlnPlot(object = CIA_i1p, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)


mito.featuresCIA_i2p <- grep(pattern="^mt-", x=rownames(x=CIA_i2p), value=T)
percent.mitoCIA_i2p <- Matrix::colSums(x = GetAssayData(object = CIA_i2p, slot = "counts")[mito.featuresCIA_i2p,]) / Matrix::colSums(x = GetAssayData(object = CIA_i2p, slot = "counts"))
CIA_i2p[["percent.mito"]] <- percent.mitoCIA_i2p
VlnPlot(object = CIA_i2p, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0)

CIA_i2p <- subset(x = CIA_i2p, subset = nFeature_RNA > 1000 & nFeature_RNA <6000 & nCount_RNA > 500 & nCount_RNA < 40000 & percent.mitoCIA_i2p < 0.1)
VlnPlot(object = CIA_i2p, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
CIA_i2p <- NormalizeData(object = CIA_i2p, verbose = F)
CIA_i2p <- FindVariableFeatures(object = CIA_i2p, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesCIA_i2p <- rownames(CIA_i2p)
CIA_i2p <- ScaleData(CIA_i2p, features = all.genesCIA_i2p)
VlnPlot(object = CIA_i2p, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)


mito.featuresCIA_i3p <- grep(pattern="^mt-", x=rownames(x=CIA_i3p), value=T)
percent.mitoCIA_i3p <- Matrix::colSums(x = GetAssayData(object = CIA_i3p, slot = "counts")[mito.featuresCIA_i3p,]) / Matrix::colSums(x = GetAssayData(object = CIA_i3p, slot = "counts"))
CIA_i3p[["percent.mito"]] <- percent.mitoCIA_i3p
VlnPlot(object = CIA_i3p, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0)

CIA_i3p <- subset(x = CIA_i3p, subset = nFeature_RNA > 1000 & nFeature_RNA <8000 & nCount_RNA > 500 & nCount_RNA < 60000 & percent.mitoCIA_i3p < 0.1)
VlnPlot(object = CIA_i3p, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
CIA_i3p <- NormalizeData(object = CIA_i3p, verbose = F)
CIA_i3p <- FindVariableFeatures(object = CIA_i3p, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesCIA_i3p <- rownames(CIA_i3p)
CIA_i3p <- ScaleData(CIA_i3p, features = all.genesCIA_i3p)
VlnPlot(object = CIA_i3p, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)


mito.featuresCIA_c1p <- grep(pattern="^mt-", x=rownames(x=CIA_c1p), value=T)
percent.mitoCIA_c1p <- Matrix::colSums(x = GetAssayData(object = CIA_c1p, slot = "counts")[mito.featuresCIA_c1p,]) / Matrix::colSums(x = GetAssayData(object = CIA_c1p, slot = "counts"))
CIA_c1p[["percent.mito"]] <- percent.mitoCIA_c1p
VlnPlot(object = CIA_c1p, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0)

CIA_c1p <- subset(x = CIA_c1p, subset = nFeature_RNA > 1000 & nFeature_RNA <7500 & nCount_RNA > 500 & nCount_RNA < 50000 & percent.mitoCIA_c1p < 0.1)
VlnPlot(object = CIA_c1p, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
CIA_c1p <- NormalizeData(object = CIA_c1p, verbose = F)
CIA_c1p <- FindVariableFeatures(object = CIA_c1p, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesCIA_c1p <- rownames(CIA_c1p)
CIA_c1p <- ScaleData(CIA_c1p, features = all.genesCIA_c1p)
VlnPlot(object = CIA_c1p, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)


mito.featuresCIA_c2p <- grep(pattern="^mt-", x=rownames(x=CIA_c2p), value=T)
percent.mitoCIA_c2p <- Matrix::colSums(x = GetAssayData(object = CIA_c2p, slot = "counts")[mito.featuresCIA_c2p,]) / Matrix::colSums(x = GetAssayData(object = CIA_c2p, slot = "counts"))
CIA_c2p[["percent.mito"]] <- percent.mitoCIA_c2p
VlnPlot(object = CIA_c2p, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0)

CIA_c2p <- subset(x = CIA_c2p, subset = nFeature_RNA > 1000 & nFeature_RNA <7500 & nCount_RNA > 500 & nCount_RNA < 50000 & percent.mitoCIA_c2p < 0.1)
VlnPlot(object = CIA_c2p, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
CIA_c2p <- NormalizeData(object = CIA_c2p, verbose = F)
CIA_c2p <- FindVariableFeatures(object = CIA_c2p, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesCIA_c2p <- rownames(CIA_c2p)
CIA_c2p <- ScaleData(CIA_c2p, features = all.genesCIA_c2p)
VlnPlot(object = CIA_c2p, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)


mito.featuresCIA_c3p <- grep(pattern="^mt-", x=rownames(x=CIA_c3p), value=T)
percent.mitoCIA_c3p <- Matrix::colSums(x = GetAssayData(object = CIA_c3p, slot = "counts")[mito.featuresCIA_c3p,]) / Matrix::colSums(x = GetAssayData(object = CIA_c3p, slot = "counts"))
CIA_c3p[["percent.mito"]] <- percent.mitoCIA_c3p
VlnPlot(object = CIA_c3p, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0)

CIA_c3p <- subset(x = CIA_c3p, subset = nFeature_RNA > 2000 & nFeature_RNA <7500 & nCount_RNA > 500 & nCount_RNA < 40000 & percent.mitoCIA_c3p < 0.1)
VlnPlot(object = CIA_c3p, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
CIA_c3p <- NormalizeData(object = CIA_c3p, verbose = F)
CIA_c3p <- FindVariableFeatures(object = CIA_c3p, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesCIA_c3p <- rownames(CIA_c3p)
CIA_c3p <- ScaleData(CIA_c3p, features = all.genesCIA_c3p)
VlnPlot(object = CIA_c3p, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)


## Identification of integration anchors
reference.list <- c(CIA_r1p, CIA_r2p, CIA_r3p, CIA_i1p, CIA_i2p, CIA_i3p, CIA_c1p, CIA_c2p, CIA_c3p)

anchors <- FindIntegrationAnchors(object.list = reference.list, anchor.features = 2000, dims = 1:30)
integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
saveRDS(integrated, file = "cd45pos.rds")


DefaultAssay(object=integrated) <- "integrated"
integrated <- ScaleData(object = integrated, verbose=F)
integrated <- RunPCA(object = integrated, verbose=F, npcs = 50)
ElbowPlot(object = integrated, ndims = 50)
integrated <- FindNeighbors(object = integrated, dims = 1:30)

integrated <- RunUMAP(object = integrated, reduction = "pca", dims = 1:30)
DimPlot(integrated, reduction = "umap", label = F, pt.size = 0.1, repel = F, split.by = "State") + NoAxes()+NoLegend()

DefaultAssay(object=integrated) <- "RNA"
FeaturePlot(integrated, features = "Col1a1")

table(integrated$sample_id)
Idents(integrated)<-"sample_id"

###down sample to no more than 5000 cells per sample#####
integrated1 <- subset(integrated, downsample = 5000)
table(integrated1$sample_id)

DefaultAssay(object=integrated1) <- "integrated"
integrated1 <- FindVariableFeatures(integrated1, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(integrated1)
integrated1 <- ScaleData(object = integrated1, verbose=F)
integrated1 <- RunPCA(object = integrated1, verbose=F, npcs = 50)
ElbowPlot(object = integrated1, ndims = 50)
integrated1 <- FindNeighbors(object = integrated1, dims = 1:30)
integrated1 <- FindClusters(integrated1, graph.name = "integrated_snn", resolution = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7))
clustree(integrated1, assay = "integrated")
integrated1 <- FindClusters(integrated1, graph.name = "integrated_snn", resolution = 0.2)
integrated1 <- RunUMAP(integrated1, reduction = "pca", dims = 1:30)
DimPlot(integrated1, reduction = "umap", label = F, pt.size = 0.1, repel = T, label.box = T, split.by = "State") + NoAxes()

####check gene lists for annotation####
c0 <- FindMarkers(integrated1, ident.1 = 0, min.pct = 0.2, only.pos = T)
c1 <- FindMarkers(integrated1, ident.1 = 1, min.pct = 0.2, only.pos = T)
c2 <- FindMarkers(integrated1, ident.1 = 2, min.pct = 0.2, only.pos = T)
c3 <- FindMarkers(integrated1, ident.1 = 3, min.pct = 0.2, only.pos = T)
c4 <- FindMarkers(integrated1, ident.1 = 4, min.pct = 0.2, only.pos = T)
c5 <- FindMarkers(integrated1, ident.1 = 5, min.pct = 0.2, only.pos = T)
c6 <- FindMarkers(integrated1, ident.1 = 6, min.pct = 0.2, only.pos = T)
c7 <- FindMarkers(integrated1, ident.1 = 7, min.pct = 0.2, only.pos = T)
c8 <- FindMarkers(integrated1, ident.1 = 8, min.pct = 0.2, only.pos = T)
c9 <- FindMarkers(integrated1, ident.1 = 9, min.pct = 0.2, only.pos = T)
c10 <- FindMarkers(integrated1, ident.1 = 10, min.pct = 0.2, only.pos = T)
c11 <- FindMarkers(integrated1, ident.1 = 11, min.pct = 0.2, only.pos = T)
c12 <- FindMarkers(integrated1, ident.1 = 12, min.pct = 0.2, only.pos = T)
c13 <- FindMarkers(integrated1, ident.1 = 13, min.pct = 0.2, only.pos = T)
c14 <- FindMarkers(integrated1, ident.1 = 14, min.pct = 0.2, only.pos = T)


#####annotate#####
current.sample.ids <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14")
new.sample.ids <- c("Macrophages", "Mono_Macrophages", "Mono_Macrophages", "Mono_Macrophages", "Macrophages", "Cycling", "Cycling",
                    "T_Lymphocytes", "B_Lymphocytes", "Macrophages", "Con", "Mast_Cells", "Con", "Con", "Cycling")


integrated1$clusters <- integrated1@active.ident
integrated1@meta.data[["clusters"]] <- plyr::mapvalues(x = integrated1@meta.data[["clusters"]], from = current.sample.ids, to = new.sample.ids)
Idents(integrated1)<-"clusters"


####remove contaminants####
Intergated2<-subset(integrated1, idents=c("Macrophages", "Mono_Macrophages", "T_Lymphocytes", "B_Lymphocytes", "Mast_Cells"))
DimPlot(Intergated2, reduction = "umap", label = F, pt.size = 0.1, repel = F, label.box = T)+ NoAxes()#+scale_color_avatar(palette = "FireNation")

getwd()
saveRDS(Intergated2, file = "cd45pos.rds")


####add CD45neg data####
cd45neg<-readRDS(file = "cd45neg.rds")
cd45pos<-readRDS(file = "cd45pos.rds")

all <- merge(cd45neg, cd45pos)
##############all stromal cell data analysis#########

all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(all)
all <- ScaleData(all, features = all.genes)
all.H <- RunHarmony(all, "orig.ident", assay.use = "RNA", plot_convergence = TRUE, max.iter.harmony = 20)
all.H <- FindNeighbors(all.H, reduction = "pca", graph.name = "integrated_snn", dims = 1:50)
DefaultAssay(object=all.H) <- "integrated"
all.H <- FindClusters(all.H, graph.name = "integrated_snn", resolution = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7))
clustree(all.H, assay = "integrated")
all.H <- FindClusters(all.H, resolution = 0.5)
all.H <- RunUMAP(all.H, reduction = "harmony", dims = 1:50)
DimPlot(all.H, reduction = "HarmonyUMAP", label = T) +NoAxes()

FindAllMarkers_all <- FindAllMarkers(all.H, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.5)

#####annotate#####
current.sample.ids <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21")
new.sample.ids <- c("Fibroblasts", "Fibroblasts","Mono_Macrophages", "Macrophages", "Fibroblasts", "Macrophages", "Mono_Macrophages", "Fibroblasts", "Osteoblast", "Mono_Macrophages", "Fibroblasts", "Mural", "Macrophages", "Vascular", "Cycling",
                    "Fibroblasts", "T_Lymphocytes", "Vadcular", "B_Lymphocytes",  "Mono_Macrophages", "Vascular", "Mast_Cells")


all.H$clusters <- all.H@active.ident
all.H@meta.data[["clusters"]] <- plyr::mapvalues(x = all.H@meta.data[["clusters"]], from = current.sample.ids, to = new.sample.ids)
Idents(all.H)<-"clusters"


####remove cycling cells####
all.H.res<-subset(integrated1, idents=c("Fibroblasts", "Osteoblast", "Mural", "Vascular", "Macrophages", "Mono_Macrophages", "T_Lymphocytes", "B_Lymphocytes", "Mast_Cells"))

########extended data figure 9a
DimPlot(all.H.res, reduction = "HarmonyUMAP", label = F, pt.size = 0.1, repel = F, label.box = T)+ NoAxes()#+scale_color_avatar(palette = "FireNation")



FindAllMarkers_all <- FindAllMarkers(all.H.res, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.5)
write.csv(FindAllMarkers_all, file = "Supplementary Table 15.csv")

allcell.genes<-c("Col1a1", "Pdpn", "Cd14", "Il1b", "Cd68", "Apoe", "Runx2", "Bglap", "Pecam1", "Cd36", "Notch3", "Rgs4", "Cd3g", "Icos", "Cd19", "Cd79a", "Mcpt14", "Cma1")

########extended data figure 9b
#plot dotplot
dotplot <- DotPlot(all.H.res, features = allcell.genes, dot.scale = 6)+ theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) +ylab(NULL) +xlab(NULL)+ scale_color_gradientn(colours = magma(20))+ FontSize(x.text = 14, y.text = 14) 
dotplot


########extended data figure 9c
test <- sc_utils(all.H.res)
prop.test <- permutation_test(test, cluster_identity = "clusters", sample_1="CAR", sample_2="Inflamed", sample_identity="State", n_permutations=10000)
permutation_plot(prop.test, FDR_threshold = 0.01, log2FD_threshold = 0.58, order_clusters = T)+ FontSize(x.text = 12, y.text = 12, x.title = 12, y.title = 12)+xlab(NULL)+NoLegend()

####mac1 subset



#pathway analysis###########################
Idents(CARcd45)<-"clusters"
CARcd45<-subset(Intergated2, idents="CAR")
Inflamedcd45<-subset(Intergated2, idents="Inflamed")
Restingcd45<-subset(Intergated2, idents="Resting")


#########pathway analysis for clusters##########change Ident and add this to seurat_obj below
Idents(mac1)<-"clusters"

########extended data figure 9d
#########pathway analysis for tissue state##########change Ident and add this to seurat_obj below
Idents(mac1)<-"State"

#######Add your object to "seurat_obj <-"
library(gsfisher)
seurat_obj <- mac1
#Edit to suitable output folder
getwd() 
GSout <- "/rds/projects/c/croftap-croftapcarcia/GOmac-state"

###saved as annotation.rds do not run again if akready done so
fetchAnnotation(species = "mm", ensembl_version = NULL,
                ensembl_host = NULL)


## run this everytime
if(!file.exists("annotation.rds"))
{
  annotation <- fetchAnnotation(species="mm",
                                ensembl_version=NULL,
                                ensembl_host = NULL)
  saveRDS(annotation, "annotation.rds")
}
annotation <- readRDS("annotation.rds")

##declare function - change name for data set
getExpressedGenesFromSeuratObjectRNA <- function(seurat_obj,
                                                 clusters,
                                                 min.pct=0.1)
{
  expressed <- c()
  for(cluster in clusters)
  {
    # get genes detected in the cluster
    cluster_cells <- names(seurat_obj@active.ident[seurat_obj@active.ident==cluster])
    clust_pcts <- apply(seurat_obj@assays$RNA@data[,cluster_cells],
                        1, function(x) sum(x>0)/length(x))
    
    detected_in_clust <- names(clust_pcts[clust_pcts>min.pct])
    
    # get genes detected in the other cells
    other_cells <- names(seurat_obj@active.ident[seurat_obj@active.ident!=cluster])
    other_pcts <- apply(seurat_obj@assays$RNA@data[,other_cells],
                        1, function(x) sum(x>0)/length(x))
    
    detected_in_other_cells <- names(other_pcts[other_pcts>min.pct])
    
    expressed <- c(expressed, detected_in_clust, detected_in_other_cells)
  }
  expressed <- unique(expressed)
  expressed
}


#Setup for gsfisher run####
expressed_genes <- getExpressedGenesFromSeuratObjectRNA(seurat_obj,levels(seurat_obj@active.ident), min.pct=0.1)


seurat_obj.res <- FindAllMarkers(seurat_obj, only.pos = T, min.pct = 0.2, logfc.threshold = 0.5) 
#This should be your features, probably easiest to run FindAllMarkers then use the table out of that here


seurat_obj.res$entrez_id <- as.character(annotation$entrez_id[
  match(seurat_obj.res$gene, annotation$gene_name)])
seurat_obj.res <- seurat_obj.res[!is.na(seurat_obj.res$entrez_id),]
background_entrez <- as.character(annotation$entrez_id[
  match(expressed_genes, annotation$gene_name)])
background_entrez <- background_entrez[!is.na(background_entrez)]


#Create output dir
if (dir.exists(GSout)){
  print("Directory already exists")
} else{
  dir.create(GSout)
}

#Run enrichment tests and filter results
for(cluster in levels(seurat_obj@active.ident)){
  cat("Running go enrichment for cluster", cluster, "...")
  Go <- runGO(foreground_ids = seurat_obj.res$entrez_id[seurat_obj.res$cluster == cluster],
              background_ids = background_entrez,
              gene_id_type = "entrez", species = "mm")
  
  cat("Running Kegg enrichment for cluster", cluster, "...")
  Kegg <- runKEGG(foreground_ids = seurat_obj.res$entrez_id[seurat_obj.res$cluster == cluster],
                  background_ids = background_entrez,
                  gene_id_type = "entrez", species = "mm")
  
  cat("Joining results tables for cluster", cluster, "...")
  results <- full_join(Go, Kegg)
  
  results$p.adj <- p.adjust(results$p.val, method = "BH")
  
  cat("Filtering genesets for ", cluster, "...")
  filtered_genesets <- filterGenesets(results,
                                      min_foreground_genes = 5,
                                      max_genes_geneset = 500,
                                      min_odds_ratio = 4,
                                      padjust_method="BH",
                                      use_adjusted_pvalues=TRUE,
                                      pvalue_threshold=0.05)
  
  filtered_genesets$cluster <- cluster
  
  cat("Writing out results table for ", cluster, "...")
  write.csv(filtered_genesets, file.path(GSout, paste("Geneset_Enrichment", 
                                                      cluster,"csv", sep=".")))
}
library(plyr)
#Read in enrichments and bind into single table
All_enrichment_path = list.files(path=GSout, pattern="Geneset_Enrichment", full.names=TRUE)
All_enrichments = ldply(All_enrichment_path, read.csv)
All_enrichments <- All_enrichments[,-1]
write.csv(All_enrichments, file.path(GSout, paste("All.enrichments", "csv", sep=".")))


##plots
all_results_top <- All_enrichments %>% group_by(cluster) %>% top_n(n=5, -p.val)

sampleEnrichmentDotplot(all_results_top, selection_col = "description", selected_genesets = unique(All_enrichments$description), fill_colors = c("yellow3", "brown", "black"), sample_id_col = "cluster", fill_var = "odds.ratio", maxl=50, rotate_sample_labels = TRUE, min_dot_size = 5, pvalue_threshold = 0.05) + FontSize(x.text = 12, y.text = 12) + ylab(NULL) +xlab(NULL)
sampleEnrichmentHeatmap(all_results_top)


###########macrophage analysis##########

#### macrophage annotation
macs<-subset(Intergated2, idents=c("Macrophages"))
DefaultAssay(object=macs) <- "integrated"
macs <- FindVariableFeatures(macs, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(macs)
macs <- ScaleData(object = macs, verbose=F)
macs <- RunPCA(object = macs, verbose=F, npcs = 50)
ElbowPlot(object = macs, ndims = 50)
macs <- FindNeighbors(object = macs, dims = 1:20)
macs <- FindClusters(macs, graph.name = "integrated_snn", resolution = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7))
clustree(macs, assay = "integrated")
macs <- FindClusters(macs, graph.name = "integrated_snn", resolution = 0.2)
macs <- RunUMAP(macs, reduction = "pca", dims = 1:20)
DimPlot(macs, reduction = "umap", label = F, pt.size = 0.1, repel = T, label.box = T, split.by = "State") + NoAxes()
DefaultAssay(object=macs) <- "RNA"

m0 <- FindMarkers(macs, ident.1 = "0", min.pct = 0.2, only.pos = T)
m1 <- FindMarkers(macs, ident.1 = "1", min.pct = 0.2, only.pos = T)
m2 <- FindMarkers(macs, ident.1 = "2", min.pct = 0.2, only.pos = T)
m3 <- FindMarkers(macs, ident.1 = "3", min.pct = 0.2, only.pos = T)
m4 <- FindMarkers(macs, ident.1 = "3", min.pct = 0.2, only.pos = T)
m5 <- FindMarkers(macs, ident.1 = "5", min.pct = 0.2, only.pos = T)
m6 <- FindMarkers(macs, ident.1 = "6", min.pct = 0.2, only.pos = T)
m7 <- FindMarkers(macs, ident.1 = "7", min.pct = 0.2, only.pos = T)

current.sample.ids <- c("0", "1", "2", "3", "4", "5", "6", "7")
new.sample.ids <- c("Cxcr2", "Aqp1", "H2DMb1", "Tnf", "Vsig4", "Gas6", "Plac8", "Itgax_DC")

macs$mac_clusters <- macs@active.ident
macs@meta.data[["mac_clusters"]] <- plyr::mapvalues(x = macs@meta.data[["mac_clusters"]], from = current.sample.ids, to = new.sample.ids)
Idents(macs)<-"mac_clusters"

#####remove Itgax_DC - dendritic cells#####
mac1 <- subset(macs, idents=c("Cxcr2", "Aqp1", "H2DMb1","Tnf", "Vsig4", "Gas6", "Plac8"))


#########figure 5a
Idents(mac1)<-'mac_clusters'
DimPlot(macs, reduction = "umap", label = F, pt.size = 0.1, repel = T, label.box = T) + NoAxes()+scale_color_simpsons(n, type = "discrete")
DimPlot(mac1, reduction = "umap", label = F, pt.size = 0.1, repel = T, label.box = T, split.by = "State") + NoAxes()+scale_color_simpsons(n, type = "discrete")


######extended data figure 9f
#####cd14 expression
VlnPlot(mac1,features = "Cd14", cols = c("#FED439FF", "#8A9197FF", "#D2AF81FF",   "#F05C3BFF", "#D5E4A2FF", "#197EC0FF", "brown"), pt.size = 0)+NoLegend()

#############figure 5b
########macrophage proportion bar chart#######
Idents(mac1)<-"mac_clusters"

pt <- table(mac1$mac_clusters, mac1$State)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)
pt$Var1 <- factor(pt$Var1,levels=c("Cx3cr1", "Gas6", "H2DMb1", "Aqp1", "Tnf", "Cxcr2", "Plac8"))

scales::show_col(tvthemes:::simpsons_palette)

ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab(NULL) +
  ylab("Proportion") +
  theme(legend.title = element_blank())+RotatedAxis() +
  scale_fill_manual(values = c("#ABC67D", "#0363C3", "#d1b271", "#424f46", "#FC0209", "#fed90f", "chocolate4"))+NoLegend()

###########figure 5c
####scproprtion test
test <- sc_utils(mac1)
prop.test <- permutation_test(test, cluster_identity = "mac_clusters", sample_1="Resting", sample_2="Inflamed", sample_identity="State", n_permutations=10000)
permutation_plot(prop.test, FDR_threshold = 0.01, log2FD_threshold = 0.58, order_clusters = T)+ FontSize(x.text = 12, y.text = 12, x.title = 12, y.title = 12)+xlab(NULL)+NoLegend()

prop.test <- permutation_test(test, cluster_identity = "clusters", sample_1="CAR", sample_2="Inflamed", sample_identity="State", n_permutations=10000)
permutation_plot(prop.test, FDR_threshold = 0.01, log2FD_threshold = 0.58, order_clusters = T)+ FontSize(x.text = 12, y.text = 12, x.title = 12, y.title = 12)+xlab(NULL)+NoLegend()


#########Gene lists###########
####vsig4 lining gene list
cx3cr1.genes<-c("Vsig4", "Sparc", "Ctsd", "Lyz2", "Srgn", "S100b", "Fn1", "Apoe", "Wfdc17", "Hoxb", "Nupr1", "Grn", "Ltc4s", "Trem2", "Cd68", "Folr2", "Cx3cr1")
#####relma macrophage gene list#####
Relm.genes<-c("Ccl8", "Ccl2", "Ccl7", "Cxcl2", "Dusp1", "Zfp36", "Junb", "Pf4", "Fos", "Atf3", "Mt1", "Ier3", "Ccl14", "Retnla", "Ccl24", "Marcksl1", "Cxcl13", "Jun", "Mt2")
#####MHCII macrophage gene list 
MHC.genes<-c("H2-Ab1", "H2-Eb1", "H2-Aa", "Cd74", "H2-DMa", "Cd52", "Mg12", "S100a11", "Gm2a", "Lsp1", "Coro1a", "Pim1", "Retnla", "Clec4a2")
#####aquaporin macrophage gene list
Aqua.genes<-c("Aqp1", "Fxyd2", "Tpp3", "Cd9", "Lyve1", "Anxa1", "Cyb5r3", "Ifi27l2a")
####immunoregulatory macrophage gene list
immreg.genes<-c("Rps4y1", "Timd4", "Cfd", "Cd99", "Mtrnr2l8", "Rbp4", "Hpgds", "Folr2", "C1qc", "Aldh1a1", "Tgfb1", "Calm2", "Ppdpf", "Lilrb5", "Fxyd5", "S100a4", "Ctsa", "Serpinf1", "Cpvl", "Grn", "Tspan15", "Lipa", "Acp5", "S100a10", "Lgmn", "Tspo", "Vsig4", "Pebp1", "Dstn", "Apoe")


#######extended data figure 9e
mac1 <- AddModuleScore(
  object = mac1,
  features = list(cx3cr1.genes),
  ctrl = 5,
  name = 'cx3cr1')

plot_density(mac1, features = "cx3cr11", combine = F, method = "wkde", size = 1, adjust = 1)+NoAxes()+ scale_color_gradientn(colours = magma(20))+facet_grid(.~mac1$State)


mac1 <- AddModuleScore(
  object = mac1,
  features = list(Relm.genes),
  ctrl = 5,
  name = 'Relm')

plot_density(mac1, features = "Relm1", combine = F, method = "wkde", size = 1, adjust = 1)+NoAxes()+ scale_color_gradientn(colours = magma(20))+facet_grid(.~mac1$State)

mac1 <- AddModuleScore(
  object = mac1,
  features = list(MHC.genes),
  ctrl = 5,
  name = 'MHC')

plot_density(mac1, features = "MHC1", combine = F, method = "wkde", size = 1, adjust = 1)+NoAxes()+ scale_color_gradientn(colours = magma(20))+facet_grid(.~mac1$State)


mac1 <- AddModuleScore(
  object = mac1,
  features = list(Aqua.genes),
  ctrl = 5,
  name = 'Aqua')

plot_density(mac1, features = "Aqua1", combine = F, method = "wkde", size = 1, adjust = 1)+NoAxes()+ scale_color_gradientn(colours = magma(20))+facet_grid(.~mac1$State)

#############figure 5e
mac1 <- AddModuleScore(
  object = mac1,
  features = list(immreg.genes),
  ctrl = 5,
  name = 'immreg')

plot_density(mac1, features = "immreg1", combine = F, method = "wkde", size = 1, adjust = 1)+NoAxes()+ scale_color_gradientn(colours = magma(20))+facet_grid(.~mac1$State)

Idents(mac1)<-"State"
VlnPlot(mac1, features = "immreg1")


############boxplot figure 5e
expr <- FetchData(mac1, vars = "immreg1")
expr$State <- mac1$State   # add metadata column
colnames(expr)[1] <- "expression"

library(ggplot2)
library(ggpubr)

expr <- expr[!is.na(expr$State), ]
expr$State <- factor(expr$State,
                     levels = c("Resting", "Inflamed", "CAR"))

my_comparisons <- list(
  c("Inflamed", "Resting"),
  c("Inflamed", "CAR")
)

ggplot(expr, aes(x = State, y = expression, fill = State)) +
  geom_boxplot(outlier.size = 0.5) 
  theme_classic() +
  ylab("Immreg1 expression") +
  xlab("Treatment") +
  theme(text = element_text(size = 14),
        axis.text = element_text(size = 14)) +
  scale_fill_manual(values = c("Resting" = "goldenrod",
                               "Inflamed" = "grey",
                               "CAR" = "firebrick")) +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "kruskal.test",
    label = "p.signif"
  ) +
  NoLegend()

######gene score across tissue state#########  
immreg.df <- FetchData(mac1, vars = c("immreg1", "State"))
colnames(immreg.df)[1] <- "immreg"


immreg.avg_scores <- immreg.df %>%
  group_by(State) %>%
  summarise(ModuleScore = mean(immreg))

imm.mat <- as.matrix(immreg.avg_scores$ModuleScore)
rownames(imm.mat) <- immreg.avg_scores$State
colnames(imm.mat) <- "ImmReg"

imm.mat <- as.matrix(imm.mat)
mode(imm.mat) <- "numeric"

pheatmap(imm.mat, scale = "row", cluster_rows = FALSE, cluster_cols = FALSE,
         main = "Immunoregulatory Module")


immreg.df <- data.frame(
  State = c("Resting", "Inflamed", "CAR"),
  Score = c(0.35668676, 0.03375275, 0.15599322)
)
# convert to matrix
imm.mat <- as.matrix(immreg.df$Score)
rownames(imm.mat) <- immreg.df$State
colnames(imm.mat) <- "Immunoregulatory"

imm.mat_h <- t(imm.mat)

# plot
pheatmap(imm.mat_h,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         main = NULL,
         fontsize = 14,
     
         
             color = colorRampPalette(c("#999999", "#E69F00", "red"))(100))


########figure 5e
##########macrophage pseudotime#############
Idents(mac1)<-"State"
Restingmac<-subset(mac1, idents="Resting")
inflmac<-subset(mac1, idents="Inflamed")
CARmac<-subset(mac1, idents="CAR")
Idents(Restingmac)<-"mac_clusters"
Idents(inflmac)<-"mac_clusters"
Idents(CARmac)<-"mac_clusters"

#######get gene lists
FindAllMarkers_rest <- FindAllMarkers(Restingmac, only.pos = F, min.pct = 0.2, logfc.threshold = 0.5)
FindAllMarkers_inflam <- FindAllMarkers(inflmac, only.pos = F, min.pct = 0.2, logfc.threshold = 0.5)
FindAllMarkers_car <- FindAllMarkers(CARmac, only.pos = F, min.pct = 0.2, logfc.threshold = 0.5)
all.markers_top50.rest<- FindAllMarkers_rest %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
all.markers_top50.inflam<- FindAllMarkers_inflam %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
all.markers_top50.car<- FindAllMarkers_car %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
####rna is my seurat obj###

########resting analysis
data <- GetAssayData(Restingmac[["RNA"]], slot="data")
pd <- new('AnnotatedDataFrame', data = Restingmac@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
data <- data[, rownames(pd)]

rest.monocle2_cds <- newCellDataSet(data,
                                                  phenoData = pd,
                                                  featureData = fd,
                                                  lowerDetectionLimit = 0.5,
                                                  expressionFamily = negbinomial.size())
rest.monocle2_cds <- estimateSizeFactors(rest.monocle2_cds)
#my FindAllMarkers df is 'all.markers'
#my meta.data slot is "res_0.15_cmcleaned", you need to chnage this to your one
rest.monocle2_cds <- setOrderingFilter(rest.monocle2_cds, all.markers_top50.rest$gene)
rest.monocle2_cds <- reduceDimension(rest.monocle2_cds, max_components = 2, method = 'DDRTree')
rest.monocle2_cds <- orderCells(rest.monocle2_cds)
rest.monocle2_cds_plotclus <- plot_cell_trajectory(rest.monocle2_cds, color_by = "mac_clusters", cell_size = 1,show_branch_points = F)
rest.monocle2_cds_plotPseu <- plot_cell_trajectory(rest.monocle2_cds, color_by = "Pseudotime", cell_size = 1,show_branch_points = F)
###plot results###
rest.monocle2_cds_plotclus + NoAxes()+scale_color_simpsons(n, type = "discrete")+NoLegend()  ##this shows cluster ID
###plot results###
rest.monocle2_cds_plotPseu + NoAxes()#+scale_color_simpsons(n, type = "discrete")+NoLegend()  ##this shows cluster ID


####inflamed macs or non-treated
data <- GetAssayData(inflmac[["RNA"]], slot="data")
pd <- new('AnnotatedDataFrame', data = inflmac@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
data <- data[, rownames(pd)]

infl.monocle2_cds <- newCellDataSet(data,
                                    phenoData = pd,
                                    featureData = fd,
                                    lowerDetectionLimit = 0.5,
                                    expressionFamily = negbinomial.size())
infl.monocle2_cds <- estimateSizeFactors(infl.monocle2_cds)
#my FindAllMarkers df is 'all.markers'
#my meta.data slot is "res_0.15_cmcleaned", you need to chnage this to your one
infl.monocle2_cds <- setOrderingFilter(infl.monocle2_cds, all.markers_top50.inflam$gene)
infl.monocle2_cds <- reduceDimension(infl.monocle2_cds, max_components = 2, method = 'DDRTree')
infl.monocle2_cds <- orderCells(infl.monocle2_cds)
infl.monocle2_cds_plotclus <- plot_cell_trajectory(infl.monocle2_cds, color_by = "mac_clusters", cell_size = 1,show_branch_points = F)
infl.monocle2_cds_plotPseu <- plot_cell_trajectory(infl.monocle2_cds, color_by = "Pseudotime", cell_size = 1,show_branch_points = F)
###plot results###
infl.monocle2_cds_plotclus + NoAxes()+scale_color_simpsons(n, type = "discrete")+NoLegend()  ##this shows cluster ID
###plot results###
infl.monocle2_cds_plotPseu + NoAxes()#+scale_color_simpsons(n, type = "discrete")+NoLegend()  ##this shows cluster ID


########CAR analysis
data <- GetAssayData(CARmac[["RNA"]], slot="data")
pd <- new('AnnotatedDataFrame', data = CARmac@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
data <- data[, rownames(pd)]

car.monocle2_cds <- newCellDataSet(data,
                                    phenoData = pd,
                                    featureData = fd,
                                    lowerDetectionLimit = 0.5,
                                    expressionFamily = negbinomial.size())
car.monocle2_cds <- estimateSizeFactors(car.monocle2_cds)
#my FindAllMarkers df is 'all.markers'
#my meta.data slot is "res_0.15_cmcleaned", you need to chnage this to your one
car.monocle2_cds <- setOrderingFilter(car.monocle2_cds, all.markers_top50.car$gene)
car.monocle2_cds <- reduceDimension(car.monocle2_cds, max_components = 2, method = 'DDRTree')
car.monocle2_cds <- orderCells(car.monocle2_cds)
car.monocle2_cds_plotclus <- plot_cell_trajectory(car.monocle2_cds, color_by = "mac_clusters", cell_size = 1,show_branch_points = F)
car.monocle2_cds_plotPseu <- plot_cell_trajectory(car.monocle2_cds, color_by = "Pseudotime", cell_size = 1,show_branch_points = F)
###plot results###
car.monocle2_cds_plotclus + NoAxes()+scale_color_simpsons(n, type = "discrete")+NoLegend()  ##this shows cluster ID
###plot results###
car.monocle2_cds_plotPseu + NoAxes()#+scale_color_simpsons(n, type = "discrete")+NoLegend()  ##this shows cluster ID




####load firboblasts
fibs <- readRDS(file = "fibs.rds")
#####subset lining layer fibroblasts
liningfib<-subset(fibs, idents="Prg4")
#####subset lining macrophages
liningmac<-subset(mac1, idents="Vsig4")

##############figure 5f
######match metadata cluster names to "clusters"
current.sample.ids <- c("Prg4")
new.sample.ids <- c("Prg4")

liningfib$clusters <- liningfib@active.ident
liningfib@meta.data[["clusters"]] <- plyr::mapvalues(x = liningfib@meta.data[["clusters"]], from = current.sample.ids, to = new.sample.ids)
Idents(liningfib)<-"clusters"

######merge fibroblasts and macrophages
liningall <- merge(liningfib,liningmac)

#####subset tissue states
Idents(liningall)<-"State"
rest.lining<-subset(liningall, idents="Resting")
inflamed.lining<-subset(liningall, idents="Inflamed")
liningCAR<-subset(liningall, idents="CAR")


#####resting analysis
data.input <- GetAssayData(rest.lining, assay = "RNA", slot = "data") 
labels <- Idents(rest.lining)
meta <- data.frame(group = labels, row.names = names(labels))
cellchatrest <- createCellChat(object = data.input, meta = meta, group.by = "group")

cellchatrest <- addMeta(cellchatrest, meta = meta, meta.name = "labels")
cellchatrest <- setIdent(cellchatrest, ident.use = "labels") # set "labels" as default cell identity
levels(cellchatrest@idents) # show factor levels of the cell labels

CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB
cellchatrest@DB <- CellChatDB.use
cellchatrest <- subsetData(cellchatrest)
cellchatrest <- identifyOverExpressedGenes(cellchatrest)
cellchatrest <- identifyOverExpressedInteractions(cellchatrest)
cellchatrest <- projectData(cellchatrest, PPI.mouse)
cellchatrest <- computeCommunProb(cellchatrest)
cellchatrest <- filterCommunication(cellchatrest, min.cells = 10)
df.net <- subsetCommunication(cellchatrest)
cellchatrest <- computeCommunProbPathway(cellchatrest)
cellchatrest <- aggregateNet(cellchatrest)
groupSize <- as.numeric(table(cellchatrest@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchatrest@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchatrest@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
mat <- cellchatrest@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, arrow.width = 4, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
mat <- cellchatrest@net$weight
par(mfrow = c(2,5), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}


cellchatrest@netP$pathways

pathways.show <- c("CSF", "SPP1")

par(mfrow=c(1,1))


netVisual_aggregate(cellchatrest, signaling = pathways.show, layout = "chord", color.use = c("firebrick", "darkolivegreen3"), show.legend = T, vertex.label.cex = 1)

netAnalysis_contribution(cellchatrest, signaling = pathways.show)
netVisual_chord_gene(cellchatrest, sources.use = c(1, 2), targets.use = c(2, 1), signaling = pathways.show, lab.cex = 2,legend.pos.y = 0, small.gap = 2.5, color.use = c("firebrick", "darkolivegreen3"))+NoLegend()
netVisual_chord_gene(cellchatinfl, sources.use = c(1, 2), targets.use = c(2, 1), signaling = pathways.show, lab.cex = 2,legend.pos.y = 30, small.gap = 2.5, color.use = c("darkolivegreen3", "firebrick"))+NoLegend()


#####inflamed or non-treated analysis
data.input <- GetAssayData(inflamed.lining, assay = "RNA", slot = "data") 
labels <- Idents(inflamed.lining)
meta <- data.frame(group = labels, row.names = names(labels))
cellchatinfl <- createCellChat(object = data.input, meta = meta, group.by = "group")

cellchatinfl <- addMeta(cellchatinfl, meta = meta, meta.name = "labels")
cellchatinfl <- setIdent(cellchatinfl, ident.use = "labels") # set "labels" as default cell identity
levels(cellchatinfl@idents) # show factor levels of the cell labels

CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB
cellchatinfl@DB <- CellChatDB.use
cellchatinfl <- subsetData(cellchatinfl)
cellchatinfl <- identifyOverExpressedGenes(cellchatinfl)
cellchatinfl <- identifyOverExpressedInteractions(cellchatinfl)
cellchatinfl <- projectData(cellchatinfl, PPI.mouse)
cellchatinfl <- computeCommunProb(cellchatinfl)
cellchatinfl <- filterCommunication(cellchatinfl, min.cells = 10)
df.net <- subsetCommunication(cellchatinfl)
cellchatinfl <- computeCommunProbPathway(cellchatinfl)
cellchatinfl <- aggregateNet(cellchatinfl)
groupSize <- as.numeric(table(cellchatinfl@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchatinfl@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchatinfl@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
mat <- cellchatinfl@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, arrow.width = 4, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
mat <- cellchatinfl@net$weight
par(mfrow = c(2,5), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}


cellchaticar@netP$pathways

pathways.show <- c("CSF", "MIF", "SPP1", "TGFb", "GAS")

par(mfrow=c(1,1))


netVisual_aggregate(cellchatinfl, signaling = pathways.show, layout = "chord", color.use = c("firebrick", "darkolivegreen3"), show.legend = T, vertex.label.cex = 1)

netAnalysis_contribution(cellchatinfl, signaling = pathways.show)
netVisual_chord_gene(cellchatinfl, sources.use = c(1, 2), targets.use = c(2, 1), signaling = pathways.show, lab.cex = 0,legend.pos.y = 0, small.gap = 2.5, color.use = c("firebrick", "darkolivegreen3"))+NoLegend()


#####CAR analysis
data.input <- GetAssayData(liningCAR, assay = "RNA", slot = "data") 
labels <- Idents(liningCAR)
meta <- data.frame(group = labels, row.names = names(labels))
cellchatcar <- createCellChat(object = data.input, meta = meta, group.by = "group")

cellchatcar <- addMeta(cellchatcar, meta = meta, meta.name = "labels")
cellchatcar <- setIdent(cellchatcar, ident.use = "labels") # set "labels" as default cell identity
levels(cellchatcar@idents) # show factor levels of the cell labels

CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB
cellchatcar@DB <- CellChatDB.use
cellchatcar <- subsetData(cellchatcar)
cellchatcar <- identifyOverExpressedGenes(cellchatcar)
cellchatcar <- identifyOverExpressedInteractions(cellchatcar)
cellchatcar <- projectData(cellchatcar, PPI.mouse)
cellchatcar <- computeCommunProb(cellchatcar)
cellchatcar <- filterCommunication(cellchatcar, min.cells = 10)
df.net <- subsetCommunication(cellchatcar)
cellchatcar <- computeCommunProbPathway(cellchatcar)
cellchatcar <- aggregateNet(cellchatcar)
groupSize <- as.numeric(table(cellchatcar@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchatcar@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchatcar@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
mat <- cellchatcar@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, arrow.width = 4, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
mat <- cellchatcar@net$weight
par(mfrow = c(2,5), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}


cellchaticar@netP$pathways

pathways.show <- c("CSF", "MIF", "SPP1", "TGFb", "GAS")

par(mfrow=c(1,1))


netVisual_aggregate(cellchatcar, signaling = pathways.show, layout = "chord", color.use = c("firebrick", "darkolivegreen3"), show.legend = T, vertex.label.cex = 1)

netAnalysis_contribution(cellchatcar, signaling = pathways.show)
netVisual_chord_gene(cellchatcar, sources.use = c(1, 2), targets.use = c(2, 1), signaling = pathways.show, lab.cex = 0,legend.pos.y = 0, small.gap = 2.5, color.use = c("firebrick", "darkolivegreen3"))+NoLegend()


