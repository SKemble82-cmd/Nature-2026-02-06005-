############lines 0-293 sample intergation and processing. lines 294 onwards analysis for figures.

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
library(scater)
library(patchwork)
library(SingleCellExperiment)
library(monocle)
library(DDRTree)
library(SeuratWrappers)
library(devtools)
library(ggthemes)
library(tvthemes)
library(SeuratObject)
library(viridis)
library(Nebulosa)

#######data download, metadata, processing#######

CIA_r1 <- Read10X(data.dir = "/rds/projects/c/croftap-ktr-cia-sc-01/CIA_CD45neg/Con2_CD45neg/outs/filtered_feature_bc_matrix")
CIA_r1 <- CreateSeuratObject(counts = CIA_r1, min.cells=3, min.features=100, project="CIA_r1")
CIA_r1 <- RenameCells(CIA_r1, add.cell.id = "CIA_r1")
CIA_r1$sample_id <- "Resting1"
CIA_r1$CD45 <- "neg"
CIA_r1$State <- "Resting"
head(CIA_r1)

CIA_r2 <- Read10X(data.dir = "/rds/projects/c/croftap-ktr-cia-sc-01/CIA_CD45neg/Con3_CD45neg/outs/filtered_feature_bc_matrix")
CIA_r2 <- CreateSeuratObject(counts = CIA_r2, min.cells=3, min.features=100, project="CIA_r2")
CIA_r2 <- RenameCells(CIA_r2, add.cell.id = "CIA_r2")
CIA_r2$sample_id <- "Resting2"
CIA_r2$CD45 <- "neg"
CIA_r2$State <- "Resting"
head(CIA_r2)

CIA_r3 <- Read10X(data.dir = "/rds/projects/c/croftap-ktr-cia-sc-01/CIA_CD45_pos/Con4_CD45neg/outs/filtered_feature_bc_matrix")
CIA_r3 <- CreateSeuratObject(counts = CIA_r3, min.cells=3, min.features=100, project="CIA_r3")
CIA_r3 <- RenameCells(CIA_r3, add.cell.id = "CIA_r3")
CIA_r3$sample_id <- "Resting3"
CIA_r3$CD45 <- "neg"
CIA_r3$State <- "Resting"
head(CIA_r3)

####inflamed = non-treated
CIA_i1 <- Read10X(data.dir = "/rds/projects/c/croftap-sgcar/X204SC24074958-Z01-F001/count/CON1NEG/outs/per_sample_outs/CON1NEG/count/sample_filtered_feature_bc_matrix")
CIA_i1 <- CreateSeuratObject(counts = CIA_i1, min.cells=3, min.features=100, project="CIA_i1")
CIA_i1 <- RenameCells(CIA_i1, add.cell.id = "CIA_i1")
CIA_i1$sample_id <- "Inflamed1"
CIA_i1$CD45 <- "neg"
CIA_i1$State <- "Inflamed"
head(CIA_i1)

CIA_i2 <- Read10X(data.dir = "/rds/projects/c/croftap-sgcar/X204SC24074958-Z01-F001/count/CON2NEG/outs/per_sample_outs/CON2NEG/count/sample_filtered_feature_bc_matrix")
CIA_i2 <- CreateSeuratObject(counts = CIA_i2, min.cells=3, min.features=100, project="CIA_i2")
CIA_i2 <- RenameCells(CIA_i2, add.cell.id = "CIA_i2")
CIA_i2$sample_id <- "Inflamed2"
CIA_i2$CD45 <- "neg"
CIA_i2$State <- "Inflamed"
head(CIA_i2)

CIA_i3 <- Read10X(data.dir = "/rds/projects/c/croftap-sgcar/X204SC24074958-Z01-F001/count/CON3NEG/outs/per_sample_outs/CON3NEG/count/sample_filtered_feature_bc_matrix")
CIA_i3 <- CreateSeuratObject(counts = CIA_i3, min.cells=3, min.features=100, project="CIA_i3")
CIA_i3 <- RenameCells(CIA_i3, add.cell.id = "CIA_i3")
CIA_i3$sample_id <- "Inflamed3"
CIA_i3$CD45 <- "neg"
CIA_i3$State <- "Inflamed"
head(CIA_i3)

####CAR = FAPaCAR-T treated
CIA_c1 <- Read10X(data.dir = "/rds/projects/c/croftap-sgcar/X204SC24074958-Z01-F001/count/CAR1NEG/outs/per_sample_outs/CAR1NEG/count/sample_filtered_feature_bc_matrix")
CIA_c1 <- CreateSeuratObject(counts = CIA_c1, min.cells=3, min.features=100, project="CIA_c1")
CIA_c1 <- RenameCells(CIA_c1, add.cell.id = "CIA_c1")
CIA_c1$sample_id <- "CAR1"
CIA_c1$CD45 <- "neg"
CIA_c1$State <- "CAR"
head(CIA_c1)

CIA_c2 <- Read10X(data.dir = "/rds/projects/c/croftap-sgcar/X204SC24074958-Z01-F001/count/CAR2NEG/outs/per_sample_outs/CAR2NEG/count/sample_filtered_feature_bc_matrix")
CIA_c2 <- CreateSeuratObject(counts = CIA_c2, min.cells=3, min.features=100, project="CIA_c2")
CIA_c2 <- RenameCells(CIA_c2, add.cell.id = "CIA_c2")
CIA_c2$sample_id <- "CAR2"
CIA_c2$CD45 <- "neg"
CIA_c2$State <- "CAR"
head(CIA_c2)

CIA_c3 <- Read10X(data.dir = "/rds/projects/c/croftap-sgcar/X204SC24074958-Z01-F001/count/CAR3NEG/outs/per_sample_outs/CAR3NEG/count/sample_filtered_feature_bc_matrix")
CIA_c3 <- CreateSeuratObject(counts = CIA_c3, min.cells=3, min.features=100, project="CIA_c3")
CIA_c3 <- RenameCells(CIA_c3, add.cell.id = "CIA_c3")
CIA_c3$sample_id <- "CAR3"
CIA_c3$CD45 <- "neg"
CIA_c3$State <- "CAR"
head(CIA_c3)

#####QC and cell selection#######
mito.featuresCIA_r1 <- grep(pattern="^mt-", x=rownames(x=CIA_r1), value=T)
percent.mitoCIA_r1 <- Matrix::colSums(x = GetAssayData(object = CIA_r1, slot = "counts")[mito.featuresCIA_r1,]) / Matrix::colSums(x = GetAssayData(object = CIA_r1, slot = "counts"))
CIA_r1[["percent.mito"]] <- percent.mitoCIA_r1
VlnPlot(object = CIA_r1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

CIA_r1 <- subset(x = CIA_r1, subset = nFeature_RNA > 200 & nFeature_RNA <4000 & nCount_RNA > 200 & nCount_RNA < 25000 & percent.mitoCIA_r1 < 0.1)
VlnPlot(object = CIA_r1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
CIA_r1 <- NormalizeData(object = CIA_r1, verbose = F)
CIA_r1 <- FindVariableFeatures(object = CIA_r1, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesCIA_r1 <- rownames(CIA_r1)
CIA_r1 <- ScaleData(CIA_r1, features = all.genesCIA_r1)
VlnPlot(object = CIA_r1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)


mito.featuresCIA_r2 <- grep(pattern="^mt-", x=rownames(x=CIA_r2), value=T)
percent.mitoCIA_r2 <- Matrix::colSums(x = GetAssayData(object = CIA_r2, slot = "counts")[mito.featuresCIA_r2,]) / Matrix::colSums(x = GetAssayData(object = CIA_r2, slot = "counts"))
CIA_r2[["percent.mito"]] <- percent.mitoCIA_r2
VlnPlot(object = CIA_r2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

CIA_r2 <- subset(x = CIA_r2, subset = nFeature_RNA > 500 & nFeature_RNA <6000 & nCount_RNA > 200 & nCount_RNA < 30000 & percent.mitoCIA_r2 < 0.1)
VlnPlot(object = CIA_r2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
CIA_r2 <- NormalizeData(object = CIA_r2, verbose = F)
CIA_r2 <- FindVariableFeatures(object = CIA_r2, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesCIA_r2 <- rownames(CIA_r2)
CIA_r2 <- ScaleData(CIA_r2, features = all.genesCIA_r2)
VlnPlot(object = CIA_r2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)


mito.featuresCIA_r3 <- grep(pattern="^mt-", x=rownames(x=CIA_r3), value=T)
percent.mitoCIA_r3 <- Matrix::colSums(x = GetAssayData(object = CIA_r3, slot = "counts")[mito.featuresCIA_r3,]) / Matrix::colSums(x = GetAssayData(object = CIA_r3, slot = "counts"))
CIA_r3[["percent.mito"]] <- percent.mitoCIA_r3
VlnPlot(object = CIA_r3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

CIA_r3 <- subset(x = CIA_r3, subset = nFeature_RNA > 500 & nFeature_RNA <5000 & nCount_RNA > 200 & nCount_RNA < 25000 & percent.mitoCIA_r3 < 0.1)
VlnPlot(object = CIA_r3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
CIA_r3 <- NormalizeData(object = CIA_r3, verbose = F)
CIA_r3 <- FindVariableFeatures(object = CIA_r3, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesCIA_r3 <- rownames(CIA_r3)
CIA_r3 <- ScaleData(CIA_r3, features = all.genesCIA_r3)
VlnPlot(object = CIA_r3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)


mito.featuresCIA_i1 <- grep(pattern="^mt-", x=rownames(x=CIA_i1), value=T)
percent.mitoCIA_i1 <- Matrix::colSums(x = GetAssayData(object = CIA_i1, slot = "counts")[mito.featuresCIA_i1,]) / Matrix::colSums(x = GetAssayData(object = CIA_i1, slot = "counts"))
CIA_i1[["percent.mito"]] <- percent.mitoCIA_i1
VlnPlot(object = CIA_i1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0)

CIA_i1 <- subset(x = CIA_i1, subset = nFeature_RNA > 2500 & nFeature_RNA <75000 & nCount_RNA > 500 & nCount_RNA < 35000 & percent.mitoCIA_i1 < 0.1)
VlnPlot(object = CIA_i1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
CIA_i1 <- NormalizeData(object = CIA_i1, verbose = F)
CIA_i1 <- FindVariableFeatures(object = CIA_i1, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesCIA_i1 <- rownames(CIA_i1)
CIA_i1 <- ScaleData(CIA_i1, features = all.genesCIA_i1)
VlnPlot(object = CIA_i1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)


mito.featuresCIA_i2 <- grep(pattern="^mt-", x=rownames(x=CIA_i2), value=T)
percent.mitoCIA_i2 <- Matrix::colSums(x = GetAssayData(object = CIA_i2, slot = "counts")[mito.featuresCIA_i2,]) / Matrix::colSums(x = GetAssayData(object = CIA_i2, slot = "counts"))
CIA_i2[["percent.mito"]] <- percent.mitoCIA_i2
VlnPlot(object = CIA_i2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0)

CIA_i2 <- subset(x = CIA_i2, subset = nFeature_RNA > 2000 & nFeature_RNA <7000 & nCount_RNA > 500 & nCount_RNA < 25000 & percent.mitoCIA_i2 < 0.1)
VlnPlot(object = CIA_i2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
CIA_i2 <- NormalizeData(object = CIA_i2, verbose = F)
CIA_i2 <- FindVariableFeatures(object = CIA_i2, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesCIA_i2 <- rownames(CIA_i2)
CIA_i2 <- ScaleData(CIA_i2, features = all.genesCIA_i2)
VlnPlot(object = CIA_i2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)


mito.featuresCIA_i3 <- grep(pattern="^mt-", x=rownames(x=CIA_i3), value=T)
percent.mitoCIA_i3 <- Matrix::colSums(x = GetAssayData(object = CIA_i3, slot = "counts")[mito.featuresCIA_i3,]) / Matrix::colSums(x = GetAssayData(object = CIA_i3, slot = "counts"))
CIA_i3[["percent.mito"]] <- percent.mitoCIA_i3
VlnPlot(object = CIA_i3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0)

CIA_i3 <- subset(x = CIA_i3, subset = nFeature_RNA > 2500 & nFeature_RNA <7500 & nCount_RNA > 500 & nCount_RNA < 40000 & percent.mitoCIA_i3 < 0.1)
VlnPlot(object = CIA_i3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
CIA_i3 <- NormalizeData(object = CIA_i3, verbose = F)
CIA_i3 <- FindVariableFeatures(object = CIA_i3, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesCIA_i3 <- rownames(CIA_i3)
CIA_i3 <- ScaleData(CIA_i3, features = all.genesCIA_i3)
VlnPlot(object = CIA_i3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)


mito.featuresCIA_c1 <- grep(pattern="^mt-", x=rownames(x=CIA_c1), value=T)
percent.mitoCIA_c1 <- Matrix::colSums(x = GetAssayData(object = CIA_c1, slot = "counts")[mito.featuresCIA_c1,]) / Matrix::colSums(x = GetAssayData(object = CIA_c1, slot = "counts"))
CIA_c1[["percent.mito"]] <- percent.mitoCIA_c1
VlnPlot(object = CIA_c1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0)

CIA_c1 <- subset(x = CIA_c1, subset = nFeature_RNA > 2500 & nFeature_RNA <7500 & nCount_RNA > 500 & nCount_RNA < 40000 & percent.mitoCIA_c1 < 0.1)
VlnPlot(object = CIA_c1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
CIA_c1 <- NormalizeData(object = CIA_c1, verbose = F)
CIA_c1 <- FindVariableFeatures(object = CIA_c1, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesCIA_c1 <- rownames(CIA_c1)
CIA_c1 <- ScaleData(CIA_c1, features = all.genesCIA_c1)
VlnPlot(object = CIA_c1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)


mito.featuresCIA_c2 <- grep(pattern="^mt-", x=rownames(x=CIA_c2), value=T)
percent.mitoCIA_c2 <- Matrix::colSums(x = GetAssayData(object = CIA_c2, slot = "counts")[mito.featuresCIA_c2,]) / Matrix::colSums(x = GetAssayData(object = CIA_c2, slot = "counts"))
CIA_c2[["percent.mito"]] <- percent.mitoCIA_c2
VlnPlot(object = CIA_c2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0)

CIA_c2 <- subset(x = CIA_c2, subset = nFeature_RNA > 2500 & nFeature_RNA <7500 & nCount_RNA > 500 & nCount_RNA < 40000 & percent.mitoCIA_c2 < 0.1)
VlnPlot(object = CIA_c2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
CIA_c2 <- NormalizeData(object = CIA_c2, verbose = F)
CIA_c2 <- FindVariableFeatures(object = CIA_c2, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesCIA_c2 <- rownames(CIA_c2)
CIA_c2 <- ScaleData(CIA_c2, features = all.genesCIA_c2)
VlnPlot(object = CIA_c2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)


mito.featuresCIA_c3 <- grep(pattern="^mt-", x=rownames(x=CIA_c3), value=T)
percent.mitoCIA_c3 <- Matrix::colSums(x = GetAssayData(object = CIA_c3, slot = "counts")[mito.featuresCIA_c3,]) / Matrix::colSums(x = GetAssayData(object = CIA_c3, slot = "counts"))
CIA_c3[["percent.mito"]] <- percent.mitoCIA_c3
VlnPlot(object = CIA_c3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0)

CIA_c3 <- subset(x = CIA_c3, subset = nFeature_RNA > 2500 & nFeature_RNA <7500 & nCount_RNA > 500 & nCount_RNA < 40000 & percent.mitoCIA_c3 < 0.1)
VlnPlot(object = CIA_c3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
CIA_c3 <- NormalizeData(object = CIA_c3, verbose = F)
CIA_c3 <- FindVariableFeatures(object = CIA_c3, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesCIA_c3 <- rownames(CIA_c3)
CIA_c3 <- ScaleData(CIA_c3, features = all.genesCIA_c3)
VlnPlot(object = CIA_c3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)


## Identification of integration anchors
reference.list <- c(CIA_r1, CIA_r2, CIA_r3, CIA_i1, CIA_i2, CIA_i3, CIA_c1, CIA_c2, CIA_c3)

anchors <- FindIntegrationAnchors(object.list = reference.list, anchor.features = 2000, dims = 1:30)
integrated <- IntegrateData(anchorset = anchors, dims = 1:30)

## Identification of integration anchors
reference.list <- c(CIA_i1, CIA_i2, CIA_i3)

anchors.infl <- FindIntegrationAnchors(object.list = reference.list, anchor.features = 2000, dims = 1:30)
integrated.infl <- IntegrateData(anchorset = anchors, dims = 1:30)


#####save rds object#####
getwd()
saveRDS(integrated, file = "cd45neg.rds")


#####sample crossover after integration and stromal cell annotation#####
DefaultAssay(object=integrated) <- "integrated"
integrated <- FindVariableFeatures(integrated, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(integrated)
integrated <- ScaleData(integrated, features = all.genes)
integrated <- RunPCA(object = integrated, verbose=F, npcs = 50)
ElbowPlot(object = integrated, ndims = 50)
integrated <- FindNeighbors(object = integrated, dims = 1:30)

integrated <- RunUMAP(object = integrated, reduction = "pca", dims = 1:30)
DimPlot(integrated, reduction = "umap", label = TRUE, pt.size = 0.1, repel = F, split.by = "State") + NoAxes()

integrated <- FindClusters(integrated, graph.name = "integrated_snn", resolution = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7))
clustree(integrated, assay = "integrated")
integrated <- FindClusters(integrated, graph.name = "integrated_snn", resolution = 0.1)
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:30)
DimPlot(integrated, reduction = "umap", label = TRUE, pt.size = 0.1, repel = T, label.box = T) + NoAxes()

current.sample.ids <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11")
new.sample.ids <- c("Fibroblastic", "Fibroblastic", "Fibroblastic", "Con", "Endothelial", "Osteoblastic", "Mural",
                    "Con", "Con", "Con", "Glial", "Lymphatic")
#####con = contaminants 

integrated$clusters <- integrated@active.ident
integrated@meta.data[["clusters"]] <- plyr::mapvalues(x = integrated@meta.data[["clusters"]], from = current.sample.ids, to = new.sample.ids)
Idents(integrated)<-"clusters"

######extended data fig 4b&c
####subset stromal cells only######
int.cells <- subset(integrated, idents = c("Fibroblastic", "Endothelial", "Osteoblastic", "Mural", "Lymphatic"))
int.cells$State <- factor(int.cells$State,levels=c("Resting", "Inflamed", "CAR"))
DimPlot(int.cells, reduction = "umap", label = F, pt.size = 0.1, repel = F, label.box = T, group.by = "sample_id")+ NoAxes()+scale_color_avatar(palette = "EarthKingdom")
DimPlot(int.cells, reduction = "umap", label = F, pt.size = 0.1, repel = F, label.box = T)+ NoAxes()+scale_color_westeros(palette = "Greyjoy")
DimPlot(int.cells, reduction = "umap", label = F, pt.size = 0.1, repel = F, label.box = T, split.by = "State") + NoAxes()+scale_color_westeros(palette = "Greyjoy")

FindAllMarkers_stromal <- FindAllMarkers(int.cells, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.5)
write.csv(FindAllMarkers_stromal, file = "Supplementary Table 3")
stromal.genes<-c("Pdgfra", "Prg4", "Pdpn", "Pecam1", "Ptprb", "Cd36", "Bglap", "Omd", "Bglap2", "Notch3", "Rgs5", "Des", "Prox1", "Lyve1", "Reln")
#plot dotplot

dotplot <- DotPlot(int.cells, features = stromal.genes, dot.scale = 6)+ theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) +ylab(NULL) +xlab(NULL)+ scale_color_gradientn(colours = magma(20))+ FontSize(x.text = 14, y.text = 14) 
dotplot

#######subset and process fibroblasts######
fibs<-subset(int.cells, idents="Fibroblastic")

DefaultAssay(object=fibs) <- "integrated"
fibs <- FindVariableFeatures(fibs, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(fibs)
fibs <- ScaleData(fibs, features = all.genes)
DefaultAssay(object=fibs) <- "integrated"
fibs <- RunPCA(fibs, npcs = 50, features = VariableFeatures(object = fibs))

ElbowPlot(object = fibs, ndims = 50)
fibs <- FindNeighbors(object = fibs, dims = 1:30)
fibs <- FindClusters(fibs, graph.name = "integrated_snn", resolution = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7))
clustree(fibs, assay = "integrated")
fibs <- FindClusters(fibs, graph.name = "integrated_snn", resolution = 0.2)
fibs <- RunUMAP(fibs, reduction = "pca", dims = 1:30)
DimPlot(fibs, reduction = "umap", label = F, pt.size = 0.1, repel = T, label.box = F) + NoAxes()

fib_clusterssaveRDS(fibs.res, file = "fibs.rds")


F0 <- FindMarkers(fibs, ident.1 = 0, min.pct = 0.2, only.pos = T)
F1 <- FindMarkers(fibs, ident.1 = 1, min.pct = 0.2, only.pos = T)
F2 <- FindMarkers(fibs, ident.1 = 2, min.pct = 0.2, only.pos = T)
F3 <- FindMarkers(fibs, ident.1 = 3, min.pct = 0.2, only.pos = T)
F4 <- FindMarkers(fibs, ident.1 = 4, min.pct = 0.2, only.pos = T)
F5 <- FindMarkers(fibs, ident.1 = 5, min.pct = 0.2, only.pos = T)
F6 <- FindMarkers(fibs, ident.1 = 6, min.pct = 0.2, only.pos = T)
F7 <- FindMarkers(fibs, ident.1 = 7, min.pct = 0.2, only.pos = T)
F8 <- FindMarkers(fibs, ident.1 = 8, min.pct = 0.2, only.pos = T)
F9 <- FindMarkers(fibs, ident.1 = 9, min.pct = 0.2, only.pos = T)

####remove contaminant cluster 7, 8, 9#######
fibs.res<-subset(fibs, idents=c("0", "1", "2", "3", "4", "5", "6"))

DefaultAssay(object=fibs.res) <- "integrated"
fibs.res <- FindVariableFeatures(fibs.res, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(fibs.res)
fibs.res <- ScaleData(fibs.res, features = all.genes)
DefaultAssay(object=fibs.res) <- "integrated"
fibs.res <- RunPCA(fibs.res, npcs = 50, features = VariableFeatures(object = fibs))

ElbowPlot(object = fibs.res, ndims = 50)
fibs.res <- FindNeighbors(object = fibs.res, dims = 1:30)
fibs.res <- FindClusters(fibs.res, graph.name = "integrated_snn", resolution = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7))
clustree(fibs.res, assay = "integrated")
fibs.res <- FindClusters(fibs.res, graph.name = "integrated_snn", resolution = 0.2)
fibs.res <- RunUMAP(fibs.res, reduction = "pca", dims = 1:30)

#########figure 2f and Extended data figure 4d
DimPlot(fibs.res, reduction = "umap", label = F, pt.size = 0.1, repel = T, label.box = F) + NoAxes()
DimPlot(fibs.res, reduction = "umap", label = F, pt.size = 0.1, repel = T, label.box = F) + NoAxes()+scale_color_avatar(palette = "FireNation")
DimPlot(fibs.res, reduction = "umap", label = F, pt.size = 0.1, repel = T, label.box = F, split.by = "State") + NoAxes()+scale_color_avatar(palette = "FireNation")


FindAllMarkers_FBs <- FindAllMarkers(fibs.res, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.5)
write.csv(FindAllMarkers_FBs, file = "Supplementary Table 5")

###rename clusters#####
current.sample.ids <- c("0", "1", "2", "3", "4", "5", "6")
new.sample.ids <- c("Cxcl5", "Postn", "Ifi204", "Pi16", "Prg4", "Fmod")
fibs.res$fib_clusters <- fibs.res@active.ident
fibs.res@meta.data[["fib_clusters"]] <- plyr::mapvalues(x = fibs.res@meta.data[["fib_clusters"]], from = current.sample.ids, to = new.sample.ids)
Idents(fibs.res)<-"fib_clusters"

#######Extended data figure e
fib.genes<-c("Cthrc1", "Postn", "Cxcl5", "Tnn", "Ifi204", "Cd34", "Pi16", "Ly6c1", "Prg4", "Cd55", "Fmod", "Cd200", "Top2a", "Mki67")

#plot dotplot
dotplot <- DotPlot(fibs.res, features = fib.genes, dot.scale = 6)+ theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) +ylab(NULL) +xlab(NULL)+ scale_color_gradientn(colours = magma(20))+ FontSize(x.text = 14, y.text = 14) 
dotplot


~###check gene expression for each cluster#### eg below using violin plot 
VlnPlot(fibs.res, features = "Cxcl5", pt.size = 0, cols = c("#ecb100", "#a10000","#7E605E", "#FF4500", "#994823", "#4B4C4E","#572530"))+xlab(NULL)

#####figure 2f
#####fap expression density plot########
plot_density(fibs.res, features = "Fap", combine = F, method = "wkde", size = 1, adjust = 1)+NoAxes()+ scale_color_gradientn(colours = viridis(5))+facet_grid(.~fibs.res1$State)

save.image("/rds/projects/c/croftap-croftapcarcia/Cd45neg.RData")

library(tvthemes)

#######extended data figure 4g
########bar chart percentage fibroblasts
pt <- table(fibs.res$fib_clusters, fibs.res$State)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)
pt$Var1 <- factor(pt$Var1,levels=c("Postn", "Cxcl5", "Ifi204", "Pi16", "Prg4","Fmod", "Cycling"))

ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab(NULL) +
  ylab("Proportion") +
  theme(legend.title = element_blank())+RotatedAxis() +
  scale_fill_avatar(palette = "FireNation")+NoLegend()


#######figure 2h
#########proportional change in fibroblasts######
test <- sc_utils(fibs.res)

###car v non-treated####
prop.test <- permutation_test(test, cluster_identity = "fib_clusters", sample_1="CAR", sample_2="Inflamed", sample_identity="State", n_permutations=10000)
permutation_plot(prop.test, FDR_threshold = 0.01, log2FD_threshold = 0.58, order_clusters = T)+ FontSize(x.text = 12, y.text = 12, x.title = 12, y.title = 12)+xlab(NULL)+NoLegend()

###########extended data figure 4h
###resting v non-treated####
prop.test <- permutation_test(test, cluster_identity = "fib_clusters", sample_1="Resting", sample_2="Inflamed", sample_identity="State", n_permutations=10000)
permutation_plot(prop.test, FDR_threshold = 0.01, log2FD_threshold = 0.58, order_clusters = T)+ FontSize(x.text = 12, y.text = 12, x.title = 12, y.title = 12)+xlab(NULL)+NoLegend()

##########figure 2i extended data figure 4f
#########GO-gsfisher##########
library(gsfisher)
seurat_obj <- fibs.res
#Edit to suitable output folder
getwd() 
GSout <- "/rds/projects/c/croftap-croftapcarcia/GOfibsNEW"

###saved as annotation.rds do not run again if already done so
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
expressed_genes <- getExpressedGenesFromSeuratObjectRNA(seurat_obj,levels(seurat_obj@active.ident))


seurat_obj.res <- FindAllMarkers(seurat_obj, only.pos = T, min.pct = 0.1) 
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

####All_enrichments is Supplementary Table 6.


##plots
all_results_top1 <- All_enrichments %>% group_by(cluster) %>% top_n(n=5, -p.val) 
all_results_top2 <- All_enrichments %>%
  filter(cluster %in% c("Cxcl5", "Fmod", "Ifi204")) %>% #####repeat for clusters c("Cycling", "Postn", "Pi16", "Prg4")
  group_by(cluster) 

sampleEnrichmentDotplot(all_results_top1, selection_col = "description", selected_genesets = c("regulation of lymphocyte migration", "negative regulation of lymphocyte migration", "regulation of T cell migration", "mesenchymal cell proliferation", "regulation of mesenchymal cell proliferation", "positive regulation of mesenchymal cell proliferation", "mitochondrial translation", "endopeptidase complex", "proteasome complex"), fill_colors = c("yellow3", "brown", "black"), sample_id_col = "cluster", fill_var = "odds.ratio", maxl=50, rotate_sample_labels = TRUE, min_dot_size = 5, pvalue_threshold = 0.05) + FontSize(x.text = 12, y.text = 12) + ylab(NULL) +xlab(NULL)

####select description for appropriate clusters#####
"cytokine activity", "cellular response to chemokine", "cell chemotaxis", "regulation of metallopeptidase activity", "metalloendopeptidase inhibitor activity", "negative regulation of mononuclear cell migration", "platelet-derived growth factor beta-receptor activity", "stem cell factor receptor activity", "macrophage colony-stimulating factor receptor activity", 
"regulation of lymphocyte migration", "negative regulation of lymphocyte migration", "regulation of T cell migration", "mesenchymal cell proliferation", "regulation of mesenchymal cell proliferation", "positive regulation of mesenchymal cell proliferation", "mitochondrial translation", "endopeptidase complex", "proteasome complex"


######FAP####hi#####low#####analysis

#######inflam is the object

######subset fap>0.5 & <0.5

fibs.res1<-subset(fibs.res, idents = c("Cxcl5", "Postn", "Ifi204", "Pi16", "Prg4", "Fmod"))
Idents(fibs.res1)<-"State"
inflam <- subset(fibs.res1, idents = "Inflamed")
Fap_hi <- subset(inflam, subset = Fap >0.5)
Fap_low <- subset(inflam, subset = Fap <0.5)

Fap_hi$Fap <- "High"
Fap_no$Fap <- "Low"

Idents(Fap_hi)<-"Fap"
Idents(Fap_no)<-"Fap"

fap_all <- merge(Fap_hi, Fap_no)
F_DEGs <- FindAllMarkers(fap_all, only.pos = T, min.pct = 0.1, logfc.threshold = 0.05)
write.csv(F_DEGs, file = "Supplementary Table 4")


library(AnnotationDbi)
library(org.Mm.eg.db)
library(clusterProfiler)

degs_hi <- rownames(subset(F_DEGs, p_val_adj < 0.05 & abs(avg_log2FC) > 0.1 & cluster == "High"))
degs_low <- rownames(subset(F_DEGs, p_val_adj < 0.05 & abs(avg_log2FC) > 0.1 & cluster == "Low"))

degs_low <- gsub("\\..*", "", degs_low)

#########figure 2d
ego.hi <- enrichGO(gene = degs_hi,
                  OrgDb = org.Mm.eg.db,
                  keyType = "SYMBOL",
                  ont = "BP",
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.05,
                  readable = TRUE)
dotplot(ego.hi, showCategory=c("tissue remodeling", "ERK1 and ERK2 cascade", "complement activation", "Wnt signaling pathway", "cell-cell signaling by wnt", "canonical Wnt signaling pathway", "bone remodeling", "positive regulation of inflammatory response", "cell chemotaxis", "learning or memory"))

##########figure 2e
ego.low <- enrichGO(gene = degs_low,
                OrgDb = org.Mm.eg.db,
                keyType = "SYMBOL",
                ont = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = TRUE)
dotplot(ego.low, showCategory = c("wound healing", "regulation of peptidase activity", "regulation of wound healing", "glycolytic process through glucose-6-phosphate", "regulation of endopeptidase activity", "regulation of response to wounding", "positive regulation of cell-matrix adhesion", "myeloid cell homeostasis", "regulation of fibroblast proliferation", "hemostasis"))

############figure 2
#####priming gene module######
prime.genes<-c("Ccl2", "Ccl7", "Cxcl1", "Cxcl2", "Ptges", "Ptgs2", "Il1b", "Nlrp3", "Tnfrsf11b", "Csf1", "Mmp9", "C3", "Hc", "Ctsl", "Slc2a1")

fibs.res1 <- AddModuleScore(
  object = fibs.res1,
  features = list(prime.genes),
  ctrl = 5,
  name = 'prime')

DotPlot(
  fibs.res1,
  features = c("prime1"),
  assay = NULL,
  dot.scale = 6
) +
  scale_color_gradientn(
    colors = c("skyblue", "lightgrey", "red"),
    values = scales::rescale(c(-1.5, 0, 1.5))
  )+xlab(NULL)+ylab(NULL)

Idents(fibs.res1)<-"State"
per.fibs <- subset(fibs.res1, idents = c("Inflamed", "CAR"))
#repeat dotpot above#

#######extended figure 4i&j
######Cd200 gene module#####
Idents(fibs.res1)<-"fib_clusters"

cd200.genes<-c("Cd200", "Ogn", "Comp", "Ndnf", "Hhip", "Prelp", "Dkk3", "Sfrp1", "Cilp")

fibs.res1 <- AddModuleScore(
  object = fibs.res1,
  features = list(cd200.genes),
  ctrl = 5,
  name = 'cd200')

DotPlot(
  fibs.res1,
  features = c("cd2001"),
  assay = NULL,
  dot.scale = 6
) +
  scale_color_gradientn(
    colors = c("skyblue", "lightgrey", "red"),
    values = scales::rescale(c(-1.5, 0, 1.5))
  )+xlab(NULL)+ylab(NULL)

Idents(fibs.res1)<-"State"
#repeat dotplot above#

