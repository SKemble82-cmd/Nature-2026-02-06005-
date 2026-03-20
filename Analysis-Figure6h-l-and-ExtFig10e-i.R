###############lines 0-351 samples intergation and processing. lines 353 onwards analysis for figures

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
#library(velocyto.R)
library(SeuratWrappers)
library(devtools)
library(SCpubr)
library(dichromat)
library(ggpubr)
library(Nebulosa)
        

###########non-treated
p1 <- Read10X(data.dir = "/rds/projects/c/croftap-sgcar/single_cell_feb_2024/count/SGP1/outs/filtered_feature_bc_matrix")
p1 <- CreateSeuratObject(counts = p1, min.cells=3, min.features=100, project="p1")
p1 <- RenameCells(p1, add.cell.id = "p1")
p1$sample_id <- "PBS1"
p1$Treatment <- "Control"
p1$Inflammation <- "Inflamed"
head(p1)

p2 <- Read10X(data.dir = "/rds/projects/c/croftap-sgcar/single_cell_feb_2024/count/SGP2/outs/filtered_feature_bc_matrix")
p2 <- CreateSeuratObject(counts = p2, min.cells=3, min.features=100, project="p2")
p2 <- RenameCells(p2, add.cell.id = "p2")
p2$sample_id <- "PBS2"
p2$Treatment <- "Control"
p2$Inflammation <- "Inflamed"
head(p2)

p3 <- Read10X(data.dir = "/rds/projects/c/croftap-sgcar/single_cell_feb_2024/count/SGP3/outs/filtered_feature_bc_matrix")
p3 <- CreateSeuratObject(counts = p3, min.cells=3, min.features=100, project="p3")
p3 <- RenameCells(p3, add.cell.id = "p3")
p3$sample_id <- "PBS3"
p3$Treatment <- "Control"
p3$Inflammation <- "Inflamed"
head(p3)

############CAR-T treated
f1 <- Read10X(data.dir = "/rds/projects/c/croftap-sgcar/count/SGC1/outs/per_sample_outs/SGC1/count/sample_filtered_feature_bc_matrix")
f1 <- CreateSeuratObject(counts = f1, min.cells=3, min.features=100, project="f1")
f1 <- RenameCells(f1, add.cell.id = "f1")
f1$sample_id <- "FAP1"
f1$Treatment <- "CAR"
f1$Inflammation <- "CAR_Inflamed"
head(f1)

f2 <- Read10X(data.dir = "/rds/projects/c/croftap-sgcar/count/SGC2/outs/per_sample_outs/SGC2/count/sample_filtered_feature_bc_matrix")
f2 <- CreateSeuratObject(counts = f2, min.cells=3, min.features=100, project="f2")
f2 <- RenameCells(f2, add.cell.id = "f2")
f2$sample_id <- "FAP2"
f2$Treatment <- "CAR"
f2$Inflammation <- "CAR_Inflamed"
head(f2)

f3 <- Read10X(data.dir = "/rds/projects/c/croftap-sgcar/count/SGC3/outs/per_sample_outs/SGC3/count/sample_filtered_feature_bc_matrix")
f3 <- CreateSeuratObject(counts = f3, min.cells=3, min.features=100, project="f3")
f3 <- RenameCells(f3, add.cell.id = "f3")
f3$sample_id <- "FAP3"
f3$Treatment <- "CAR"
f3$Inflammation <- "CAR_Inflamed"
head(f3)

################resting
h1 <- Read10X(data.dir = "/rds/projects/c/croftap-sgcar/count/SGH1/outs/per_sample_outs/SGH1/count/sample_filtered_feature_bc_matrix")
h1 <- CreateSeuratObject(counts = h1, min.cells=3, min.features=100, project="h1")
h1 <- RenameCells(h1, add.cell.id = "h1")
h1$sample_id <- "REST1"
h1$Treatment <- "Naive"
h1$Inflammation <- "Resting"
head(h1)

h2 <- Read10X(data.dir = "/rds/projects/c/croftap-sgcar/count/SGH2/outs/per_sample_outs/SGH2/count/sample_filtered_feature_bc_matrix")
h2 <- CreateSeuratObject(counts = h2, min.cells=3, min.features=100, project="h2")
h2 <- RenameCells(h2, add.cell.id = "h2")
h2$sample_id <- "REST2"
h2$Treatment <- "Naive"
h2$Inflammation <- "Resting"
head(h2)

h3 <- Read10X(data.dir = "/rds/projects/c/croftap-sgcar/count/SGH3/outs/per_sample_outs/SGH3/count/sample_filtered_feature_bc_matrix")
h3 <- CreateSeuratObject(counts = h3, min.cells=3, min.features=100, project="h3")
h3 <- RenameCells(h3, add.cell.id = "h3")
h3$sample_id <- "REST3"
h3$Treatment <- "Naive"
h3$Inflammation <- "Resting"
head(h3)

###QC
mito.featuresbsa_p1 <- grep(pattern="^mt-", x=rownames(x=p1), value=T)
percent.mitobsa_p1 <- Matrix::colSums(x = GetAssayData(object = p1, slot = "counts")[mito.featuresbsa_p1,]) / Matrix::colSums(x = GetAssayData(object = p1, slot = "counts"))
p1[["percent.mito"]] <- percent.mitobsa_p1
VlnPlot(object = p1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

p1 <- subset(x = p1, subset = nFeature_RNA > 1000 & nFeature_RNA <6000 & nCount_RNA > 500 & nCount_RNA < 25000 & percent.mitobsa_p1 < 0.05)
VlnPlot(object = p1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
p1 <- NormalizeData(object = p1, verbose = F)
p1 <- FindVariableFeatures(object = p1, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesbsa_p1 <- rownames(p1)
p1 <- ScaleData(p1, features = all.genesbsa_p1)
VlnPlot(object = p1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0)

mito.featuresbsa_p2 <- grep(pattern="^mt-", x=rownames(x=p2), value=T)
percent.mitobsa_p2 <- Matrix::colSums(x = GetAssayData(object = p2, slot = "counts")[mito.featuresbsa_p2,]) / Matrix::colSums(x = GetAssayData(object = p2, slot = "counts"))
p2[["percent.mito"]] <- percent.mitobsa_p2
VlnPlot(object = p2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0)

p2 <- subset(x = p2, subset = nFeature_RNA > 1000 & nFeature_RNA <6000 & nCount_RNA > 500 & nCount_RNA < 20000 & percent.mitobsa_p2 < 0.1)
VlnPlot(object = p2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
p2 <- NormalizeData(object = p2, verbose = F)
p2 <- FindVariableFeatures(object = p2, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesbsa_p2 <- rownames(p2)
p2 <- ScaleData(p2, features = all.genesbsa_p2)
VlnPlot(object = p2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0)

mito.featuresbsa_p3 <- grep(pattern="^mt-", x=rownames(x=p3), value=T)
percent.mitobsa_p3 <- Matrix::colSums(x = GetAssayData(object = p3, slot = "counts")[mito.featuresbsa_p3,]) / Matrix::colSums(x = GetAssayData(object = p3, slot = "counts"))
p3[["percent.mito"]] <- percent.mitobsa_p3
VlnPlot(object = p3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3,pt.size = 0)

p3 <- subset(x = p3, subset = nFeature_RNA > 1000 & nFeature_RNA <7000 & nCount_RNA > 500 & nCount_RNA < 25000 & percent.mitobsa_p3 < 0.1)
VlnPlot(object = p3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0)
p3 <- NormalizeData(object = p3, verbose = F)
p3 <- FindVariableFeatures(object = p3, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesbsa_p3 <- rownames(p3)
p3 <- ScaleData(p3, features = all.genesbsa_p3)
VlnPlot(object = p3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3,pt.size = 0)


mito.featuresbsa_f1 <- grep(pattern="^mt-", x=rownames(x=f1), value=T)
percent.mitobsa_f1 <- Matrix::colSums(x = GetAssayData(object = f1, slot = "counts")[mito.featuresbsa_f1,]) / Matrix::colSums(x = GetAssayData(object = f1, slot = "counts"))
f1[["percent.mito"]] <- percent.mitobsa_f1
VlnPlot(object = f1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3,pt.size = 0)

f1 <- subset(x = f1, subset = nFeature_RNA > 1000 & nFeature_RNA <5000 & nCount_RNA > 500 & nCount_RNA < 20000 & percent.mitobsa_f1 < 0.1)
VlnPlot(object = f1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
f1 <- NormalizeData(object = f1, verbose = F)
f1 <- FindVariableFeatures(object = f1, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesbsa_f1 <- rownames(f1)
f1 <- ScaleData(f1, features = all.genesbsa_f1)
VlnPlot(object = f1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

mito.featuresbsa_f2 <- grep(pattern="^mt-", x=rownames(x=f2), value=T)
percent.mitobsa_f2 <- Matrix::colSums(x = GetAssayData(object = f2, slot = "counts")[mito.featuresbsa_f2,]) / Matrix::colSums(x = GetAssayData(object = f2, slot = "counts"))
f2[["percent.mito"]] <- percent.mitobsa_f2
VlnPlot(object = f2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

f2 <- subset(x = f2, subset = nFeature_RNA > 1000 & nFeature_RNA <5000 & nCount_RNA > 500 & nCount_RNA < 20000 & percent.mitobsa_f2 < 0.1)
VlnPlot(object = f2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
f2 <- NormalizeData(object = f2, verbose = F)
f2 <- FindVariableFeatures(object = f2, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesbsa_f2 <- rownames(f2)
f2 <- ScaleData(f2, features = all.genesbsa_f2)
VlnPlot(object = f2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

mito.featuresbsa_f3 <- grep(pattern="^mt-", x=rownames(x=f3), value=T)
percent.mitobsa_f3 <- Matrix::colSums(x = GetAssayData(object = f3, slot = "counts")[mito.featuresbsa_f3,]) / Matrix::colSums(x = GetAssayData(object = f3, slot = "counts"))
f3[["percent.mito"]] <- percent.mitobsa_f3
VlnPlot(object = f3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

f3 <- subset(x = f3, subset = nFeature_RNA > 1000 & nFeature_RNA <5000 & nCount_RNA > 500 & nCount_RNA < 25000 & percent.mitobsa_f3 < 0.1)
VlnPlot(object = f3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
f3 <- NormalizeData(object = f3, verbose = F)
f3 <- FindVariableFeatures(object = f3, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesbsa_f3 <- rownames(f3)
f3 <- ScaleData(f3, features = all.genesbsa_f3)
VlnPlot(object = f3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)


mito.featuresbsa_h1 <- grep(pattern="^mt-", x=rownames(x=h1), value=T)
percent.mitobsa_h1 <- Matrix::colSums(x = GetAssayData(object = h1, slot = "counts")[mito.featuresbsa_h1,]) / Matrix::colSums(x = GetAssayData(object = h1, slot = "counts"))
h1[["percent.mito"]] <- percent.mitobsa_h1
VlnPlot(object = h1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

h1 <- subset(x = h1, subset = nFeature_RNA > 1000 & nFeature_RNA <6000 & nCount_RNA > 500 & nCount_RNA < 25000 & percent.mitobsa_h1 < 0.1)
VlnPlot(object = h1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
h1 <- NormalizeData(object = h1, verbose = F)
h1 <- FindVariableFeatures(object = h1, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesbsa_h1 <- rownames(h1)
h1 <- ScaleData(h1, features = all.genesbsa_h1)
VlnPlot(object = h1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0)

mito.featuresbsa_h2 <- grep(pattern="^mt-", x=rownames(x=h2), value=T)
percent.mitobsa_h2 <- Matrix::colSums(x = GetAssayData(object = h2, slot = "counts")[mito.featuresbsa_h2,]) / Matrix::colSums(x = GetAssayData(object = h2, slot = "counts"))
h2[["percent.mito"]] <- percent.mitobsa_h2
VlnPlot(object = h2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0)

h2 <- subset(x = h2, subset = nFeature_RNA > 1000 & nFeature_RNA <5000 & nCount_RNA > 500 & nCount_RNA < 20000 & percent.mitobsa_h2 < 0.1)
VlnPlot(object = h2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3,pt.size = 0)
h2 <- NormalizeData(object = h2, verbose = F)
h2 <- FindVariableFeatures(object = h2, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesbsa_h2 <- rownames(h2)
h2 <- ScaleData(h2, features = all.genesbsa_h2)
VlnPlot(object = h2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0)

mito.featuresbsa_h3 <- grep(pattern="^mt-", x=rownames(x=h3), value=T)
percent.mitobsa_h3 <- Matrix::colSums(x = GetAssayData(object = h3, slot = "counts")[mito.featuresbsa_h3,]) / Matrix::colSums(x = GetAssayData(object = h3, slot = "counts"))
h3[["percent.mito"]] <- percent.mitobsa_h3
VlnPlot(object = h3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3,pt.size = 0)

h3 <- subset(x = h3, subset = nFeature_RNA > 1000 & nFeature_RNA <5000 & nCount_RNA > 500 & nCount_RNA < 20000 & percent.mitobsa_h3 < 0.1)
VlnPlot(object = h3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0)
h3 <- NormalizeData(object = h3, verbose = F)
h3 <- FindVariableFeatures(object = h3, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesbsa_h3 <- rownames(h3)
h3 <- ScaleData(h3, features = all.genesbsa_h3)
VlnPlot(object = h3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3,pt.size = 0)


total.genes_int <- list(rownames(p1),
                        rownames(p2),
                        rownames(p3),
                        rownames(f1),
                        rownames(f2),
                        rownames(f3),
                        rownames(h1),
                        rownames(h2),
                        rownames(h3))

#find common genes (genes appearing in only one object not used)
common.genes_int <- Reduce(f = intersect, x = total.genes_int)
#check numbers
length(common.genes_int)
nrow(p1)
nrow(p2)
nrow(p3)
nrow(f1)
nrow(f2)
nrow(f3)
nrow(h1)
nrow(h2)
nrow(h3)


#integrate all
#rm(int_all.res)
reference.list_all <- c(p1,p2,p3,f1,f2,f3,h1,h2,h3)
anchors_all <- FindIntegrationAnchors(object.list = reference.list_all, anchor.features = common.genes_int, dims = 1:30)
int_all <- IntegrateData(anchorset = anchors_all, dims = 1:30)

Idents(p1)<-"sample_id"
p1<-subset(p1, downsample = 5000)

Idents(p2)<-"sample_id"
p2<-subset(p2, downsample = 5000)

Idents(p3)<-"sample_id"
p3<-subset(p3, downsample = 5000)

Idents(f1)<-"sample_id"
f1<-subset(f1, downsample = 5000)

Idents(f2)<-"sample_id"
f2<-subset(f2, downsample = 5000)

Idents(f3)<-"sample_id"
f3<-subset(f3, downsample = 5000)

Idents(h1)<-"sample_id"
h1<-subset(h1, downsample = 5000)

Idents(h2)<-"sample_id"
h2<-subset(h2, downsample = 5000)

Idents(h3)<-"sample_id"
h3<-subset(h3, downsample = 5000)


saveRDS(int_all, file = "int_all.rds")

Idents(int_all)<-"sample_id"
int_all<-subset(int_all, downsample = 5000)

DefaultAssay(object=int_all) <- "integrated"
int_all <- ScaleData(object = int_all, verbose=F)
int_all <- RunPCA(object = int_all, verbose=F, npcs = 50)
ElbowPlot(object = int_all, ndims = 50)
#int_all <- SCTransform(int_all)
#DefaultAssay(object=int_all) <- "integrated"
int_all <- FindNeighbors(object = int_all, dims = 1:30)

int_all <- RunUMAP(object = int_all, reduction = "pca", dims = 1:30)
DimPlot(int_all, reduction = "umap", label = T, pt.size = 0.01, repel = F, group.by = "sample_id", raster = F) + NoAxes()


DefaultAssay(object=int_all) <- "integrated"
int_all <- FindVariableFeatures(int_all, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(int_all)
int_all <- ScaleData(int_all, features = all.genes)
DefaultAssay(object=int_all) <- "integrated"
int_all <- RunPCA(int_all, npcs = 50, features = VariableFeatures(object = int_all))
ElbowPlot(int_all,  ndims = 50)
#int_all.H <- RunHarmony(int_all, "orig.ident", reduction = "pca", assay.use = "integrated", plot_convergence = TRUE, max.iter.harmony = 30)

all.cells <- FindNeighbors(all.cells, reduction = "pca", dims = 1:30)
DefaultAssay(object=all.cells) <- "integrated"
all.cells <- FindClusters(all.cells, graph.name = "integrated_snn", resolution = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7))
clustree(all.cells, assay = "integrated")
all.cells <- RunUMAP(all.cells, dims = 1:30)
all.cells <- FindClusters(all.cells, graph.name = "integrated_snn", resolution = 0.1)

DimPlot(all.cells, label = F, pt.size = 0.01, repel =T, label.box = F, raster = F, cols = "Set2") +NoAxes()
DefaultAssay(object=all.cells) <- "RNA"
FeaturePlot(all.cells ,features = "Mpz", raster = F)
VlnPlot(all.cells, features = "Lyve1", pt.size=0)


###remove contaminants


new.cluster.ids <- c("Fibroblasts", "Vascular", "Fibroblasts", "Fibroblasts", "Mural", "Neural")
names(new.cluster.ids) <- levels(all.cells)
all.cells <- RenameIdents(all.cells, new.cluster.ids)

all.cells <-subset(int_all, idents = c("Fibroblasts", "Vascular", "Mural", "Neural"))


####fibroblasts only
fibs <-subset(all.cells, idents = c("Fibroblasts"))
table(fibs$sample_id)
Idents(fibs.1)<-"clusters"

DefaultAssay(object=fibs.1) <- "integrated"
fibs.1 <- FindVariableFeatures(fibs.1, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(fibs.1)
fibs.1 <- ScaleData(fibs.1, features = all.genes)
fibs.1 <- FindClusters(fibs.1, graph.name = "integrated_snn", resolution = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7))
clustree(fibs.1, assay = "integrated")
fibs.1 <- RunUMAP(fibs.1, reduction = "pca", dims = 1:30)
fibs.1 <- FindClusters(fibs.1, graph.name = "integrated_snn", resolution = 0.3)

###Figure 6h
DimPlot(fibs.1, label = F, pt.size = 0.01, repel =T, label.box = F, raster = F) +NoAxes()+scale_color_manual(values = c("burlywood3", "gray20", "violetred4","olivedrab", "darkslategray3", "mistyrose2"))

fibs.1$Treatment <- factor(fibs.1$Treatment,levels=c("Naive", "Control", "CAR"))

DefaultAssay(object=fibs.1) <- "RNA"
comp <- list(c("Naive", "Control"), c("Naive", "CAR"))
(fibs.1, features = "Fap", pt.size=0.1, group.by = "Treatment")+NoLegend()
FeaturePlot(fibs.1, features = "Dpt")

cluster0 <- FindMarkers(fibs.1, ident.1 = 0, min.pct = 0.2, only.pos = T)
cluster1 <- FindMarkers(fibs.1, ident.1 = 1, min.pct = 0.2, only.pos = T)
cluster2 <- FindMarkers(fibs.1, ident.1 = 2, min.pct = 0.2, only.pos = T)
cluster3 <- FindMarkers(fibs.1, ident.1 = 3, min.pct = 0.2, only.pos = T)
cluster4 <- FindMarkers(fibs.1, ident.1 = 4, min.pct = 0.2, only.pos = T)
cluster5 <- FindMarkers(fibs.1, ident.1 = 5, min.pct = 0.2, only.pos = T)
cluster6 <- FindMarkers(fibs.1, ident.1 = 6, min.pct = 0.2, only.pos = T)


new.cluster.ids <- c("Klf4", "Pappa", "Lama2", "Lama2", "Cxcl10", "Gas6", "Pi16")
names(new.cluster.ids) <- levels(fibs.1)
fibs.1 <- RenameIdents(fibs.1, new.cluster.ids)

current.sample.ids <- c("Klf4", "Pappa", "Lama2", "Cxcl10", "Gas6", "Pi16")
new.sample.ids <- c("Bmper", "Rspo3", "Mgp", "Cxcl10", "Gas6", "Pi16")
fibs.1$clusters <- fibs.1@active.ident
fibs.1@meta.data[["clusters"]] <- plyr::mapvalues(x = fibs.1@meta.data[["clusters"]], from = current.sample.ids, to = new.sample.ids)
Idents(fibs.1)<-"clusters"


saveRDS(fibs.1, file = "fibs.rds")

#####fibroblast dotplot
fib.genes<-c("Bmper", "Cdh11", "Rspo3", "Pappa", "Mgp", "Bgn", "Cxcl10", "Cxcl9", "Gas6", "Sfrp4", "Pi16", "Creb5")

dotplot <- DotPlot(fibs.1, features = fib.genes, dot.scale = 5)+ylab(NULL) +xlab(NULL)+ scale_color_gradientn(colours = magma(20))+ FontSize(x.text = 14, y.text = 14) 
dotplot

#####Figure 6h
Idents(fib.pt)<-"clusters"

pt <- table(fibs.1$clusters, fibs.1$Treatment)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)
pt$Var1 <- factor(pt$Var1,levels=c("Bmper", "Mgp", "Rspo3", "Cxcl10", "Gas6", "Pi16"))

ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab(NULL) +
  ylab("Proportion") +
  theme(legend.title = element_blank())+RotatedAxis() +
  scale_fill_manual(values = c("burlywood3", "gray20", "violetred4","olivedrab", "darkslategray3", "mistyrose2"))+NoLegend()


test <- sc_utils(fibs.1)

######Extended data figure 10f
###resting v non-treated
prop.test <- permutation_test(test, cluster_identity = "clusters", sample_1="Naive", sample_2="Control", sample_identity="Treatment", n_permutations=10000)
permutation_plot(prop.test, FDR_threshold = 0.05, log2FD_threshold = 0.5, order_clusters = T)+NoLegend()+xlab(NULL)

######Figure 6i
###car v non-treated
prop.test <- permutation_test(test, cluster_identity = "clusters", sample_1="CAR", sample_2="Control", sample_identity="Treatment", n_permutations=10000)
permutation_plot(prop.test, FDR_threshold = 0.05, log2FD_threshold = 0.5, order_clusters = T)+NoLegend()+xlab(NULL)



#fibs1<-subset(fibs.H, idents = "Fib")



##############Extended Data Figure 10g
library(gsfisher)
seurat_obj <- fibs.1
#Edit to suitable output folder
getwd() 
GSout <- "/rds/projects/c/croftap-sgcar/fibs" 

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
expressed_genes <- getExpressedGenesFromSeuratObjectRNA(seurat_obj,levels(seurat_obj@active.ident), min.pct=0.2)


FindAllMarkers_cells <- FindAllMarkers(seurat_obj, only.pos = T, min.pct = 0.2, logfc.threshold = 0.5) 
#This should be your features, probably easiest to run FindAllMarkers then use the table out of that here


FindAllMarkers_cells$entrez_id <- as.character(annotation$entrez_id[
  match(FindAllMarkers_cells$gene, annotation$gene_name)])
FindAllMarkers_cells <- FindAllMarkers_cells[!is.na(FindAllMarkers_cells$entrez_id),]
background_entrez <- as.character(annotation$entrez_id[
  match(expressed_genes, annotation$gene_name)])
background_entrez <- background_entrez[!is.na(background_entrez)]


#Create output dir
if (dir.exists(GSout)){
  print("Directory already exists")
} else{
  dir.create(GSout)
}

###run go only
#go.results <- runGO.all(results=All_marker,
#background_ids = background_entrez, gene_id_col="entrez_id", gene_id_type="entrez", sample_col="cluster", p_col="p_val_adj", p_threshold=0.05,
#species = "mm")
#go.results <- filterGenesets(go.results)
#go.results.top <- go.results %>% group_by(cluster) %>% top_n(n=5, -p.val)
#sampleEnrichmentDotplot(go.results.top, selection_col = "description", selected_genesets = unique(go.results.top$description), sample_id_col = "cluster", fill_var = "odds.ratio", maxl=50, rotate_sample_labels = TRUE)


#Run enrichment tests and filter results
for(cluster in levels(seurat_obj@active.ident)){
  cat("Running go enrichment for cluster", cluster, "...")
  Go <- runGO(foreground_ids = FindAllMarkers_cells$entrez_id[FindAllMarkers_cells$cluster == cluster],
              background_ids = background_entrez,
              gene_id_type = "entrez", species = "mm")
  
  cat("Running Kegg enrichment for cluster", cluster, "...")
  Kegg <- runKEGG(foreground_ids = FindAllMarkers_cells$entrez_id[FindAllMarkers_cells$cluster == cluster],
                  background_ids = background_entrez,
                  gene_id_type = "entrez", species = "mm")
  
  cat("Joining results tables for cluster", cluster, "...")
  results <- full_join(Go, Kegg)
  
  results$p.adj <- p.adjust(results$p.val, method = "BH")
  
  cat("Filtering genesets for ", cluster, "...")
  filtered_genesets <- filterGenesets(results,
                                      min_foreground_genes = 3,
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

sampleEnrichmentDotplot(All_enrichments, selection_col = "description", selected_genesets = c("negative regulation of cell migration", "negative regulation of cell motility", "fat cell differentiation", "chemokine receptor binding", "chemokine activity", "cytokine activity", "extracellular matrix binding", "extracellular matrix organization", "extracellular matrix structural constituent", "cell chemotaxis", "cell adhesion molecule binding", "regulation of chemotaxis", "BMP signaling pathway", "response to BMP", "regulation of BMP signaling pathway", "regulation of canonical Wnt signaling pathway", "regulation of Wnt signaling", "tissue migration"), sample_id_col = "cluster", fill_var = "odds.ratio", maxl=50, rotate_sample_labels = TRUE,fill_colors = c("yellow3", "red", "black")) + FontSize(x.text = 12, y.text = 12) + ylab(NULL) +xlab(NULL)
sampleEnrichmentHeatmap(All_enrichments, sample_id_col = "cluster", max_rows = 30, min_odds_ratio = 1.5)




                        
########pi16 and immunofibroblast gene modules plotted by density

############Figure 6j                        
Imm.genes <- c("Icam1", "Vcam1", "Pdgfra", "Pdgfrb", "Cxcl13", "Ccl19", "Cxcl10", "Cxcl9", "Cd82", "Tnfsf13b")

fibs.1 <- AddModuleScore(
  object = fibs.1,
  features = list(Imm.genes),
  ctrl = 5,
  name = 'Imm'
)

plot_density(fibs.1, features = c("Imm1"), combine = F, method = "wkde", size = 1, adjust = 1)+NoAxes()#+ scale_color_gradientn(colours = magma(20))#+facet_grid(.~fibs.1$Treatment)


############Figure 6k
pi16.genes <- c("Anxa3", "Ly6c1",	"Fn1",	"Ly6a",	"Dpp4",	"Cd248", 	"Igfbp6",	"Timp2",	"Sema3c",	"Anxa1",	"Pi16", "Emilin2",	"Efhd1",	"Metrnl",	"Pcolce2",	"Efemp1",	"Smpd3",	"Ugdh",	"Cadm3",
                "Sfrp4",	"bn1",	"Ackr3",	"Lrrn4cl",	"Pla1a",	"Fndc1",	"Il1r2",	"Mustn1",	"Plac8",	"Tmem100",	"Basp1", "Marcks",	"Scara3",	"Igfbp5",	"Scara5",	"Akr1c18",	"Clec3b",	"Adgrd1",
                "Fstl1",	"Osr",	"Axl",	"Procr",	"Opcml",	"Tgfbr2",	"Il18",	"Prrx1",	"Creb5",	"Serpinb6a",	"Mfap5",	"Car8",	"Cmah",	"Has1",	"Lurap1l",	"Prss23",	"Gfpt2",	"Tmem158",	"Dpt",	"Errfi1",	"Col14a1",	"S100a13", "Gas7")
                        
fibs.1 <- AddModuleScore(
  object = fibs.1,
  features = list(pi16.genes),
  ctrl = 5,
  name = 'univ'
)
plot_density(fibs.1, features = c("univ1"), combine = F, method = "wkde", size = 1, adjust = 1)+NoAxes()#+ scale_color_gradientn(colours = magma(20))#+facet_grid(.~fibs.1$Treatment)
plot_density(fibs.1, features = c("Pi16"), combine = F, method = "wkde", size = 1, adjust = 1)+NoAxes()+ scale_color_gradientn(colours = viridis(20))#+facet_grid(.~resting.fb$InflammationStatus)

########Extended data figure 10h
VlnPlot(fibs.1, features = c("Icam1"), pt.size = 0, cols = c("burlywood3", "violetred4", "gray20", "olivedrab", "darkslategray3", "mistyrose2"))+xlab(NULL)+NoLegend()
VlnPlot(fibs.1, features = c("Vcam1"), pt.size = 0, cols = c("burlywood3", "violetred4", "gray20", "olivedrab", "darkslategray3", "mistyrose2"))+xlab(NULL)+NoLegend()
VlnPlot(fibs.1, features = c("Cxcl13"), pt.size = 0, cols = c("burlywood3", "violetred4", "gray20", "olivedrab", "darkslategray3", "mistyrose2"))+xlab(NULL)+NoLegend()
VlnPlot(fibs.1, features = c("Ccl19"), pt.size = 0, cols = c("burlywood3", "violetred4", "gray20", "olivedrab", "darkslategray3", "mistyrose2"))+xlab(NULL)+NoLegend()
VlnPlot(fibs.1, features = c("Cxcl10"), pt.size = 0, cols = c("burlywood3", "violetred4", "gray20", "olivedrab", "darkslategray3", "mistyrose2"))+xlab(NULL)+NoLegend()
VlnPlot(fibs.1, features = c("Tnfsf13b"), pt.size = 0, cols = c("burlywood3", "violetred4", "gray20", "olivedrab", "darkslategray3", "mistyrose2"))+xlab(NULL)+NoLegend()


############Extended data figure 10i
###pseudotime
Idents(fibs.1)<-"Treatment"
Naive.fibs<-subset(fibs.1, idents="Naive")
control.fibs<-subset(fibs.1, idents="Control")
CAR.fibs<-subset(fibs.1, idents="CAR")

FindAllMarkers_control <- FindAllMarkers(control.fibs, only.pos = F, min.pct = 0.2, logfc.threshold = 0.5) 
FindAllMarkers_car <- FindAllMarkers(CAR.fibs, only.pos = F, min.pct = 0.2, logfc.threshold = 0.5) 
FindAllMarkers_naive<-FindAllMarkers(Naive.fibs, only.pos = F, min.pct = 0.2, logfc.threshold = 0.5)


all.markers_top50.control<- FindAllMarkers_control %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
all.markers_top50.car<- FindAllMarkers_car %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
all.markers_top50.naive<- FindAllMarkers_naive %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)


######resting analysis####rna is my seurat obj###
data <- GetAssayData(Naive.fibs[["RNA"]], slot="data")
pd <- new('AnnotatedDataFrame', data = Naive.fibs@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
data <- data[, rownames(pd)]

naive.monocle2_cds <- newCellDataSet(data,
                                       phenoData = pd,
                                       featureData = fd,
                                       lowerDetectionLimit = 0.5,
                                       expressionFamily = negbinomial.size())

naive.monocle2_cds <- estimateSizeFactors(naive.monocle2_cds)



#my FindAllMarkers df is 'all.markers'
#my meta.data slot is "res_0.15_cmcleaned", you need to chnage this to your one
naive.monocle2_cds <- setOrderingFilter(naive.monocle2_cds, all.markers_top50.naive$gene)
naive.monocle2_cds <- reduceDimension(naive.monocle2_cds, max_components = 2, method = 'DDRTree')
naive.monocle2_cds <- orderCells(naive.monocle2_cds)
naive.monocle2_cds_plotclus <- plot_cell_trajectory(naive.monocle2_cds, color_by = "clusters", cell_size = 1,show_branch_points = F)
naive.monocle2_cds_plotPseu <- plot_cell_trajectory(naive.monocle2_cds, color_by = "Pseudotime", cell_size = 1,show_branch_points = F)

naive.monocle2_cds_plotclus <- plot_cell_trajectory(naive.monocle2_cds, color_by = "clusters", show_branch_points = F) +facet_wrap(~clusters, nrow = 1)


###plot results###
naive.monocle2_cds_plotclus + NoAxes()+scale_color_manual(values = c("burlywood3", "gray20", "violetred4","olivedrab", "darkslategray3", "mistyrose2"))+NoLegend()  ##this shows cluster ID
naive.monocle2_cds_plotPseu  + NoAxes() + NoLegend() ##this show pseudotime

                        
imm.genes<-c("Cxcl10", "Cxcl9", "Cxcl13", "Icam1", "Tnfsf13b")
naive.monocle2_cds <- naive.monocle2_cds["imm.genes,"]
plot_genes_in_pseudotime(naive.monocle2_cds, color_by = "clusters")+scale_color_manual(values = c("burlywood3", "gray20", "violetred4","olivedrab", "darkslategray3", "mistyrose2"))


####non-treated####rna is my seurat obj###
data <- GetAssayData(control.fibs[["RNA"]], slot="data")
pd <- new('AnnotatedDataFrame', data = control.fibs@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
data <- data[, rownames(pd)]

control.monocle2_cds <- newCellDataSet(data,
                                     phenoData = pd,
                                     featureData = fd,
                                     lowerDetectionLimit = 0.5,
                                     expressionFamily = negbinomial.size())

control.monocle2_cds <- estimateSizeFactors(control.monocle2_cds)



#my FindAllMarkers df is 'all.markers'
#my meta.data slot is "res_0.15_cmcleaned", you need to chnage this to your one
control.monocle2_cds <- setOrderingFilter(control.monocle2_cds, all.markers_top50.naive$gene)
control.monocle2_cds <- reduceDimension(control.monocle2_cds, max_components = 2, method = 'DDRTree')
control.monocle2_cds <- orderCells(control.monocle2_cds)
control.monocle2_cds_plotclus <- plot_cell_trajectory(control.monocle2_cds, color_by = "clusters", cell_size = 1,show_branch_points = F)
control.monocle2_cds_plotPseu <- plot_cell_trajectory(control.monocle2_cds, color_by = "Pseudotime", cell_size = 1,show_branch_points = F)

control.monocle2_cds_plotclus <- plot_cell_trajectory(control.monocle2_cds, color_by = "clusters", show_branch_points = F) +facet_wrap(~clusters, nrow = 1)


###plot results###
control.monocle2_cds_plotclus + NoAxes()+scale_color_manual(values = c("burlywood3", "gray20", "violetred4","olivedrab", "darkslategray3", "mistyrose2"))+NoLegend()  ##this shows cluster ID
control.monocle2_cds_plotPseu  + NoAxes() + NoLegend() ##this show pseudotime

imm.genes<-c("Cxcl10", "Cxcl9", "Cxcl13", "Icam1", "Tnfsf13b")
control.monocle2_cds <- control.monocle2_cds["imm.genes,"]
plot_genes_in_pseudotime(control.monocle2_cds, color_by = "clusters")+scale_color_manual(values = c("burlywood3", "gray20", "violetred4","olivedrab", "darkslategray3", "mistyrose2"))


######CAR analysis####rna is my seurat obj###
data <- GetAssayData(CAR.fibs[["RNA"]], slot="data")
pd <- new('AnnotatedDataFrame', data = CAR.fibs@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
data <- data[, rownames(pd)]

CAR.monocle2_cds <- newCellDataSet(data,
                                       phenoData = pd,
                                       featureData = fd,
                                       lowerDetectionLimit = 0.5,
                                       expressionFamily = negbinomial.size())

CAR.monocle2_cds <- estimateSizeFactors(CAR.monocle2_cds)



#my FindAllMarkers df is 'all.markers'
#my meta.data slot is "res_0.15_cmcleaned", you need to chnage this to your one
CAR.monocle2_cds <- setOrderingFilter(CAR.monocle2_cds, all.markers_top50.naive$gene)
CAR.monocle2_cds <- reduceDimension(CAR.monocle2_cds, max_components = 2, method = 'DDRTree')
CAR.monocle2_cds <- orderCells(CAR.monocle2_cds)
CAR.monocle2_cds_plotclus <- plot_cell_trajectory(CAR.monocle2_cds, color_by = "clusters", cell_size = 1,show_branch_points = F)
CAR.monocle2_cds_plotPseu <- plot_cell_trajectory(CAR.monocle2_cds, color_by = "Pseudotime", cell_size = 1,show_branch_points = F)

CAR.monocle2_cds_plotclus <- plot_cell_trajectory(CAR.monocle2_cds, color_by = "clusters", show_branch_points = F) +facet_wrap(~clusters, nrow = 1)


###plot results###
CAR.monocle2_cds_plotclus + NoAxes()+scale_color_manual(values = c("burlywood3", "gray20", "violetred4","olivedrab", "darkslategray3", "mistyrose2"))+NoLegend()  ##this shows cluster ID
CAR.monocle2_cds_plotPseu  + NoAxes() + NoLegend() ##this show pseudotime

imm.genes<-c("Cxcl10", "Cxcl9", "Cxcl13", "Icam1", "Tnfsf13b")
CAR.monocle2_cds <- CAR.monocle2_cds["imm.genes,"]
plot_genes_in_pseudotime(CAR.monocle2_cds, color_by = "clusters")+scale_color_manual(values = c("burlywood3", "gray20", "violetred4","olivedrab", "darkslategray3", "mistyrose2"))


#############creb5 boxplot
gene <- "Creb5"   # your gene/feature of interest

expr <- FetchData(fib.pt, vars = gene)
expr$Treatment <- fib.pt$Treatment   # add metadata column
colnames(expr)[1] <- "expression"



library(ggplot2)
library(ggpubr)

ggplot(expr, aes(x = Treatment, y = expression, fill = Treatment)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_classic() +
  ylab(paste(gene, "expression")) +
  xlab("Treatment") + theme(text = element_text(size = 14),
                            axis.text.x = element_text(size = 14),
                            axis.text.y = element_text(size = 14)) + scale_fill_manual(values = c("CAR" = "firebrick", 
                                                                                                  "Control" = "grey"))+stat_compare_means(
                                                                                                    method = "wilcox.test",
                                                                                                    label = "p.signif",
                                                                                                    paired = F,
                                                                                                    ref.group = "Control") +NoLegend()

########Extended data figure 10d
                        
Idents(control.fibs)<-"clusters"
fap.fib<-subset(x = control.fibs, subset = Fap > 1)
fap.lo<-subset(x = control.fibs, subset = Fap < 1)

fap.fib$Fap <- "High"
fap.lo$Fap <- "Low"

Idents(fap.fib)<-"Fap"
Idents(fap.lo)<-"Fap"


fap_all <- merge(fap.fib, fap.lo)
F_DEGs <- FindAllMarkers(fap_all, only.pos = T, min.pct = 0.1, logfc.threshold = 0.05)

library(AnnotationDbi)
library(org.Mm.eg.db)
library(clusterProfiler)


degs_hi <- rownames(subset(F_DEGs, p_val_adj < 0.05 & abs(avg_log2FC) > 0.1 & cluster == "High"))
degs_low <- rownames(subset(F_DEGs, p_val_adj < 0.05 & abs(avg_log2FC) > 0.1 & cluster == "Low"))

degs_low <- gsub("\\..*", "", degs_low)

ego.lo <- enrichGO(gene = degs_low,
                   OrgDb = org.Mm.eg.db,
                   keyType = "SYMBOL",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.25,
                   readable = TRUE)
dotplot(ego.lo, showCategory=c("negative regulation of catalytic activity", "negative regulation of peptidase activity", "response to interleukin-4", "cellular response to interleukin-4", "negative regulation of proteolysis", "membraneless organelle assembly", "translation at presynapse", "translation at synapse"))

"mesenchymal cell differentiation", "regulation of canonical Wnt signaling pathway", "sprouting angiogenesis", "canonical Wnt signaling pathway", "vascular endothelial growth factor signaling pathway", "mesenchymal cell proliferation", "cell migration involved in sprouting angiogenesis", "response to fibroblast growth factor", "learning or memory", "tissue migration", "Wnt signaling pathway", "blood vessel endothelial cell migration", "epithelial cell migration", "regulation of endothelial cell migration", "fibroblast growth factor receptor signaling pathway", "regulation of epithelial cell proliferation"

"regulation of response to wounding", "negative regulation of canonical NF-kappaB signal transduction", "myeloid cell homeostasis", "negative regulation of endopeptidase activity", "negative regulation of immune effector process", "response to type II interferon", "negative regulation of interleukin-8 production", "negative regulation of peptidase activity", "response to interleukin-4", "cellular response to interleukin-4"



hi <- head(degs_hi, 20)
lo <- head(degs_low, 20)


fibs.1 <- AddModuleScore(
  object = fibs.1,
  features = list(hi),
  ctrl = 5,
  name = 'hi')

fibs.1 <- AddModuleScore(
  object = fibs.1,
  features = list(lo),
  ctrl = 5,
  name = 'low')

DotPlot(
  fibs.1,
  features = c("hi1", "low1"),
  assay = NULL,
  dot.scale = 6
) +
  scale_color_gradientn(
    colors = c("skyblue", "lightgrey", "red"),
    values = scales::rescale(c(-1.5, 0, 1.5))
  )+xlab(NULL)+ylab(NULL)
