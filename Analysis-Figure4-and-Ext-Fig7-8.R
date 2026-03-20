############lines 0-292 sample intergation and processing. line 293 onwards analysis for figures

library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(reticulate)
library(RColorBrewer)
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
library(scater)
library(patchwork)
library(SingleCellExperiment)
library(monocle)
library(DDRTree)
library(SeuratWrappers)
library(devtools)
library(ggthemes)
library(tvthemes)




############data integration

#########non-treated
p1 <- Read10X(data.dir = "/rds/projects/c/croftap-cartcell/count/P1/outs/filtered_feature_bc_matrix")
p1 <- CreateSeuratObject(counts = p1, min.cells=3, min.features=100, project="p1")
p1 <- RenameCells(p1, add.cell.id = "p1")
p1$sample_id <- "PBS1"
head(p1)

p2 <- Read10X(data.dir = "/rds/projects/c/croftap-cartcell/count/P2/outs/filtered_feature_bc_matrix")
p2 <- CreateSeuratObject(counts = p2, min.cells=3, min.features=100, project="p2")
p2 <- RenameCells(p2, add.cell.id = "p2")
p2$sample_id <- "PBS2"
head(p2)

p3 <- Read10X(data.dir = "/rds/projects/c/croftap-cartcell/count/P3/outs/filtered_feature_bc_matrix")
p3 <- CreateSeuratObject(counts = p3, min.cells=3, min.features=100, project="p3")
p3 <- RenameCells(p3, add.cell.id = "p3")
p3$sample_id <- "PBS3"
head(p3)

#############CAR-T treated
f1 <- Read10X(data.dir = "/rds/projects/c/croftap-cartcell/count/F1/outs/filtered_feature_bc_matrix")
f1 <- CreateSeuratObject(counts = f1, min.cells=3, min.features=100, project="f1")
f1 <- RenameCells(f1, add.cell.id = "f1")
f1$sample_id <- "FAP1"
head(f1)

f2 <- Read10X(data.dir = "/rds/projects/c/croftap-cartcell/count/F2/outs/filtered_feature_bc_matrix")
f2 <- CreateSeuratObject(counts = f2, min.cells=3, min.features=100, project="f2")
f2 <- RenameCells(f2, add.cell.id = "f2")
f2$sample_id <- "FAP2"
head(f2)

f3 <- Read10X(data.dir = "/rds/projects/c/croftap-cartcell/count/F3/outs/filtered_feature_bc_matrix")
f3 <- CreateSeuratObject(counts = f3, min.cells=3, min.features=100, project="f3")
f3 <- RenameCells(f3, add.cell.id = "f3")
f3$sample_id <- "FAP3"
head(f3)

###QC
mito.featuresbsa_p1 <- grep(pattern="^mt-", x=rownames(x=p1), value=T)
percent.mitobsa_p1 <- Matrix::colSums(x = GetAssayData(object = p1, slot = "counts")[mito.featuresbsa_p1,]) / Matrix::colSums(x = GetAssayData(object = p1, slot = "counts"))
p1[["percent.mito"]] <- percent.mitobsa_p1
VlnPlot(object = p1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

p1 <- subset(x = p1, subset = nFeature_RNA > 1000 & nFeature_RNA <8000 & nCount_RNA > 500 & nCount_RNA < 60000 & percent.mitobsa_p1 < 0.1)
VlnPlot(object = p1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
p1 <- NormalizeData(object = p1, verbose = F)
p1 <- FindVariableFeatures(object = p1, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesbsa_p1 <- rownames(p1)
p1 <- ScaleData(p1, features = all.genesbsa_p1)
VlnPlot(object = p1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

mito.featuresbsa_p2 <- grep(pattern="^mt-", x=rownames(x=p2), value=T)
percent.mitobsa_p2 <- Matrix::colSums(x = GetAssayData(object = p2, slot = "counts")[mito.featuresbsa_p2,]) / Matrix::colSums(x = GetAssayData(object = p2, slot = "counts"))
p2[["percent.mito"]] <- percent.mitobsa_p2
VlnPlot(object = p2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

p2 <- subset(x = p2, subset = nFeature_RNA > 1000 & nFeature_RNA <8000 & nCount_RNA > 500 & nCount_RNA < 50000 & percent.mitobsa_p2 < 0.1)
VlnPlot(object = p2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
p2 <- NormalizeData(object = p2, verbose = F)
p2 <- FindVariableFeatures(object = p2, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesbsa_p2 <- rownames(p2)
p2 <- ScaleData(p2, features = all.genesbsa_p2)
VlnPlot(object = p2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

mito.featuresbsa_p3 <- grep(pattern="^mt-", x=rownames(x=p3), value=T)
percent.mitobsa_p3 <- Matrix::colSums(x = GetAssayData(object = p3, slot = "counts")[mito.featuresbsa_p3,]) / Matrix::colSums(x = GetAssayData(object = p3, slot = "counts"))
p3[["percent.mito"]] <- percent.mitobsa_p3
VlnPlot(object = p3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

p3 <- subset(x = p3, subset = nFeature_RNA > 1000 & nFeature_RNA <8000 & nCount_RNA > 500 & nCount_RNA < 60000 & percent.mitobsa_p3 < 0.1)
VlnPlot(object = p3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
p3 <- NormalizeData(object = p3, verbose = F)
p3 <- FindVariableFeatures(object = p3, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesbsa_p3 <- rownames(p3)
p3 <- ScaleData(p3, features = all.genesbsa_p3)
VlnPlot(object = p3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

##########subset high quality cells from each sample
mito.featuresbsa_f1 <- grep(pattern="^mt-", x=rownames(x=f1), value=T)
percent.mitobsa_f1 <- Matrix::colSums(x = GetAssayData(object = f1, slot = "counts")[mito.featuresbsa_f1,]) / Matrix::colSums(x = GetAssayData(object = f1, slot = "counts"))
f1[["percent.mito"]] <- percent.mitobsa_f1
VlnPlot(object = f1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

f1 <- subset(x = f1, subset = nFeature_RNA > 1000 & nFeature_RNA <8000 & nCount_RNA > 500 & nCount_RNA < 60000 & percent.mitobsa_f1 < 0.1)
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

f2 <- subset(x = f2, subset = nFeature_RNA > 1000 & nFeature_RNA <7500 & nCount_RNA > 500 & nCount_RNA < 50000 & percent.mitobsa_f2 < 0.1)
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

f3 <- subset(x = f3, subset = nFeature_RNA > 1000 & nFeature_RNA <8000 & nCount_RNA > 500 & nCount_RNA < 60000 & percent.mitobsa_f3 < 0.1)
VlnPlot(object = f3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
f3 <- NormalizeData(object = f3, verbose = F)
f3 <- FindVariableFeatures(object = f3, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesbsa_f3 <- rownames(f3)
f3 <- ScaleData(f3, features = all.genesbsa_f3)
VlnPlot(object = f3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

total.genes_int <- list(rownames(p1),
                        rownames(p2),
                        rownames(p3),
                        rownames(f1),
                        rownames(f2),
                        rownames(f3))

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


#integrate all
#rm(int_all.res)
reference.list_all <- c(p1,p2,p3,f1,f2,f3)
anchors_all <- FindIntegrationAnchors(object.list = reference.list_all, anchor.features = common.genes_int, dims = 1:30)
int_all <- IntegrateData(anchorset = anchors_all, dims = 1:30)


DefaultAssay(object=int_all) <- "integrated"
int_all <- ScaleData(object = int_all, verbose=F)
int_all <- RunPCA(object = int_all, verbose=F, npcs = 50)
ElbowPlot(object = int_all, ndims = 50)
int_all <- SCTransform(int_all)
DefaultAssay(object=int_all) <- "integrated"
int_all <- FindNeighbors(object = int_all, dims = 1:30)

int_all <- RunUMAP(object = int_all, reduction = "pca", dims = 1:30)
DimPlot(int_all, reduction = "umap", label = F, pt.size = 0.01, repel = F, group.by = "sample_id") + NoLegend() + NoAxes()

##add model meta data
Idents(int_all)<-'orig.ident'
int_all$model<-int_all@active.ident


current.sample.ids <- c("FAP1", "FAP2", "FAP3", "PBS1", "PBS2", "PBS3")

new.sample.ids<-c("CART", "CART", "CART", "PBS", "PBS", "PBS")

int_all@meta.data[["model"]] <- plyr::mapvalues(x = int_all@meta.data[["model"]], from = current.sample.ids, to = new.sample.ids)

int_all <- FindVariableFeatures(int_all, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(int_all)
int_all <- ScaleData(int_all, features = all.genes)
int_all <- RunPCA(int_all, npcs = 50, features = VariableFeatures(object = int_all))
ElbowPlot(int_all.res,  ndims = 50)
int_all <- FindNeighbors(int_all, reduction = "pca", dims = 1:30)
# to find out graph.name integrated.H@graphs
int_all <- FindClusters(int_all, graph.name = "integrated_snn", resolution = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7))
clustree(int_all, assay = "integrated")
int_all <- FindClusters(int_all, graph.name = "integrated_snn", resolution = 0.3)
int_all <- RunUMAP(int_all, reduction = "pca", dims = 1:30)
DimPlot(int_all, reduction = "umap", label = T, pt.size = 0.01, repel = T, label.box = F) + NoAxes()
DimPlot(int_all, reduction = "umap", label = F, pt.size = 0.01, repel = T, label.box = F, split.by = "model") + NoLegend() + NoAxes()

cluster0 <- FindMarkers(int_all, ident.1 = 0, min.pct = 0.7, only.pos = T)
cluster1 <- FindMarkers(int_all, ident.1 = 1, min.pct = 0.7, only.pos = T)
cluster2 <- FindMarkers(int_all, ident.1 = 2, min.pct = 0.7, only.pos = T)
cluster3 <- FindMarkers(int_all, ident.1 = 3, min.pct = 0.7, only.pos = T)
cluster4 <- FindMarkers(int_all, ident.1 = 4, min.pct = 0.7, only.pos = T)
cluster5 <- FindMarkers(int_all, ident.1 = 5, min.pct = 0.7, only.pos = T)
cluster6 <- FindMarkers(int_all, ident.1 = 6, min.pct = 0.7, only.pos = T)
cluster7 <- FindMarkers(int_all, ident.1 = 7, min.pct = 0.7, only.pos = T)
cluster8 <- FindMarkers(int_all, ident.1 = 8, min.pct = 0.7, only.pos = T)
cluster9 <- FindMarkers(int_all, ident.1 = 9, min.pct = 0.7, only.pos = T)
cluster10 <- FindMarkers(int_all, ident.1 = 10, min.pct = 0.7, only.pos = T)
cluster11 <- FindMarkers(int_all, ident.1 = 11, min.pct = 0.7, only.pos = T)
cluster12 <- FindMarkers(int_all, ident.1 = 12, min.pct = 0.7, only.pos = T)
cluster13 <- FindMarkers(int_all, ident.1 = 13, min.pct = 0.7, only.pos = T)
cluster14 <- FindMarkers(int_all, ident.1 = 14, min.pct = 0.7, only.pos = T)
cluster15 <- FindMarkers(int_all, ident.1 = 15, min.pct = 0.7, only.pos = T)
cluster16 <- FindMarkers(int_all, ident.1 = 16, min.pct = 0.7, only.pos = T)
cluster17 <- FindMarkers(int_all, ident.1 = 17, min.pct = 0.7, only.pos = T)
cluster18 <- FindMarkers(int_all, ident.1 = 18, min.pct = 0.7, only.pos = T)
cluster19 <- FindMarkers(int_all, ident.1 = 19, min.pct = 0.7, only.pos = T)

DefaultAssay(object=int_all) <- "RNA"
FeaturePlot(int_all, features = "Cd34", split.by = "model")
DimPlot(int_all.res, reduction = "umap", label = F, pt.size = 0.01, repel = T, label.box = T,split.by = "model") + NoAxes()+scale_fill_jco()+scale_color_jco()


new.cluster.ids <- c("Fibroblasts", "Osteoblasts", "Fibroblasts", "Fibroblasts", "Vascular", "Fibroblasts", "Fibroblasts", "Con", "Cycling", "Fibroblasts", "Cycling Fibroblasts"
                     , "Con", "Mural", "Chondrocytes", "Con", "Fibroblasts", "Lymphatics", "Muscle", "Glial", "Con")
names(new.cluster.ids) <- levels(int_all)
int_all <- RenameIdents(int_all, new.cluster.ids)

int_all.res <- subset(int_all, idents = grep(c("Fibroblasts|Osteoblasts|Vascular|Cycling Fibroblasts|Mural|Chondrocytes|Lymphatics"), 
                                                       levels(int_all@active.ident), value = T))

DefaultAssay(object=int_all.res) <- "integrated"
int_all.res <- FindVariableFeatures(int_all.res, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(int_all.res)
int_all.res <- ScaleData(int_all.res, features = all.genes)
int_all.res <- RunPCA(int_all.res, npcs = 50, features = VariableFeatures(object = int_all.res))
ElbowPlot(int_all.res,  ndims = 50)
int_all.res <- FindNeighbors(int_all.res, reduction = "pca", dims = 1:30)
# to find out graph.name integrated.H@graphs
int_all.res <- FindClusters(int_all.res, graph.name = "integrated_snn", resolution = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7))
clustree(int_all.res, assay = "integrated")
int_all.res <- FindClusters(int_all.res, graph.name = "integrated_snn", resolution = 0.3)
int_all.res <- RunUMAP(int_all.res, reduction = "pca", dims = 1:50)
DimPlot(int_all.res, reduction = "umap", label = TRUE, pt.size = 0.1, repel = T, label.box = T) + NoAxes()


DefaultAssay(object=int_all.res) <- "RNA"
FeaturePlot(int_all.res, features = "Acan", split.by = "model")
DimPlot(int_all.res.1, reduction = "umap", label = F, pt.size = 0.01, repel = T, label.box = T) + NoAxes()+scale_fill_npg()+scale_color_npg()


new.cluster.ids <- c("Fibroblasts", "Osteoblasts", "Fibroblasts", "Fibroblasts", "Vascular", "Fibroblasts", "Fibroblasts", "Fibroblasts", "Cycling"
                     , "Mural", "Lymphatics", "Chondrocytes", "Fibroblasts")
names(new.cluster.ids) <- levels(int_all.res)
int_all.res <- RenameIdents(int_all.res, new.cluster.ids)

int_all.res.1 <- subset(int_all.res, idents = grep(c("Fibroblasts"), 
                                             levels(int_all.res@active.ident), value = T))

#########save fibroblasts
saveRDS(int_all.res.1, file = "Fbs.rds")

#############extended data figure 6a
###########fibroblast tissue state pathway analysis
Idents(int_all.res.1)<-"model"
library(gsfisher)
seurat_obj <- int_all.res.1
#Edit to suitable output folder
getwd() 
GSout <- "/rds/projects/c/croftap-cartcell/count/analysis/model" 

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
expressed_genes <- getExpressedGenesFromSeuratObjectRNA(seurat_obj,levels(seurat_obj@active.ident), min.pct=0.25)


FindAllMarkers_pFbs <- FindAllMarkers(seurat_obj, only.pos = T, min.pct = 0.25, logfc.threshold = 0.05) 
#This should be your features, probably easiest to run FindAllMarkers then use the table out of that here


FindAllMarkers_pFbs$entrez_id <- as.character(annotation$entrez_id[
  match(FindAllMarkers_pFbs$gene, annotation$gene_name)])
FindAllMarkers_pFbs <- FindAllMarkers_pFbs[!is.na(FindAllMarkers_pFbs$entrez_id),]
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
  Go <- runGO(foreground_ids = FindAllMarkers_pFbs$entrez_id[FindAllMarkers_pFbs$cluster == cluster],
              background_ids = background_entrez,
              gene_id_type = "entrez", species = "mm")
  
  cat("Running Kegg enrichment for cluster", cluster, "...")
  Kegg <- runKEGG(foreground_ids = FindAllMarkers_pFbs$entrez_id[FindAllMarkers_pFbs$cluster == cluster],
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

sampleEnrichmentDotplot(all_results_top, selection_col = "description", selected_genesets = unique(All_enrichments$description), sample_id_col = "cluster", fill_var = "odds.ratio", maxl=50, rotate_sample_labels = TRUE,fill_colors = c("yellow3", "red", "black")) + FontSize(x.text = 12, y.text = 12) + ylab(NULL) +xlab(NULL)


#########reload fibroblasts - annotation
int_all.res.1 <- readRDS("Fbs.rds")

DefaultAssay(object=int_all.res.1) <- "integrated"
int_all.res.1 <- FindVariableFeatures(int_all.res.1, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(int_all.res.1)
int_all.res.1 <- ScaleData(int_all.res.1, features = all.genes)
int_all.res.1 <- RunPCA(int_all.res.1, npcs = 50, features = VariableFeatures(object = int_all.res.1))
ElbowPlot(int_all.res.1,  ndims = 50)
int_all.res.1 <- FindNeighbors(int_all.res.1, reduction = "pca", dims = 1:30)
# to find out graph.name integrated.H@graphs
int_all.res.1 <- FindClusters(int_all.res.1, graph.name = "integrated_snn", resolution = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7))
clustree(int_all.res.1, assay = "integrated")
int_all.res.1 <- FindClusters(int_all.res.1, graph.name = "integrated_snn", resolution = 0.6)
int_all.res.1 <- RunUMAP(int_all.res.1, reduction = "pca", dims = 1:50)
DimPlot(int_all.res.1, reduction = "umap", label = F, pt.size = 0.01, repel = T, label.box = F)+NoAxes()+scale_fill_igv()+scale_color_igv()

cluster0 <- FindMarkers(int_all.res.1, ident.1 = 0, min.pct = 0.4, only.pos = T)
cluster1 <- FindMarkers(int_all.res.1, ident.1 = 1, min.pct = 0.4, only.pos = T)
cluster2 <- FindMarkers(int_all.res.1, ident.1 = 2, min.pct = 0.4, only.pos = T)
cluster3 <- FindMarkers(int_all.res.1, ident.1 = 3, min.pct = 0.4, only.pos = T)
cluster4 <- FindMarkers(int_all.res.1, ident.1 = 4, min.pct = 0.4, only.pos = T)
cluster5 <- FindMarkers(int_all.res.1, ident.1 = 5, min.pct = 0.4, only.pos = T)
cluster6 <- FindMarkers(int_all.res.1, ident.1 = 6, min.pct = 0.4, only.pos = T)
cluster7 <- FindMarkers(int_all.res.1, ident.1 = 7, min.pct = 0.4, only.pos = T)
cluster8 <- FindMarkers(int_all.res.1, ident.1 = 8, min.pct = 0.4, only.pos = T)
cluster9 <- FindMarkers(int_all.res.1, ident.1 = 8, min.pct = 0.4, only.pos = T)
cluster10 <- FindMarkers(int_all.res.1, ident.1 = 8, min.pct = 0.4, only.pos = T)
cluster11 <- FindMarkers(int_all.res.1, ident.1 = 8, min.pct = 0.4, only.pos = T)

current.sample.ids <- c("0", "1", "2", "3", "4","5", "6", "7", "8", "9", "10", "11")
new.sample.ids <- c("Postn", "Postn", "Prg4", "Mdk", "Mdk", "Fbn1", "Spp1", "Angptl7", "Spp1", "Chodl","Pi16", "con")

FindAllMarkers_FBs <- FindAllMarkers(int_all.res.1, only.pos = TRUE, min.pct = 0.7, logfc.threshold = 0.01) 

int_all.res.1b <- subset(int_all.res.1, idents = grep(c("Postn|Prg4|Mdk|Fbn1|Spp1|Angptl7|Chodl|Pi16"), 
                                                      levels(int_all.res.1@active.ident), value = T))

current.sample.ids <- c("Postn", "Prg4", "Col8a1", "Cd34", "Sox9", "Fmod", "Chodl","Pi16")
new.sample.ids <- c("Postn", "Prg4", "Col8a1", "Cd34", "Sox9", "Fmod", "Col8a1","Pi16")

int_all.res.1b$FB_clusters <- int_all.res.1b@active.ident
int_all.res.1b@meta.data[["FB_clusters"]] <- plyr::mapvalues(x = int_all.res.1b@meta.data[["FB_clusters"]], from = current.sample.ids, to = new.sample.ids)


#######extended data figure 6b
fib.genes<-c("Postn", "Cthrc1", "Cxcl5", "Prg4", "F13a1", "Sox5", "Col8a1", "Chad", "Sfrp2", "Ifi204", "Bmper", "Clec3b", "Spp1", "Sox9", "Mmp3", "Fmod", "Angptl7", "Thbs4", "Pi16", "Cd34", "Cd55")

#plot dotplot
dotplot <- DotPlot(int_all.res.1b, features = fib.genes, dot.scale = 5, cols = "RdGy") + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) +ylab("Cluster") +xlab("Gene")
dotplot
head(dotplot)
top5 <- top5[-c(30), ]

########figure 4a
DimPlot(int_all.res.1b, reduction = "umap", label = F, pt.size = 0.01, repel = T, label.box = F)+NoAxes()+scale_fill_igv()+scale_color_igv()

##########barplot

pt <- table(int_all.res.1b$FB_clusters, int_all.res.1b$model)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)
pt$Var1 <- factor(pt$Var1,levels=c("Postn", "Prg4", "Col8a1", "Cd34", "Sox9", "Fmod", "Col8a1","Pi16"))

ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab(NULL) +
  ylab("Proportion") +
  theme(legend.title = element_blank())+RotatedAxis() +
  scale_fill_igv()+NoLegend()

###########figure 4b
#########scproportiontest
test <- sc_utils(int_all.res.1b)
prop.test <- permutation_test(test, cluster_identity = "FB_clusters", sample_1="CART", sample_2="PBS", sample_identity="Model", n_permutations=10000)
permutation_plot(prop.test, FDR_threshold = 0.01, log2FD_threshold = 0.58, order_clusters = T)+ FontSize(x.text = 12, y.text = 12, x.title = 12, y.title = 12)+xlab(NULL)+NoLegend()

#########extended data figure 6b
###########fibroblast pathway analysis
library(gsfisher)
seurat_obj <- int_all.res.1b
#Edit to suitable output folder
getwd() 
GSout <- "/rds/projects/c/croftap-cartcell/count/analysis" 

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
expressed_genes <- getExpressedGenesFromSeuratObjectRNA(seurat_obj,levels(seurat_obj@active.ident), min.pct=0.25)


FindAllMarkers_pFbs <- FindAllMarkers(seurat_obj, only.pos = T, min.pct = 0.25, logfc.threshold = 0.05) 
#This should be your features, probably easiest to run FindAllMarkers then use the table out of that here


FindAllMarkers_pFbs$entrez_id <- as.character(annotation$entrez_id[
  match(FindAllMarkers_pFbs$gene, annotation$gene_name)])
FindAllMarkers_pFbs <- FindAllMarkers_pFbs[!is.na(FindAllMarkers_pFbs$entrez_id),]
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
  Go <- runGO(foreground_ids = FindAllMarkers_pFbs$entrez_id[FindAllMarkers_pFbs$cluster == cluster],
              background_ids = background_entrez,
              gene_id_type = "entrez", species = "mm")
  
  cat("Running Kegg enrichment for cluster", cluster, "...")
  Kegg <- runKEGG(foreground_ids = FindAllMarkers_pFbs$entrez_id[FindAllMarkers_pFbs$cluster == cluster],
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

sampleEnrichmentDotplot(all_results_top, selection_col = "description", selected_genesets = unique(All_enrichments$description), sample_id_col = "cluster", fill_var = "odds.ratio", maxl=50, rotate_sample_labels = TRUE,fill_colors = c("yellow3", "red", "black")) + FontSize(x.text = 12, y.text = 12) + ylab(NULL) +xlab(NULL)


########Pi16 gene expression and module score plotted as density
pi16.genes <- c("Anxa3", "Ly6c1",	"Fn1",	"Ly6a",	"Dpp4",	"Cd248", 	"Igfbp6",	"Timp2",	"Sema3c",	"Anxa1",	"Pi16", "Emilin2",	"Efhd1",	"Metrnl",	"Pcolce2",	"Efemp1",	"Smpd3",	"Ugdh",	"Cadm3",
                "Sfrp4",	"bn1",	"Ackr3",	"Lrrn4cl",	"Pla1a",	"Fndc1",	"Il1r2",	"Mustn1",	"Plac8",	"Tmem100",	"Basp1", "Marcks",	"Scara3",	"Igfbp5",	"Scara5",	"Akr1c18",	"Clec3b",	"Adgrd1",
                "Fstl1",	"Osr",	"Axl",	"Procr",	"Opcml",	"Tgfbr2",	"Il18",	"Prrx1",	"Creb5",	"Serpinb6a",	"Mfap5",	"Car8",	"Cmah",	"Has1",	"Lurap1l",	"Prss23",	"Gfpt2",	"Tmem158",	"Dpt",	"Errfi1",	"Col14a1",	"S100a13", "Gas7")

int_all.res.1b <- AddModuleScore(
  object = int_all.res.1b,
  features = list(pi16.genes),
  ctrl = 5,
  name = 'univ'
)

#####figure 4c
plot_density(int_all.res.1b, features = c("Pi16"), combine = F, method = "wkde", size = 1, adjust = 1)+NoAxes()+ scale_color_gradientn(colours = viridis(20))#+facet_grid(.~resting.fb$InflammationStatus)
########figure 4d
plot_density(int_all.res.1b, features = c("univ1"), combine = F, method = "wkde", size = 1, adjust = 1)+NoAxes()+ scale_color_gradientn(colours = viridis(20))#+facet_grid(.~resting.fb$InflammationStatus)

##########figure 4e
##########boxplot klf2 klf4
gene <- "Klf2"   # your gene/feature of interest-also run "klf4"

expr <- FetchData(int_all.res.1b, vars = gene)
expr$model <- int_all.res.1b$model   # add metadata column
colnames(expr)[1] <- "expression"

library(ggplot2)
library(ggpubr)

ggplot(expr, aes(x = model, y = expression, fill = model)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_classic() +
  ylab(paste(gene, "expression")) +
  xlab("Treatment") + theme(text = element_text(size = 14),
                            axis.text.x = element_text(size = 14),
                            axis.text.y = element_text(size = 14)) + scale_fill_manual(values = c("PBS" = "grey", 
                                                                                                  "CART" = "firebrick"))+stat_compare_means(
                                                                                                    method = "wilcox.test",
                                                                                                    label = "p.signif",
                                                                                                    ref.group = "PBS") +NoLegend()

##########Pi16 DEGS and volcano plot
DEG_Pi16p <- FindMarkers(int_all.res.1b, ident.1 = "Pi16.PBS", ident.2 = "Pi16.CART", show.ress = F, only.pos =T)
DEG_Pi16c <- FindMarkers(int_all.res.1b, ident.1 = "Pi16.CART", ident.2 = "Pi16.PBS", show.ress = F, only.pos =T)

write.csv(DEG_Pi16p, file = "Pi16p.csv")
write.csv(DEG_Pi16c, file = "Pi16c.csv")
write.csv(DEG_Cd55.PBS, file = "Cd55.pbs.csv")
write.csv(DEG_Cd55.CART, file = "Cd55.cart.csv")

DEG_Pi16_T <-FindMarkers(int_all.res.1b, ident.1 = "Pi16.PBS", ident.2 = "Pi16.CART", show.ress = F, only.pos =F) 

keyvals <- ifelse(
  DEG_Pi16_T$avg_log2FC < -0.5, 'firebrick',
  ifelse(DEG_Pi16_T$avg_log2FC > 0.5, 'grey',
         'skyblue2'))

keyvals[is.na(keyvals)] <- 'skyblue2'
names(keyvals)[keyvals == 'grey'] <- 'Control'
names(keyvals)[keyvals == 'skyblue2'] <- 'NA'
names(keyvals)[keyvals == 'firebrick'] <- 'CAR'

celltype1 <- c('Dpp4','Clec3b', 'Ly6c', 'Csf1', 'Pi16', 'Adgrd1')
celltype2 <- c('S100a9', 'S100a8', 'Cxcl5', 'Postn', 'Axl', 'Lox')

library(ggplot2)
library(gridExtra)
remotes::install_github("hrbrmstr/ggalt", ref = "noproj")
library(ggalt)
library(scales)
has_ggalt <- ! is(try(find.package("ggalt")), "try-error")


##########figure 4f
EnhancedVolcano(DEG_Pi16_T,
                lab = rownames(DEG_Pi16_T),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlim = c(-1.5, 1.5),
                ylim = c(0, 100),
                title = 'Lining Layer',
                subtitle = "(p<0.05 & FC > 0.5)",
                pCutoff = 0.05,
                FCcutoff = 0.4,
                pointSize = 1.0,
                labSize = 4,
                cutoffLineType = 'blank',
                cutoffLineWidth = 0.3,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                legendPosition = 'none',
                legendLabSize = 4,
                legendIconSize =3.0,
                legendLabels=c('NS','p<0.05 & FC > 0.25'),
                colAlpha = 1,
                colCustom = keyvals,
                drawConnectors = T,
                widthConnectors = 0,
                xlab = bquote(~Log[2]~ 'fold change'),
                boxedLabels = F,
                borderWidth = 0.5,
                max.overlaps =Inf,
                selectLab = c("Ly6c1",	"Fn1",	"Ly6a",	"Dpp4",	"Cd248", 	"Igfbp6",	"Timp2", "Pi16","Sfrp4",	"Il1r2",	"Clec3b",	"Adgrd1", "Osr",	"Axl",	"Procr",	"Opcml",	"Tgfbr2",	"Il18",	"Mfap5",	"Has1",	"Dpt", "Col14a1", "S100a9", "S100a8", "Col12a1", "Cthrc1", "Mif", "Apod", "Cxcl5", "Postn")
)
                
selectLab = c("Itih5", "Sash1", "Col4a1", "Clic5", "Col22a1", "Lbp", "Grid2", "Hbegf", "Cxcl5", "Tnn", "Cthrc1", "Postn", "Mmp13", "Mmp3", "S100a9", "S100a8", "Csrp2", "Gm42418", "Fbln2", "Ptgs2", "Mif", "Cxcl10")



##############Extended data figure 7a

FAPdegs <- FindMarkers(FAP, ident.1 = "PBS", ident.2 = "CART", show.ress = F)


#########FAP fibroblast analysis

#####Extended data figure 8a
######add sublining and lining layer annotation
Idents(int_all.res.1b)<-"arch"
current.sample.ids <- c("Postn", "Prg4", "Mdk", "Fbn1", "Spp1", "Angptl7", "Pi16")
new.sample.ids <- c("Sub_Lining", "Lining", "Sub_Lining", "Sub_Lining", "Sub_Lining", "Sub_Lining", "Sub_Lining")

int_all.res.1b$arch <- int_all.res.1b@active.ident
int_all.res.1b@meta.data[["arch"]] <- plyr::mapvalues(x = int_all.res.1b@meta.data[["arch"]], from = current.sample.ids, to = new.sample.ids)

DimPlot(int_all.res.1b, cols = c("red3", "blue4"))+NoAxes()

#####subset FAP fibroblasts
FAP <-subset(x = int_all.res.1b, subset = Fap > 1)

#####Extended data figure 8a
DimPlot(FAP, order = "PBS", cols = c("brown", "steelblue"), pt.size = 0.01, split.by = "model")+NoLegend()+NoAxes()

Idents(FAP)<-"model"
FAPdegs <- FindMarkers(FAP, ident.1 = "PBS", ident.2 = "CART", show.ress = F)

#####Extended data figure 8b
VlnPlot(int_all.res.1b, features = "Cdh11", pt.size = 0, cols = c("grey", "firebrick"))+ylim(0, 5)+NoLegend()
VlnPlot(int_all.res.1b, features = "Cxcl5", pt.size = 0, cols = c("grey", "firebrick"))+ylim(0, 5)+NoLegend()
VlnPlot(int_all.res.1b, features = "Cthrc1", pt.size = 0, cols = c("grey", "firebrick"))+ylim(0, 5)+NoLegend()
VlnPlot(int_all.res.1b, features = "Postn", pt.size = 0, cols = c("grey", "firebrick"))+ylim(0, 5)+NoLegend()
VlnPlot(int_all.res.1b, features = "Il6", pt.size = 0, cols = c("grey", "firebrick"))+ylim(0, 5)+NoLegend()
VlnPlot(int_all.res.1b, features = "S100a9", pt.size = 0, cols = c("grey", "firebrick"))+ylim(0, 5)+NoLegend()
VlnPlot(int_all.res.1b, features = "S100a8", pt.size = 0, cols = c("grey", "firebrick"))+ylim(0, 5)+NoLegend()
VlnPlot(int_all.res.1b, features = "Mmp3", pt.size = 0, cols = c("grey", "firebrick"))+ylim(0, 5)+NoLegend()
VlnPlot(int_all.res.1b, features = "Mmp13", pt.size = 0, cols = c("grey", "firebrick"))+ylim(0, 5)+NoLegend()
VlnPlot(int_all.res.1b, features = "Tnfsf11", pt.size = 0, cols = c("grey", "firebrick"))+ylim(0, 5)+NoLegend()
VlnPlot(int_all.res.1b, features = "Pi16", pt.size = 0, cols = c("grey", "firebrick"))+ylim(0, 5)+NoLegend()
VlnPlot(int_all.res.1b, features = "Creb5", pt.size = 0, cols = c("grey", "firebrick"))+ylim(0, 5)+NoLegend()

Idents(FAP)<-"model"
pbs1<-subset(FAP, idents="PBS")
cart1<-subset(FAP, idents="CART")

DefaultAssay(object=pbs1) <- "integrated"
pbs1 <- FindVariableFeatures(pbs1, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbs1)
pbs1 <- ScaleData(pbs1, features = all.genes)
pbs1 <- RunPCA(pbs1, npcs = 50, features = VariableFeatures(object = pbs1))
ElbowPlot(pbs1,  ndims = 50)
pbs1 <- FindNeighbors(pbs1, reduction = "pca", dims = 1:30)
# to find out graph.name integrated.H@graphs
pbs1 <- FindClusters(pbs1, graph.name = "integrated_snn", resolution = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7))
clustree(pbs1, assay = "integrated")
pbs1 <- FindClusters(pbs1, graph.name = "integrated_snn", resolution = 0.1)
pbs1 <- RunUMAP(pbs1, reduction = "pca", dims = 1:50)

#########figure 4h
DimPlot(pbs1, reduction = "umap", label = F, pt.size = 0.01, repel = T, label.box = F)+NoAxes()+scale_fill_igv()+scale_color_igv()


cluster0p <- FindMarkers(pbs1, ident.1 = 0, min.pct = 0.4, only.pos = T)
cluster1p <- FindMarkers(pbs1, ident.1 = 1, min.pct = 0.4, only.pos = T)
cluster2p <- FindMarkers(pbs1, ident.1 = 2, min.pct = 0.4, only.pos = T)

DefaultAssay(object=cart1) <- "integrated"
cart1 <- FindVariableFeatures(cart1, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(cart1)
cart1 <- ScaleData(cart1, features = all.genes)
cart1 <- RunPCA(cart1, npcs = 50, features = VariableFeatures(object = cart1))
ElbowPlot(cart1,  ndims = 50)
cart1 <- FindNeighbors(cart1, reduction = "pca", dims = 1:30)
# to find out graph.name integrated.H@graphs
cart1 <- FindClusters(cart1, graph.name = "integrated_snn", resolution = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7))
clustree(cart1, assay = "integrated")
cart1 <- FindClusters(cart1, graph.name = "integrated_snn", resolution = 0.2)
cart1 <- RunUMAP(cart1, reduction = "pca", dims = 1:50)

#########figure 4h
DimPlot(cart1, reduction = "umap", label = F, pt.size = 0.01, repel = T, label.box = F)+NoAxes()+scale_fill_igv()+scale_color_igv()


cluster0 <- FindMarkers(cart1, ident.1 = 0, min.pct = 0.4, only.pos = T)
cluster1 <- FindMarkers(cart1, ident.1 = 1, min.pct = 0.4, only.pos = T)
cluster2 <- FindMarkers(cart1, ident.1 = 2, min.pct = 0.4, only.pos = T)
cluster3 <- FindMarkers(cart1, ident.1 = 3, min.pct = 0.4, only.pos = T)

write.csv(cluster3, file = "Fibroblast.csv")

gene.fap <- c("Prg4", "Cd55", "Has1", "Clic5", "Tspan15", "Hbegf", "F13a1")

#########Extended data figure 8d
dotplot <- DotPlot(pbs1, features = gene.fap, dot.scale = 5, cols = "BrBG") + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) +ylab(NULL) +xlab(NULL)
dotplot
#########Extended data figure 8d
dotplot <- DotPlot(cart1, features = gene.fap, dot.scale = 5, cols = "BrBG") + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) +ylab(NULL) +xlab(NULL)
dotplot


#######creb5 progenitor analysis
prog.genes<-c("Thy1", "Prg4", "Cd55", "Clic5", "Tspan15", "Itga6", "Creb5", "Megf10", "Enpp5", "Rspo2", "Cdkn2c", "Pla1a", "Sparcl1")

#########figure 4i
####build pheatmap from dotplot FAP+ non-treated fibroblasts
dotplot <- DotPlot(pbs1, features = prog.genes, dot.scale = 5) + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))+scale_colour_viridis(option="magma")
dotplot

fbdotplot.data <- dotplot[["data"]]
#remove non essesntial info
fbdotplot.data_edit = subset(fbdotplot.data, select = -c(1:3) )

fbdotplot.data_edit_0 = subset(fbdotplot.data_edit, id == "0")
fbdotplot.data_edit_0 = subset(fbdotplot.data_edit_0, select = -c(1) )
names(fbdotplot.data_edit_0)[names(fbdotplot.data_edit_0) == 'avg.exp.scaled'] <- 'Cluster0'

fbdotplot.data_edit_1 = subset(fbdotplot.data_edit, id == "1")
fbdotplot.data_edit_1 = subset(fbdotplot.data_edit_1, select = -c(1) )
names(fbdotplot.data_edit_1)[names(fbdotplot.data_edit_1) == 'avg.exp.scaled'] <- 'Cluster1'

fbdotplot.data_edit_2 = subset(fbdotplot.data_edit, id == "2")
fbdotplot.data_edit_2 = subset(fbdotplot.data_edit_2, select = -c(1) )
names(fbdotplot.data_edit_2)[names(fbdotplot.data_edit_2) == 'avg.exp.scaled'] <- 'Cluster2'


fbdotplot.data_edit_0$Cluster1 <- fbdotplot.data_edit_1$Cluster1
fbdotplot.data_edit_0$Cluster2 <- fbdotplot.data_edit_2$Cluster2

mat_fb <- fbdotplot.data_edit_0
mat_fb <- as.matrix(mat_fb)

cat_macs=data.frame("seurat_clusters"=c("0", "1", "2"))
rownames(cat_macs)=colnames(mat_fb)

anncols <- list("Fibroblasts" = c(Cluster0="#800000FF", Cluster1="#767676FF", Cluster2="#CC8214FF"))
pheatmap(mat_fb, main="Prog Genes", cluster_rows=F , cluster_cols = F, annotation_col=cat_macs, cellwidth = 20, show_colnames = F, color = brewer.pal(5, "Blues"), annotation_colors = anncols) 

#########figure 4i
####build pheatmap from dotplot FAP+ CAR fibroblasts
dotplot <- DotPlot(cart1, features = prog.genes, dot.scale = 5) + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))+scale_colour_viridis(option="magma")
dotplot

fbdotplot.data <- dotplot[["data"]]
#remove non essesntial info
fbdotplot.data_edit = subset(fbdotplot.data, select = -c(1:3) )

fbdotplot.data_edit_0 = subset(fbdotplot.data_edit, id == "0")
fbdotplot.data_edit_0 = subset(fbdotplot.data_edit_0, select = -c(1) )
names(fbdotplot.data_edit_0)[names(fbdotplot.data_edit_0) == 'avg.exp.scaled'] <- 'Cluster0'

fbdotplot.data_edit_1 = subset(fbdotplot.data_edit, id == "1")
fbdotplot.data_edit_1 = subset(fbdotplot.data_edit_1, select = -c(1) )
names(fbdotplot.data_edit_1)[names(fbdotplot.data_edit_1) == 'avg.exp.scaled'] <- 'Cluster1'

fbdotplot.data_edit_2 = subset(fbdotplot.data_edit, id == "2")
fbdotplot.data_edit_2 = subset(fbdotplot.data_edit_2, select = -c(1) )
names(fbdotplot.data_edit_2)[names(fbdotplot.data_edit_2) == 'avg.exp.scaled'] <- 'Cluster2'

fbdotplot.data_edit_3 = subset(fbdotplot.data_edit, id == "3")
fbdotplot.data_edit_3 = subset(fbdotplot.data_edit_3, select = -c(1) )
names(fbdotplot.data_edit_3)[names(fbdotplot.data_edit_3) == 'avg.exp.scaled'] <- 'Cluster3'

fbdotplot.data_edit_0$Cluster1 <- fbdotplot.data_edit_1$Cluster1
fbdotplot.data_edit_0$Cluster2 <- fbdotplot.data_edit_2$Cluster2
fbdotplot.data_edit_0$Cluster3 <- fbdotplot.data_edit_3$Cluster3

mat_fb <- fbdotplot.data_edit_0
mat_fb <- as.matrix(mat_fb)

cat_macs=data.frame("seurat_clusters"=c("0", "1", "2", "3"))
rownames(cat_macs)=colnames(mat_fb)

anncols <- list("Fibroblasts" = c(Cluster0="#800000FF", Cluster1="#767676FF", Cluster2="#CC8214FF", Cluster3="#616530FF"))
pheatmap(mat_fb, main="Prog Genes", cluster_rows=F , cluster_cols = F, annotation_col=cat_macs, cellwidth = 20, show_colnames = F, color = brewer.pal(5, "Blues"), annotation_colors = anncols) 
