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


getwd()

###########resting
bsa_c1 <- Read10X(data.dir = "/rds/projects/c/croftap-aia-seq-data/ALL_data/all_neg/D0_cd45neg_1_b5")
bsa_c1 <- CreateSeuratObject(counts = bsa_c1, min.cells=3, min.features=100, project="bsa_c1")
bsa_c1 <- RenameCells(bsa_c1, add.cell.id = "bsa_c1")
bsa_c1$sample_id <- "Rest1"
bsa_c1$type <- "CD45neg"
bsa_c1$InflammationStatus <- "Resting"
head(bsa_c1)

bsa_c2 <- Read10X(data.dir = "/rds/projects/c/croftap-aia-seq-data/ALL_data/all_neg/D0_cd45neg_2_b3")
bsa_c2 <- CreateSeuratObject(counts = bsa_c2, min.cells=3, min.features=100, project="bsa_c2")
bsa_c2 <- RenameCells(bsa_c2, add.cell.id = "bsa_c2")
bsa_c2$sample_id <- "Rest2"
bsa_c2$type <- "CD45neg"
bsa_c2$InflammationStatus <- "Resting"
head(bsa_c2)

bsa_c3 <- Read10X(data.dir = "/rds/projects/c/croftap-aia-seq-data/ALL_data/all_neg/D0_cd45neg_3_b6")
bsa_c3 <- CreateSeuratObject(counts = bsa_c3, min.cells=3, min.features=100, project="bsa_c3")
bsa_c3 <- RenameCells(bsa_c3, add.cell.id = "bsa_c3")
bsa_c3$sample_id <- "Rest3"
bsa_c3$type <- "CD45neg"
bsa_c3$InflammationStatus <- "Resting"
head(bsa_c3)


###########AIA
##############peak of inflammation
bsa_In1 <- Read10X(data.dir = "/rds/projects/c/croftap-aia-seq-data/ALL_data/all_neg/D2_cd45neg_1_b5")
bsa_In1 <- CreateSeuratObject(counts = bsa_In1, min.cells=3, min.features=100, project="bsa_In1")
bsa_In1 <- RenameCells(bsa_In1, add.cell.id = "bsa_In1")
bsa_In1$sample_id <- "Peak1"
bsa_In1$type <- "CD45neg"
bsa_In1$InflammationStatus <- "Peak"
head(bsa_In1)

bsa_In2 <- Read10X(data.dir = "/rds/projects/c/croftap-aia-seq-data/ALL_data/all_neg/D2_cd45neg_2_b5")
bsa_In2 <- CreateSeuratObject(counts = bsa_In2, min.cells=3, min.features=100, project="bsa_In2")
bsa_In2 <- RenameCells(bsa_In2, add.cell.id = "bsa_In2")
bsa_In2$sample_id <- "Peak2"
bsa_In2$type <- "CD45neg"
bsa_In2$InflammationStatus <- "Peak"
head(bsa_In2)

bsa_In3 <- Read10X(data.dir = "/rds/projects/c/croftap-aia-seq-data/ALL_data/all_neg/D2_cd45neg_3_b1")
bsa_In3 <- CreateSeuratObject(counts = bsa_In3, min.cells=3, min.features=100, project="bsa_In3")
bsa_In3 <- RenameCells(bsa_In3, add.cell.id = "bsa_In3")
bsa_In3$sample_id <- "Peak3"
bsa_In3$type <- "CD45neg"
bsa_In3$InflammationStatus <- "Peak"
head(bsa_In3)

############early resolving
bsa_Ereslvn1 <- Read10X(data.dir = "/rds/projects/c/croftap-aia-seq-data/ALL_data/all_neg/D4_cd45neg_1_b1")
bsa_Ereslvn1 <- CreateSeuratObject(counts = bsa_Ereslvn1, min.cells=3, min.features=100, project="bsa_Ereslvn1")
bsa_Ereslvn1 <- RenameCells(bsa_Ereslvn1, add.cell.id = "bsa_Ereslvn1")
bsa_Ereslvn1$sample_id <- "Ereslvn1"
bsa_Ereslvn1$type <- "CD45neg"
bsa_Ereslvn1$InflammationStatus <- "Ereslvn"
head(bsa_Ereslvn1)

bsa_Ereslvn2 <- Read10X(data.dir = "/rds/projects/c/croftap-aia-seq-data/ALL_data/all_neg/D4_cd45neg_2_b1")
bsa_Ereslvn2 <- CreateSeuratObject(counts = bsa_Ereslvn2, min.cells=3, min.features=100, project="bsa_Ereslvn2")
bsa_Ereslvn2 <- RenameCells(bsa_Ereslvn2, add.cell.id = "bsa_Ereslvn2")
bsa_Ereslvn2$sample_id <- "Ereslvn2"
bsa_Ereslvn2$type <- "CD45neg"
bsa_Ereslvn2$InflammationStatus <- "Ereslvn"
head(bsa_Ereslvn2)

bsa_Ereslvn3 <- Read10X(data.dir = "/rds/projects/c/croftap-aia-seq-data/ALL_data/all_neg/D4_cd45neg_3_b1")
bsa_Ereslvn3 <- CreateSeuratObject(counts = bsa_Ereslvn3, min.cells=3, min.features=100, project="bsa_Ereslvn3")
bsa_Ereslvn3 <- RenameCells(bsa_Ereslvn3, add.cell.id = "bsa_Ereslvn3")
bsa_Ereslvn3$sample_id <- "Ereslvn3"
bsa_Ereslvn3$type <- "CD45neg"
bsa_Ereslvn3$InflammationStatus <- "Ereslvn"
head(bsa_Ereslvn3)

############Resolving
bsa_Reslvn1 <- Read10X(data.dir = "/rds/projects/c/croftap-aia-seq-data/ALL_data/all_neg/D7_cd45neg_1_b2")
bsa_Reslvn1 <- CreateSeuratObject(counts = bsa_Reslvn1, min.cells=3, min.features=100, project="bsa_Reslvn1")
bsa_Reslvn1 <- RenameCells(bsa_Reslvn1, add.cell.id = "bsa_Reslvn1")
bsa_Reslvn1$sample_id <- "Reslvn1"
bsa_Reslvn1$type <- "CD45neg"
bsa_Reslvn1$InflammationStatus <- "Reslvn"
head(bsa_Reslvn1)

bsa_Reslvn2 <- Read10X(data.dir = "/rds/projects/c/croftap-aia-seq-data/ALL_data/all_neg/D7_cd45neg_2_b2")
bsa_Reslvn2 <- CreateSeuratObject(counts = bsa_Reslvn2, min.cells=3, min.features=100, project="bsa_Reslvn2")
bsa_Reslvn2 <- RenameCells(bsa_Reslvn2, add.cell.id = "bsa_Reslvn2")
bsa_Reslvn2$sample_id <- "Reslvn2"
bsa_Reslvn2$type <- "CD45neg"
bsa_Reslvn2$InflammationStatus <- "Reslvn"
head(bsa_Reslvn2)

bsa_Reslvn3 <- Read10X(data.dir = "/rds/projects/c/croftap-aia-seq-data/ALL_data/all_neg/D7_cd45neg_3_b2")
bsa_Reslvn3 <- CreateSeuratObject(counts = bsa_Reslvn3, min.cells=3, min.features=100, project="bsa_Reslvn3")
bsa_Reslvn3 <- RenameCells(bsa_Reslvn3, add.cell.id = "bsa_Reslvn3")
bsa_Reslvn3$sample_id <- "Reslvn3"
bsa_Reslvn3$type <- "CD45neg"
bsa_Reslvn3$InflammationStatus <- "Reslvn"
head(bsa_Reslvn3)

############resolved
bsa_Reslvd1 <- Read10X(data.dir = "/rds/projects/c/croftap-aia-seq-data/ALL_data/all_neg/D14_cd45neg_1_b6")
bsa_Reslvd1 <- CreateSeuratObject(counts = bsa_Reslvd1, min.cells=3, min.features=100, project="bsa_Reslvd1")
bsa_Reslvd1 <- RenameCells(bsa_Reslvd1, add.cell.id = "bsa_Reslvd1")
bsa_Reslvd1$sample_id <- "Reslvd1"
bsa_Reslvd1$type <- "CD45neg"
bsa_Reslvd1$InflammationStatus <- "Reslvd"
head(bsa_Reslvd1)

bsa_Reslvd2 <- Read10X(data.dir = "/rds/projects/c/croftap-aia-seq-data/ALL_data/all_neg/D14_cd45neg_2_b6")
bsa_Reslvd2 <- CreateSeuratObject(counts = bsa_Reslvd2, min.cells=3, min.features=100, project="bsa_Reslvd2")
bsa_Reslvd2 <- RenameCells(bsa_Reslvd2, add.cell.id = "bsa_Reslvd2")
bsa_Reslvd2$sample_id <- "Reslvd2"
bsa_Reslvd2$type <- "CD45neg"
bsa_Reslvd2$InflammationStatus <- "Reslvd"
head(bsa_Reslvd2)

bsa_Reslvd3 <- Read10X(data.dir = "/rds/projects/c/croftap-aia-seq-data/ALL_data/all_neg/D14_cd45neg_3_b6")
bsa_Reslvd3 <- CreateSeuratObject(counts = bsa_Reslvd3, min.cells=3, min.features=100, project="bsa_Reslvd3")
bsa_Reslvd3 <- RenameCells(bsa_Reslvd3, add.cell.id = "bsa_Reslvd3")
bsa_Reslvd3$sample_id <- "Reslvd3"
bsa_Reslvd3$type <- "CD45neg"
bsa_Reslvd3$InflammationStatus <- "Reslvd"
head(bsa_Reslvd3)

#############STIA
###########resting
STIA_c1 <- Read10X(data.dir = "/rds/projects/2018/croftap-time-course-rna-ss/STIA_data_timecourse/FASTQ/CON_A/Align/CON_A/outs/filtered_feature_bc_matrix")
STIA_c1 <- CreateSeuratObject(counts = STIA_c1, min.cells=3, min.features=100, project="STIA_c1")
STIA_c1 <- RenameCells(STIA_c1, add.cell.id = "STIA_c1")
STIA_c1$sample_id <- "STIAControl1"
STIA_c1$model <- "STIA"
STIA_c1$InflammationStatus <- "Control"
head(STIA_c1)

STIA_c2 <- Read10X(data.dir = "/rds/projects/2018/croftap-time-course-rna-ss/STIA_data_timecourse/FASTQ/CON_B/Align/CON_B/outs/filtered_feature_bc_matrix")
STIA_c2 <- CreateSeuratObject(counts = STIA_c2, min.cells=3, min.features=100, project="STIA_c2")
STIA_c2 <- RenameCells(STIA_c2, add.cell.id = "STIA_c2")
STIA_c2$sample_id <- "STIAControl2"
STIA_c2$model <- "STIA"
STIA_c2$InflammationStatus <- "Control"
head(STIA_c2)

STIA_c3 <- Read10X(data.dir = "/rds/projects/2018/croftap-time-course-rna-ss/STIA_data_timecourse/FASTQ/CON_C/Align/CON_C/outs/filtered_feature_bc_matrix")
STIA_c3 <- CreateSeuratObject(counts = STIA_c3, min.cells=3, min.features=100, project="STIA_c3")
STIA_c3 <- RenameCells(STIA_c3, add.cell.id = "STIA_c3")
STIA_c3$sample_id <- "STIAControl3"
STIA_c3$model <- "STIA"
STIA_c3$InflammationStatus <- "Control"
head(STIA_c3)

##############peak of inflammation
STIA_in1 <- Read10X(data.dir = "/rds/projects/2018/croftap-time-course-rna-ss/STIA_data_timecourse/FASTQ/PEAK_A/Align/PEAK_A/outs/filtered_feature_bc_matrix")
STIA_in1 <- CreateSeuratObject(counts = STIA_in1, min.cells=3, min.features=100, project="STIA_in1")
STIA_in1 <- RenameCells(STIA_in1, add.cell.id = "STIA_in1")
STIA_in1$sample_id <- "STIAInflamed1"
STIA_in1$model <- "STIA"
STIA_in1$InflammationStatus <- "Inflamed"
head(STIA_in1)

STIA_in2 <- Read10X(data.dir = "/rds/projects/2018/croftap-time-course-rna-ss/STIA_data_timecourse/FASTQ/PEAK_B/Align/PEAK_B/outs/filtered_feature_bc_matrix")
STIA_in2 <- CreateSeuratObject(counts = STIA_in2, min.cells=3, min.features=100, project="STIA_in2")
STIA_in2 <- RenameCells(STIA_in2, add.cell.id = "STIA_in2")
STIA_in2$sample_id <- "STIAInflamed2"
STIA_in2$model <- "STIA"
STIA_in2$InflammationStatus <- "Inflamed"
head(STIA_in2)

STIA_in3 <- Read10X(data.dir = "/rds/projects/2018/croftap-time-course-rna-ss/STIA_data_timecourse/FASTQ/PEAK_C/Align/PEAK_C/outs/filtered_feature_bc_matrix")
STIA_in3 <- CreateSeuratObject(counts = STIA_in3, min.cells=3, min.features=100, project="STIA_in3")
STIA_in3 <- RenameCells(STIA_in3, add.cell.id = "STIA_in3")
STIA_in3$sample_id <- "STIAInflamed3"
STIA_in3$model <- "STIA"
STIA_in3$InflammationStatus <- "Inflamed"
head(STIA_in3)

###########resolving
STIA_Reslvn1 <- Read10X(data.dir = "/rds/projects/2018/croftap-time-course-rna-ss/STIA_data_timecourse/FASTQ/RESING_A/Align/RESING_A/outs/filtered_feature_bc_matrix")
STIA_Reslvn1 <- CreateSeuratObject(counts = STIA_Reslvn1, min.cells=3, min.features=100, project="STIA_Reslvn1")
STIA_Reslvn1 <- RenameCells(STIA_Reslvn1, add.cell.id = "STIA_Reslvn1")
STIA_Reslvn1$sample_id <- "STIAReslvn1"
STIA_Reslvn1$model <- "STIA"
STIA_Reslvn1$InflammationStatus <- "Resolving"
head(STIA_Reslvn1)

STIA_Reslvn2 <- Read10X(data.dir = "/rds/projects/2018/croftap-time-course-rna-ss/STIA_data_timecourse/FASTQ/RESING_B/Align/RESING_B/outs/filtered_feature_bc_matrix")
STIA_Reslvn2 <- CreateSeuratObject(counts = STIA_Reslvn2, min.cells=3, min.features=100, project="STIA_Reslvn2")
STIA_Reslvn2 <- RenameCells(STIA_Reslvn2, add.cell.id = "STIA_Reslvn2")
STIA_Reslvn2$sample_id <- "STIAReslvn2"
STIA_Reslvn2$model <- "STIA"
STIA_Reslvn2$InflammationStatus <- "Resolving"
head(STIA_Reslvn2)

STIA_Reslvn3 <- Read10X(data.dir = "/rds/projects/2018/croftap-time-course-rna-ss/STIA_data_timecourse/FASTQ/RESING_C/Align/RESING_C/outs/filtered_feature_bc_matrix")
STIA_Reslvn3 <- CreateSeuratObject(counts = STIA_Reslvn3, min.cells=3, min.features=100, project="STIA_Reslvn3")
STIA_Reslvn3 <- RenameCells(STIA_Reslvn3, add.cell.id = "STIA_Reslvn3")
STIA_Reslvn3$sample_id <- "STIAReslvn3"
STIA_Reslvn3$model <- "STIA"
STIA_Reslvn3$InflammationStatus <- "Resolving"
head(STIA_Reslvn3)

############resolved
STIA_Reslvd1 <- Read10X(data.dir = "/rds/projects/2018/croftap-time-course-rna-ss/STIA_data_timecourse/FASTQ/RESED_A/Align/RESED_A/outs/filtered_feature_bc_matrix")
STIA_Reslvd1 <- CreateSeuratObject(counts = STIA_Reslvd1, min.cells=3, min.features=100, project="STIA_Reslvd1")
STIA_Reslvd1 <- RenameCells(STIA_Reslvd1, add.cell.id = "STIA_Reslvd1")
STIA_Reslvd1$sample_id <- "STIAReslvd1"
STIA_Reslvd1$model <- "STIA"
STIA_Reslvd1$InflammationStatus <- "Resolved"
head(STIA_Reslvd1)

STIA_Reslvd2 <- Read10X(data.dir = "/rds/projects/2018/croftap-time-course-rna-ss/STIA_data_timecourse/FASTQ/RESED_B/Align/RESED_B/outs/filtered_feature_bc_matrix")
STIA_Reslvd2 <- CreateSeuratObject(counts = STIA_Reslvd2, min.cells=3, min.features=100, project="STIA_Reslvd2")
STIA_Reslvd2 <- RenameCells(STIA_Reslvd2, add.cell.id = "STIA_Reslvd2")
STIA_Reslvd2$sample_id <- "STIAReslvd2"
STIA_Reslvd2$model <- "STIA"
STIA_Reslvd2$InflammationStatus <- "Resolved"
head(STIA_Reslvd2)

STIA_Reslvd3 <- Read10X(data.dir = "/rds/projects/2018/croftap-time-course-rna-ss/STIA_data_timecourse/FASTQ/RESED_C/Align/RESED_C/outs/filtered_feature_bc_matrix")
STIA_Reslvd3 <- CreateSeuratObject(counts = STIA_Reslvd3, min.cells=3, min.features=100, project="STIA_Reslvd3")
STIA_Reslvd3 <- RenameCells(STIA_Reslvd3, add.cell.id = "STIA_Reslvd3")
STIA_Reslvd3$sample_id <- "STIAReslvd3"
STIA_Reslvd3$model <- "STIA"
STIA_Reslvd3$InflammationStatus <- "Resolved"
head(STIA_Reslvd3)

#####################CIA
############resting
CIA_c2 <- Read10X(data.dir = "/rds/projects/c/croftap-ktr-cia-sc-01/CIA_CD45neg/Con2_CD45neg/outs/filtered_feature_bc_matrix")
CIA_c2 <- CreateSeuratObject(counts = CIA_c2, min.cells=3, min.features=100, project="CIA_c2")
CIA_c2 <- RenameCells(CIA_c2, add.cell.id = "CIA_c2")
CIA_c2$sample_id <- "CIAControl2"
CIA_c2$model <- "CIA"
CIA_c2$InflammationStatus <- "Control"
head(CIA_c2)

CIA_c3 <- Read10X(data.dir = "/rds/projects/c/croftap-ktr-cia-sc-01/CIA_CD45neg/Con3_CD45neg/outs/filtered_feature_bc_matrix")
CIA_c3 <- CreateSeuratObject(counts = CIA_c3, min.cells=3, min.features=100, project="CIA_c3")
CIA_c3 <- RenameCells(CIA_c3, add.cell.id = "CIA_c3")
CIA_c3$sample_id <- "CIAControl3"
CIA_c3$model <- "CIA"
CIA_c3$InflammationStatus <- "Control"
head(CIA_c3)

CIA_c4 <- Read10X(data.dir = "/rds/projects/c/croftap-ktr-cia-sc-01/CIA_CD45_pos/Con4_CD45neg/outs/filtered_feature_bc_matrix")
CIA_c4 <- CreateSeuratObject(counts = CIA_c4, min.cells=3, min.features=100, project="CIA_c4")
CIA_c4 <- RenameCells(CIA_c4, add.cell.id = "CIA_c4")
CIA_c4$sample_id <- "CIAControl4"
CIA_c4$model <- "CIA"
CIA_c4$InflammationStatus <- "Control"
head(CIA_c4)

##############inflamed
CIA_in1 <- Read10X(data.dir = "/rds/projects/c/croftap-ktr-cia-sc-01/CIA_CD45neg/Infla1_CD45neg/outs/filtered_feature_bc_matrix")
CIA_in1 <- CreateSeuratObject(counts = CIA_in1, min.cells=3, min.features=100, project="CIA_in1")
CIA_in1 <- RenameCells(CIA_in1, add.cell.id = "CIA_in1")
CIA_in1$sample_id <- "CIAInflamed1"
CIA_in1$model <- "CIA"
CIA_in1$InflammationStatus <- "Inflamed"
head(CIA_in1)

CIA_in2 <- Read10X(data.dir = "/rds/projects/c/croftap-ktr-cia-sc-01/CIA_CD45neg/Infla2_CD45neg/outs/filtered_feature_bc_matrix")
CIA_in2 <- CreateSeuratObject(counts = CIA_in2, min.cells=3, min.features=100, project="CIA_in2")
CIA_in2 <- RenameCells(CIA_in2, add.cell.id = "CIA_in2")
CIA_in2$sample_id <- "CIAInflamed2"
CIA_in2$model <- "CIA"
CIA_in2$InflammationStatus <- "Inflamed"
head(CIA_in2)

CIA_in3 <- Read10X(data.dir = "/rds/projects/c/croftap-ktr-cia-sc-01/CIA_CD45neg/Infla3_CD45neg/outs/filtered_feature_bc_matrix")
CIA_in3 <- CreateSeuratObject(counts = CIA_in3, min.cells=3, min.features=100, project="CIA_in3")
CIA_in3 <- RenameCells(CIA_in3, add.cell.id = "CIA_in3")
CIA_in3$sample_id <- "CIAInflamed3"
CIA_in3$model <- "CIA"
CIA_in3$InflammationStatus <- "Inflamed"
head(CIA_in3)


###QC
mito.featuresbsa_c1 <- grep(pattern="^mt-", x=rownames(x=bsa_c1), value=T)
percent.mitobsa_c1 <- Matrix::colSums(x = GetAssayData(object = bsa_c1, slot = "counts")[mito.featuresbsa_c1,]) / Matrix::colSums(x = GetAssayData(object = bsa_c1, slot = "counts"))
bsa_c1[["percent.mito"]] <- percent.mitobsa_c1
VlnPlot(object = bsa_c1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

bsa_c1 <- subset(x = bsa_c1, subset = nFeature_RNA > 500 & nFeature_RNA <7000 & nCount_RNA > 500 & nCount_RNA < 50000 & percent.mitobsa_c1 < 0.1)
VlnPlot(object = bsa_c1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
bsa_c1 <- NormalizeData(object = bsa_c1, verbose = F)
bsa_c1 <- FindVariableFeatures(object = bsa_c1, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesbsa_c1 <- rownames(bsa_c1)
bsa_c1 <- ScaleData(bsa_c1, features = all.genesbsa_c1)
VlnPlot(object = bsa_c1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)



mito.featuresbsa_c2 <- grep(pattern="^mt-", x=rownames(x=bsa_c2), value=T)
percent.mitobsa_c2 <- Matrix::colSums(x = GetAssayData(object = bsa_c2, slot = "counts")[mito.featuresbsa_c2,]) / Matrix::colSums(x = GetAssayData(object = bsa_c2, slot = "counts"))
bsa_c2[["percent.mito"]] <- percent.mitobsa_c2
VlnPlot(object = bsa_c2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

bsa_c2 <- subset(x = bsa_c2, subset = nFeature_RNA > 500 & nFeature_RNA <7000 & nCount_RNA > 500 & nCount_RNA < 40000 & percent.mitobsa_c2 < 0.1)
VlnPlot(object = bsa_c2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
bsa_c2 <- NormalizeData(object = bsa_c2, verbose = F)
bsa_c2 <- FindVariableFeatures(object = bsa_c2, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesbsa_c2 <- rownames(bsa_c2)
bsa_c2 <- ScaleData(bsa_c2, features = all.genesbsa_c2)
VlnPlot(object = bsa_c2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)



mito.featuresbsa_c3 <- grep(pattern="^mt-", x=rownames(x=bsa_c3), value=T)
percent.mitobsa_c3 <- Matrix::colSums(x = GetAssayData(object = bsa_c3, slot = "counts")[mito.featuresbsa_c3,]) / Matrix::colSums(x = GetAssayData(object = bsa_c3, slot = "counts"))
bsa_c3[["percent.mito"]] <- percent.mitobsa_c3
VlnPlot(object = bsa_c3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

bsa_c3 <- subset(x = bsa_c3, subset = nFeature_RNA > 500 & nFeature_RNA <6000 & nCount_RNA > 500 & nCount_RNA < 40000 & percent.mitobsa_c3 < 0.1)
VlnPlot(object = bsa_c3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
bsa_c3 <- NormalizeData(object = bsa_c3, verbose = F)
bsa_c3 <- FindVariableFeatures(object = bsa_c3, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesbsa_c3 <- rownames(bsa_c3)
bsa_c3 <- ScaleData(bsa_c3, features = all.genesbsa_c3)
VlnPlot(object = bsa_c3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)



mito.featuresbsa_In1 <- grep(pattern="^mt-", x=rownames(x=bsa_In1), value=T)
percent.mitobsa_In1 <- Matrix::colSums(x = GetAssayData(object = bsa_In1, slot = "counts")[mito.featuresbsa_In1,]) / Matrix::colSums(x = GetAssayData(object = bsa_In1, slot = "counts"))
bsa_In1[["percent.mito"]] <- percent.mitobsa_In1
VlnPlot(object = bsa_In1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

bsa_In1 <- subset(x = bsa_In1, subset = nFeature_RNA > 500 & nFeature_RNA <6000 & nCount_RNA > 500 & nCount_RNA < 40000 & percent.mitobsa_In1 < 0.1)
VlnPlot(object = bsa_In1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
bsa_In1 <- NormalizeData(object = bsa_In1, verbose = F)
bsa_In1 <- FindVariableFeatures(object = bsa_In1, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesbsa_In1 <- rownames(bsa_In1)
bsa_In1 <- ScaleData(bsa_In1, features = all.genesbsa_In1)
VlnPlot(object = bsa_In1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

mito.featuresbsa_In2 <- grep(pattern="^mt-", x=rownames(x=bsa_In2), value=T)
percent.mitobsa_In2 <- Matrix::colSums(x = GetAssayData(object = bsa_In2, slot = "counts")[mito.featuresbsa_In2,]) / Matrix::colSums(x = GetAssayData(object = bsa_In2, slot = "counts"))
bsa_In2[["percent.mito"]] <- percent.mitobsa_In2
VlnPlot(object = bsa_In2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

bsa_In2 <- subset(x = bsa_In2, subset = nFeature_RNA > 500 & nFeature_RNA <6500 & nCount_RNA > 500 & nCount_RNA < 50000 & percent.mitobsa_In2 < 0.1)
VlnPlot(object = bsa_In2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
bsa_In2 <- NormalizeData(object = bsa_In2, verbose = F)
bsa_In2 <- FindVariableFeatures(object = bsa_In2, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesbsa_In2 <- rownames(bsa_In2)
bsa_In2 <- ScaleData(bsa_In2, features = all.genesbsa_In2)
VlnPlot(object = bsa_In2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)



mito.featuresbsa_In3 <- grep(pattern="^mt-", x=rownames(x=bsa_In3), value=T)
percent.mitobsa_In3 <- Matrix::colSums(x = GetAssayData(object = bsa_In3, slot = "counts")[mito.featuresbsa_In3,]) / Matrix::colSums(x = GetAssayData(object = bsa_In3, slot = "counts"))
bsa_In3[["percent.mito"]] <- percent.mitobsa_In3
VlnPlot(object = bsa_In3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

bsa_In3 <- subset(x = bsa_In3, subset = nFeature_RNA > 500 & nFeature_RNA <7000 & nCount_RNA > 500 & nCount_RNA < 60000 & percent.mitobsa_In3 < 0.1)
VlnPlot(object = bsa_In3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
bsa_In3 <- NormalizeData(object = bsa_In3, verbose = F)
bsa_In3 <- FindVariableFeatures(object = bsa_In3, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesbsa_In3 <- rownames(bsa_In3)
bsa_In3 <- ScaleData(bsa_In3, features = all.genesbsa_In3)
VlnPlot(object = bsa_In3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)



mito.featuresbsa_Ereslvn1 <- grep(pattern="^mt-", x=rownames(x=bsa_Ereslvn1), value=T)
percent.mitobsa_Ereslvn1 <- Matrix::colSums(x = GetAssayData(object = bsa_Ereslvn1, slot = "counts")[mito.featuresbsa_Ereslvn1,]) / Matrix::colSums(x = GetAssayData(object = bsa_Ereslvn1, slot = "counts"))
bsa_Ereslvn1[["percent.mito"]] <- percent.mitobsa_Ereslvn1
VlnPlot(object = bsa_Ereslvn1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

bsa_Ereslvn1 <- subset(x = bsa_Ereslvn1, subset = nFeature_RNA > 500 & nFeature_RNA <7000 & nCount_RNA > 500 & nCount_RNA < 50000 & percent.mitobsa_Ereslvn1 < 0.1)
VlnPlot(object = bsa_Ereslvn1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
bsa_Ereslvn1 <- NormalizeData(object = bsa_Ereslvn1, verbose = F)
bsa_Ereslvn1 <- FindVariableFeatures(object = bsa_Ereslvn1, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesbsa_Ereslvn1 <- rownames(bsa_Ereslvn1)
bsa_Ereslvn1 <- ScaleData(bsa_Ereslvn1, features = all.genesbsa_Ereslvn1)
VlnPlot(object = bsa_Ereslvn1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

mito.featuresbsa_Ereslvn2 <- grep(pattern="^mt-", x=rownames(x=bsa_Ereslvn2), value=T)
percent.mitobsa_Ereslvn2 <- Matrix::colSums(x = GetAssayData(object = bsa_Ereslvn2, slot = "counts")[mito.featuresbsa_Ereslvn2,]) / Matrix::colSums(x = GetAssayData(object = bsa_Ereslvn2, slot = "counts"))
bsa_Ereslvn2[["percent.mito"]] <- percent.mitobsa_Ereslvn2
VlnPlot(object = bsa_Ereslvn2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

bsa_Ereslvn2 <- subset(x = bsa_Ereslvn2, subset = nFeature_RNA > 500 & nFeature_RNA <7000 & nCount_RNA > 500 & nCount_RNA < 60000 & percent.mitobsa_Ereslvn2 < 0.1)
VlnPlot(object = bsa_Ereslvn2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
bsa_Ereslvn2 <- NormalizeData(object = bsa_Ereslvn2, verbose = F)
bsa_Ereslvn2 <- FindVariableFeatures(object = bsa_Ereslvn2, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesbsa_Ereslvn2 <- rownames(bsa_Ereslvn2)
bsa_Ereslvn2 <- ScaleData(bsa_Ereslvn2, features = all.genesbsa_Ereslvn2)
VlnPlot(object = bsa_Ereslvn2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)



mito.featuresbsa_Ereslvn3 <- grep(pattern="^mt-", x=rownames(x=bsa_Ereslvn3), value=T)
percent.mitobsa_Ereslvn3 <- Matrix::colSums(x = GetAssayData(object = bsa_Ereslvn3, slot = "counts")[mito.featuresbsa_Ereslvn3,]) / Matrix::colSums(x = GetAssayData(object = bsa_Ereslvn3, slot = "counts"))
bsa_Ereslvn3[["percent.mito"]] <- percent.mitobsa_Ereslvn3
VlnPlot(object = bsa_Ereslvn3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

bsa_Ereslvn3 <- subset(x = bsa_Ereslvn3, subset = nFeature_RNA > 500 & nFeature_RNA <7000 & nCount_RNA > 500 & nCount_RNA < 60000 & percent.mitobsa_Ereslvn3 < 0.1)
VlnPlot(object = bsa_Ereslvn3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
bsa_Ereslvn3 <- NormalizeData(object = bsa_Ereslvn3, verbose = F)
bsa_Ereslvn3 <- FindVariableFeatures(object = bsa_Ereslvn3, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesbsa_Ereslvn3 <- rownames(bsa_Ereslvn3)
bsa_Ereslvn3 <- ScaleData(bsa_Ereslvn3, features = all.genesbsa_Ereslvn3)
VlnPlot(object = bsa_Ereslvn3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)



mito.featuresbsa_Reslvd1 <- grep(pattern="^mt-", x=rownames(x=bsa_Reslvd1), value=T)
percent.mitobsa_Reslvd1 <- Matrix::colSums(x = GetAssayData(object = bsa_Reslvd1, slot = "counts")[mito.featuresbsa_Reslvd1,]) / Matrix::colSums(x = GetAssayData(object = bsa_Reslvd1, slot = "counts"))
bsa_Reslvd1[["percent.mito"]] <- percent.mitobsa_Reslvd1
VlnPlot(object = bsa_Reslvd1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

bsa_Reslvd1 <- subset(x = bsa_Reslvd1, subset = nFeature_RNA > 500 & nFeature_RNA <7000 & nCount_RNA > 500 & nCount_RNA < 50000 & percent.mitobsa_Reslvd1 < 0.1)
VlnPlot(object = bsa_Reslvd1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
bsa_Reslvd1 <- NormalizeData(object = bsa_Reslvd1, verbose = F)
bsa_Reslvd1 <- FindVariableFeatures(object = bsa_Reslvd1, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesbsa_Reslvd1 <- rownames(bsa_Reslvd1)
bsa_Reslvd1 <- ScaleData(bsa_Reslvd1, features = all.genesbsa_Reslvd1)
VlnPlot(object = bsa_Reslvd1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)



mito.featuresbsa_Reslvd2 <- grep(pattern="^mt-", x=rownames(x=bsa_Reslvd2), value=T)
percent.mitobsa_Reslvd2 <- Matrix::colSums(x = GetAssayData(object = bsa_Reslvd2, slot = "counts")[mito.featuresbsa_Reslvd2,]) / Matrix::colSums(x = GetAssayData(object = bsa_Reslvd2, slot = "counts"))
bsa_Reslvd2[["percent.mito"]] <- percent.mitobsa_Reslvd2
VlnPlot(object = bsa_Reslvd2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

bsa_Reslvd2 <- subset(x = bsa_Reslvd2, subset = nFeature_RNA > 500 & nFeature_RNA <7000 & nCount_RNA > 500 & nCount_RNA < 50000 & percent.mitobsa_Reslvd2 < 0.1)
VlnPlot(object = bsa_Reslvd2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
bsa_Reslvd2 <- NormalizeData(object = bsa_Reslvd2, verbose = F)
bsa_Reslvd2 <- FindVariableFeatures(object = bsa_Reslvd2, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesbsa_Reslvd2 <- rownames(bsa_Reslvd2)
bsa_Reslvd2 <- ScaleData(bsa_Reslvd2, features = all.genesbsa_Reslvd2)
VlnPlot(object = bsa_Reslvd2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)



mito.featuresbsa_Reslvd3 <- grep(pattern="^mt-", x=rownames(x=bsa_Reslvd3), value=T)
percent.mitobsa_Reslvd3 <- Matrix::colSums(x = GetAssayData(object = bsa_Reslvd3, slot = "counts")[mito.featuresbsa_Reslvd3,]) / Matrix::colSums(x = GetAssayData(object = bsa_Reslvd3, slot = "counts"))
bsa_Reslvd3[["percent.mito"]] <- percent.mitobsa_Reslvd3
VlnPlot(object = bsa_Reslvd3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

bsa_Reslvd3 <- subset(x = bsa_Reslvd3, subset = nFeature_RNA > 500 & nFeature_RNA <7000 & nCount_RNA > 500 & nCount_RNA < 50000 & percent.mitobsa_Reslvd3 < 0.1)
VlnPlot(object = bsa_Reslvd3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
bsa_Reslvd3 <- NormalizeData(object = bsa_Reslvd3, verbose = F)
bsa_Reslvd3 <- FindVariableFeatures(object = bsa_Reslvd3, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesbsa_Reslvd3 <- rownames(bsa_Reslvd3)
bsa_Reslvd3 <- ScaleData(bsa_Reslvd3, features = all.genesbsa_Reslvd3)
VlnPlot(object = bsa_Reslvd3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)



mito.featuresSTIA_c1 <- grep(pattern="^mt-", x=rownames(x=STIA_c1), value=T)
percent.mitoSTIA_c1 <- Matrix::colSums(x = GetAssayData(object = STIA_c1, slot = "counts")[mito.featuresSTIA_c1,]) / Matrix::colSums(x = GetAssayData(object = STIA_c1, slot = "counts"))
STIA_c1[["percent.mito"]] <- percent.mitoSTIA_c1
VlnPlot(object = STIA_c1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

STIA_c1 <- subset(x = STIA_c1, subset = nFeature_RNA > 500 & nFeature_RNA <5000 & nCount_RNA > 500 & nCount_RNA < 20000 & percent.mitoSTIA_c1 < 0.1)
VlnPlot(object = STIA_c1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
STIA_c1 <- NormalizeData(object = STIA_c1, verbose = F)
STIA_c1 <- FindVariableFeatures(object = STIA_c1, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesSTIA_c1 <- rownames(STIA_c1)
STIA_c1 <- ScaleData(STIA_c1, features = all.genesSTIA_c1)
VlnPlot(object = STIA_c1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)



mito.featuresSTIA_c2 <- grep(pattern="^mt-", x=rownames(x=STIA_c2), value=T)
percent.mitoSTIA_c2 <- Matrix::colSums(x = GetAssayData(object = STIA_c2, slot = "counts")[mito.featuresSTIA_c2,]) / Matrix::colSums(x = GetAssayData(object = STIA_c2, slot = "counts"))
STIA_c2[["percent.mito"]] <- percent.mitoSTIA_c2
VlnPlot(object = STIA_c2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

STIA_c2 <- subset(x = STIA_c2, subset = nFeature_RNA > 500 & nFeature_RNA <5000 & nCount_RNA > 500 & nCount_RNA < 15000 & percent.mitoSTIA_c2 < 0.1)
VlnPlot(object = STIA_c2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
STIA_c2 <- NormalizeData(object = STIA_c2, verbose = F)
STIA_c2 <- FindVariableFeatures(object = STIA_c2, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesSTIA_c2 <- rownames(STIA_c2)
STIA_c2 <- ScaleData(STIA_c2, features = all.genesSTIA_c2)
VlnPlot(object = STIA_c2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)



mito.featuresSTIA_c3 <- grep(pattern="^mt-", x=rownames(x=STIA_c3), value=T)
percent.mitoSTIA_c3 <- Matrix::colSums(x = GetAssayData(object = STIA_c3, slot = "counts")[mito.featuresSTIA_c3,]) / Matrix::colSums(x = GetAssayData(object = STIA_c3, slot = "counts"))
STIA_c3[["percent.mito"]] <- percent.mitoSTIA_c3
VlnPlot(object = STIA_c3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

STIA_c3 <- subset(x = STIA_c3, subset = nFeature_RNA > 500 & nFeature_RNA <4000 & nCount_RNA > 500 & nCount_RNA < 15000 & percent.mitoSTIA_c3 < 0.1)
VlnPlot(object = STIA_c3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
STIA_c3 <- NormalizeData(object = STIA_c3, verbose = F)
STIA_c3 <- FindVariableFeatures(object = STIA_c3, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesSTIA_c3 <- rownames(STIA_c3)
STIA_c3 <- ScaleData(STIA_c3, features = all.genesSTIA_c3)
VlnPlot(object = STIA_c3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)



mito.featuresSTIA_in1 <- grep(pattern="^mt-", x=rownames(x=STIA_in1), value=T)
percent.mitoSTIA_in1 <- Matrix::colSums(x = GetAssayData(object = STIA_in1, slot = "counts")[mito.featuresSTIA_in1,]) / Matrix::colSums(x = GetAssayData(object = STIA_in1, slot = "counts"))
STIA_in1[["percent.mito"]] <- percent.mitoSTIA_in1
VlnPlot(object = STIA_in1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

STIA_in1 <- subset(x = STIA_in1, subset = nFeature_RNA > 500 & nFeature_RNA <6000 & nCount_RNA > 500 & nCount_RNA < 40000 & percent.mitoSTIA_in1 < 0.1)
VlnPlot(object = STIA_in1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
STIA_in1 <- NormalizeData(object = STIA_in1, verbose = F)
STIA_in1 <- FindVariableFeatures(object = STIA_in1, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesSTIA_in1 <- rownames(STIA_in1)
STIA_in1 <- ScaleData(STIA_in1, features = all.genesSTIA_in1)
VlnPlot(object = STIA_in1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)



mito.featuresSTIA_in2 <- grep(pattern="^mt-", x=rownames(x=STIA_in2), value=T)
percent.mitoSTIA_in2 <- Matrix::colSums(x = GetAssayData(object = STIA_in2, slot = "counts")[mito.featuresSTIA_in2,]) / Matrix::colSums(x = GetAssayData(object = STIA_in2, slot = "counts"))
STIA_in2[["percent.mito"]] <- percent.mitoSTIA_in2
VlnPlot(object = STIA_in2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

STIA_in2 <- subset(x = STIA_in2, subset = nFeature_RNA > 500 & nFeature_RNA <5000 & nCount_RNA > 500 & nCount_RNA < 30000 & percent.mitoSTIA_in2 < 0.1)
VlnPlot(object = STIA_in2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
STIA_in2 <- NormalizeData(object = STIA_in2, verbose = F)
STIA_in2 <- FindVariableFeatures(object = STIA_in2, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesSTIA_in2 <- rownames(STIA_in2)
STIA_in2 <- ScaleData(STIA_in2, features = all.genesSTIA_in2)
VlnPlot(object = STIA_in2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)



mito.featuresSTIA_in3 <- grep(pattern="^mt-", x=rownames(x=STIA_in3), value=T)
percent.mitoSTIA_in3 <- Matrix::colSums(x = GetAssayData(object = STIA_in3, slot = "counts")[mito.featuresSTIA_in3,]) / Matrix::colSums(x = GetAssayData(object = STIA_in3, slot = "counts"))
STIA_in3[["percent.mito"]] <- percent.mitoSTIA_in3
VlnPlot(object = STIA_in3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

STIA_in3 <- subset(x = STIA_in3, subset = nFeature_RNA > 500 & nFeature_RNA <6000 & nCount_RNA > 500 & nCount_RNA < 35000 & percent.mitoSTIA_in3 < 0.1)
VlnPlot(object = STIA_in3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
STIA_in3 <- NormalizeData(object = STIA_in3, verbose = F)
STIA_in3 <- FindVariableFeatures(object = STIA_in3, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesSTIA_in3 <- rownames(STIA_in3)
STIA_in3 <- ScaleData(STIA_in3, features = all.genesSTIA_in3)
VlnPlot(object = STIA_in3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)



mito.featuresSTIA_Reslvn1 <- grep(pattern="^mt-", x=rownames(x=STIA_Reslvn1), value=T)
percent.mitoSTIA_Reslvn1 <- Matrix::colSums(x = GetAssayData(object = STIA_Reslvn1, slot = "counts")[mito.featuresSTIA_Reslvn1,]) / Matrix::colSums(x = GetAssayData(object = STIA_Reslvn1, slot = "counts"))
STIA_Reslvn1[["percent.mito"]] <- percent.mitoSTIA_Reslvn1
VlnPlot(object = STIA_Reslvn1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

STIA_Reslvn1 <- subset(x = STIA_Reslvn1, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & nCount_RNA > 500 & nCount_RNA < 25000 & percent.mitoSTIA_Reslvn1 < 0.1)
VlnPlot(object = STIA_Reslvn1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
STIA_Reslvn1 <- NormalizeData(object = STIA_Reslvn1, verbose = F)
STIA_Reslvn1 <- FindVariableFeatures(object = STIA_Reslvn1, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesSTIA_Reslvn1 <- rownames(STIA_Reslvn1)
STIA_Reslvn1 <- ScaleData(STIA_Reslvn1, features = all.genesSTIA_Reslvn1)
VlnPlot(object = STIA_Reslvn1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)



mito.featuresSTIA_Reslvn2 <- grep(pattern="^mt-", x=rownames(x=STIA_Reslvn2), value=T)
percent.mitoSTIA_Reslvn2 <- Matrix::colSums(x = GetAssayData(object = STIA_Reslvn2, slot = "counts")[mito.featuresSTIA_Reslvn2,]) / Matrix::colSums(x = GetAssayData(object = STIA_Reslvn2, slot = "counts"))
STIA_Reslvn2[["percent.mito"]] <- percent.mitoSTIA_Reslvn2
VlnPlot(object = STIA_Reslvn2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

STIA_Reslvn2 <- subset(x = STIA_Reslvn2, subset = nFeature_RNA > 500 & nFeature_RNA <5000 & nCount_RNA > 500 & nCount_RNA < 25000 & percent.mitoSTIA_Reslvn2 < 0.1)
VlnPlot(object = STIA_Reslvn2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
STIA_Reslvn2 <- NormalizeData(object = STIA_Reslvn2, verbose = F)
STIA_Reslvn2 <- FindVariableFeatures(object = STIA_Reslvn2, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesSTIA_Reslvn2 <- rownames(STIA_Reslvn2)
STIA_Reslvn2 <- ScaleData(STIA_Reslvn2, features = all.genesSTIA_Reslvn2)
VlnPlot(object = STIA_Reslvn2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)



mito.featuresSTIA_Reslvn3 <- grep(pattern="^mt-", x=rownames(x=STIA_Reslvn3), value=T)
percent.mitoSTIA_Reslvn3 <- Matrix::colSums(x = GetAssayData(object = STIA_Reslvn3, slot = "counts")[mito.featuresSTIA_Reslvn3,]) / Matrix::colSums(x = GetAssayData(object = STIA_Reslvn3, slot = "counts"))
STIA_Reslvn3[["percent.mito"]] <- percent.mitoSTIA_Reslvn3
VlnPlot(object = STIA_Reslvn3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

STIA_Reslvn3 <- subset(x = STIA_Reslvn3, subset = nFeature_RNA > 500 & nFeature_RNA <5000 & nCount_RNA > 200 & nCount_RNA < 25000 & percent.mitoSTIA_Reslvn3 < 0.1)
VlnPlot(object = STIA_Reslvn3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
STIA_Reslvn3 <- NormalizeData(object = STIA_Reslvn3, verbose = F)
STIA_Reslvn3 <- FindVariableFeatures(object = STIA_Reslvn3, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesSTIA_Reslvn3 <- rownames(STIA_Reslvn3)
STIA_Reslvn3 <- ScaleData(STIA_Reslvn3, features = all.genesSTIA_Reslvn3)
VlnPlot(object = STIA_Reslvn3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)



mito.featuresSTIA_Reslvd1 <- grep(pattern="^mt-", x=rownames(x=STIA_Reslvd1), value=T)
percent.mitoSTIA_Reslvd1 <- Matrix::colSums(x = GetAssayData(object = STIA_Reslvd1, slot = "counts")[mito.featuresSTIA_Reslvd1,]) / Matrix::colSums(x = GetAssayData(object = STIA_Reslvd1, slot = "counts"))
STIA_Reslvd1[["percent.mito"]] <- percent.mitoSTIA_Reslvd1
VlnPlot(object = STIA_Reslvd1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

STIA_Reslvd1 <- subset(x = STIA_Reslvd1, subset = nFeature_RNA > 500 & nFeature_RNA <3500 & nCount_RNA > 500 & nCount_RNA < 15000 & percent.mitoSTIA_Reslvd1 < 0.1)
VlnPlot(object = STIA_Reslvd1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
STIA_Reslvd1 <- NormalizeData(object = STIA_Reslvd1, verbose = F)
STIA_Reslvd1 <- FindVariableFeatures(object = STIA_Reslvd1, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesSTIA_Reslvd1 <- rownames(STIA_Reslvd1)
STIA_Reslvd1 <- ScaleData(STIA_Reslvd1, features = all.genesSTIA_Reslvd1)
VlnPlot(object = STIA_Reslvd1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)



mito.featuresSTIA_Reslvd2 <- grep(pattern="^mt-", x=rownames(x=STIA_Reslvd2), value=T)
percent.mitoSTIA_Reslvd2 <- Matrix::colSums(x = GetAssayData(object = STIA_Reslvd2, slot = "counts")[mito.featuresSTIA_Reslvd2,]) / Matrix::colSums(x = GetAssayData(object = STIA_Reslvd2, slot = "counts"))
STIA_Reslvd2[["percent.mito"]] <- percent.mitoSTIA_Reslvd2
VlnPlot(object = STIA_Reslvd2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

STIA_Reslvd2 <- subset(x = STIA_Reslvd2, subset = nFeature_RNA > 500 & nFeature_RNA <4000 & nCount_RNA > 500 & nCount_RNA < 15000 & percent.mitoSTIA_Reslvd2 < 0.1)
VlnPlot(object = STIA_Reslvd2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
STIA_Reslvd2 <- NormalizeData(object = STIA_Reslvd2, verbose = F)
STIA_Reslvd2 <- FindVariableFeatures(object = STIA_Reslvd2, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesSTIA_Reslvd2 <- rownames(STIA_Reslvd2)
STIA_Reslvd2 <- ScaleData(STIA_Reslvd2, features = all.genesSTIA_Reslvd2)
VlnPlot(object = STIA_Reslvd2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)



mito.featuresSTIA_Reslvd3 <- grep(pattern="^mt-", x=rownames(x=STIA_Reslvd3), value=T)
percent.mitoSTIA_Reslvd3 <- Matrix::colSums(x = GetAssayData(object = STIA_Reslvd3, slot = "counts")[mito.featuresSTIA_Reslvd3,]) / Matrix::colSums(x = GetAssayData(object = STIA_Reslvd3, slot = "counts"))
STIA_Reslvd3[["percent.mito"]] <- percent.mitoSTIA_Reslvd3
VlnPlot(object = STIA_Reslvd3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

STIA_Reslvd3 <- subset(x = STIA_Reslvd3, subset = nFeature_RNA > 500 & nFeature_RNA <3500 & nCount_RNA > 500 & nCount_RNA < 15000 & percent.mitoSTIA_Reslvd3 < 0.1)
VlnPlot(object = STIA_Reslvd3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
STIA_Reslvd3 <- NormalizeData(object = STIA_Reslvd3, verbose = F)
STIA_Reslvd3 <- FindVariableFeatures(object = STIA_Reslvd3, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesSTIA_Reslvd3 <- rownames(STIA_Reslvd3)
STIA_Reslvd3 <- ScaleData(STIA_Reslvd3, features = all.genesSTIA_Reslvd3)
VlnPlot(object = STIA_Reslvd3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)



mito.featuresCIA_c2 <- grep(pattern="^mt-", x=rownames(x=CIA_c2), value=T)
percent.mitoCIA_c2 <- Matrix::colSums(x = GetAssayData(object = CIA_c2, slot = "counts")[mito.featuresCIA_c2,]) / Matrix::colSums(x = GetAssayData(object = CIA_c2, slot = "counts"))
CIA_c2[["percent.mito"]] <- percent.mitoCIA_c2
VlnPlot(object = CIA_c2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

CIA_c2 <- subset(x = CIA_c2, subset = nFeature_RNA > 200 & nFeature_RNA <4000 & nCount_RNA > 200 & nCount_RNA < 25000 & percent.mitoCIA_c2 < 0.1)
VlnPlot(object = CIA_c2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
CIA_c2 <- NormalizeData(object = CIA_c2, verbose = F)
CIA_c2 <- FindVariableFeatures(object = CIA_c2, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesCIA_c2 <- rownames(CIA_c2)
CIA_c2 <- ScaleData(CIA_c2, features = all.genesCIA_c2)
VlnPlot(object = CIA_c2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)



mito.featuresCIA_c3 <- grep(pattern="^mt-", x=rownames(x=CIA_c3), value=T)
percent.mitoCIA_c3 <- Matrix::colSums(x = GetAssayData(object = CIA_c3, slot = "counts")[mito.featuresCIA_c3,]) / Matrix::colSums(x = GetAssayData(object = CIA_c3, slot = "counts"))
CIA_c3[["percent.mito"]] <- percent.mitoCIA_c3
VlnPlot(object = CIA_c3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

CIA_c3 <- subset(x = CIA_c3, subset = nFeature_RNA > 500 & nFeature_RNA <5000 & nCount_RNA > 500 & nCount_RNA < 25000 & percent.mitoCIA_c3 < 0.1)
VlnPlot(object = CIA_c3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
CIA_c3 <- NormalizeData(object = CIA_c3, verbose = F)
CIA_c3 <- FindVariableFeatures(object = CIA_c3, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesCIA_c3 <- rownames(CIA_c3)
CIA_c3 <- ScaleData(CIA_c3, features = all.genesCIA_c3)
VlnPlot(object = CIA_c3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)



mito.featuresCIA_c4 <- grep(pattern="^mt-", x=rownames(x=CIA_c4), value=T)
percent.mitoCIA_c4 <- Matrix::colSums(x = GetAssayData(object = CIA_c4, slot = "counts")[mito.featuresCIA_c4,]) / Matrix::colSums(x = GetAssayData(object = CIA_c4, slot = "counts"))
CIA_c4[["percent.mito"]] <- percent.mitoCIA_c4
VlnPlot(object = CIA_c4, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

CIA_c4 <- subset(x = CIA_c4, subset = nFeature_RNA > 500 & nFeature_RNA <4500 & nCount_RNA > 500 & nCount_RNA < 15000 & percent.mitoCIA_c4 < 0.1)
VlnPlot(object = CIA_c4, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
CIA_c4 <- NormalizeData(object = CIA_c4, verbose = F)
CIA_c4 <- FindVariableFeatures(object = CIA_c4, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesCIA_c4 <- rownames(CIA_c4)
CIA_c4 <- ScaleData(CIA_c4, features = all.genesCIA_c4)
VlnPlot(object = CIA_c4, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)



mito.featuresCIA_in1 <- grep(pattern="^mt-", x=rownames(x=CIA_in1), value=T)
percent.mitoCIA_in1 <- Matrix::colSums(x = GetAssayData(object = CIA_in1, slot = "counts")[mito.featuresCIA_in1,]) / Matrix::colSums(x = GetAssayData(object = CIA_in1, slot = "counts"))
CIA_in1[["percent.mito"]] <- percent.mitoCIA_in1
VlnPlot(object = CIA_in1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

CIA_in1 <- subset(x = CIA_in1, subset = nFeature_RNA > 500 & nFeature_RNA <7000 & nCount_RNA > 500 & nCount_RNA < 45000 & percent.mitoCIA_in1 < 0.1)
VlnPlot(object = CIA_in1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
CIA_in1 <- NormalizeData(object = CIA_in1, verbose = F)
CIA_in1 <- FindVariableFeatures(object = CIA_in1, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesCIA_in1 <- rownames(CIA_in1)
CIA_in1 <- ScaleData(CIA_in1, features = all.genesCIA_in1)
VlnPlot(object = CIA_in1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)



mito.featuresCIA_in2 <- grep(pattern="^mt-", x=rownames(x=CIA_in2), value=T)
percent.mitoCIA_in2 <- Matrix::colSums(x = GetAssayData(object = CIA_in2, slot = "counts")[mito.featuresCIA_in2,]) / Matrix::colSums(x = GetAssayData(object = CIA_in2, slot = "counts"))
CIA_in2[["percent.mito"]] <- percent.mitoCIA_in2
VlnPlot(object = CIA_in2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

CIA_in2 <- subset(x = CIA_in2, subset = nFeature_RNA > 500 & nFeature_RNA <6000 & nCount_RNA > 500 & nCount_RNA < 35000 & percent.mitoCIA_in2 < 0.1)
VlnPlot(object = CIA_in2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
CIA_in2 <- NormalizeData(object = CIA_in2, verbose = F)
CIA_in2 <- FindVariableFeatures(object = CIA_in2, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesCIA_in2 <- rownames(CIA_in2)
CIA_in2 <- ScaleData(CIA_in2, features = all.genesCIA_in2)
VlnPlot(object = CIA_in2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)



mito.featuresCIA_in3 <- grep(pattern="^mt-", x=rownames(x=CIA_in3), value=T)
percent.mitoCIA_in3 <- Matrix::colSums(x = GetAssayData(object = CIA_in3, slot = "counts")[mito.featuresCIA_in3,]) / Matrix::colSums(x = GetAssayData(object = CIA_in3, slot = "counts"))
CIA_in3[["percent.mito"]] <- percent.mitoCIA_in3
VlnPlot(object = CIA_in3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

CIA_in3 <- subset(x = CIA_in3, subset = nFeature_RNA > 500 & nFeature_RNA <5000 & nCount_RNA > 500 & nCount_RNA < 30000 & percent.mitoCIA_in3 < 0.1)
VlnPlot(object = CIA_in3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
CIA_in3 <- NormalizeData(object = CIA_in3, verbose = F)
CIA_in3 <- FindVariableFeatures(object = CIA_in3, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genesCIA_in3 <- rownames(CIA_in3)
CIA_in3 <- ScaleData(CIA_in3, features = all.genesCIA_in3)
VlnPlot(object = CIA_in3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

## Identification of integration anchors
c.reference.list <- c(bsa_c1, bsa_c2, bsa_c3, STIA_c1, STIA_c2, STIA_c3, CIA_c2, CIA_c3, CIA_c4)
#get list of all genes from objects
c.total.genes <- list(rownames(bsa_c1),
                    rownames(bsa_c2),
                    rownames(bsa_c3),
                    rownames(STIA_c1),
                    rownames(STIA_c2),
                    rownames(STIA_c3),
                    rownames(CIA_c2),
                    rownames(CIA_c3),
                    rownames(CIA_c4))

#find common genes (genes appearing in only one object not used)
c.common.genes <- Reduce(f = intersect, x = c.total.genes)
#check numbers
length(c.common.genes)
nrow(bsa_c1)
nrow(bsa_c2)
nrow(bsa_c3)
nrow(STIA_c1)
nrow(STIA_c2)
nrow(STIA_c3)
nrow(CIA_c2)
nrow(CIA_c3)
nrow(CIA_c4)

c.anchors <- FindIntegrationAnchors(object.list = c.reference.list, anchor.features = c.common.genes, dims = 1:30)
c.integrated <- IntegrateData(anchorset = c.anchors, dims = 1:30, features.to.integrate = c.common.genes)


in.reference.list <- c(bsa_In1, bsa_In2, bsa_In3, CIA_in1,CIA_in2,CIA_in3,STIA_in1, STIA_in2,STIA_in3)
#get list of all genes from objects
in.total.genes <- list(rownames(bsa_In1),
                       rownames(bsa_In2),
                       rownames(bsa_In3),
                       rownames(CIA_in1),
                       rownames(CIA_in2),
                       rownames(CIA_in3),
                       rownames(STIA_in1),
                       rownames(STIA_in2),
                       rownames(STIA_in3))

#find common genes (genes appearing in only one object not used)
in.common.genes <- Reduce(f = intersect, x = in.total.genes)
#check numbers
length(in.common.genes)
nrow(bsa_In1)
nrow(bsa_In2)
nrow(bsa_In3)
nrow(CIA_in1)
nrow(CIA_in2)
nrow(CIA_in3)
nrow(STIA_in1)
nrow(STIA_in2)
nrow(STIA_in3)

in.anchors <- FindIntegrationAnchors(object.list = in.reference.list, anchor.features = in.common.genes, dims = 1:30)
in.integrated <- IntegrateData(anchorset = in.anchors, dims = 1:30, features.to.integrate = in.common.genes)


rsln.reference.list <- c(bsa_Ereslvn1, bsa_Ereslvn2, bsa_Ereslvn3, bsa_Reslvn1,bsa_Reslvn2,bsa_Reslvn3,STIA_Reslvn1, STIA_Reslvn2,STIA_Reslvn3)
#get list of all genes from objects
rsln.total.genes <- list(rownames(bsa_Ereslvn1),
                       rownames(bsa_Ereslvn2),
                       rownames(bsa_Ereslvn3),
                       rownames(bsa_Reslvn1),
                       rownames(bsa_Reslvn2),
                       rownames(bsa_Reslvn3),
                       rownames(STIA_Reslvn1),
                       rownames(STIA_Reslvn2),
                       rownames(STIA_Reslvn3))

#find common genes (genes appearing in only one object not used)
rsln.common.genes <- Reduce(f = intersect, x = rsln.total.genes)
#check numbers
length(rsln.common.genes)
nrow(bsa_Ereslvn1)
nrow(bsa_Ereslvn2)
nrow(bsa_Ereslvn3)
nrow(bsa_Reslvn1)
nrow(bsa_Reslvn2)
nrow(bsa_Reslvn3)
nrow(STIA_Reslvn1)
nrow(STIA_Reslvn2)
nrow(STIA_Reslvn3)

rsln.anchors <- FindIntegrationAnchors(object.list = rsln.reference.list, anchor.features = rsln.common.genes, dims = 1:30)
rsln.integrated <- IntegrateData(anchorset = rsln.anchors, dims = 1:30, features.to.integrate = rsln.common.genes)


rsld.reference.list <- c(bsa_Reslvd1,bsa_Reslvd2,bsa_Reslvd3,STIA_Reslvd1, STIA_Reslvd2,STIA_Reslvd3)
#get list of all genes from objects
rsld.total.genes <- list(rownames(bsa_Reslvd1),
                         rownames(bsa_Reslvd2),
                         rownames(bsa_Reslvd3),
                         rownames(STIA_Reslvd1),
                         rownames(STIA_Reslvd2),
                         rownames(STIA_Reslvd3))

#find common genes (genes appearing in only one object not used)
rsld.common.genes <- Reduce(f = intersect, x = rsld.total.genes)
#check numbers
length(rsld.common.genes)
nrow(bsa_Reslvd1)
nrow(bsa_Reslvd2)
nrow(bsa_Reslvd3)
nrow(STIA_Reslvd1)
nrow(STIA_Reslvd2)
nrow(STIA_Reslvd3)

rsld.anchors <- FindIntegrationAnchors(object.list = rsld.reference.list, anchor.features = rsld.common.genes, dims = 1:30)
rsld.integrated <- IntegrateData(anchorset = rsld.anchors, dims = 1:30, features.to.integrate = rsld.common.genes)



total.genes_int <- list(rownames(c.integrated),
                        rownames(in.integrated),
                        rownames(rsln.integrated),
                        rownames(rsld.integrated))

#find common genes (genes appearing in only one object not used)
common.genes_int <- Reduce(f = intersect, x = total.genes_int)
#check numbers
length(common.genes_int)
nrow(c.integrated)
nrow(in.integrated)
nrow(rsln.integrated)
nrow(rsld.integrated)


#integrate all
#rm(integrated_all)
reference.list_all <- c(c.integrated,in.integrated,rsln.integrated, rsld.integrated)
anchors_all <- FindIntegrationAnchors(object.list = reference.list_all, anchor.features = 2000, dims = 1:30)
integrated_all <- IntegrateData(anchorset = anchors_all, dims = 1:30)


DefaultAssay(object=integrated_all) <- "integrated"
integrated_all <- ScaleData(object = integrated_all, verbose=F)
integrated_all <- RunPCA(object = integrated_all, verbose=F, npcs = 50)
ElbowPlot(object = integrated_all, ndims = 50)
integrated_all <- SCTransform(integrated_all)
DefaultAssay(object=integrated_all) <- "integrated"
integrated_all <- FindNeighbors(object = integrated_all, dims = 1:30)

integrated_all <- RunUMAP(object = integrated_all, reduction = "pca", dims = 1:40)
DimPlot(integrated_all, reduction = "umap", label = TRUE, pt.size = 0.1, repel = F, group.by = "sample_id") + NoLegend() + NoAxes()

DefaultAssay(integrated_all) <- "RNA"

FeaturePlot(integrated_all, features = "Ptprc", cols = c("grey", "yellow", "red")) + NoAxes()

DefaultAssay(object=integrated_all) <- "integrated"
integrated_all <- FindVariableFeatures(integrated_all, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(integrated_all)
integrated_all <- ScaleData(integrated_all, features = all.genes)
integrated_all <- SCTransform(integrated_all)
DefaultAssay(object=integrated_all) <- "integrated"
integrated_all <- RunPCA(integrated_all, npcs = 50, features = VariableFeatures(object = integrated_all))
ElbowPlot(integrated_all,  ndims = 50)
integrated_all <- FindNeighbors(integrated_all, reduction = "pca", dims = 1:30)
# to find out graph.name integrated.H@graphs
integrated_all <- FindClusters(integrated_all, graph.name = "integrated_snn", resolution = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7))
clustree(integrated_all, assay = "integrated")
integrated_all <- FindClusters(integrated_all, graph.name = "integrated_snn", resolution = 0.1)
integrated_all <- RunUMAP(integrated_all, reduction = "pca", dims = 1:50)
DimPlot(integrated_all, reduction = "umap", label = TRUE, pt.size = 0.1, repel = T, label.box = T) + NoAxes()
DimPlot(integrated_all, reduction = "umap", label = F, pt.size = 0.5, repel = T, label.box = F, split.by = "model") + NoLegend() + NoAxes()

###downsample  re-preprocess
head(integrated_all)
Idents(integrated_all) <- "orig.ident"
integrated.5k <- subset(integrated_all, downsample = 5000)
DimPlot(integrated.5k, reduction = "umap", label = F, pt.size = 0.5, repel = T, label.box = F, split.by = "model") + NoLegend() + NoAxes()

DefaultAssay(object=integrated.5k) <- "integrated"
integrated.5k <- FindVariableFeatures(integrated.5k, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(integrated.5k)
integrated.5k <- ScaleData(integrated.5k, features = all.genes)
integrated.5k <- SCTransform(integrated.5k)
DefaultAssay(object=integrated.5k) <- "integrated"
integrated.5k <- RunPCA(integrated.5k, npcs = 50, features = VariableFeatures(object = integrated.5k))
ElbowPlot(integrated.5k,  ndims = 50)
integrated.5k <- FindNeighbors(integrated.5k, reduction = "pca", dims = 1:30)
# to find out graph.name integrated.H@graphs
integrated.5k <- FindClusters(integrated.5k, graph.name = "integrated_snn", resolution = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7))
clustree(integrated.5k, assay = "integrated")
integrated.5k <- FindClusters(integrated.5k, graph.name = "integrated_snn", resolution = 0.1)
integrated.5k <- RunUMAP(integrated.5k, reduction = "pca", dims = 1:30)
DimPlot(integrated.5k, reduction = "umap", label = TRUE, pt.size = 0.1, repel = T, label.box = T) + NoAxes()
DimPlot(integrated.5k, reduction = "umap", label = F, pt.size = 0.1, repel = T, label.box = F, split.by = "model") + NoLegend() + NoAxes()

getwd()
saveRDS(integrated.5k, file = "integrated.5k.rds")
integrated.5k <- readRDS("integrated.5k.rds")


cluster0 <- FindMarkers(integrated.5k, ident.1 = 0, min.pct = 0.5, only.pos = T)
cluster1 <- FindMarkers(integrated.5k, ident.1 = 1, min.pct = 0.5, only.pos = T)
cluster2 <- FindMarkers(integrated.5k, ident.1 = 2, min.pct = 0.5, only.pos = T)
cluster3 <- FindMarkers(integrated.5k, ident.1 = 3, min.pct = 0.5, only.pos = T)
cluster4 <- FindMarkers(integrated.5k, ident.1 = 4, min.pct = 0.5, only.pos = T)
cluster5 <- FindMarkers(integrated.5k, ident.1 = 5, min.pct = 0.5, only.pos = T)
cluster6 <- FindMarkers(integrated.5k, ident.1 = 6, min.pct = 0.5, only.pos = T)
cluster7 <- FindMarkers(integrated.5k, ident.1 = 7, min.pct = 0.5, only.pos = T)
cluster8 <- FindMarkers(integrated.5k, ident.1 = 8, min.pct = 0.5, only.pos = T)
cluster9 <- FindMarkers(integrated.5k, ident.1 = 9, min.pct = 0.5, only.pos = T)
cluster10 <- FindMarkers(integrated.5k, ident.1 = 10, min.pct = 0.5, only.pos = T)
cluster11 <- FindMarkers(integrated.5k, ident.1 = 11, min.pct = 0.5, only.pos = T)
cluster12 <- FindMarkers(integrated.5k, ident.1 = 12, min.pct = 0.5, only.pos = T)
cluster13 <- FindMarkers(integrated.5k, ident.1 = 13, min.pct = 0.5, only.pos = T)

new.cluster.ids <- c("Fibroblasts", "Fibroblasts", "Fibroblasts", "Fibroblasts", "EC", "Mural", "Fibroblasts", "Cycling", "Fibroblasts", "Myocytes",
                     "Chondrocytes", "Lymphatics", "Glial", "Fibroblasts") 
names(new.cluster.ids) <- levels(integrated.5k)
integrated.5k <- RenameIdents(integrated.5k, new.cluster.ids)
DimPlot(integrated.5k, reduction = "umap", label = TRUE, pt.size = 0.1, repel = T, label.box = T, cols = "Accent") + NoAxes() +NoLegend()

integrated.5k <- subset(integrated.5k, idents = grep(c("Fibroblasts|Cycling|Mural|EC|Lymphatics"),
                                                           levels(integrated.5k@active.ident), value = T))

DefaultAssay(object=integrated.5k) <- "integrated"

integrated.5k <- FindNeighbors(integrated.5k, reduction = "pca", dims = 1:30)
# to find out graph.name integrated.H@graphs
integrated.5k <- FindClusters(integrated.5k, graph.name = "integrated_snn", resolution = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7))
clustree(integrated.5k, assay = "integrated")

integrated.5k <- FindClusters(integrated.5k, graph.name = "integrated_snn", resolution = 0.2)
integrated.5k <- RunUMAP(integrated.5k, reduction = "pca", dims = 1:30)
DimPlot(integrated.5k, reduction = "umap", label = TRUE, pt.size = 0.1, repel = T, label.box = T) + NoAxes()

#######check stormal cell gene expression
FindAllMarkers_stromal <- FindAllMarkers(integrated.5k, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.5)

#####subset Fibroblasts

FB.integrated.5k.time<-subset(integrated.5k, idents="Fibroblasts")

DefaultAssay(object=FB.integrated.5k.time) <- "integrated"
FB.integrated.5k.time <- FindVariableFeatures(FB.integrated.5k.time, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(FB.integrated.5k.time)
FB.integrated.5k.time <- ScaleData(FB.integrated.5k.time, features = all.genes)
FB.integrated.5k.time <- SCTransform(FB.integrated.5k.time)
DefaultAssay(object=FB.integrated.5k.time) <- "integrated"
FB.integrated.5k.time <- RunPCA(FB.integrated.5k.time, npcs = 50, features = VariableFeatures(object = FB.integrated.5k.time))
ElbowPlot(FB.integrated.5k.time,  ndims = 50)
FB.integrated.5k.time <- FindNeighbors(FB.integrated.5k.time, reduction = "pca", dims = 1:30)
# to find out graph.name integrated.H@graphs
FB.integrated.5k.time <- FindClusters(FB.integrated.5k.time, graph.name = "integrated_snn", resolution = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7))
clustree(FB.integrated.5k.time, assay = "integrated")
FB.integrated.5k.time <- FindClusters(FB.integrated.5k.time, graph.name = "integrated_snn", resolution = 0.1)
FB.integrated.5k.time <- RunUMAP(FB.integrated.5k.time, reduction = "pca", dims = 1:30)
DimPlot(FB.integrated.5k.time, reduction = "umap", label = F, pt.size = 0.1, repel = T, label.box = T, cols = "Set2", split.by = "InflammationStatus") + NoAxes() +NoLegend()

FB.integrated.5k.time <- subset(FB.integrated.5k.time, subset = InflammationStatus == c("Control", "Inflamed", "Resolving", "Resolved"))
new.cluster.ids <- c("0", "1", "2", "3", "4", "1")
names(new.cluster.ids) <- levels(FB.integrated.5k.time)
FB.integrated.5k.time <- RenameIdents(FB.integrated.5k.time, new.cluster.ids)

########Extended data figure 6d
DimPlot(FB.integrated.5k.time, reduction = "umap", label = T, pt.size = 0.1, repel = T, label.box = T, cols = "Set2", split.by = "InflammationStatus") + NoLegend() + NoAxes()



#######pseudobulk fibroblasts
FB.integrated.5k.time <- FindClusters(FB.integrated.5k.time, graph.name = "integrated_snn", resolution = 0)
DimPlot(FB.integrated.5k, reduction = "umap", label = T, pt.size = 0.1, repel = T, label.box = T, cols = "Set2") + NoLegend() + NoAxes()


###Pi16

pi16p<-read.csv(file = "Pi16p.csv")
pi16c<-read.csv(file = "Pi16c.csv")
pi16p.genes <- list(pi16p$X)
pi16c.genes <- list(pi16c$X)

FB.integrated.5k.time <- AddModuleScore(
  object = FB.integrated.5k.time,
  features = pi16p.genes,
  ctrl = 5,
  name = 'Pi16_PBS'
)

FB.integrated.5k.time <- AddModuleScore(
  object = FB.integrated.5k.time,
  features = pi16c.genes,
  ctrl = 5,
  name = 'Pi16_CAR'
)

head(x = FB.integrated.5k.time[])

comp <- list(c("Control", "Inflamed"), c("Control", "Resolving"), c("Control", "Resolved"))
comp <- list(c("Control", "Inflamed", "Resolving", "Resolved"))
comp <- list(c("Inflamed", "Control"), c("Inflamed", "Resolving"), c("Inflamed", "Resolved"))

############figure 4g
VlnPlot(FB.integrated.5k.time, features = "Pi16_PBS1", adjust = 1, pt.size = 0,split.by = "InflammationStatus")+NoLegend()+xlab(NULL)+stat_compare_means(label = "p.signif", comparisons = comp, method = "t.test")+ylim(-0.5, 1.2)+scale_fill_manual(values = c("grey", "grey", "grey", "grey"))
VlnPlot(FB.integrated.5k.time, features = "Pi16_CAR1", adjust = 1, pt.size = 0,split.by = "InflammationStatus")+NoLegend()+xlab(NULL)+stat_compare_means(label = "p.signif", comparisons = comp, method = "t.test")+ylim(-0.5, 1.2)+scale_fill_manual(values = c("firebrick3", "firebrick3", "firebrick3", "firebrick3"))

############extended data figure 7e
plot_density(FB.integrated.5k.time, features = "Pi16_PBS1", combine = F, method = "wkde", size = 1, adjust = 1)+NoAxes()+ scale_color_gradientn(colours = magma(20))+facet_grid(.~FB.integrated.5k.time$InflammationStatus)
plot_density(FB.integrated.5k.time, features = "Pi16_CAR1", combine = F, method = "wkde", size = 1, adjust = 1)+NoAxes()+ scale_color_gradientn(colours = magma(20))+facet_grid(.~FB.integrated.5k.time$InflammationStatus)

############extended data figure 7f
plot_density(FB.integrated.5k.time, features = "Pi16", combine = F, method = "wkde", size = 1, adjust = 1)+NoAxes()+ scale_color_gradientn(colours = magma(20))+facet_grid(.~FB.integrated.5k.time$InflammationStatus)

############extended data figure 7g
plot_density(FB.integrated.5k.time, features = "Cxcl5", combine = F, method = "wkde", size = 1, adjust = 1)+NoAxes()+ scale_color_gradientn(colours = magma(20))+facet_grid(.~FB.integrated.5k.time$InflammationStatus)



DotPlot(
  FB.integrated.5k.time,
  features = c("Pi16_CAR1", "Pi16_PBS1"),
  assay = NULL,
  dot.scale = 6,
  scale.min = 0
)+
  scale_color_gradientn(
    colors = c("lightgrey", "skyblue", "red"),
    values = scales::rescale(c(-1, 0, 1))
  )+RotatedAxis()

DefaultAssay(Inflamed)<-"RNA"


##########fibroblast3
fib3c<-read.csv(file = "fib3c.csv")
fib3c.genes <- list(fib3c$X)

FB.integrated.5k.time <- AddModuleScore(
  object = FB.integrated.5k.time,
  features = fib3c.genes,
  ctrl = 5,
  name = 'fib3'
)

############extended data figure 8d
VlnPlot(FB.integrated.5k.time, features = "fib31", adjust = 1, pt.size = 0,split.by = "InflammationStatus")+NoLegend()+xlab(NULL)+stat_compare_means(label = "p.signif", comparisons = comp, method = "t.test")+ylim(-0.5, 1.2)+scale_fill_manual(values = c("olivedrab", "olivedrab", "olivedrab", "olivedrab"))


