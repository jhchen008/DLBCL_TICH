
###### load required packages ######
library(ggplot2)
library(ggridges)
library(ggpubr)
library(ggthemes)
library(cowplot)
library(dplyr)
library(plyr)
library(reshape2)
library(dittoSeq)
library(RColorBrewer)
library(grDevices)
library(viridis)
library(scDataviz)
library(GSEABase)
library(GGally)

library(scran)
library(Seurat)
library(SeuratWrappers)
library(scCustomize)
library(umap)
library(harmony)
library(SingleR)
library(celldex)

library(UCell)
library(AUCell)
library(CellChat)


##############################################################################################################################
# Human CH dataset by Ben-Crentsil et al.
##############################################################################################################################

###### seurat object and QC ######
scech_all <- readRDS(file = "GSE210433_seurat_object.rds") # obtain from GEO database under accession number GSE210433

# subset samples, TET2 vs CH-
scech = subset(scech_all, group %in% c("COVID (no CHIP)","TET2 CHIP + COVID"))

# filter
scech = subset(scech, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA < 35000 & nCount_RNA > 1000)

# normalization 
scech <- NormalizeData(scech)

# dimplot, use umap and cell type annotations from the original study
DimPlot(scech, reduction = "umap", label = TRUE,  group.by = "celltype")



###### GSEA analysis ######

counts <- LayerData(scech, assay = "RNA", layer = "counts")
expressed_gene = rownames(counts)
rm(counts)

gset_cur <- getGmt("gene_sets.gmt")
gset_cur <- subsetGeneSets(gset_cur, expressed_gene) 
scech <- AddModuleScore_UCell(scech, features=geneIds(gset_cur), name=NULL)



###### differential expression analysis ######

# TET2 vs CH- CD14+ Mono
marker_CD14_TETvsWT = FindMarkers(subset(scech,celltype %in% c("CD14 Mono")), ident.1 = "TET2 CHIP + COVID", ident.2 = "COVID (no CHIP)", group.by = 'group')

# TET2 vs CH- Int Mono
marker_IntMono_TETvsWT = FindMarkers(subset(scech,celltype %in% c("Int Mono")), ident.1 = "TET2 CHIP + COVID", ident.2 = "COVID (no CHIP)", group.by = 'group')



###### FeaturePlot ######
# define a cut-off for low expressor with na_cutoff
FeaturePlot_scCustom(scech, colors_use = viridis_magma_dark_high, features = c("DLBCL_CH_Up","DLBCL_CH_Dn","DLBCL_CH_Up_CAPS"),na_cutoff = 0.01, num_columns = 3)
FeaturePlot_scCustom(scech, colors_use = viridis_magma_dark_high, features = c("CD163","CCL8","S100A12","LILRA5"),na_cutoff = 0.01, num_columns = 3)



###### boxplot for gene set scores ######

metadata = scech@meta.data
factor(metadata$celltype)


var_list <- c("group","DLBCL_CH_Up","DLBCL_CH_Dn","DLBCL_CH_Up_CAPS")

## CD14 Mono
selected_var <- subset(metadata,celltype %in% c("CD14 Mono"))[var_list]
long_data <- melt(as.data.frame(selected_var), id.vars="group", variable.name = "Variable", value.name = "Score")

ggplot(subset(long_data,group != "NA"), aes(x = group, y = Score, fill = group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(shape = 16, position = position_jitter(0.2), size = 0, alpha = 0.6) +
    facet_wrap(~Variable, scales = "free") +
    facet_wrap(~Variable, scales = "free", ncol = 5) +
    theme_pubr() +
    labs(title = "",
         x = "",
         y = "Z-score") +
    scale_fill_manual(values = c("#96969680","#B2DF8A70","#5757F970"))


## Int Mono
selected_var <- subset(metadata,celltype %in% c("Int Mono"))[var_list]
long_data <- melt(as.data.frame(selected_var), id.vars="group", variable.name = "Variable", value.name = "Score")

ggplot(subset(long_data,group != "NA"), aes(x = group, y = Score, fill = group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(shape = 16, position = position_jitter(0.2), size = 0, alpha = 0.6) +
    facet_wrap(~Variable, scales = "free") +
    facet_wrap(~Variable, scales = "free", ncol = 5) +
    theme_pubr() +
    labs(title = "",
         x = "",
         y = "Z-score") +
    scale_fill_manual(values = c("#96969680","#B2DF8A70","#5757F970"))

t.test(subset(metadata, celltype %in% c("CD14 Mono") & group == "TET2 CHIP + COVID")$DLBCL_CH_Up_CAPS, subset(metadata, celltype %in% c("CD14 Mono") & group == "COVID (no CHIP)")$DLBCL_CH_Up_CAPS)
# p-value < 2.2e-16
t.test(subset(metadata, celltype %in% c("CD14 Mono") & group == "TET2 CHIP + COVID")$DLBCL_CH_Up, subset(metadata, celltype %in% c("CD14 Mono") & group == "COVID (no CHIP)")$DLBCL_CH_Up)
# p-value = 3.642e-10
t.test(subset(metadata, celltype %in% c("Int Mono") & group == "TET2 CHIP + COVID")$DLBCL_CH_Up_CAPS, subset(metadata, celltype %in% c("Int Mono") & group == "COVID (no CHIP)")$DLBCL_CH_Up_CAPS)
# p-value = 0.001355
t.test(subset(metadata, celltype %in% c("Int Mono") & group == "TET2 CHIP + COVID")$DLBCL_CH_Up, subset(metadata, celltype %in% c("Int Mono") & group == "COVID (no CHIP)")$DLBCL_CH_Up)
# p-value = 5.649e-05

###### cellChat differential analysis ######

## cellChat object
sce_tet2ch = subset(scech,  group == "TET2 CHIP + COVID" & celltype %in% c("CD4 T","Treg","CD8 T","gdT","MAIT","NK","B","CD14 Mono","CD16 Mono","Int Mono","cDC","HSPC","Neutrophils"))
sce_control = subset(scech,  group == "COVID (no CHIP)" & celltype %in% c("CD4 T","Treg","CD8 T","gdT","MAIT","NK","B","CD14 Mono","CD16 Mono","Int Mono","cDC","HSPC","Neutrophils"))

sce_tet2ch$celltype = factor(sce_tet2ch$celltype, levels=c("CD4 T","Treg","CD8 T","gdT","MAIT","NK","B","CD14 Mono","CD16 Mono","Int Mono","cDC","HSPC","Neutrophils"))
sce_control$celltype = factor(sce_control$celltype, levels=c("CD4 T","Treg","CD8 T","gdT","MAIT","NK","B","CD14 Mono","CD16 Mono","Int Mono","cDC","HSPC","Neutrophils"))

cellchat_tet2ch <- createCellChat(sce_tet2ch, group.by = "celltype", assay = "RNA")
cellchat_control <- createCellChat(sce_control, group.by = "celltype", assay = "RNA")

# cellchat Db
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation")

# run cellchat for each sample, cellchat_tet2ch
cellchat_tet2ch@DB <- CellChatDB.use
cellchat_tet2ch <- subsetData(cellchat_tet2ch) 
future::plan("multisession", workers = 4)
cellchat_tet2ch <- identifyOverExpressedGenes(cellchat_tet2ch)
cellchat_tet2ch <- identifyOverExpressedInteractions(cellchat_tet2ch)
cellchat_tet2ch <- computeCommunProb(cellchat_tet2ch, type = "triMean")
cellchat_tet2ch <- filterCommunication(cellchat_tet2ch, min.cells = 10)
cellchat_tet2ch <- computeCommunProbPathway(cellchat_tet2ch)
cellchat_tet2ch <- aggregateNet(cellchat_tet2ch)

# significant interaction
cellchat_tet2ch@netP$pathways

# run cellchat for each sample, cellchat_control

# set the used database in the object
cellchat_control@DB <- CellChatDB.use
cellchat_control <- subsetData(cellchat_control)
future::plan("multisession", workers = 4)
cellchat_control <- identifyOverExpressedGenes(cellchat_control)
cellchat_control <- identifyOverExpressedInteractions(cellchat_control)
cellchat_control <- computeCommunProb(cellchat_control, type = "triMean")
cellchat_control <- filterCommunication(cellchat_control, min.cells = 10)
cellchat_control <- computeCommunProbPathway(cellchat_control)
cellchat_control <- aggregateNet(cellchat_control)

# significant interaction
cellchat_control@netP$pathways


# Merge analysis
object.list <- list(CH_neg = cellchat_control, CH_TET2 = cellchat_tet2ch)
cellchat_com <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = T)

# circle plot 
netVisual_diffInteraction(cellchat_com, weight.scale = T, measure = "weight")

# rankNet plot
gg1 <- rankNet(cellchat_com, mode = "comparison", stacked = T, do.stat = TRUE, color.use = c("#96969680","#B2DF8A70"))
gg2 <- rankNet(cellchat_com, mode = "comparison", stacked = F, do.stat = TRUE, color.use = c("#96969680","#B2DF8A70"))
gg1$data ## values in the plots



##############################################################################################################################
# Mouse CH dataset by Rauch et al.
##############################################################################################################################

###### seurat object and QC ######

# load cellranger files downloaded from GEO, and create seurat object
d3a_wt = Read10X(data.dir="./GSM7056033_30week_dnmt3a_45_1")
d3a_ko = Read10X(data.dir="./GSM7056034_30week_dnmt3a_45_2")
tet_wt = Read10X(data.dir="./GSM7056035_30week_tet2_45_1")
tet_ko = Read10X(data.dir="./GSM7056036_30week_tet2_45_2")

# create seurat object
count_list <- list(
  d3a_wt = d3a_wt,
  d3a_ko = d3a_ko,
  tet_wt = tet_wt,
  tet_ko = tet_ko

)

seurat_objects <- list()
seurat_objects <- lapply(names(count_list), function(sample_name) {
    CreateSeuratObject(
        counts = count_list[[sample_name]],
        project = sample_name,
        min.cells = 5,
        min.features = 500
    )
})

names(seurat_objects) <- names(count_list)
original_names <- names(seurat_objects)


# QC
## add MT percentage
seurat_objects <- lapply(seurat_objects, function(object) {
  object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
  return(object)
})

# filter
seurat_objects <- lapply(seurat_objects, function(object) {
  object = subset(object, subset = nFeature_RNA > 200 & nFeature_RNA < 4400 & percent.mt < 15)
  return(object)
})


# add meta.data
seurat_objects <- lapply(seurat_objects, function(object) {
  sample_name = object@project.name
  object$sample <- sample_name
  object$group <- ifelse(grepl("wt", object@project.name, fixed = FALSE),"wt",sample_name)
  return(object)
})


# normalization
seurat_objects <- lapply(seurat_objects, function(object) {
  object = NormalizeData(object)
  return(object)
})


# merge objects of different samples

# add sample identity
seurat_objects <- lapply(names(seurat_objects), function(i) {
  RenameCells(seurat_objects[[i]], add.cell.id = i)
})
names(seurat_objects) <- original_names

smerge <- Reduce(function(x, y) merge(x, y), seurat_objects)
rm(seurat_objects)

###### UMAP clustering ######

# clustering, with harmony correction
smerge = FindVariableFeatures(smerge, selection.method = "vst")
smerge = ScaleData(smerge, verbose = FALSE)
smerge = RunPCA(smerge, verbose = FALSE)
ElbowPlot(object = smerge, ndims = 40)

smerge = RunHarmony(smerge, "sample")
smerge = FindNeighbors(smerge, reduction = "harmony", dims = 1:15)
smerge = FindClusters(smerge)
smerge = FindClusters(smerge,resolution =0.6)
smerge = RunUMAP(smerge, reduction = "harmony", dims = 1:15)



###### SingleR annotation ######

# v5 is not supported, needs to be converted to v3 before SingleR
smerge[["RNA"]] <- as(smerge[["RNA"]], Class="Assay")

# main annotation
ref <- fetchReference("immgen", "2024-02-26")
pred.immgen <- SingleR(test = smerge@assays$RNA@data, ref = ref, labels = ref$label.main)
smerge@meta.data$CellType.immgen <- pred.immgen$labels


# cell types across seurat_clusters
DimPlot(smerge, reduction = "umap", pt.size = 0.5, group.by="CellType.immgen", label = TRUE, label.size = 4)
table(smerge$seurat_clusters, smerge$CellType.immgen)

# add cell types information based on SingleR results
Idents(smerge) = smerge$seurat_clusters
smerge <- RenameIdents(object = smerge,
'0'='Inflam_Mac',
'1'='B',
'2'='T',
'3'='Resident-like_Mac',
'4'='Mono',
'5'='T',
'6'='B',
'7'='TREM2hi_Mac',
'8'='DC',
'9'='NKT',
'10'='Mixed_lymphocytes',
'11'='Resident-like_Mac',
'12'='DC',
'13'='Mac',
'14'='T',
'15'='B',
'16'='Fibroblasts',
'17'='DC'
)
smerge$celltype = Idents(smerge)
Idents(smerge) = smerge$seurat_clusters

# UMAP by cell types
DimPlot(smerge, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 4, group.by = "celltype")



###### GSEA analysis ######

counts <- LayerData(smerge, assay = "RNA", layer = "counts")
expressed_gene = rownames(counts)
rm(counts)

gset_cur <- getGmt("gene_sets.mus.gmt")
gset_cur <- subsetGeneSets(gset_cur, expressed_gene) 
smerge <- AddModuleScore_UCell(smerge, features=geneIds(gset_cur), name=NULL)



###### differential expression analysis ######

# D3A-KO vs WT Resident-like_Mac
marker_TRM_D3AvsWT = FindMarkers(subset(smerge,celltype %in% c("Resident-like_Mac")), ident.1 = "d3a_ko", ident.2 = "wt", group.by = 'group')

# TET2-KO vs WT Resident-like_Mac
marker_TRM_TETvsWT = FindMarkers(subset(smerge,celltype %in% c("Resident-like_Mac")), ident.1 = "tet_ko", ident.2 = "wt", group.by = 'group')


###### FeaturePlot ######
FeaturePlot_scCustom(scech, colors_use = viridis_magma_dark_high, features = c("DLBCL_CH_Up","DLBCL_CH_Dn","DLBCL_CH_Up_CAPS"),na_cutoff = 0.01, num_columns = 3)
FeaturePlot_scCustom(scech, colors_use = viridis_magma_dark_high, features = c("Cd163","Ccl8","Lilra5"),na_cutoff = 0.01, num_columns = 3)

###### boxplot for gene set scores ######
metadata = smerge_hm@meta.data
metadata$group = factor(metadata$group, levels=c("wt","d3a_ko","tet_ko"))

var_list <- c("group","DLBCL_CH_Up","DLBCL_CH_Dn","DLBCL_CH_Up_CAPS") 

selected_var <- subset(metadata,celltype %in% c("Resident-like_Mac") & sample != "wt")[var_list]
long_data <- melt(as.data.frame(selected_var), id.vars="group", variable.name = "Variable", value.name = "Zscore")
ggplot(subset(long_data,group != "NA"), aes(x = group, y = Zscore, fill = group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(shape = 16, position = position_jitter(0.2), size = 0, alpha = 0.6) +
    facet_wrap(~Variable, scales = "free") +
    facet_wrap(~Variable, scales = "free", ncol = 5) +
    theme_pubr() +
    labs(title = "",
         x = "",
         y = "Z-score") +
    scale_fill_manual(values = c("#96969680", "#1B9E7770","#B2DF8A70"))

t.test(subset(metadata, celltype %in% c("Resident-like_Mac") & group == "d3a_ko")$DLBCL_CH_Up, subset(metadata, celltype %in% c("Resident-like_Mac") & group == "wt")$DLBCL_CH_Up)
# p-value = 1.912e-10
t.test(subset(metadata, celltype %in% c("Resident-like_Mac") & group == "tet_ko")$DLBCL_CH_Up, subset(metadata, celltype %in% c("Resident-like_Mac") & group == "wt")$DLBCL_CH_Up)
# p-value = 2.563e-12

t.test(subset(metadata, celltype %in% c("Resident-like_Mac") & group == "d3a_ko")$DLBCL_CH_Up_CAPS, subset(metadata, celltype %in% c("Resident-like_Mac") & group == "wt")$DLBCL_CH_Up_CAPS)
# p-value < 2.2e-16
t.test(subset(metadata, celltype %in% c("Resident-like_Mac") & group == "tet_ko")$DLBCL_CH_Up_CAPS, subset(metadata, celltype %in% c("Resident-like_Mac") & group == "wt")$DLBCL_CH_Up_CAPS)
# p-value < 2.2e-16



###### cellChat differential analysis ######

## cellChat object
sce_d3a_ko = subset(smerge, group == "d3a_ko")
sce_tet_ko = subset(smerge, group == "tet_ko")
sce_wt = subset(smerge, group %in% c("wt"))

cellchat_d3a_ko <- createCellChat(object = sce_d3a_ko, group.by = "celltype", assay = "RNA")
cellchat_tet_ko <- createCellChat(object = sce_tet_ko, group.by = "celltype", assay = "RNA")
cellchat_wt <- createCellChat(object = sce_wt, group.by = "celltype", assay = "RNA")


# cellchat Db
CellChatDB <- CellChatDB.mouse 
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") 

# run cellchat for each sample, cellchat_wt
# set the used database in the object
cellchat_wt@DB <- CellChatDB.use
cellchat_wt <- subsetData(cellchat_wt)
future::plan("multisession", workers = 4) # do parallel
cellchat_wt <- identifyOverExpressedGenes(cellchat_wt)
cellchat_wt <- identifyOverExpressedInteractions(cellchat_wt)
cellchat_wt <- computeCommunProb(cellchat_wt, type = "triMean")
cellchat_wt <- filterCommunication(cellchat_wt, min.cells = 10)
cellchat_wt <- computeCommunProbPathway(cellchat_wt)
cellchat_wt <- aggregateNet(cellchat_wt)

# significant interaction
cellchat_wt@netP$pathways


# run cellchat for each sample, cellchat_d3a_ko
# set the used database in the object
cellchat_d3a_ko@DB <- CellChatDB.use
cellchat_d3a_ko <- subsetData(cellchat_d3a_ko) 
future::plan("multisession", workers = 4)
cellchat_d3a_ko <- identifyOverExpressedGenes(cellchat_d3a_ko)
cellchat_d3a_ko <- identifyOverExpressedInteractions(cellchat_d3a_ko)
cellchat_d3a_ko <- computeCommunProb(cellchat_d3a_ko, type = "triMean")
cellchat_d3a_ko <- filterCommunication(cellchat_d3a_ko, min.cells = 10)
cellchat_d3a_ko <- computeCommunProbPathway(cellchat_d3a_ko)
cellchat_d3a_ko <- aggregateNet(cellchat_d3a_ko)

# significant interaction
cellchat_d3a_ko@netP$pathways


# run cellchat for each sample, cellchat_tet_ko
# set the used database in the object
cellchat_tet_ko@DB <- CellChatDB.use
cellchat_tet_ko <- subsetData(cellchat_tet_ko)
future::plan("multisession", workers = 4)
cellchat_tet_ko <- identifyOverExpressedGenes(cellchat_tet_ko)
cellchat_tet_ko <- identifyOverExpressedInteractions(cellchat_tet_ko)
cellchat_tet_ko <- computeCommunProb(cellchat_tet_ko, type = "triMean")
cellchat_tet_ko <- filterCommunication(cellchat_tet_ko, min.cells = 10)
cellchat_tet_ko <- computeCommunProbPathway(cellchat_tet_ko)
cellchat_tet_ko <- aggregateNet(cellchat_tet_ko)

# significant interaction
cellchat_tet_ko@netP$pathways



# Merge analysis, Dnmt3a-ko vs WT
cell_color = c("#984EA3","#E41A1C","#54B0E4","#1B9E77","#BC9DCC","#222F75","#377EB8","#A65628","#F781BF","#F29403","#4DAF4A")

object_list_d3a <- list(WT = cellchat_wt, D3A_KO = cellchat_d3a_ko)
cellchat_d3a_wt <- mergeCellChat(object_list_d3a, add.names = names(object_list_d3a), cell.prefix = T)


# circle plot 
netVisual_diffInteraction(cellchat_d3a_wt, weight.scale = T, measure = "weight")

# heatmap
gg1 <- netVisual_heatmap(cellchat_d3a_wt,row.show = c("B","DC","Fibroblasts","Inflam_Mac","Mac","Mixed_lymphocytes","Mono","NKT","T","TREM2hi_Mac","Resident-like_Mac"),col.show = c("B","DC","Fibroblasts","Inflam_Mac","Mac","Mixed_lymphocytes","Mono","NKT","T","TREM2hi_Mac","Resident-like_Mac"),color.use = cell_color)
gg2 <- netVisual_heatmap(cellchat_d3a_wt, measure = "weight",row.show = c("B","DC","Fibroblasts","Inflam_Mac","Mac","Mixed_lymphocytes","Mono","NKT","T","TREM2hi_Mac","Resident-like_Mac"),col.show = c("B","DC","Fibroblasts","Inflam_Mac","Mac","Mixed_lymphocytes","Mono","NKT","T","TREM2hi_Mac","Resident-like_Mac"),color.use = cell_color)
gg1 + gg2

# rankNet plot
gg1 <- rankNet(cellchat_d3a_wt, mode = "comparison", stacked = T, do.stat = TRUE,color.use = c("#969696","#1B9E77"))
gg2 <- rankNet(cellchat_d3a_wt, mode = "comparison", stacked = F, do.stat = TRUE,color.use = c("#969696","#1B9E77"))
gg1$data ## values in the plots


# Merge analysis, Tet2-ko vs WT
object_list_tet <- list(WT = cellchat_wt, TET2_KO = cellchat_tet_ko)
cellchat_tet_wt <- mergeCellChat(object_list_tet, add.names = names(object_list_tet), cell.prefix = T)


# circle plot 
netVisual_diffInteraction(cellchat_tet_wt, weight.scale = T, measure = "weight")

# heatmap
gg1 <- netVisual_heatmap(cellchat_tet_wt,row.show = c("B","DC","Fibroblasts","Inflam_Mac","Mac","Mixed_lymphocytes","Mono","NKT","T","TREM2hi_Mac","Resident-like_Mac"),col.show = c("B","DC","Fibroblasts","Inflam_Mac","Mac","Mixed_lymphocytes","Mono","NKT","T","TREM2hi_Mac","Resident-like_Mac"),color.use = cell_color)
gg2 <- netVisual_heatmap(cellchat_tet_wt, measure = "weight",row.show = c("B","DC","Fibroblasts","Inflam_Mac","Mac","Mixed_lymphocytes","Mono","NKT","T","TREM2hi_Mac","Resident-like_Mac"),col.show = c("B","DC","Fibroblasts","Inflam_Mac","Mac","Mixed_lymphocytes","Mono","NKT","T","TREM2hi_Mac","Resident-like_Mac"),color.use = cell_color)
gg1 + gg2

# rankNet plot
gg1 <- rankNet(cellchat_tet_wt, mode = "comparison", stacked = T, do.stat = TRUE,color.use = c("#969696","#B2DF8A"))
gg2 <- rankNet(cellchat_tet_wt, mode = "comparison", stacked = F, do.stat = TRUE,color.use = c("#969696","#B2DF8A"))
gg1$data ## values in the plots


##############################################################################################################################
# DLBCL dataset by Wang et al.
##############################################################################################################################

###### seurat object and QC ######

# raw gene expression results, cellranger outputs, from HRA007235
s01 = Read10X(data.dir="/S01")
s02 = Read10X(data.dir="/S02")
s08 = Read10X(data.dir="/S08")
s16 = Read10X(data.dir="/S16")


# create seurat object
count_list <- list(
  s02 = s02,
  s01 = s01,
  s08 = s08,
  s16 = s16

)

seurat_objects <- list()
seurat_objects <- lapply(names(count_list), function(sample_name) {
    CreateSeuratObject(
        counts = count_list[[sample_name]],
        project = sample_name,
        min.cells = 5,
        min.features = 500
    )
})

names(seurat_objects) <- names(count_list)
original_names <- names(seurat_objects)


# QC
## add MT percentage
seurat_objects <- lapply(seurat_objects, function(object) {
  object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
  return(object)
})

# filter
seurat_objects <- lapply(seurat_objects, function(object) {
  object = subset(object, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA < 30000 & nCount_RNA > 1000)
  return(object)
})


# add meta.data
seurat_objects <- lapply(seurat_objects, function(object) {
  sample_name = object@project.name

  object$disease <- "DBLCL"
  object$sample <- sample_name

  return(object)
})


# normalization
seurat_objects <- lapply(seurat_objects, function(object) {
  object = NormalizeData(object)
  return(object)
})


# merge objects of different samples

# add sample identity
seurat_objects <- lapply(names(seurat_objects), function(i) {
  RenameCells(seurat_objects[[i]], add.cell.id = i)
})
names(seurat_objects) <- original_names

smerge <- Reduce(function(x, y) merge(x, y), seurat_objects)
rm(seurat_objects)



###### UMAP clustering ######

# clustering, with harmony correction
smerge = FindVariableFeatures(smerge, selection.method = "vst")
smerge = ScaleData(smerge, verbose = FALSE)
smerge = RunPCA(smerge, verbose = FALSE, npcs = 50)
ElbowPlot(object = smerge, ndims = 40)

smerge = RunHarmony(smerge, "sample")
smerge = FindNeighbors(smerge, reduction = "harmony", dims = 1:20)
smerge = FindClusters(smerge)
smerge = FindClusters(smerge,resolution =0.5)
smerge = RunUMAP(smerge, reduction = "harmony", dims = 1:20)

# UMAP plots
DimPlot(smerge, reduction = "umap", pt.size = 0.5, group.by="seurat_clusters", label = TRUE, label.size = 4) + theme_few()
DimPlot(smerge, reduction = "umap", pt.size = 0.5, group.by="sample", label = TRUE, label.size = 4) + theme_few()



###### SingleR annotation ######

# v5 is not supported, needs to be converted to v3 before SingleR
smerge[["RNA"]] <- as(smerge[["RNA"]], Class="Assay")

# main annotation
ref <- fetchReference("blueprint_encode", "2024-02-26")
pred.BlueprintEncodeData <- SingleR(test = smerge@assays$RNA@data, ref = ref, labels = ref$label.main)
smerge@meta.data$CellType_BlueprintEncodeData <- pred.BlueprintEncodeData$labels


# cell types across seurat_clusters
DimPlot(smerge, reduction = "umap", pt.size = 0.5, group.by="CellType_BlueprintEncodeData", label = TRUE, label.size = 4) + theme_few()
table(smerge$seurat_clusters, smerge$CellType_BlueprintEncodeData)

# add cell types information based on SingleR results
smerge <- RenameIdents(object = smerge,
'0'='T CD8',
'1'='B',
'2'='T CD4',
'3'='B',
'4'='B',
'5'='T CD8',
'6'='B',
'7'='T CD4',
'8'='B',
'9'='T CD4',
'10'='T CD8',
'11'='B',
'12'='MonoMac',
'13'='B',
'14'='NK',
'15'='B',
'16'='FIB',
'17'='T CD4',
'18'='DC',
'19'='B',
'20'='pDC',
'21'='B'
)
smerge$major_cell_name = Idents(smerge)
Idents(smerge) = smerge$RNA_snn_res.0.5

# UMAP by cell types
DimPlot(smerge, reduction = "umap", pt.size = 0.5, group.by="major_cell_name", label = TRUE, label.size = 4) + theme_few()



###### geneset analysis ######

counts <- LayerData(smerge, assay = "RNA", layer = "counts")
expressed_gene = rownames(counts)
rm(counts)

gset_cur <- getGmt("gene_sets.gmt")
gset_cur <- subsetGeneSets(gset_cur, expressed_gene) 
smerge <- AddModuleScore_UCell(smerge, features=geneIds(gset_cur), name=NULL)



###### FeaturePlot ######

FeaturePlot_scCustom(smerge, colors_use = viridis_magma_dark_high, features = c("DLBCL_CH_Up","DLBCL_CH_Dn","DLBCL_CH_Up_CAPS"),na_cutoff = 0.01, num_columns = 3)

FeaturePlot_scCustom(smerge, colors_use = viridis_magma_dark_high, features = c("CD163","CCL8","S100A12","LILRA5"),na_cutoff = 0.01, num_columns = 3)




###### ecotyper ######

# export expression matrics for ecotyper analysis
df = LayerData(subset(smerge, major_cell_name == "MonoMac"), assay = "RNA", layer = "data")
df = as.data.frame(df)
write.table(df,"./ecotyper/normalizedata.MonoMac.txt", sep="\t", quote=F)
write.table(subset(smerge, major_cell_name == "MonoMac")@meta.data,"./ecotyper/metadata.MonoMac.txt", sep="\t", quote=F)


cellstat = read.table("./ecotyper/Monocytes.and.Macrophages/state_assignment.txt",header = T)
rownames(cellstat) = cellstat$ID
smerge = AddMetaData(object =smerge, cellstat)

metadata = smerge@meta.data

var_list <- c("State","DLBCL_CH_Up_CAPS")  

long_data <- melt(as.data.frame(selected_var), id.vars="State", variable.name = "Variable", value.name = "Score")

ggplot(subset(long_data,State != "NA"), aes(x = State, y = Score, fill = State)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(shape = 16, position = position_jitter(0.2), size = 0.8, alpha = 0.6) +
    facet_wrap(~Variable, scales = "free") +
    facet_wrap(~Variable, scales = "free", ncol = 5) +
    theme_pubr() + ylim(0,0.18) +
    labs(title = "",
         x = "",
         y = "Z-score") +
    scale_fill_manual(values = c("#55A0FB80", "#FFA04070","#05BE7870"))



###### cellChat differential analysis ######

# define CAPS+ vs CAPS- MonoMac
genelist = c("CD163","CCL8","S100A12","DLBCL_CH_Up_CAPS")
ft_gene = FetchData(smerge, vars = c(genelist,"major_cell_name")) 

ft_gene <- ft_gene %>%
  mutate(CellType_CAPS = case_when(
    major_cell_name == "MonoMac" & DLBCL_CH_Up_CAPS > 0 ~ "CAPS_pos",
    major_cell_name == "MonoMac" & DLBCL_CH_Up_CAPS <= 0 ~ "CAPS_neg",
    TRUE ~ as.character(major_cell_name)
  ))

smerge = AddMetaData(object =smerge, ft_gene[,"CellType_CAPS")])

## cellChat object
sce_caps_pos = subset(smerge_sub, CellType3 %in% c("B","CAPS_pos","DC","FIB","NK","pDC","T CD4","T CD8"))
sce_caps_neg = subset(smerge_sub, CellType3 %in% c("B","CAPS_neg","DC","FIB","NK","pDC","T CD4","T CD8"))

cellchat_pos <- createCellChat(object = sce_caps_pos, group.by = "major_cell_name", assay = "RNA")
cellchat_neg <- createCellChat(object = sce_caps_neg, group.by = "major_cell_name", assay = "RNA")

# cellchat Db
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# use Secreted Signaling of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling

# run cellchat for each sample, positive

# set the used database in the object
cellchat_pos@DB <- CellChatDB.use

cellchat_pos <- subsetData(cellchat_pos) 
future::plan("multisession", workers = 4) 
cellchat_pos <- identifyOverExpressedGenes(cellchat_pos)
cellchat_pos <- identifyOverExpressedInteractions(cellchat_pos)
cellchat_pos <- computeCommunProb(cellchat_pos, type = "triMean")
cellchat_pos <- filterCommunication(cellchat_pos, min.cells = 10)
cellchat_pos <- computeCommunProbPathway(cellchat_pos)
cellchat_pos <- aggregateNet(cellchat_pos)

# significant interaction
cellchat_pos@netP$pathways

# run cellchat for each sample, negative

# set the used database in the object
cellchat_neg@DB <- CellChatDB.use

cellchat_neg <- subsetData(cellchat_neg) 
future::plan("multisession", workers = 4) 
cellchat_neg <- identifyOverExpressedGenes(cellchat_neg)
cellchat_neg <- identifyOverExpressedInteractions(cellchat_neg)
cellchat_neg <- computeCommunProb(cellchat_neg, type = "triMean")
cellchat_neg <- filterCommunication(cellchat_neg, min.cells = 10)
cellchat_neg <- computeCommunProbPathway(cellchat_neg)
cellchat_neg <- aggregateNet(cellchat_neg)

# significant interaction
cellchat_neg@netP$pathways

# Merge analysis

object.list <- list(CAPS_neg = cellchat_neg, CAPS_pos = cellchat_pos)
cellchat_com <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = T)

# circle plot 
netVisual_diffInteraction(cellchat_com, weight.scale = T, measure = "weight")

# rankNet plot
gg1 <- rankNet(cellchat_com, mode = "comparison", stacked = T, do.stat = TRUE, color.use = c("#377EB8","#4DAF4A"))
gg2 <- rankNet(cellchat_com, mode = "comparison", stacked = F, do.stat = TRUE, color.use = c("#377EB8","#4DAF4A"))
gg1 + gg2

gg1$data  # values from rankNet plots



##############################################################################################################################
# DLBCL dataset by Steen et al.
##############################################################################################################################

###### seurat object and QC ######
# raw matrix
rawmt = as.matrix(read.table(gzfile("GSE182434_raw_count_matrix.txt.gz"),header=TRUE, row.names=1)) # obtain from GEO database under accession number GSE182434
# meta.data
metd = read.table(gzfile("GSE182434_cell_annotation.txt.gz"),header=TRUE, row.names=1)

# seurat obj
scedlbc = CreateSeuratObject(counts = rawmt, project = "dlbc", meta.data =metd, min.cells = 5, min.features =200)
# only include DLBCL samples in the analysis
scedlbc = subset(scedlbc, Patient %in% c("DLBCL002","DLBCL007","DLBCL008","DLBCL111"))


# fitler
scedlbc[["percent.mt"]] = PercentageFeatureSet(scedlbc, pattern = "^MT-")
scedlbc = subset(scedlbc, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 10 & nCount_RNA < 20000 & nCount_RNA > 1000)

# normalization 
scedlbc <- NormalizeData(scedlbc)



###### UMAP clustering ######
# clustering, with harmony correction
scedlbc = FindVariableFeatures(scedlbc, selection.method = "vst")
scedlbc = ScaleData(scedlbc, verbose = FALSE)
scedlbc = RunPCA(scedlbc, npcs = 30, verbose = FALSE)
ElbowPlot(object = scedlbc)

scedlbc = RunHarmony(scedlbc, "Patient")
scedlbc = FindNeighbors(scedlbc, reduction = "harmony", dims = 1:20)
scedlbc = FindClusters(scedlbc)
scedlbc = RunUMAP(scedlbc, reduction = "harmony", dims = 1:20)

# dimplot, use cell type annotations from the original study
DimPlot(scedlbc, reduction = "umap", label = TRUE,  group.by = "CellType",label.size = 5)



###### GSEA analysis ######

counts <- LayerData(scedlbc, assay = "RNA", layer = "counts")
expressed_gene = rownames(counts)
rm(counts)

gset_cur <- getGmt("gene_sets.gmt")
gset_cur <- subsetGeneSets(gset_cur, expressed_gene) 
scedlbc <- AddModuleScore_UCell(scedlbc, features=geneIds(gset_cur), name=NULL)



###### FeaturePlot ######
# define a cut-off for low expressor with na_cutoff
FeaturePlot_scCustom(scedlbc, colors_use = viridis_magma_dark_high, features = c("DLBCL_CH_Up","DLBCL_CH_Dn","DLBCL_CH_Up_CAPS"),na_cutoff = 0.01, num_columns = 3)
FeaturePlot_scCustom(scedlbc, colors_use = viridis_magma_dark_high, features = c("CD163","CCL8","S100A12","LILRA5"),na_cutoff = 0.01, num_columns = 3)



###### ecotyper ######

# export expression matrics for ecotyper analysis
df = LayerData(subset(scedlbc, CellType %in% c("Monocytes and Macrophages")), assay = "RNA", layer = "data")
df = as.data.frame(df)
write.table(df,"./ecotyper/scedlbc.normalizedata.Mye.txt", sep="\t", quote=F)
write.table(subset(scedlbc, CellType %in% c("Monocytes and Macrophages"))@meta.data,"./ecotyper/scedlbc.metadata.Mye.txt", sep="\t", quote=F)

# import cell states from ecotyper output resuls
cellstat = read.table("./ecotyper/Monocytes.and.Macrophages/state_assignment.txt",header = T)
rownames(cellstat) = cellstat$ID
scedlbc = AddMetaData(object =scedlbc, cellstat)


# boxplot for CAPS score
metadata = scedlbc@meta.data

var_list <- c("State","DLBCL_CH_Up_CAPS")  
selected_var <- subset(metadata,CellType %in% c("Monocytes and Macrophages"))[var_list]

long_data <- melt(as.data.frame(selected_var), id.vars="State", variable.name = "Variable", value.name = "Score")
ggplot(subset(long_data,State != "NA"), aes(x = State, y = Score, fill = State)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(shape = 16, position = position_jitter(0.2), size = 0.8, alpha = 0.6) +
    facet_wrap(~Variable, scales = "free") +
    facet_wrap(~Variable, scales = "free", ncol = 5) +
    theme_pubr() + ylim(0,0.18) +
    labs(title = "",
         x = "",
         y = "Z-score") +
    scale_fill_manual(values = c("#55A0FB80", "#FFA04070","#05BE7870"))


t.test(subset(metadata, CellType %in% c("Monocytes and Macrophages") & State == "S03")$DLBCL_CH_Up_CAPS, subset(metadata, CellType %in% c("Monocytes and Macrophages") & State == "S01")$DLBCL_CH_Up_CAPS)
# p-value = 0.003553
t.test(subset(metadata, CellType %in% c("Monocytes and Macrophages") & State == "S03")$DLBCL_CH_Up_CAPS, subset(metadata, CellType %in% c("Monocytes and Macrophages") & State == "S02")$DLBCL_CH_Up_CAPS)
# p-value = 0.007388



###### cellChat differential analysis ######

# define CAPS+ vs CAPS- Monocytes and Macrophages
genelist = c("CD163","CCL8","S100A12","DLBCL_CH_Up_CAPS")
ft_gene = FetchData(scedlbc, vars = c(genelist,"CellType")) 

ft_gene <- ft_gene %>%
  mutate(CellType_CAPS = case_when(
    CellType == "Monocytes and Macrophages" & DLBCL_CH_Up_CAPS > 0 ~ "CAPS_pos",
    CellType == "Monocytes and Macrophages" & DLBCL_CH_Up_CAPS <= 0 ~ "CAPS_neg",
    TRUE ~ as.character(CellType)
  ))

scedlbc = AddMetaData(object =scedlbc, ft_gene[,"CellType_CAPS")])

## cellChat object
sce_caps_pos = subset(scedlbc, CellType_CAPS %in% c("B cells","CH_Up_pos","NK cells","Plasma cells","T cells CD4","T cells CD8","TFH","Tregs"))
sce_caps_neg = subset(scedlbc, CellType_CAPS %in% c("B cells","CH_Up_neg","NK cells","Plasma cells","T cells CD4","T cells CD8","TFH","Tregs"))

sce_caps_pos$CellType = factor(sce_caps_pos$CellType, levels=c("T cells CD8","B cells","T cells CD4","Monocytes and Macrophages","NK cells","Plasma cells","TFH","Tregs"))
sce_caps_neg$CellType = factor(sce_caps_neg$CellType, levels=c("T cells CD8","B cells","T cells CD4","Monocytes and Macrophages","NK cells","Plasma cells","TFH","Tregs"))

cellchat_pos <- createCellChat(object = sce_caps_pos, group.by = "CellType", assay = "RNA")
cellchat_neg <- createCellChat(object = sce_caps_neg, group.by = "CellType", assay = "RNA")


# cellchat Db
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# use Secreted Signaling of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling

# run cellchat for each sample, positive
# set the used database in the object
cellchat_pos@DB <- CellChatDB.use
cellchat_pos <- subsetData(cellchat_pos) 
future::plan("multisession", workers = 4) 
cellchat_pos <- identifyOverExpressedGenes(cellchat_pos)
cellchat_pos <- identifyOverExpressedInteractions(cellchat_pos)
cellchat_pos <- computeCommunProb(cellchat_pos, type = "triMean")
cellchat_pos <- filterCommunication(cellchat_pos)
cellchat_pos <- computeCommunProbPathway(cellchat_pos)
cellchat_pos <- aggregateNet(cellchat_pos)



# run cellchat for each sample, negative
# set the used database in the object
cellchat_neg@DB <- CellChatDB.use
cellchat_neg <- subsetData(cellchat_neg) 
future::plan("multisession", workers = 4) 
cellchat_neg <- identifyOverExpressedGenes(cellchat_neg)
cellchat_neg <- identifyOverExpressedInteractions(cellchat_neg)
cellchat_neg <- computeCommunProb(cellchat_neg, type = "triMean")
cellchat_neg <- filterCommunication(cellchat_neg)
cellchat_neg <- computeCommunProbPathway(cellchat_neg)
cellchat_neg <- aggregateNet(cellchat_neg)


# Merge analysis

object.list <- list(CAPS_neg = cellchat_neg, CAPS_pos = cellchat_pos)
cellchat_com <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = T)

# circle plot 
netVisual_diffInteraction(cellchat_com, weight.scale = T, measure = "weight")

# rankNet plot
gg1 <- rankNet(cellchat_com, mode = "comparison", stacked = T, do.stat = TRUE, color.use = c("#377EB8","#4DAF4A"))
gg2 <- rankNet(cellchat_com, mode = "comparison", stacked = F, do.stat = TRUE, color.use = c("#377EB8","#4DAF4A"))
gg1 + gg2

gg1$data  # values from rankNet plots
