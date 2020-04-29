rm(list=ls())
library("reticulate")
use_condaenv(condaenv="Renv", conda="/home/xiao/anaconda3/bin/conda")

library(Seurat)
library(ggplot2)
library(pheatmap)

in.dir <- "merge_seurat"
out.dir <- "cluster_tumor"
if(!dir.exists(out.dir)){
  dir.create(out.dir)
}

#read the seurat obj
tumor.obj <- readRDS(file.path(in.dir,"merged_seurat_obj_tumor.rds"))

#normalization
tumor.obj <- NormalizeData(tumor.obj)
#identification of highly variable features (feature selection)
tumor.obj<- FindVariableFeatures(tumor.obj)
top10 <- head(VariableFeatures(tumor.obj))
plot1 <- VariableFeaturePlot(tumor.obj)
plot2 <- LabelPoints(plot=plot1,points=top10,repel=TRUE)
CombinePlots(plots=list(plot1,plot2))
dev.off()

#scaling the data for PCA
all.genes <- rownames(tumor.obj)
tumor.obj<- ScaleData(tumor.obj,features = all.genes)

#PCA
tumor.obj <- RunPCA(tumor.obj, features = VariableFeatures(object = tumor.obj))
VizDimLoadings(tumor.obj,dims = 1:2, reduction="pca")
DimPlot(tumor.obj,reduction="pca")
DimHeatmap(tumor.obj, dims =1:20, cells = 500, balanced = TRUE)


#determine the dimensionality of the data
tumor.obj <- JackStraw(tumor.obj, num.replicate = 100)
tumor.obj <- ScoreJackStraw(tumor.obj, dims = 1:20)
JackStrawPlot(tumor.obj, dims =1:20)
ElbowPlot(tumor.obj,ndims=50)


#cluster the cells
tumor.obj <- FindNeighbors(tumor.obj, dims = 1:20)
tumor.obj <- FindClusters(tumor.obj, resolution = 0.5)

tumor.obj <- RunTSNE(tumor.obj,dims=1:20)
p <- DimPlot(tumor.obj,reduction="tsne",group.by = "orig.ident")
ggsave(file.path(out.dir,"tumor_tsne.pdf"),p,width=5,height=4,device="pdf",useDingbats=FALSE)

saveRDS(tumor.obj,file.path(out.dir,"tumor_obj.rds"))
